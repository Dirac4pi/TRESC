"""
Compute overlapping Gaussian cube domains and approximate bounding boxes.
author: Dirac4pi
env: vis2c
"""

from dataclasses import dataclass
from itertools import product
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull, HalfspaceIntersection, QhullError


@dataclass
class GridBox:
  """A cube grid-center domain; positions and steps are in Bohr."""

  origin: np.ndarray
  axes: np.ndarray
  shape: tuple[int, int, int]

  def __post_init__(self):
    self.origin = np.asarray(self.origin, dtype=float)
    self.axes = np.asarray(self.axes, dtype=float)
    self.shape = tuple(int(n) for n in self.shape)
    if self.origin.shape != (3,) or self.axes.shape != (3, 3):
      raise ValueError("origin must be (3,) and axes must be (3, 3).")
    if len(self.shape) != 3 or any(n < 2 for n in self.shape):
      raise ValueError("Each grid dimension needs at least two points.")
    if not np.all(np.isfinite(self.origin)) or not np.all(np.isfinite(self.axes)):
      raise ValueError("Grid origin and axes must be finite.")
    if np.linalg.cond(self.axes) > 1e12:
      raise ValueError("Grid axes are singular or poorly conditioned.")

  def vertices(self):
    """Return eight grid-center domain corners in physical coordinates."""
    limits = np.asarray(self.shape, dtype=float) - 1.0
    corners = np.asarray(list(product((0.0, 1.0), repeat=3))) * limits
    return self.origin + corners @ self.axes

  def halfspaces(self):
    """Return six inequalities [normal, offset] with normal @ x + offset <= 0."""
    inverse = np.linalg.inv(self.axes)
    normals = inverse.T
    upper = np.asarray(self.shape, dtype=float) - 1.0
    rows = []
    for axis in range(3):
      normal = normals[axis]
      rows.append(np.r_[-normal, normal @ self.origin])
      rows.append(np.r_[normal, -normal @ self.origin - upper[axis]])
    return np.asarray(rows)

  def contains(self, points, tolerance=1e-7):
    """Test physical points against continuous grid index bounds."""
    points = np.asarray(points, dtype=float)
    if points.shape[-1] != 3:
      raise ValueError("points must have a final dimension of 3.")
    indices = (points - self.origin) @ np.linalg.inv(self.axes)
    return np.all((indices >= -tolerance) &
                  (indices <= np.asarray(self.shape) - 1 + tolerance), axis=-1)


@dataclass
class IntersectionGeometry:
  """Continuous full-dimensional intersection in the reference frame, Bohr."""

  halfspaces: np.ndarray
  vertices: np.ndarray
  ref_box: GridBox
  aligned_probe_box: GridBox
  interior_point: np.ndarray
  inscribed_radius: float

  def contains(self, points, tolerance=1e-7):
    """Check whether physical points lie within both sampling domains."""
    return (self.ref_box.contains(points, tolerance) &
            self.aligned_probe_box.contains(points, tolerance))


@dataclass
class BoundingBox:
  """A conservative, approximately minimum-volume orthogonal bounding box."""

  center: np.ndarray
  rotation: np.ndarray
  extent: np.ndarray

  def vertices(self):
    """Return eight corners in reference Bohr coordinates."""
    offsets = np.asarray(list(product((-0.5, 0.5), repeat=3))) * self.extent
    return self.center + offsets @ self.rotation.T

  @property
  def volume(self):
    """Return the volume in Bohr cubed."""
    return float(np.prod(self.extent))


def cube_grid_box(cube_file):
  """Read only origin, three grid steps and counts from a Gaussian cube header.

  Parameters:
    cube_file: Path to a Gaussian cube with positive axis counts; its grid
      coordinates are assumed to be in Bohr.

  Returns:
    GridBox describing the convex hull of the sampled grid-point centers.
    This does not read orbital identifiers or volumetric values.
  """
  with open(cube_file, "r", encoding="utf-8") as source:
    if not source.readline() or not source.readline():
      raise ValueError("Cube file is missing its comment lines.")
    fields = source.readline().split()
    if len(fields) < 4:
      raise ValueError("Invalid cube atom-count/origin line.")
    try:
      int(fields[0])
      origin = np.array([float(v.replace("D", "E").replace("d", "E"))
                         for v in fields[1:4]])
      counts = []
      axes = []
      for _ in range(3):
        fields = source.readline().split()
        if len(fields) < 4:
          raise ValueError("Incomplete cube grid header.")
        count = int(fields[0])
        if count < 2:
          raise ValueError("Only positive grid counts >= 2 are supported.")
        counts.append(count)
        axes.append([float(v.replace("D", "E").replace("d", "E"))
                     for v in fields[1:4]])
    except (TypeError, ValueError) as exc:
      raise ValueError(f"Invalid Gaussian cube grid header: {cube_file}") from exc
  return GridBox(origin, np.asarray(axes), tuple(counts))


def aligned_probe_grid_box(probe_box, transform, bohr_to_angstrom=0.529177210903):
  """Move a probe sampling domain to the reference coordinate frame.

  Parameters:
    probe_box: Original probe GridBox, in Bohr.
    transform: RDKit probe-to-reference 4x4 rigid transform; its translation
      is in Angstrom because fitting XYZ coordinates were in Angstrom.
    bohr_to_angstrom: Length conversion factor in Angstrom per Bohr.

  Returns:
    GridBox with transformed origin and grid-step vectors in reference Bohr.
  """
  transform = np.asarray(transform, dtype=float)
  if transform.shape != (4, 4) or not np.all(np.isfinite(transform)):
    raise ValueError("transform must be a finite 4x4 matrix.")
  rotation = transform[:3, :3]
  if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6) or not np.isclose(
      np.linalg.det(rotation), 1.0, atol=1e-6):
    raise ValueError("transform must contain a proper rigid rotation.")
  if not np.allclose(transform[3], [0, 0, 0, 1], atol=1e-6):
    raise ValueError("Invalid homogeneous transform last row.")
  if not np.isfinite(bohr_to_angstrom) or bohr_to_angstrom <= 0:
    raise ValueError("bohr_to_angstrom must be positive.")
  origin = probe_box.origin @ rotation.T + transform[:3, 3] / bohr_to_angstrom
  axes = probe_box.axes @ rotation.T
  return GridBox(origin, axes, probe_box.shape)


def intersect_grid_boxes(ref_box, aligned_probe_box, radius_tolerance=1e-7):
  """Intersect two grid-center domains in the same Bohr coordinate frame.

  Parameters:
    ref_box: Unchanged reference GridBox.
    aligned_probe_box: Probe GridBox already transformed into the ref frame.
    radius_tolerance: Minimum Chebyshev radius in Bohr accepted as a
      full-dimensional intersection; touching/thin intersections return None.

  Returns:
    IntersectionGeometry with halfspaces, vertices, boxes and an interior
    point, or None if no sufficiently thick 3D intersection exists.
  """
  if radius_tolerance <= 0:
    raise ValueError("radius_tolerance must be positive.")
  halfspaces = np.vstack((ref_box.halfspaces(), aligned_probe_box.halfspaces()))
  normals = halfspaces[:, :3]
  offsets = halfspaces[:, 3]
  lengths = np.linalg.norm(normals, axis=1)
  unit_normals = normals / lengths[:, None]
  unit_offsets = offsets / lengths
  constraints = np.column_stack((unit_normals, np.ones(len(unit_normals))))
  solution = linprog([0.0, 0.0, 0.0, -1.0], A_ub=constraints,
                     b_ub=-unit_offsets,
                     bounds=[(None, None)] * 3 + [(0, None)], method="highs")
  if not solution.success or solution.x[3] <= radius_tolerance:
    return None
  center = solution.x[:3]
  if np.max(unit_normals @ center + unit_offsets) >= -radius_tolerance:
    return None
  try:
    vertices = HalfspaceIntersection(halfspaces, center).intersections
  except QhullError as exc:
    raise RuntimeError("Failed to compute intersection vertices.") from exc
  if not np.all(ref_box.contains(vertices, 1e-5) &
                aligned_probe_box.contains(vertices, 1e-5)):
    raise RuntimeError("Intersection vertices failed grid-domain validation.")
  return IntersectionGeometry(halfspaces, vertices, ref_box,
                              aligned_probe_box, center, float(solution.x[3]))


def cube_intersection(ref_cube, probe_cube, transform, radius_tolerance=1e-7):
  """Find the 3D overlap of reference and aligned probe cube domains.

  Parameters:
    ref_cube: Reference Gaussian cube path; its original Bohr frame is retained.
    probe_cube: Probe Gaussian cube path; its geometry is moved by transform.
    transform: RDKit probe-to-reference 4x4 transform obtained from XYZ files
      in Angstrom, for example from align_xyz_fragment().
    radius_tolerance: Minimum accepted interior-ball radius, in Bohr.

  Returns:
    IntersectionGeometry in reference Bohr coordinates, or None if the
    common sampled volume is absent or too thin. Does not write a new cube.
  """
  ref_box = cube_grid_box(ref_cube)
  probe_box = cube_grid_box(probe_cube)
  aligned_box = aligned_probe_grid_box(probe_box, transform)
  return intersect_grid_boxes(ref_box, aligned_box, radius_tolerance)


def _minimum_area_rectangle(points, tolerance):
  """Fit a two-dimensional minimum-area rectangle to projected 3D points."""
  hull = ConvexHull(points)
  polygon = points[hull.vertices]
  best_area = np.inf
  best = None
  for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
    edge = end - start
    length = np.linalg.norm(edge)
    if length <= tolerance:
      continue
    u = edge / length
    v = np.array([-u[1], u[0]])
    pu = points @ u
    pv = points @ v
    low_u, high_u = pu.min(), pu.max()
    low_v, high_v = pv.min(), pv.max()
    width = high_u - low_u
    height = high_v - low_v
    area = width * height
    if area < best_area:
      best_area = area
      best = (u, v, width, height,
              0.5 * (low_u + high_u), 0.5 * (low_v + high_v))
  if best is None:
    raise ValueError("Projected points have no valid 2D bounding rectangle.")
  return best


def approximate_mvbb(intersection, padding=1e-7):
  """Fit a conservative approximate minimum-volume orthogonal bounding box.

  Parameters:
    intersection: Full-dimensional IntersectionGeometry in reference Bohr
      coordinates, returned by intersect_grid_boxes() or cube_intersection().
    padding: Absolute margin in Bohr added on each side along every box axis.
      This accommodates roundoff; it is not the later grid-spacing margin.

  Returns:
    BoundingBox(center, rotation, extent), where rotation columns form a
    right-handed orthonormal frame, extent lists the three side lengths in
    Bohr, and all original intersection vertices lie inside the box.

  Notes:
    Candidate first-axis directions come from all 3D convex-hull edges;
    for each direction a 2D minimum-area rectangle determines the other axes.
    This is an approximate OBB, not a proven global 3D volume minimum.
  """
  if not np.isfinite(padding) or padding < 0:
    raise ValueError("padding must be finite and nonnegative.")
  vertices = np.asarray(intersection.vertices, dtype=float)
  if vertices.ndim != 2 or vertices.shape[1] != 3 or len(vertices) < 4:
    raise ValueError("At least four three-dimensional vertices are required.")
  span = np.linalg.norm(np.ptp(vertices, axis=0))
  if span <= 0 or not np.isfinite(span):
    raise ValueError("Intersection vertices have no spatial extent.")
  tolerance = max(1e-12 * span, 1e-12)
  try:
    hull = ConvexHull(vertices)
  except QhullError as exc:
    raise ValueError("Intersection vertices are not full-dimensional.") from exc
  edges = set()
  for triangle in hull.simplices:
    for i in range(3):
      a, b = sorted((int(triangle[i]), int(triangle[(i + 1) % 3])))
      if np.linalg.norm(vertices[a] - vertices[b]) > tolerance:
        edges.add((a, b))
  if not edges:
    raise ValueError("No non-degenerate hull edges found.")

  best_volume = np.inf
  best = None
  for i, j in sorted(edges):
    u = (vertices[j] - vertices[i]) / np.linalg.norm(vertices[j] - vertices[i])
    tmp = np.eye(3)[np.argmin(np.abs(u))]
    a = tmp - (tmp @ u) * u
    a /= np.linalg.norm(a)
    b = np.cross(u, a)
    projected = np.column_stack((vertices @ a, vertices @ b))
    try:
      v2, w2, width, height, _, _ = _minimum_area_rectangle(projected, tolerance)
    except (QhullError, ValueError):
      continue
    v = v2[0] * a + v2[1] * b
    w = w2[0] * a + w2[1] * b
    rotation = np.column_stack((u, v, w))
    if np.linalg.det(rotation) < 0:
      rotation[:, 2] *= -1
    extent_u = np.ptp(vertices @ u)
    volume = extent_u * width * height
    if volume < best_volume:
      best_volume = volume
      best = rotation
  if best is None:
    raise ValueError("Could not fit an oriented box to the intersection.")

  # Recompute bounds in the selected basis to contain unjoggled input points.
  projected = vertices @ best
  lower = projected.min(axis=0) - padding
  upper = projected.max(axis=0) + padding
  center = (0.5 * (lower + upper)) @ best.T
  extent = upper - lower
  box = BoundingBox(center, best, extent)
  reconstructed = (vertices - box.center) @ box.rotation
  if not np.all(np.abs(reconstructed) <= box.extent / 2 + tolerance):
    raise RuntimeError("Computed box does not enclose the intersection vertices.")
  return box
