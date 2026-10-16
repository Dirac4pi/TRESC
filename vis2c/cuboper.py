#!/home/lky/miniconda3/envs/vis2c/bin/python

'''
Perform mathematical operations on the scalar real-space functions of
the two systems.
author: Dirac4pi
env: vis2c
'''

import viskit as vk
import numpy as np
from rdkit import Chem
from pathlib import Path
from scipy.ndimage import map_coordinates
from alignmol import align_xyz_fragment
from mvbb import GridBox, aligned_probe_grid_box, intersect_grid_boxes, approximate_mvbb
import argparse
import sys

BOHR_TO_ANGSTROM = 0.529177210903

#------------------------------------------------------------------------------
def cub2xyz(cube_file: str, xyz_file: str | None = None) -> None:
  """
  Extract cube atom coordinates to an XYZ file in angstroms.
  cube_file: name of Gaussian cube format file
  xyz_file: name of generated xyz file
  Returns: None
  """

  cube_path = Path(cube_file)
  xyz_path = (
    Path(xyz_file)
    if xyz_file is not None
    else cube_path.with_suffix(".xyz")
  )

  periodic_table = Chem.GetPeriodicTable()
  with cube_path.open("r", encoding="utf-8") as source:
    for _ in range(2):
      if not source.readline():
        raise ValueError("Cube file is missing a comment line.")
    header = source.readline().split()
    if len(header) < 4:
      raise ValueError("Invalid cube atom-count/origin line.")
    try:
      atom_count = abs(int(header[0]))
    except ValueError as exc:
      raise ValueError("Invalid cube atom count.") from exc
    if atom_count == 0:
      raise ValueError("Cube file contains no atoms.")
    # Skip the three grid vectors.
    for _ in range(3):
      if len(source.readline().split()) < 4:
        raise ValueError("Invalid cube grid header.")
    atoms = []
    for atom_index in range(1, atom_count + 1):
      fields = source.readline().split()
      if len(fields) < 5:
        raise ValueError(f"Invalid cube atom record: {atom_index}.")
      try:
        atomic_number = int(fields[0])
        if not 1 <= atomic_number <= 118:
          raise ValueError(f"Unsupported atomic number: {atomic_number}.")
        symbol = periodic_table.GetElementSymbol(atomic_number)
        coordinates = [
          float(value.replace("D", "E")) * BOHR_TO_ANGSTROM
          for value in fields[2:5]
        ]
      except ValueError as exc:
        raise ValueError(f"Invalid cube atom record: {atom_index}.") from exc
      atoms.append((symbol, coordinates))
    # export XYZ
  with xyz_path.open("w", encoding="utf-8") as destination:
    destination.write(f"{atom_count}\n")
    destination.write(
      f"Extracted from {cube_path.name}; coordinates in angstroms\n"
    )
    for symbol, (x, y, z) in atoms:
      destination.write(
        f"{symbol:<2s} {x: .10f} {y: .10f} {z: .10f}\n"
      )
  return None

#------------------------------------------------------------------------------
def _load_cube_metadata(cube_file: str):
  """Read grid metadata not returned by vk.load_cube().

  Parameters:
    cube_file: Input Gaussian cube path; positive grid counts and one scalar
      value per grid point are required.

  Returns:
    Tuple (origin, axes, shape). Origin and the three row-wise grid-step
    vectors are in Bohr; shape is (n1, n2, n3).
  """
  with Path(cube_file).open("r", encoding="utf-8") as source:
    source.readline()
    source.readline()
    header = source.readline().split()
    if len(header) != 4 or int(header[0]) <= 0:
      raise ValueError("Only standard single-scalar Gaussian cubes are supported.")
    origin = np.array([float(v.replace("D", "E")) for v in header[1:4]])
    shape = []
    axes = []
    for _ in range(3):
      fields = source.readline().split()
      if len(fields) != 4 or int(fields[0]) < 2:
        raise ValueError("Each Gaussian cube grid dimension must have at least 2 points.")
      shape.append(int(fields[0]))
      axes.append([float(v.replace("D", "E")) for v in fields[1:4]])
    axes = np.asarray(axes, dtype=float)
    if not np.all(np.isfinite(origin)) or not np.all(np.isfinite(axes)):
      raise ValueError("Cube grid metadata must be finite.")
    if abs(np.linalg.det(axes)) < 1e-12:
      raise ValueError("Cube grid-step vectors are linearly dependent.")
  return origin, axes, tuple(shape)

#------------------------------------------------------------------------------
def _write_cube(output_cube: str, atoms, origin, axes, values,
                operation: str, ref_cube: str, probe_cube: str) -> str:
  """Write a Gaussian cube using reference geometry and grid.

  Parameters:
    output_cube: Destination path; the caller must prevent overwriting inputs.
    atoms: Atomic records returned by vk.load_cube(ref_cube), in Bohr.
    origin: Reference grid origin (3,), in Bohr.
    axes: Reference grid-step vectors (3, 3), rows in i, j, k order, in Bohr.
    values: Computed scalar values, shaped (n1, n2, n3).
    operation: Name of the grid operation for the first cube comment line.
    ref_cube: Original reference path for output provenance.
    probe_cube: Original probe path for output provenance.

  Returns:
    Path to the written cube file as a string.
  """
  path = Path(output_cube)
  with path.open("w", encoding="utf-8") as dest:
    dest.write(f"{operation}: {Path(probe_cube).name} and {Path(ref_cube).name}\n")
    dest.write("Probe aligned to reference; coordinates in Bohr\n")
    dest.write(f"{len(atoms):5d}" + "".join(f" {v: .8f}" for v in origin) + "\n")
    for count, vector in zip(values.shape, axes):
      dest.write(f"{count:5d}" + "".join(f" {v: .8f}" for v in vector) + "\n")
    for atom in atoms:
      dest.write(f"{int(atom['index']):5d} {atom['charge']: .8f}"
                 f" {atom['x']: .8f} {atom['y']: .8f} {atom['z']: .8f}\n")
    for i in range(values.shape[0]):
      for j in range(values.shape[1]):
        for start in range(0, values.shape[2], 6):
          line = values[i, j, start:start + 6]
          dest.write("".join(f" {v: .8E}" for v in line) + "\n")
  return str(path)

#------------------------------------------------------------------------------
def _build_output_grid(box, ref_axes, probe_axes):
  """Build an MVBB-oriented grid using the coarsest input step length.

  Parameters:
    box: BoundingBox from mvbb.approximate_mvbb(), in reference Bohr.
    ref_axes: Three reference grid-step vectors (rows), in Bohr.
    probe_axes: Three original probe grid-step vectors (rows), in Bohr.

  Returns:
    Output origin (3,), row-wise step vectors (3, 3), grid shape, and the
    requested maximum input step length. Integer counts are chosen so the
    outermost grid-point centers cover every MVBB vertex. The actual step
    length along each axis is no larger than the requested length.
  """
  spacing = max(np.linalg.norm(ref_axes, axis=1).max(),
                np.linalg.norm(probe_axes, axis=1).max())
  if not np.isfinite(spacing) or spacing <= 0:
    raise ValueError("Input grid spacing must be positive and finite.")
  extent = np.asarray(box.extent, dtype=float)
  if np.any(extent <= 0) or not np.all(np.isfinite(extent)):
    raise ValueError("MVBB extent must be positive and finite.")
  intervals = np.maximum(1, np.ceil(extent / spacing).astype(int))
  steps = extent / intervals
  axes = box.rotation.T * steps[:, None]
  origin = box.center - 0.5 * np.sum(extent[:, None] * box.rotation.T, axis=0)
  return origin, axes, tuple((intervals + 1).tolist()), spacing

#------------------------------------------------------------------------------
def operate_cubes(ref_cube: str, probe_cube: str, operation: str,
                  output_cube: str, chunk_size: int = 200000) -> str:
  """Align two scalar cubes, operate on the overlap, and export an MVBB cube.

  Parameters:
    ref_cube: Reference Gaussian cube path. Its geometry and Bohr coordinate
      frame define the output atomic records and spatial reference frame.
    probe_cube: Probe Gaussian cube path. Its geometry is aligned to ref_cube
      via interactive, 1-based fragment selections; coordinates are in Bohr.
    operation: One of difference, sum, average, positive_difference,
      negative_difference, absolute_difference. Difference means aligned
      probe minus reference; both input scalar fields must have matching units.
    output_cube: Destination cube path. Cannot equal either input; existing
      unrelated files at this path will be replaced.
    chunk_size: Maximum number of output grid points evaluated in one batch.
      Lower this to reduce peak temporary memory use.

  Returns:
    Path to the generated cube as a string. Its box is an approximate MVBB of
    the continuous sampling-domain intersection in the ref Bohr frame; grid
    steps are at most the largest original step. Points outside the overlap
    are 0.0 for visualization, not physically evaluated field values.

  Raises:
    ValueError: Invalid data, no meaningful 3D overlap, or no valid output
      grid points. No output cube is written in these cases.
  """
  operations = {"difference", "sum", "average", "positive_difference",
                "negative_difference", "absolute_difference"}
  if operation not in operations:
    raise ValueError(f"Unknown operation: {operation}; choose from {sorted(operations)}")
  if chunk_size < 1:
    raise ValueError("chunk_size must be positive.")
  ref_path = Path(ref_cube)
  probe_path = Path(probe_cube)
  output_path = Path(output_cube)
  if output_path.resolve() in (ref_path.resolve(), probe_path.resolve()):
    raise ValueError("Output cube must not overwrite an input cube.")
  if not ref_path.is_file() or not probe_path.is_file():
    raise FileNotFoundError("Reference or probe cube file does not exist.")

  ref_origin, ref_axes, ref_shape = _load_cube_metadata(str(ref_path))
  probe_origin, probe_axes, probe_shape = _load_cube_metadata(str(probe_path))
  ref_atoms, ref_grid = vk.load_cube(str(ref_path))
  probe_atoms, probe_grid = vk.load_cube(str(probe_path))
  if ref_grid.shape != ref_shape or probe_grid.shape != probe_shape:
    raise ValueError("vk.load_cube shape differs from cube header counts.")
  if len(ref_atoms) == 0 or len(probe_atoms) == 0:
    raise ValueError("Cube file contains no atomic geometry.")
  ref_values = ref_grid['value']
  probe_values = probe_grid['value']
  if not np.all(np.isfinite(ref_values)) or not np.all(np.isfinite(probe_values)):
    raise ValueError("Cube values must be finite.")
  ref_box = GridBox(ref_origin, ref_axes, ref_shape)
  probe_box = GridBox(probe_origin, probe_axes, probe_shape)

  ref_xyz = ref_path.with_suffix(".xyz")
  probe_xyz = probe_path.with_suffix(".xyz")
  if ref_xyz.resolve() == probe_xyz.resolve():
    raise ValueError(
      "Reference and probe cubes would produce the same XYZ path."
    )
  cub2xyz(str(ref_path), str(ref_xyz))
  cub2xyz(str(probe_path), str(probe_xyz))
  print(f"Reference cube: {ref_path}")
  print(f"Probe cube: {probe_path}")
  print(f"Reference XYZ: {ref_xyz}")
  print(f"Probe XYZ: {probe_xyz}")
  transform = align_xyz_fragment(str(ref_xyz), str(probe_xyz))

  rotation = np.asarray(transform[:3, :3], dtype=float)
  translation_bohr = np.asarray(transform[:3, 3], dtype=float) / BOHR_TO_ANGSTROM
  aligned_probe_box = aligned_probe_grid_box(probe_box, transform,
                                             bohr_to_angstrom=BOHR_TO_ANGSTROM)
  intersection = intersect_grid_boxes(ref_box, aligned_probe_box)
  if intersection is None:
    raise ValueError("Cube sampling domains have no full-dimensional overlap.")
  box = approximate_mvbb(intersection)
  origin, axes, shape, spacing = _build_output_grid(box, ref_axes, probe_axes)
  count = int(np.prod(shape))
  result = np.zeros(shape, dtype=float)
  valid_count = 0
  max_abs_boundary_value = 0.0
  ref_inverse = np.linalg.inv(ref_axes)
  probe_inverse = np.linalg.inv(probe_axes)
  ref_upper = np.array(ref_shape, dtype=float) - 1
  probe_upper = np.array(probe_shape, dtype=float) - 1
  tolerance = 1e-7

  for start in range(0, count, chunk_size):
    end = min(start + chunk_size, count)
    indices = np.column_stack(np.unravel_index(np.arange(start, end), shape))
    output_points = origin + indices @ axes
    original_probe_points = (output_points - translation_bohr) @ rotation
    ref_indices_float = (output_points - ref_origin) @ ref_inverse
    probe_indices_float = (original_probe_points - probe_origin) @ probe_inverse
    valid = (np.all((ref_indices_float >= -tolerance) &
                    (ref_indices_float <= ref_upper + tolerance), axis=1) &
             np.all((probe_indices_float >= -tolerance) &
                    (probe_indices_float <= probe_upper + tolerance), axis=1))
    if not np.any(valid):
      continue
    ref_positions = np.clip(ref_indices_float[valid], 0, ref_upper)
    probe_positions = np.clip(probe_indices_float[valid], 0, probe_upper)
    ref_sample = map_coordinates(ref_values, ref_positions.T, order=1,
                                 mode="nearest", prefilter=False)
    probe_sample = map_coordinates(probe_values, probe_positions.T, order=1,
                                   mode="nearest", prefilter=False)
    delta = probe_sample - ref_sample
    if operation == "difference":
      sampled = delta
    elif operation == "sum":
      sampled = probe_sample + ref_sample
    elif operation == "average":
      sampled = 0.5 * (probe_sample + ref_sample)
    elif operation == "positive_difference":
      sampled = np.maximum(delta, 0.0)
    elif operation == "negative_difference":
      sampled = np.maximum(-delta, 0.0)
    else:
      sampled = np.abs(delta)
    if not np.all(np.isfinite(sampled)):
      raise ValueError("Interpolation produced non-finite values.")
    values = np.zeros(end - start, dtype=float)
    values[valid] = sampled
    result.flat[start:end] = values
    valid_count += int(valid.sum())
    boundary = np.any((indices == 0) | (indices == np.array(shape) - 1), axis=1)
    if np.any(valid & boundary):
      max_abs_boundary_value = max(
        max_abs_boundary_value, float(np.max(np.abs(values[valid & boundary]))))

  if valid_count == 0:
    raise ValueError("No output grid point lies in both cube sampling domains; "
                     "choose a finer grid or regenerate the cubes.")
  voxel_volume = abs(np.linalg.det(axes))
  print(f"Intersection vertices: {len(intersection.vertices)}")
  print(f"Approximate MVBB extent (Bohr): {box.extent}")
  print(f"Requested maximum spacing: {spacing:.6g} Bohr")
  print(f"Output grid shape: {shape}; actual spacings: {np.linalg.norm(axes, axis=1)} Bohr")
  print(f"Valid output points: {valid_count}/{count} ({valid_count / count:.2%})")
  print("Outside-overlap points are visualization-only zeros.")
  print(f"Result integral on valid grid (truncated domain): {result.sum() * voxel_volume:.8g}")
  print(f"Result range: {result.min():.8g} to {result.max():.8g}")
  print(f"Max |value| on valid outer-grid boundary: {max_abs_boundary_value:.8g}")
  return _write_cube(str(output_path), ref_atoms, origin, axes, result,
                     operation, str(ref_path), str(probe_path))


#==============================================================================
def main(argv=None) -> None:
  """Parse one operation and two Gaussian cube paths, then write a result cube.

  Parameters:
    argv: Optional command-line argument list for testing. None reads sys.argv.

  Returns:
    None. Writes <probe-stem>_<operation>_<ref-stem>.cub in the current
    directory unless --output is supplied; prints the written file path.
  """
  class OperationParser(argparse.ArgumentParser):
    def error(self, message):
      self.print_usage(sys.stderr)
      self.exit(2, f"{self.prog}: error: {message}\n"
                   "Operations: --add, --diff, --avg, --posdiff, "
                   "--negdiff, --absdiff.\n"
                   "Use --help for details.\n")

  parser = OperationParser(
    description=(
      "Align the probe cube to the reference cube and operate "
      "on their overlapping region. "
      "All differences are defined as aligned probe minus reference."
    ),
    allow_abbrev=False,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=(
      "Definitions:\n"
      "  R = reference field\n"
      "  P = probe field after alignment to the reference\n"
      "  Both fields are sampled on the output grid before the operation.\n\n"
      "Operations (choose exactly one):\n"
      "  --add      P + R\n"
      "  --diff     P - R\n"
      "  --avg      (P + R) / 2\n"
      "  --posdiff  max(P - R, 0)\n"
      "  --negdiff  max(R - P, 0)\n"
      "             Nonnegative magnitude of the negative part of P - R.\n"
      "  --absdiff  abs(P - R)\n\n"
      "Examples:\n"
      "  cuboper.py --add ref.cub probe.cub\n"
      "  cuboper.py --diff ref.cub probe.cub "
      "-o probe_minus_ref.cub\n\n"
      "The first input is the reference; the second input is the probe.\n"
      "The inputs must be Gaussian .cub or .cube files, not XYZ files.\n"
      "The input fields must have compatible physical meanings and units.\n"
      "Fitting atom indices are requested interactively; enter mcs for "
      "automatic matching."
    )
  )
  operations = parser.add_mutually_exclusive_group(required=True)
  operations.add_argument("--add", dest="operation", action="store_const",
                          const="sum", help="Add aligned probe to reference.")
  operations.add_argument("--diff", dest="operation", action="store_const",
                          const="difference", help="Aligned probe minus reference.")
  operations.add_argument("--avg", dest="operation", action="store_const",
                          const="average", help="Average both aligned fields.")
  operations.add_argument("--posdiff", dest="operation", action="store_const",
                          const="positive_difference", help="Positive part of probe - ref.")
  operations.add_argument("--negdiff", dest="operation", action="store_const",
                          const="negative_difference", help="Positive part of ref - probe.")
  operations.add_argument("--absdiff", dest="operation", action="store_const",
                          const="absolute_difference", help="Absolute ref - probe.")
  parser.add_argument("ref_cube", type=Path, help="Reference .cub or .cube file.")
  parser.add_argument("probe_cube", type=Path, help="Probe .cub or .cube file.")
  parser.add_argument("-o", "--output", type=Path, help="Output cube path.")
  parser.add_argument("--chunk-size", type=int, default=200000,
                      help="Output-grid points processed per batch (default: 200000).")
  args = parser.parse_args(argv)

  for label, path in (("Reference", args.ref_cube), ("Probe", args.probe_cube)):
    if path.suffix.lower() not in (".cub", ".cube"):
      parser.error(f"{label} input must be a .cub or .cube file, not {path.name!r}.")
    if not path.is_file():
      parser.error(f"{label} cube does not exist: {path}.")
  if args.chunk_size < 1:
    parser.error("--chunk-size must be a positive integer.")
  if args.ref_cube.resolve() == args.probe_cube.resolve():
    parser.error("Reference and probe paths must differ for CLI use.")
  operation_labels = {
    "sum": "add",
    "difference": "diff",
    "average": "avg",
    "positive_difference": "posdiff",
    "negative_difference": "negdiff",
    "absolute_difference": "absdiff",
  }
  output = args.output or Path.cwd() / (
    f"{args.probe_cube.stem}_{operation_labels[args.operation]}_"
    f"{args.ref_cube.stem}.cub"
  )
  if output.suffix.lower() not in (".cub", ".cube"):
    parser.error("Output must have a .cub or .cube extension.")
  if output.resolve() in (args.ref_cube.resolve(), args.probe_cube.resolve()):
    parser.error("Output cube must not overwrite an input cube.")
  if output.exists():
    parser.error(f"Output already exists: {output}. Choose another -o path.")
  if not output.parent.is_dir():
    parser.error(f"Output directory does not exist: {output.parent}.")
  if not sys.stdin.isatty():
    parser.error("Interactive atom selection requires a terminal on stdin.")

  try:
    written = operate_cubes(str(args.ref_cube), str(args.probe_cube),
                            args.operation, str(output), args.chunk_size)
  except (OSError, ValueError, RuntimeError, EOFError) as exc:
    parser.exit(1, f"{parser.prog}: {type(exc).__name__}: {exc}\n"
                   "No valid result was produced. Use --help for usage.\n")
  print(f"Output cube: {written}")
  print(f"If you want to visualize it, use real3damp.py {written}")
  return None


if __name__ == "__main__":
  main()


