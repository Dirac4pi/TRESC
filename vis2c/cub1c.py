#!/home/lky/miniconda3/envs/vis2c/bin/python
"""
Visualize a real scalar orbital or export a fixed-material OBJ bundle.
author: Dirac4pi
env: vis2c
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import uuid


# -----------------------------------------------------------------------------
# Input validation and orbital cube generation

def _hash(path):
  """Return a file's SHA-256 digest with bounded read buffers."""
  digest = hashlib.sha256()
  with Path(path).open('rb') as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b''):
      digest.update(block)
  return digest.hexdigest()


def _level(value):
  """Require a positive finite orbital-amplitude contour magnitude."""
  value = float(value)
  if not math.isfinite(value) or value <= 0:
    raise ValueError('isovalue must be positive and finite.')
  return value


def _index(value):
  """Canonicalize an explicitly supplied 1-based orbital/spin index."""
  value = str(value)
  if re.fullmatch(r'[1-9]\d*[AB]', value, re.IGNORECASE) is None:
    raise ValueError('index must include a spin label, for example 24A or 25B.')
  return value.upper()


def _molden_path(filename):
  """Resolve a full Molden filename, or append a supported suffix safely."""
  path = Path(filename).expanduser()
  if path.name.endswith(('.molden.input', '.molden')):
    candidates = [path]
  else:
    candidates = [Path(str(path) + '.molden.input'), Path(str(path) + '.molden')]
  for candidate in candidates:
    if candidate.is_file():
      return candidate.resolve()
  raise FileNotFoundError('No Molden file found: ' + ', '.join(map(str, candidates)))


def _prepare(molden, index, isovalue, *, progress=None):
  """Generate and validate the selected signed orbital cube in the CWD."""
  import numpy as np
  import viskit as vk

  progress = progress if progress is not None else sys.stdout
  path, index, level = _molden_path(molden), _index(index), _level(isovalue)
  cube = Path.cwd() / f'{index}-real.cub'
  before = None
  if cube.exists():
    stat = cube.stat()
    before = (stat.st_mtime_ns, stat.st_ctime_ns, stat.st_size)
  command = ['tshell.sh', '-cub1c', str(path), index]
  source_hash = _hash(path)
  print('Calling TRESC: ' + ' '.join(command), file=progress, flush=True)
  subprocess.run(command, check=True, stdout=progress, stderr=progress)
  if not cube.is_file() or cube.stat().st_size == 0:
    raise RuntimeError(f'TRESC did not produce a nonempty cube in the CWD: {cube}')
  stat = cube.stat()
  if before == (stat.st_mtime_ns, stat.st_ctime_ns, stat.st_size):
    raise RuntimeError(f'The existing cube was not updated; refusing stale data: {cube}')
  if _hash(path) != source_hash:
    raise RuntimeError('The Molden source changed during cube generation.')
  print(f'Loading signed orbital field from {cube}...', file=progress, flush=True)
  atoms, grid = vk.load_cube(str(cube))
  shape = np.asarray(grid['value']).shape
  if len(shape) != 3 or min(shape) < 2:
    raise ValueError('Cube values must have three dimensions, each >= 2.')
  for key in ('x', 'y', 'z', 'value'):
    data = np.asarray(grid[key])
    if data.shape != shape or not np.all(np.isfinite(data)):
      raise ValueError(f'Invalid shape or nonfinite cube {key} values.')
  return path, index, level, cube, atoms, grid, command, source_hash


# -----------------------------------------------------------------------------
def cub1c(molden: str, index: str, isovalue: float = 0.05,
          slice: bool = False) -> None:
  """Generate a real orbital cube and open the original isovalue controller.

  The existing positional arguments and optional slice behavior are preserved.
  The selected orbital must be supplied explicitly; no HOMO inference is added.
  """
  _, index, level, _, atoms, grid, _, _ = _prepare(molden, index, isovalue)
  from mayavi import mlab
  import plot3d as p3
  from traits.api import HasTraits, Float, observe
  from traitsui.api import View, Item

  x, y, z, values = (grid[key] for key in ('x', 'y', 'z', 'value'))
  print('Plotting...', flush=True)
  if slice:
    p3.slice_cub(
      'Orbit Isosurfase Slice', atoms, x, y, z, values,
      'viridis', values.max(), values.min(),
    )
  used = p3.real_orb_plot_cub(index, atoms, values, x, y, z, level)

  class IsoValueController(HasTraits):
    isov = Float(used, desc='isov', auto_set=False, enter_set=True)

    @observe('isov')
    def update_isov(self, event):
      try:
        new_level = _level(event.new)
      except ValueError as exc:
        print(f'Cannot update isovalue: {exc}')
        return
      mlab.close(index)
      p3.real_orb_plot_cub(index, atoms, values, x, y, z, new_level)
      mlab.draw()
      mlab.view()

    view = View(
      Item('isov', label='isovalue', show_label=True),
      width=300, height=200, resizable=True,
    )

  controller = IsoValueController()
  controller.configure_traits()
  return None


# -----------------------------------------------------------------------------
# Output paths and material extraction

def _output_path(filename):
  """Refuse existing bundle destinations and nonportable material basenames."""
  target = Path(filename).expanduser().absolute()
  if target.suffix != '.obj':
    raise ValueError('Export filename must end with lowercase .obj.')
  if any(c.isspace() for c in target.stem) or '#' in target.stem:
    raise ValueError('Use an OBJ basename without whitespace or #.')
  if not target.parent.is_dir():
    raise FileNotFoundError(f'Output parent directory does not exist: {target.parent}')
  for candidate in (
    target, target.with_suffix('.mtl'), target.with_suffix('.scene.json'),
  ):
    if os.path.lexists(candidate):
      raise FileExistsError(f'Refusing existing output: {candidate}')
  return target


def _actor_key(actor):
  """Identify a VTK actor by its underlying object address, not wrapper identity."""
  return actor.GetAddressAsString('')


def _orbital_materials(figure, used):
  """Read the two contour modules and their actual constant LUT RGBA values.

  Export-only material conversion leaves the Mayavi actors and LUTs unchanged.
  Never infer positive/negative roles from color or renderer actor order.
  """
  import numpy as np
  from tvtk.api import tvtk

  styles, visited = {}, set()
  stack = list(figure.children)
  while stack:
    node = stack.pop()
    if id(node) in visited:
      continue
    visited.add(id(node))
    stack.extend(list(getattr(node, 'children', ())))
    if type(node).__name__ != 'IsoSurface':
      continue
    contours = np.asarray(node.contour.contours, dtype=float)
    if contours.size != 1:
      continue
    contour = float(contours.ravel()[0])
    if math.isclose(contour, used, rel_tol=1e-8, abs_tol=0):
      role = 'orbital_positive'
    elif math.isclose(contour, -used, rel_tol=1e-8, abs_tol=0):
      role = 'orbital_negative'
    else:
      continue
    rgba = np.asarray(
      node.module_manager.scalar_lut_manager.lut.table.to_array(), dtype=float,
    )
    if (rgba.ndim != 2 or rgba.shape[1] != 4 or len(rgba) == 0
        or not np.all(np.isfinite(rgba)) or np.any(rgba < 0)
        or np.any(rgba > 255) or not np.all(rgba == rgba[0])):
      raise RuntimeError('Orbital LUT must contain one constant valid RGBA color.')
    actor = tvtk.to_vtk(node.actor.actor)
    opacity = float(actor.GetProperty().GetOpacity())
    if not 0 <= opacity <= 1:
      raise RuntimeError('Actor opacity is outside [0, 1].')
    styles[_actor_key(actor)] = {
      'role': role, 'contour': contour,
      'lut_rgba_255': rgba[0].tolist(),
      'color': (rgba[0, :3] / 255).tolist(),
      'opacity': float(opacity * rgba[0, 3] / 255),
      'actor_opacity': opacity,
    }
  roles = [style['role'] for style in styles.values()]
  if sorted(roles) != ['orbital_negative', 'orbital_positive']:
    raise RuntimeError(
      'Cannot uniquely identify positive/negative IsoSurface modules. '
      'The local plot3d pipeline needs a version-specific adapter.'
    )
  return styles


def _collect_actors(renderer, styles):
  """Collect fixed-material geometry, rejecting unhandled scalar/texture colors."""
  collection = renderer.GetActors()
  collection.InitTraversal()
  actors, nonempty_roles = [], []
  for index in range(collection.GetNumberOfItems()):
    actor = collection.GetNextActor()
    if not actor.GetVisibility():
      continue
    mapper = actor.GetMapper()
    if mapper is None:
      continue
    mapper.Update()
    data = mapper.GetInput()
    if data is None or data.GetNumberOfPoints() == 0:
      continue
    if not data.IsA('vtkPolyData'):
      raise RuntimeError('The selected renderer contains non-polydata geometry.')
    prop = actor.GetProperty()
    style = styles.get(_actor_key(actor))
    if actor.GetTexture() is not None:
      raise RuntimeError('Textured actors are not supported by this OBJ branch.')
    if style is None:
      if mapper.GetScalarVisibility():
        raise RuntimeError('An unclassified scalar-colored actor cannot be exported safely.')
      style = {
        'role': 'molecular_geometry', 'color': list(prop.GetColor()),
        'opacity': float(prop.GetOpacity()),
      }
    if style['role'].startswith('orbital_'):
      if data.GetNumberOfPolys() + data.GetNumberOfStrips() == 0:
        continue
      nonempty_roles.append(style['role'])
    material = dict(style)
    material.update({
      'renderer_actor_index': index,
      'ambient': float(prop.GetAmbient()),
      'diffuse': float(prop.GetDiffuse()),
      'specular': float(prop.GetSpecular()),
      'specular_power': float(prop.GetSpecularPower()),
      'specular_color': list(prop.GetSpecularColor()),
      'ambient_color': list(prop.GetAmbientColor()),
    })
    actors.append((actor, material))
  if not nonempty_roles:
    raise RuntimeError('No nonempty signed orbital surface; refusing skeleton-only output.')
  return actors, sorted(set(nonempty_roles))


# -----------------------------------------------------------------------------
# Explicit OBJ/MTL serialization

def _write_bundle(stage, stem, actors):
  """Write transformed meshes with stable groups and standard MTL dissolve values.

  A constant LUT becomes one material per actor. MTL uses d=opacity, never Tr.
  Source geometry is not smoothed or resampled here. Existing normals are kept;
  missing normals are generated only for mesh shading, not field smoothing.
  """
  import numpy as np
  from vtkmodules.util.numpy_support import vtk_to_numpy
  from vtkmodules.vtkCommonCore import vtkIdList
  from vtkmodules.vtkCommonTransforms import vtkTransform
  from vtkmodules.vtkFiltersCore import vtkPolyDataNormals, vtkTriangleFilter
  from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter

  obj_path, mtl_path = stage / f'{stem}.obj', stage / f'{stem}.mtl'
  records, offset, vertex_count, face_count = [], 0, 0, 0
  with obj_path.open('w', encoding='utf-8') as obj, \
       mtl_path.open('w', encoding='utf-8') as mtl:
    obj.write('# cub1c signed orbital scene; OBJ coordinates preserve source units\n')
    obj.write(f'mtllib {stem}.mtl\n')
    mtl.write('# Standard dissolve: d=0 is transparent; d=1 is opaque\n')
    for number, (actor, style) in enumerate(actors):
      transform = vtkTransform()
      transform.SetMatrix(actor.GetMatrix())
      transformed = vtkTransformPolyDataFilter()
      transformed.SetInputData(actor.GetMapper().GetInput())
      transformed.SetTransform(transform)
      triangulate = vtkTriangleFilter()
      triangulate.SetInputConnection(transformed.GetOutputPort())
      triangulate.PassLinesOff()
      triangulate.PassVertsOff()
      triangulate.Update()
      mesh = triangulate.GetOutput()
      if mesh.GetNumberOfPolys() == 0:
        if style['role'].startswith('orbital_'):
          raise RuntimeError('An orbital surface has no exportable polygon faces.')
        continue
      generated_normals = mesh.GetPointData().GetNormals() is None
      if generated_normals:
        normals_filter = vtkPolyDataNormals()
        normals_filter.SetInputData(mesh)
        normals_filter.SplittingOff()
        normals_filter.ConsistencyOn()
        normals_filter.Update()
        mesh = normals_filter.GetOutput()
      points = vtk_to_numpy(mesh.GetPoints().GetData())
      normals = vtk_to_numpy(mesh.GetPointData().GetNormals())
      if (points.shape != normals.shape or not np.all(np.isfinite(points))
          or not np.all(np.isfinite(normals))):
        raise ValueError('Mesh coordinates/normals are invalid.')
      name = f'{style["role"]}_{number:05d}'
      obj.write(f'o {name}\ng {name}\nusemtl {name}\n')
      for point in points:
        obj.write('v ' + ' '.join(f'{float(v):.12g}' for v in point) + '\n')
      for normal in normals:
        obj.write('vn ' + ' '.join(f'{float(v):.12g}' for v in normal) + '\n')
      cells, ids, count = mesh.GetPolys(), vtkIdList(), 0
      cells.InitTraversal()
      while cells.GetNextCell(ids):
        if ids.GetNumberOfIds() != 3:
          raise RuntimeError('Triangulation produced a non-triangle face.')
        indices = [int(ids.GetId(i)) + offset + 1 for i in range(3)]
        obj.write('f ' + ' '.join(f'{v}//{v}' for v in indices) + '\n')
        count += 1
      rgb = style['color']
      ambient = [c * style['ambient'] for c in style['ambient_color']]
      specular = [c * style['specular'] for c in style['specular_color']]
      mtl.write(f'newmtl {name}\n')
      for label, values in (('Ka', ambient), ('Kd', rgb), ('Ks', specular)):
        mtl.write(label + ' ' + ' '.join(f'{float(v):.9g}' for v in values) + '\n')
      mtl.write(f'Ns {max(0.0, min(1000.0, style["specular_power"])):.9g}\n')
      mtl.write(f'd {style["opacity"]:.9g}\nillum 2\n\n')
      record = dict(style)
      record.update({
        'obj_group': name, 'mtl_material': name,
        'vertices': len(points), 'triangles': count,
        'normals_generated': generated_normals,
      })
      records.append(record)
      offset += len(points)
      vertex_count += len(points)
      face_count += count
  if not vertex_count or not face_count:
    raise RuntimeError('No polygon geometry was exported.')
  return records, {'vertices': vertex_count, 'faces': face_count}


def _check_bundle(stage, stem):
  """Check geometry counts and exact OBJ-to-MTL references before publication."""
  obj, mtl = stage / f'{stem}.obj', stage / f'{stem}.mtl'
  for path in (obj, mtl):
    if not path.is_file() or path.stat().st_size == 0:
      raise RuntimeError(f'Missing/empty exported file: {path.name}')
  vertices, normals, faces, libraries, used, max_index = 0, 0, 0, [], set(), 0
  with obj.open(encoding='utf-8') as stream:
    for line in stream:
      fields = line.split()
      if not fields:
        continue
      if fields[0] == 'v':
        vertices += 1
      elif fields[0] == 'vn':
        normals += 1
      elif fields[0] == 'f':
        faces += 1
        for item in fields[1:]:
          vertex, normal = item.split('//')
          if vertex != normal or int(vertex) < 1:
            raise RuntimeError('Invalid vertex/normal indexing in OBJ.')
          max_index = max(max_index, int(vertex))
      elif fields[0] == 'mtllib':
        libraries.extend(fields[1:])
      elif fields[0] == 'usemtl':
        used.add(fields[1])
  definitions = set()
  with mtl.open(encoding='utf-8') as stream:
    for line in stream:
      if line.startswith('newmtl '):
        definitions.add(line.split()[1])
  if (not vertices or not faces or vertices != normals or max_index > vertices
      or libraries != [mtl.name] or not used or not used.issubset(definitions)):
    raise RuntimeError('OBJ/MTL bundle validation failed.')
  return {'vertices': vertices, 'normals': normals, 'faces': faces}


def _publish(stage, parent):
  """Publish without overwrite; roll back our new files on a collision/error."""
  published = []
  try:
    for source in sorted(stage.iterdir()):
      target = parent / source.name
      os.link(source, target)
      published.append(target)
  except Exception:
    for target in reversed(published):
      target.unlink()
    raise
  return published


# -----------------------------------------------------------------------------
# OBJ export API

def export_obj(molden: str, index: str, isovalue: float = 0.05, *, filename):
  """Export signed orbital surfaces and molecular geometry without a GUI loop.

  The normal backend is intentional: offscreen scene initialization crashed in
  the user's tested environment. A window may briefly appear; no input or GUI
  controller is required. Blender rendering remains a separate task.
  """
  import numpy as np
  target = _output_path(filename)
  _level(isovalue)
  _index(index)
  _molden_path(molden)
  from mayavi import mlab
  previous_offscreen, figure = mlab.options.offscreen, None
  title = f'{_index(index)} [obj-{uuid.uuid4().hex}]'
  try:
    mlab.options.offscreen = False
    path, index, level, cube, atoms, grid, command, source_hash = _prepare(
      molden, index, isovalue, progress=sys.stderr,
    )
    cube_hash = _hash(cube)
    import plot3d as p3
    from tvtk.api import tvtk
    x, y, z, values = (grid[key] for key in ('x', 'y', 'z', 'value'))
    print('Constructing signed orbital scene...', file=sys.stderr, flush=True)
    used = _level(p3.real_orb_plot_cub(title, atoms, values, x, y, z, level))
    figure = mlab.figure(title)
    figure.scene.show_axes = False
    figure.scene.render()
    styles = _orbital_materials(figure, used)
    actors, nonempty_roles = _collect_actors(tvtk.to_vtk(figure.scene.renderer), styles)
    with tempfile.TemporaryDirectory(prefix='.cub1c-obj-', dir=target.parent) as temporary:
      stage = Path(temporary)
      records, written = _write_bundle(stage, target.stem, actors)
      checked = _check_bundle(stage, target.stem)
      if written['vertices'] != checked['vertices'] or written['faces'] != checked['faces']:
        raise RuntimeError('Exported mesh counts do not match the validated OBJ.')
      if _hash(cube) != cube_hash or _hash(path) != source_hash:
        raise RuntimeError('Scientific input changed during export.')
      manifest = {
        'schema_version': 1,
        'scene_type': 'cub1c_signed_orbital',
        'input_molden': str(path), 'molden_sha256': source_hash,
        'orbital_index_with_spin': index,
        'orbital_number': int(index[:-1]), 'spin': index[-1],
        'orbital_selection': 'Explicit user-supplied index; HOMO/LUMO not inferred',
        'generator_command': command, 'generator_working_directory': str(Path.cwd()),
        'generated_cube': str(cube), 'cube_sha256': cube_hash,
        'script_sha256': _hash(__file__),
        'plot3d_path': str(Path(p3.__file__).resolve()),
        'plot3d_sha256': _hash(p3.__file__),
        'python_executable': sys.executable, 'python_version': sys.version,
        'requested_isovalue': level, 'effective_isovalue': used,
        'signed_contours': [used, -used], 'nonempty_surface_roles': nonempty_roles,
        'quantity': 'Real scalar orbital amplitude', 'field_unit': None,
        'geometry_length_unit': 'Bohr per reviewed loader convention; verify viskit',
        'grid_shape': list(values.shape),
        'raw_field_range': [float(np.min(values)), float(np.max(values))],
        'coordinate_bounds': {
          key: [float(np.min(grid[key])), float(np.max(grid[key]))]
          for key in ('x', 'y', 'z')
        },
        'atoms': [
          {'atomic_number': int(a['index']),
           'position': [float(a[key]) for key in ('x', 'y', 'z')]}
          for a in atoms
        ],
        'construction': {
          'function': 'plot3d.real_orb_plot_cub',
          'reviewed_smoothing': 'gaussian_filter(val, sigma=1), in grid-index units',
          'threshold_adjustment': 'Inherited from plot3d; not overridden',
          'offscreen': bool(mlab.options.offscreen),
          'orientation_axes_exported': False,
          'material_conversion': 'Constant LUT RGB; opacity = actor opacity * LUT alpha/255',
          'mesh_writer': 'Transformed VTK polydata, triangulated OBJ, per-point normals',
          'mtl_opacity_keyword': 'd',
        },
        'actors': records, 'obj_validation': checked,
        'files': [
          {'name': f.name, 'bytes': f.stat().st_size, 'sha256': _hash(f)}
          for f in sorted(stage.iterdir()) if f.is_file()
        ],
        'limitations': [
          'OBJ stores fixed contour geometry, not the original volume or editable pipeline.',
          'MTL coefficients do not guarantee the same appearance in another renderer.',
          'Per-actor groups are stable within this bundle; molecular actors are not atom/bond semantic IDs.',
          'If the local plot3d rejects a contour outside the data range, export reports that error.',
          'Field units and HOMO identity are not guessed from filenames or electron counts.',
        ],
        'validation': {'blender_import_tested': False, 'visual_review_performed': False},
      }
      (stage / f'{target.stem}.scene.json').write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
        encoding='utf-8',
      )
      published = _publish(stage, target.parent)
    print('Exported scene files:', file=sys.stderr, flush=True)
    for file in published:
      print(file, file=sys.stderr, flush=True)
    return manifest
  finally:
    try:
      if figure is not None:
        mlab.close(figure)
      else:
        for scene in list(mlab.get_engine().scenes):
          if scene.name == title:
            mlab.close(scene)
    finally:
      mlab.options.offscreen = previous_offscreen


# -----------------------------------------------------------------------------
# Command-line interface

def main(argv=None) -> None:
  """Preserve the original CLI and add an explicit OBJ export branch."""
  parser = argparse.ArgumentParser(
    description='Visualize a signed scalar orbital, or export OBJ/MTL/JSON and exit.',
    allow_abbrev=False,
  )
  parser.add_argument('molden', help='Existing .molden.input/.molden file, or basename.')
  parser.add_argument('index', help='Explicit 1-based orbital/spin index, e.g. 35A or 35b.')
  parser.add_argument('-slice', action='store_true', dest='do_slice',
                      help='Add interactive cross-sections; incompatible with --export-obj.')
  parser.add_argument('--isovalue', type=float, default=0.05,
                      help='Positive contour magnitude (default: 0.05).')
  parser.add_argument('--export-obj', metavar='FILE.obj',
                      help='Export signed surfaces and geometry without a GUI loop. No overwrite.')
  args = parser.parse_args(argv)
  if args.export_obj and args.do_slice:
    parser.error('-slice is interactive only; do not combine it with --export-obj.')
  try:
    if args.export_obj:
      export_obj(args.molden, args.index, args.isovalue, filename=args.export_obj)
    else:
      cub1c(args.molden, args.index, args.isovalue, slice=args.do_slice)
  except Exception as exc:
    parser.exit(2, f'{parser.prog}: {type(exc).__name__}: {exc}\n')
  return None


if __name__ == '__main__':
  main()
