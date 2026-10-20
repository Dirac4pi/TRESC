#!/home/lky/miniconda3/envs/vis2c/bin/python

"""
Interactive scalar-cube viewer with an optional OBJ scene-export branch.
author: Dirac4pi
env: vis2c
"""
import argparse
from contextlib import redirect_stdout
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import uuid


#------------------------------------------------------------------------------
# Input and output utilities
def _hash(path):
  """Return the SHA-256 digest of a file using bounded read buffers."""
  digest = hashlib.sha256()
  with Path(path).open('rb') as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b''):
      digest.update(block)
  return digest.hexdigest()


#------------------------------------------------------------------------------
# Cube loading and scene construction
def _read_cube(ampcube, isovalue):
  """Load a scalar cube and validate the contour level and grid arrays."""
  import numpy as np
  import viskit as vk
  level = float(isovalue)
  if not math.isfinite(level):
    raise ValueError('isovalue must be finite.')
  cube = Path(ampcube).expanduser()
  if cube.suffix.lower() not in ('.cub', '.cube'):
    raise ValueError('ampcube must be a .cub or .cube file.')
  if not cube.is_file():
    raise FileNotFoundError(cube)
  print(f'Loading real-space field from {cube}...', flush=True)
  atoms, data = vk.load_cube(str(cube))
  shape = np.asarray(data['value']).shape
  if len(shape) != 3 or min(shape) < 2:
    raise ValueError('Cube values must have three dimensions, each >= 2.')
  for key in ('x', 'y', 'z', 'value'):
    arr = np.asarray(data[key])
    if arr.shape != shape or not np.all(np.isfinite(arr)):
      raise ValueError(f'Invalid shape or nonfinite cube {key} values.')
  return cube, atoms, data, level

#------------------------------------------------------------------------------
def _make_plotter(title, atoms, data, mlab):
  """Return a redraw function using the existing plot3d implementation."""
  import plot3d as p3

  def plot(level):
    used = p3.real_func_ampplot_cub(
      title, atoms, data['value'], data['x'], data['y'], data['z'], level
    )
    return float(used), mlab.figure(title)

  return plot


#------------------------------------------------------------------------------
# Interactive entrypoint
def real3damp(ampcube: str, isovalue: float = 0.05) -> None:
  """Original interactive API and None return value; no export side effects."""
  cube, atoms, data, level = _read_cube(ampcube, isovalue)
  from mayavi import mlab
  from traits.api import HasTraits, Float, observe
  from traitsui.api import View, Item

  plot = _make_plotter(f'amplitude: {cube.stem}', atoms, data, mlab)
  print('Plotting...', flush=True)
  used, figure = plot(level)
  mlab.draw(figure=figure)

  class IsoValueController(HasTraits):
    isov = Float(used, desc='isov', auto_set=False, enter_set=True)

    @observe('isov')
    def update_isov(self, event):
      if not math.isfinite(event.new):
        print('isovalue must be finite.')
        return
      try:
        _, updated_figure = plot(event.new)
      except ValueError as exc:
        print(f'Cannot update isovalue: {exc}')
        return
      mlab.draw(figure=updated_figure)

    view = View(
      Item('isov', label='isovalue', show_label=True),
      width=300, height=200, resizable=True,
    )

  controller = IsoValueController()
  controller.configure_traits()
  return None


#------------------------------------------------------------------------------
# OBJ bundle validation and publication
def _output_path(filename):
  """Validate the OBJ destination and refuse existing bundle files."""
  target = Path(filename).expanduser().absolute()
  if target.suffix != '.obj':
    raise ValueError('Export filename must end with lowercase .obj.')
  if any(c.isspace() for c in target.stem) or '#' in target.stem:
    raise ValueError('Use an OBJ basename without whitespace or # for portable MTL references.')
  if not target.parent.is_dir():
    raise FileNotFoundError(f'Output parent directory does not exist: {target.parent}')
  for candidate in (
    target, target.with_suffix('.mtl'), target.with_suffix('.scene.json'),
  ):
    if os.path.lexists(candidate):
      raise FileExistsError(f'Refusing to overwrite output: {candidate}')
  return target

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _publish(directory, parent):
  """No-clobber same-filesystem publication; roll back files we published on failure."""
  published = []
  try:
    for source in sorted(directory.iterdir()):
      if not source.is_file():
        raise RuntimeError('Unexpected non-file exporter output.')
      target = parent / source.name
      os.link(source, target)  # Atomic refusal if another task created this name.
      published.append(target)
  except Exception:
    for file in reversed(published):
      file.unlink()
    raise
  return published



#------------------------------------------------------------------------------
# Explicit geometry/material serialization

def _collect_actors(renderer, contour):
  """Identify the actual scalar-field surface and retain its VTK actor identity.

  Role identification uses the reviewed amplitude point-data array, not actor
  order or color. The writer assigns group/material names to these exact actors.
  """
  import numpy as np

  collection = renderer.GetActors()
  collection.InitTraversal()
  actors = []
  surface_count = 0
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
      raise RuntimeError('The selected scene contains non-polydata geometry.')
    if mapper.GetScalarVisibility() or actor.GetTexture() is not None:
      raise RuntimeError(
        'Only fixed-color, untextured actors are supported. '
        'Scalar-mapped/texture colors require a separate verified exporter.'
      )
    point_data = data.GetPointData()
    arrays = [] if point_data is None else [
      point_data.GetArrayName(i)
      for i in range(point_data.GetNumberOfArrays())
    ]
    is_surface = 'amplitude' in arrays
    polygon_count = data.GetNumberOfPolys() + data.GetNumberOfStrips()
    if polygon_count == 0:
      if is_surface:
        raise RuntimeError('The amplitude surface contains no polygon geometry.')
      continue
    prop = actor.GetProperty()
    color, opacity = list(prop.GetColor()), float(prop.GetOpacity())
    coefficients = {
      'ambient': float(prop.GetAmbient()),
      'diffuse': float(prop.GetDiffuse()),
      'specular': float(prop.GetSpecular()),
      'specular_power': float(prop.GetSpecularPower()),
    }
    ambient_color, specular_color = list(prop.GetAmbientColor()), list(prop.GetSpecularColor())
    numeric = color + [opacity] + list(coefficients.values()) + ambient_color + specular_color
    if not np.all(np.isfinite(numeric)):
      raise ValueError('Actor material properties must be finite.')
    if len(color) != 3 or any(not 0 <= v <= 1 for v in color) or not 0 <= opacity <= 1:
      raise ValueError('Actor RGB and opacity must lie in [0, 1].')
    record = {
      'renderer_actor_index': index,
      'role': 'isosurface' if is_surface else 'molecular_geometry',
      'role_basis': 'amplitude point-data array' if is_surface else
                    'non-amplitude polygon actor in the reviewed real3damp scene',
      'source_points': int(data.GetNumberOfPoints()),
      'source_cells': int(data.GetNumberOfCells()),
      'point_data_arrays': arrays,
      'bounds': list(actor.GetBounds()),
      'color': color, 'opacity': opacity,
      'ambient_color': ambient_color, 'specular_color': specular_color,
      **coefficients,
    }
    if is_surface:
      record['contour'] = float(contour)
      surface_count += 1
    actors.append((actor, record))
  if surface_count != 1:
    raise RuntimeError(
      f'Expected exactly one nonempty amplitude surface, found {surface_count}. '
      'Refusing ambiguous mapping or skeleton-only output.'
    )
  return actors


def _write_bundle(stage, stem, actors):
  """Write transformed meshes with stable groups and standard MTL dissolve values.

  Each fixed-color actor becomes one explicitly named material. MTL uses d=opacity, never Tr.
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
  with (
    obj_path.open('w', encoding='utf-8') as obj,
    mtl_path.open('w', encoding='utf-8') as mtl,
  ):
    obj.write('# real3damp single-surface scene; OBJ coordinates preserve source units\n')
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
        if style['role'] == 'isosurface':
          raise RuntimeError('The scalar-field surface has no exportable polygon faces.')
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


def _validate_material_mapping(stage, stem, records):
  """Verify every serialized face group references its recorded named material."""
  expected = {}
  for record in records:
    group, material = record['obj_group'], record['mtl_material']
    if group in expected or group != material:
      raise RuntimeError('Group/material mapping must be unique and explicit.')
    expected[group] = material
  if sum(record['role'] == 'isosurface' for record in records) != 1:
    raise RuntimeError('Serialized bundle must contain exactly one isosurface material.')
  active_group = active_material = None
  seen, face_counts, definitions, dissolves = set(), {}, set(), {}
  with (stage / f'{stem}.obj').open(encoding='utf-8') as stream:
    for line in stream:
      parts = line.split()
      if not parts:
        continue
      if parts[0] == 'g':
        if len(parts) != 2:
          raise RuntimeError('Expected one named group per actor.')
        active_group = parts[1]
      elif parts[0] == 'usemtl':
        active_material = parts[1]
      elif parts[0] == 'f':
        if active_group not in expected or expected[active_group] != active_material:
          raise RuntimeError('A face references an unverified group/material pair.')
        seen.add(active_group)
        face_counts[active_group] = face_counts.get(active_group, 0) + 1
  if seen != set(expected):
    raise RuntimeError('Manifest mappings and serialized face groups differ.')
  active_material = None
  with (stage / f'{stem}.mtl').open(encoding='utf-8') as stream:
    for line in stream:
      parts = line.split()
      if not parts:
        continue
      if parts[0] == 'newmtl':
        active_material = parts[1]
        if active_material in definitions:
          raise RuntimeError('Duplicate material definition.')
        definitions.add(active_material)
      elif parts[0] == 'd':
        if active_material is None or active_material in dissolves:
          raise RuntimeError('Invalid/duplicate material dissolve entry.')
        dissolves[active_material] = float(parts[1])
      elif parts[0] == 'Tr':
        raise RuntimeError('Unexpected Tr opacity; this exporter uses standard d.')
  if definitions != set(expected.values()) or set(dissolves) != definitions:
    raise RuntimeError('MTL definitions/dissolve values do not match the manifest.')
  for record in records:
    if face_counts[record['obj_group']] != record['triangles']:
      raise RuntimeError('Group triangle count does not match the manifest.')
    if not math.isclose(dissolves[record['mtl_material']], record['opacity'],
                        rel_tol=1e-8, abs_tol=1e-10):
      raise RuntimeError('Serialized MTL opacity differs from the manifest.')
  return {
    'verified': True, 'groups': sorted(expected),
    'materials': sorted(definitions), 'opacity_keyword': 'd',
  }

#------------------------------------------------------------------------------
# Non-interactive OBJ export
def export_obj(ampcube: str, isovalue: float = 0.05, *, filename):
  """Export a single scalar-field surface with explicit OBJ material mappings.

  The normal scene backend is retained for the user's tested environment.
  No GUI controller, Blender invocation or image export is added. This still
  requires a working Mayavi display backend; a window may briefly appear.
  """
  import numpy as np

  target = _output_path(filename)
  level = float(isovalue)
  if not math.isfinite(level):
    raise ValueError('isovalue must be finite.')
  cube = Path(ampcube).expanduser().resolve()
  if cube.suffix.lower() not in ('.cub', '.cube') or not cube.is_file():
    raise ValueError('Input must be an existing .cub/.cube file.')
  cube_hash = _hash(cube)
  from mayavi import mlab
  previous_offscreen, figure = mlab.options.offscreen, None
  title = f'amplitude: {cube.stem} [obj-{uuid.uuid4().hex}]'
  try:
    mlab.options.offscreen = False
    cube, atoms, data, level = _read_cube(cube, level)
    plot = _make_plotter(title, atoms, data, mlab)
    print('Constructing single-surface scene...', flush=True)
    used, figure = plot(level)
    if not math.isfinite(used):
      raise ValueError('Renderer returned a nonfinite isovalue.')
    figure.scene.show_axes = False
    figure.scene.render()
    from tvtk.api import tvtk
    actor_pairs = _collect_actors(tvtk.to_vtk(figure.scene.renderer), used)
    with tempfile.TemporaryDirectory(
      prefix='.real3damp-obj-', dir=target.parent,
    ) as temporary:
      stage = Path(temporary)
      actor_records, counts = _write_bundle(stage, target.stem, actor_pairs)
      check = _check_bundle(stage, target.stem)
      if counts['vertices'] != check['vertices'] or counts['faces'] != check['faces']:
        raise RuntimeError('Written OBJ mesh counts differ from validated counts.')
      mapping = _validate_material_mapping(stage, target.stem, actor_records)
      if _hash(cube) != cube_hash:
        raise RuntimeError('Input cube changed during export.')
      import plot3d
      plot3d_path = Path(plot3d.__file__)
      manifest = {
        'schema_version': 2,
        'scene_type': 'real3damp_single_isosurface',
        'material_mapping_schema': 'vis2c.obj-material-map.v1',
        'input_cube': str(cube), 'input_sha256': cube_hash,
        'script_sha256': _hash(__file__),
        'plot3d_path': str(plot3d_path), 'plot3d_sha256': _hash(plot3d_path),
        'python_executable': sys.executable, 'python_version': sys.version,
        'requested_isovalue': level, 'effective_isovalue': used,
        'quantity': None, 'field_unit': None,
        'geometry_length_unit': 'Bohr, per reviewed viskit/plot3d convention; verify loader',
        'grid_shape': list(data['value'].shape),
        'raw_field_range': [float(np.min(data['value'])), float(np.max(data['value']))],
        'coordinate_bounds': {
          key: [float(np.min(data[key])), float(np.max(data[key]))]
          for key in ('x', 'y', 'z')
        },
        'atoms': [
          {
            'atomic_number': int(a['index']),
            'position': [float(a[key]) for key in ('x', 'y', 'z')],
          }
          for a in atoms
        ],
        'construction': {
          'function': 'plot3d.real_func_ampplot_cub',
          'reviewed_smoothing': 'gaussian_filter(val, sigma=1.0), in grid-index units',
          'threshold_adjustment': 'Inherited from plot3d; not overridden',
          'smoothing_note': 'Review plot3d hash/version if its implementation changes',
          'offscreen': bool(mlab.options.offscreen),
          'orientation_axes_exported': False,
          'mesh_writer': 'Transformed VTK polydata, triangulated OBJ, per-point normals',
          'material_conversion': 'Actor fixed RGB and actual actor opacity',
          'mtl_opacity_keyword': 'd',
        },
        'actors': actor_records, 'obj_validation': check,
        'material_mapping_validation': mapping,
        'files': [
          {'name': f.name, 'bytes': f.stat().st_size, 'sha256': _hash(f)}
          for f in sorted(stage.iterdir()) if f.is_file()
        ],
        'limitations': [
          'OBJ stores extracted geometry, not the full cube or editable contour pipeline.',
          'Group and material names identify roles within this bundle, not across different exports.',
          'Molecular geometry groups are not guaranteed per-atom or per-bond scientific IDs.',
          'MTL coefficients do not guarantee identical appearance in Blender and Mayavi.',
          'Quantity and field units are not inferred from the filename.',
        ],
        'validation': {'blender_import_tested': False, 'visual_review_performed': False},
      }
      meta = stage / f'{target.stem}.scene.json'
      meta.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
        encoding='utf-8',
      )
      published = _publish(stage, target.parent)
    print('Exported scene files:', flush=True)
    for file in published:
      print(file, flush=True)
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


#-----------------------------------------------------------------------------
# Command-line interface
def main(argv=None) -> None:
  """Dispatch the original interactive mode or the OBJ export mode."""
  parser = argparse.ArgumentParser(
    description='Plot a cube interactively, or export an OBJ scene and exit.',
    allow_abbrev=False,
  )
  parser.add_argument('ampcube', help='Real scalar-field cube (.cub/.cube).')
  parser.add_argument(
    '--isovalue', type=float, default=0.05,
    help='Signed contour level (default: 0.05); plot3d adjustment rules apply.',
  )
  parser.add_argument(
    '--export-obj', metavar='FILE.obj',
    help=(
      'Export OBJ, MTL and FILE.scene.json without entering the GUI loop. '
      'No overwrite.'
    ),
  )
  args = parser.parse_args(argv)
  try:
    if args.export_obj is None:
      real3damp(args.ampcube, args.isovalue)
    else:
      # Progress stays on stderr; stdout is not an image or JSON stream.
      with redirect_stdout(sys.stderr):
        export_obj(args.ampcube, args.isovalue, filename=args.export_obj)
  except Exception as exc:
    parser.exit(2, f'{parser.prog}: {type(exc).__name__}: {exc}\n')
  return None


if __name__ == '__main__':
  main()
