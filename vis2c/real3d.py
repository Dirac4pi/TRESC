#!/home/lky/miniconda3/envs/vis2c/bin/python
"""
Display a cube isosurface colored by a second cube, or export its model.
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
import struct
import sys
import tempfile
import uuid
import zlib


def _hash(path):
  digest = hashlib.sha256()
  with Path(path).open('rb') as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b''):
      digest.update(block)
  return digest.hexdigest()


def _level(value):
  value = float(value)
  if not math.isfinite(value) or value <= 0:
    raise ValueError('isovalue must be positive and finite.')
  return value


def _mapping_options(low, high):
  if (low is None) != (high is None):
    raise ValueError('--map-min and --map-max must be supplied together.')
  if low is not None and (not math.isfinite(low) or not math.isfinite(high) or low >= high):
    raise ValueError('Mapping bounds must be finite and map-min < map-max.')


def _check_mapping_grid(ampdat, mapdat, tolerance=1e-5):
  """Require the same physical grid; tolerance is in the loader's Bohr convention."""
  import numpy as np
  if np.asarray(ampdat['value']).shape != np.asarray(mapdat['value']).shape:
    raise ValueError('Mapping/amplitude shapes differ; explicit resampling is required.')
  for key in ('x', 'y', 'z'):
    delta = float(np.max(np.abs(ampdat[key] - mapdat[key])))
    if not math.isfinite(delta) or delta > tolerance:
      raise ValueError(f'Physical mapping grid differs on {key} ({delta:g} Bohr).')


def _load(ampcube, mapcube, level):
  import numpy as np
  import viskit as vk
  level = _level(level)
  paths, grids, atoms = [], [], None
  for i, filename in enumerate((ampcube, mapcube)):
    path = Path(filename).expanduser().resolve()
    if path.suffix.lower() not in ('.cub', '.cube') or not path.is_file():
      raise ValueError(f'Expected an existing .cub/.cube file: {path}')
    print(f'Loading {"amplitude" if i == 0 else "mapping"} cube: {path}', flush=True)
    current_atoms, grid = vk.load_cube(str(path))
    shape = np.asarray(grid['value']).shape
    if len(shape) != 3 or min(shape) < 2:
      raise ValueError('Cube fields must have three dimensions, each >= 2.')
    for key in ('x', 'y', 'z', 'value'):
      a = np.asarray(grid[key])
      if a.shape != shape or not np.all(np.isfinite(a)):
        raise ValueError(f'Invalid shape or nonfinite {key}: {path}')
    if i == 0:
      atoms = current_atoms
    paths.append(path)
    grids.append(grid)
  _check_mapping_grid(*grids)
  return paths, atoms, grids[0], grids[1], level


def _mapped_module(figure):
  """Find the actual mapped Surface module, not a guessed actor index."""
  from tvtk.api import tvtk
  stack, seen, found = list(figure.children), set(), []
  while stack:
    node = stack.pop()
    if id(node) in seen:
      continue
    seen.add(id(node))
    stack.extend(list(getattr(node, 'children', ())))
    if type(node).__name__ != 'Surface':
      continue
    raw = tvtk.to_vtk(node.actor.actor)
    mapper = raw.GetMapper()
    mapper.Update()
    data = mapper.GetInput()
    if data is not None and data.GetPointData().GetArray('phase') is not None:
      found.append(node)
  if len(found) != 1:
    raise RuntimeError(f'Expected exactly one mapped Surface module; found {len(found)}.')
  return found[0]


def _configure_mapping(figure, low=None, high=None, colormap=None,
                        label=None, unit=None, *, export=False):
  _mapping_options(low, high)
  module = _mapped_module(figure)
  manager = module.module_manager.scalar_lut_manager
  if colormap is not None:
    manager.lut_mode = colormap
  if low is not None:
    manager.use_default_range = False
    manager.data_range = [low, high]
  if export or label is not None or unit is not None:
    text = label if label is not None else ('Mapping value' if export else manager.scalar_bar.title.strip())
    manager.scalar_bar.title = text + (f' [{unit}]' if unit else '')
  if export:
    manager.show_scalar_bar = False
  return module


def real3d(ampcube: str, mapcube: str, isovalue: float = 0.05,
            slice: bool = False, *, map_min=None, map_max=None,
            colormap=None, map_label=None, map_unit=None) -> None:
  """Keep the original interactive API, with optional explicit mapping settings."""
  _mapping_options(map_min, map_max)
  paths, atoms, amp, mapping, level = _load(ampcube, mapcube, isovalue)
  from mayavi import mlab
  import plot3d as p3
  from traits.api import HasTraits, Float, observe
  from traitsui.api import View, Item
  title = f'amplitude: {paths[0].stem}, phase: {paths[1].stem}'
  x, y, z, values = (amp[key] for key in ('x', 'y', 'z', 'value'))
  if slice:
    p3.slice_cub(paths[0].stem, atoms, x, y, z, values,
                  'viridis', values.max(), values.min())

  def plot(contour):
    used, fig = p3.real_func_plot_cub(
      title, atoms, values, mapping['value'], x, y, z, contour,
    )
    if any(v is not None for v in (map_min, map_max, colormap, map_label, map_unit)):
      _configure_mapping(fig, map_min, map_max, colormap, map_label, map_unit)
    return used, fig

  used, figure = plot(level)
  class IsoValueController(HasTraits):
    isov = Float(used, desc='isov', auto_set=False, enter_set=True)

    @observe('isov')
    def update_isov(self, event):
      try:
        new = _level(event.new)
      except ValueError as exc:
        print(exc)
        return
      mlab.close(title)
      plot(new)
      mlab.draw()
      mlab.view()

    view = View(Item('isov', label='isovalue', show_label=True),
                width=300, height=200, resizable=True)

  IsoValueController().configure_traits()
  return None


def _output_path(filename):
  path = Path(filename).expanduser().absolute()
  if path.suffix != '.obj' or any(c.isspace() for c in path.stem) or '#' in path.stem:
    raise ValueError('Use lowercase .obj with no whitespace or # in the basename.')
  if not path.parent.is_dir():
    raise FileNotFoundError(f'Output directory does not exist: {path.parent}')
  for suffix in ('.obj', '.mtl', '.scene.json'):
    if os.path.lexists(path.with_suffix(suffix)):
      raise FileExistsError(path.with_suffix(suffix))
  return path


def _png(path, rgba):
  """Write a lossless, untagged 8-bit RGBA PNG using the standard library."""
  import numpy as np
  data = np.asarray(rgba, dtype=np.uint8)
  if data.ndim != 3 or data.shape[2] != 4:
    raise ValueError('PNG data must have shape (height, width, 4).')
  height, width = data.shape[:2]
  def chunk(kind, payload):
    return (struct.pack('>I', len(payload)) + kind + payload
            + struct.pack('>I', zlib.crc32(kind + payload) & 0xffffffff))
  raw = b''.join(b'\x00' + row.tobytes() for row in data)
  payload = (b'\x89PNG\r\n\x1a\n'
              + chunk(b'IHDR', struct.pack('>IIBBBBB', width, height, 8, 6, 0, 0, 0))
              + chunk(b'IDAT', zlib.compress(raw)) + chunk(b'IEND', b''))
  path.write_bytes(payload)


def _uv_from_values(values, low, high, size):
  """Encode the scalar linearly, with texel centers at the LUT endpoints.

  Do not clip vertex UVs: EXTEND clamps after interpolation in the renderer.
  """
  import numpy as np
  values = np.asarray(values, dtype=float)
  if not np.all(np.isfinite(values)) or not math.isfinite(low) or not math.isfinite(high) or low >= high or size < 2:
    raise ValueError('Invalid scalar-to-UV mapping.')
  return (0.5 + (size - 1) * (values - low) / (high - low)) / size


def _collect(renderer, mapped_actor, label, unit, colormap_name):
  """Capture fixed actors and one continuously mapped scalar surface."""
  import numpy as np
  key = mapped_actor.GetAddressAsString('')
  collection = renderer.GetActors()
  collection.InitTraversal()
  actors, found = [], 0
  for index in range(collection.GetNumberOfItems()):
    actor = collection.GetNextActor()
    if not actor.GetVisibility() or actor.GetMapper() is None:
      continue
    mapper = actor.GetMapper()
    mapper.Update()
    mesh = mapper.GetInput()
    if mesh is None or mesh.GetNumberOfPoints() == 0:
      continue
    if not mesh.IsA('vtkPolyData'):
      raise RuntimeError('Only polygonal scene geometry can be exported.')
    if actor.GetTexture() is not None:
      raise RuntimeError('Pre-existing actor textures need a separate adapter.')
    mapped = actor.GetAddressAsString('') == key
    if mesh.GetNumberOfPolys() + mesh.GetNumberOfStrips() == 0:
      if mapped:
        raise RuntimeError('Mapped contour is empty.')
      continue
    prop = actor.GetProperty()
    style = {
      'role': 'mapped_isosurface' if mapped else 'molecular_geometry',
      'renderer_actor_index': index,
      'color': [1., 1., 1.] if mapped else list(prop.GetColor()),
      'opacity': float(prop.GetOpacity()),
      'ambient': float(prop.GetAmbient()), 'diffuse': float(prop.GetDiffuse()),
      'specular': float(prop.GetSpecular()), 'specular_power': float(prop.GetSpecularPower()),
      'ambient_color': list(prop.GetAmbientColor()),
      'specular_color': list(prop.GetSpecularColor()),
      'backface_culling': bool(prop.GetBackfaceCulling()),
    }
    if not 0 <= style['opacity'] <= 1:
      raise ValueError('Actor opacity must lie in [0, 1].')
    if mapped:
      found += 1
      if not mapper.GetScalarVisibility() or mapper.GetScalarMode() not in (0, 1, 3):
        raise RuntimeError('Mapped surface must use point scalar coloring.')
      if mapper.GetColorMode() == 2:
        raise RuntimeError('Direct RGB scalar arrays are not this export representation.')
      selected = mapper.GetArrayName()
      if selected and selected != 'phase':
        raise RuntimeError(f'Unexpected mapped array: {selected}')
      array = mesh.GetPointData().GetArray('phase')
      if array is None or array.GetNumberOfComponents() != 1:
        raise RuntimeError('The mapped surface must carry the scalar phase array.')
      lut = mapper.GetLookupTable()
      if lut is None or not hasattr(lut, 'GetNumberOfTableValues'):
        raise RuntimeError('A discrete VTK lookup table is required.')
      if lut.GetScale() != 0 or lut.GetUseAboveRangeColor() or lut.GetUseBelowRangeColor():
        raise RuntimeError('Only linear LUTs with endpoint clamping are supported.')
      limits = lut.GetRange() if mapper.GetUseLookupTableScalarRange() else mapper.GetScalarRange()
      low, high = map(float, limits)
      count = int(lut.GetNumberOfTableValues())
      table = np.array([lut.GetTableValue(i) for i in range(count)], dtype=float)
      if (not math.isfinite(low) or not math.isfinite(high) or low >= high
          or count < 2 or count > 65536 or table.shape != (count, 4)
          or not np.all(np.isfinite(table)) or np.any(table < 0) or np.any(table > 1)):
        raise ValueError('Invalid lookup table or mapping range.')
      rgba = np.rint(table * 255).astype(np.uint8)
      style['color_mapping'] = {
        'array': 'phase', 'range': [low, high],
        'label': label, 'unit': unit, 'colormap': colormap_name,
        'lut_rgba_255': rgba.tolist(), 'texture_width': count,
        'texture_interpolation': 'Linear', 'texture_extension': 'EXTEND',
        'scalar_to_u': '(0.5 + (N-1)*(value-min)/(max-min))/N',
        'vtk_interpolate_scalars_before_mapping': bool(mapper.GetInterpolateScalarsBeforeMapping()),
        'transfer_note': 'Interpolated scalar UV followed by LUT sampling; not a pixel-identical VTK shader.',
      }
    elif mapper.GetScalarVisibility():
      raise RuntimeError('Unclassified scalar-colored molecular actor.')
    actors.append((actor, style))
  if found != 1:
    raise RuntimeError('Exactly one mapped surface must be exported.')
  return actors


def _write_bundle(stage, stem, actors):
  """Serialize actual geometry, normals, scalar UVs, named materials and LUTs."""
  import numpy as np
  from vtkmodules.util.numpy_support import vtk_to_numpy
  from vtkmodules.vtkCommonCore import vtkIdList
  from vtkmodules.vtkCommonTransforms import vtkTransform
  from vtkmodules.vtkFiltersCore import vtkPolyDataNormals, vtkTriangleFilter
  from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
  offset, face_total, records = 0, 0, []
  with (stage / f'{stem}.obj').open('w', encoding='utf-8') as obj, \
        (stage / f'{stem}.mtl').open('w', encoding='utf-8') as mtl:
    obj.write(f'# real3d mapped scalar surface\nmtllib {stem}.mtl\n')
    for number, (actor, original) in enumerate(actors):
      style = dict(original)
      transform = vtkTransform()
      transform.SetMatrix(actor.GetMatrix())
      moved = vtkTransformPolyDataFilter()
      moved.SetInputData(actor.GetMapper().GetInput())
      moved.SetTransform(transform)
      triangles = vtkTriangleFilter()
      triangles.SetInputConnection(moved.GetOutputPort())
      triangles.PassVertsOff()
      triangles.PassLinesOff()
      triangles.Update()
      mesh = triangles.GetOutput()
      if mesh.GetNumberOfPolys() == 0:
        raise RuntimeError('An exported polygon actor became empty.')
      generated = mesh.GetPointData().GetNormals() is None
      if generated:
        normal_filter = vtkPolyDataNormals()
        normal_filter.SetInputData(mesh)
        normal_filter.SplittingOff()
        normal_filter.ConsistencyOn()
        normal_filter.Update()
        mesh = normal_filter.GetOutput()
      points = vtk_to_numpy(mesh.GetPoints().GetData())
      normals = vtk_to_numpy(mesh.GetPointData().GetNormals())
      if points.shape != normals.shape or not np.all(np.isfinite(points)) or not np.all(np.isfinite(normals)):
        raise ValueError('Invalid mesh points/normals.')
      name = f'{style["role"]}_{number:05d}'
      uv = np.full(len(points), 0.5)
      if style['role'] == 'mapped_isosurface':
        mapping = style['color_mapping']
        array = mesh.GetPointData().GetArray(mapping['array'])
        if array is None:
          raise RuntimeError('Mapped scalar array was lost during mesh conversion.')
        values = np.asarray(vtk_to_numpy(array), dtype=float).reshape(-1)
        if len(values) != len(points) or not np.all(np.isfinite(values)):
          raise ValueError('Invalid mapped vertex scalar array.')
        low, high = mapping['range']
        rgba = np.asarray(mapping['lut_rgba_255'], dtype=np.uint8)
        uv = _uv_from_values(values, low, high, len(rgba))
        color = rgba.copy()
        color[:, 3] = 255
        texture = f'{stem}.{name}.color.png'
        _png(stage / texture, np.repeat(color[None, :, :], 4, axis=0))
        style['color_texture'] = texture
        if np.any(rgba[:, 3] != 255):
          alpha = np.repeat(rgba[:, 3:4], 4, axis=1)
          alpha[:, 3] = 255
          texture_alpha = f'{stem}.{name}.alpha.png'
          _png(stage / texture_alpha, np.repeat(alpha[None, :, :], 4, axis=0))
          style['alpha_texture'] = texture_alpha
        style['mapped_vertex_range'] = [float(values.min()), float(values.max())]
        style['below_range_vertex_count'] = int(np.count_nonzero(values < low))
        style['above_range_vertex_count'] = int(np.count_nonzero(values > high))
      obj.write(f'o {name}\ng {name}\nusemtl {name}\n')
      for v in points:
        obj.write('v ' + ' '.join(f'{float(x):.12g}' for x in v) + '\n')
      for u in uv:
        obj.write(f'vt {float(u):.12g} 0.5\n')
      for normal in normals:
        obj.write('vn ' + ' '.join(f'{float(x):.12g}' for x in normal) + '\n')
      cells, ids, count = mesh.GetPolys(), vtkIdList(), 0
      cells.InitTraversal()
      while cells.GetNextCell(ids):
        if ids.GetNumberOfIds() != 3:
          raise RuntimeError('Expected triangulated mesh.')
        indices = [int(ids.GetId(i)) + offset + 1 for i in range(3)]
        obj.write('f ' + ' '.join(f'{i}/{i}/{i}' for i in indices) + '\n')
        count += 1
      mtl.write(f'newmtl {name}\n')
      for label, data in (
        ('Ka', [v * style['ambient'] for v in style['ambient_color']]),
        ('Kd', style['color']),
        ('Ks', [v * style['specular'] for v in style['specular_color']]),
      ):
        mtl.write(label + ' ' + ' '.join(f'{float(v):.9g}' for v in data) + '\n')
      mtl.write(f'Ns {max(0., min(1000., style["specular_power"])):.9g}\n')
      mtl.write(f'd {style["opacity"]:.9g}\nillum 2\n')
      if style.get('color_texture'):
        mtl.write(f'map_Kd {style["color_texture"]}\n')
      mtl.write('\n')
      style.update(obj_group=name, mtl_material=name, vertices=len(points),
                    triangles=count, normals_generated=generated)
      records.append(style)
      offset += len(points)
      face_total += count
  return records, {'vertices': offset, 'faces': face_total}


def _validate_bundle(stage, stem, records):
  """Check all face/material identities, index counts, and texture references."""
  expected = {r['obj_group']: r for r in records}
  if len(expected) != len(records) or sum(r['role'] == 'mapped_isosurface' for r in records) != 1:
    raise RuntimeError('Ambiguous model/material records.')
  counts = {'v': 0, 'vt': 0, 'vn': 0, 'f': 0}
  group = material = None
  seen, group_faces, maximum = set(), {}, 0
  with (stage / f'{stem}.obj').open(encoding='utf-8') as stream:
    for line in stream:
      fields = line.split()
      if not fields:
        continue
      key = fields[0]
      if key in ('v', 'vt', 'vn'):
        counts[key] += 1
        if not all(math.isfinite(float(v)) for v in fields[1:]):
          raise ValueError('Nonfinite OBJ geometry/UV.')
      elif key == 'mtllib' and fields[1:] != [f'{stem}.mtl']:
        raise RuntimeError('Unexpected MTL library reference.')
      elif key == 'g':
        group = fields[1]
      elif key == 'usemtl':
        material = fields[1]
      elif key == 'f':
        if group not in expected or expected[group]['mtl_material'] != material:
          raise RuntimeError('Face group/material mismatch.')
        if len(fields) != 4:
          raise RuntimeError('OBJ face is not triangular.')
        for vertex in fields[1:]:
          parts = list(map(int, vertex.split('/')))
          if len(parts) != 3 or len(set(parts)) != 1 or parts[0] < 1:
            raise RuntimeError('Invalid OBJ position/UV/normal indexing.')
          maximum = max(maximum, parts[0])
        seen.add(group)
        group_faces[group] = group_faces.get(group, 0) + 1
        counts['f'] += 1
  if not counts['v'] or not counts['f'] or not counts['v'] == counts['vt'] == counts['vn'] or maximum > counts['v'] or seen != set(expected):
    raise RuntimeError('Incomplete OBJ geometry or mapping.')
  definitions, active = {}, None
  with (stage / f'{stem}.mtl').open(encoding='utf-8') as stream:
    for line in stream:
      fields = line.split()
      if not fields:
        continue
      if fields[0] == 'newmtl':
        active = fields[1]
        if active in definitions:
          raise RuntimeError('Duplicate MTL definition.')
        definitions[active] = {}
      elif active and fields[0] in ('d', 'map_Kd'):
        definitions[active][fields[0]] = fields[1]
  if set(definitions) != {r['mtl_material'] for r in records}:
    raise RuntimeError('MTL definitions do not match model records.')
  for record in records:
    data = definitions[record['mtl_material']]
    if group_faces[record['obj_group']] != record['triangles'] or not math.isclose(float(data['d']), record['opacity'], rel_tol=1e-8, abs_tol=1e-10):
      raise RuntimeError('Material opacity or triangle count mismatch.')
    if record.get('color_texture') != data.get('map_Kd'):
      raise RuntimeError('Color texture does not match the material record.')
    for key in ('color_texture', 'alpha_texture'):
      if record.get(key) and not (stage / record[key]).is_file():
        raise RuntimeError('Missing texture file.')
  return counts


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


def export_obj(ampcube, mapcube, isovalue=0.05, *, filename,
                map_min=None, map_max=None, colormap=None,
                map_label=None, map_unit=None):
  """Export geometry plus quantitative LUT/UV mapping without a GUI loop."""
  import numpy as np
  _level(isovalue)
  _mapping_options(map_min, map_max)
  target = _output_path(filename)
  source_paths = [Path(v).expanduser().resolve() for v in (ampcube, mapcube)]
  source_hashes = [_hash(v) for v in source_paths]
  from mayavi import mlab
  previous, figure = mlab.options.offscreen, None
  title = f'real3d mapped surface [obj-{uuid.uuid4().hex}]'
  try:
    mlab.options.offscreen = False
    paths, atoms, amp, mapped, level = _load(ampcube, mapcube, isovalue)
    import plot3d as p3
    from tvtk.api import tvtk
    used, figure = p3.real_func_plot_cub(
      title, atoms, amp['value'], mapped['value'], amp['x'], amp['y'], amp['z'], level,
    )
    used = _level(used)
    module = _configure_mapping(figure, map_min, map_max, colormap,
                                map_label, map_unit, export=True)
    figure.scene.show_axes = False
    figure.scene.render()
    actors = _collect(
      tvtk.to_vtk(figure.scene.renderer), tvtk.to_vtk(module.actor.actor),
      map_label or 'Mapping value', map_unit,
      str(module.module_manager.scalar_lut_manager.lut_mode),
    )
    with tempfile.TemporaryDirectory(prefix='.real3d-obj-', dir=target.parent) as name:
      stage = Path(name)
      records, counts = _write_bundle(stage, target.stem, actors)
      check = _validate_bundle(stage, target.stem, records)
      if counts['vertices'] != check['v'] or counts['faces'] != check['f']:
        raise RuntimeError('Serialized geometry counts differ.')
      if [_hash(v) for v in paths] != source_hashes:
        raise RuntimeError('Input cube changed during export.')
      mapping = next(r['color_mapping'] for r in records if r['role'] == 'mapped_isosurface')
      manifest = {
        'schema_version': 3, 'scene_type': 'real3d_mapped_isosurface',
        'material_mapping_schema': 'vis2c.obj-material-map.v1',
        'amplitude_cube': str(paths[0]), 'amplitude_sha256': source_hashes[0],
        'mapping_cube': str(paths[1]), 'mapping_sha256': source_hashes[1],
        'script_sha256': _hash(__file__), 'plot3d_sha256': _hash(p3.__file__),
        'requested_isovalue': level, 'effective_isovalue': used,
        'grid_shape': list(amp['value'].shape), 'grid_match_tolerance_bohr': 1e-5,
        'amplitude_range': [float(amp['value'].min()), float(amp['value'].max())],
        'mapping_raw_range': [float(mapped['value'].min()), float(mapped['value'].max())],
        'color_mapping': mapping, 'actors': records, 'obj_validation': check,
        'construction': {
          'function': 'plot3d.real_func_plot_cub',
          'reviewed_smoothing': 'Amplitude gaussian_filter sigma=1; mapping array not smoothed',
          'offscreen': bool(mlab.options.offscreen),
          'color_transfer': 'Vertex scalar UV with sampled 8-bit LUT texture, endpoint extension',
        },
        'files': [{'name': f.name, 'bytes': f.stat().st_size, 'sha256': _hash(f)}
                  for f in sorted(stage.iterdir()) if f.is_file()],
        'limitations': [
          'Keep OBJ, MTL, LUT textures and scene.json together.',
          'OBJ stores the extracted contour, not the original volume.',
          'Texture interpolation and lighting may differ from VTK; this is not pixel-identical export.',
          'Surface lighting changes apparent color; the legend displays the unlit LUT.',
          'Label and unit are user definitions; no ESP/IGMH identity is inferred.',
          'Grid coordinates retain the source plotting pipeline behavior; no new resampling is added.',
        ],
        'validation': {'blender_import_tested': False, 'visual_review_performed': False},
      }
      (stage / f'{target.stem}.scene.json').write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + '\n', encoding='utf-8',
      )
      files = _publish(stage, target.parent)
    for file in files:
      print(f'Exported: {file}', flush=True)
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
      mlab.options.offscreen = previous


def main(argv=None):
  parser = argparse.ArgumentParser(
    description='Display a scalar isosurface with a second field mapped to color, or export its model.',
    allow_abbrev=False,
  )
  parser.add_argument('ampcube')
  parser.add_argument('mapcube')
  parser.add_argument('-slice', action='store_true', dest='do_slice')
  parser.add_argument('--isovalue', type=float, default=0.05)
  parser.add_argument('--export-obj', metavar='FILE.obj')
  parser.add_argument('--map-min', type=float)
  parser.add_argument('--map-max', type=float)
  parser.add_argument('--colormap', help='Mayavi LUT name; default preserves the plotting function.')
  parser.add_argument('--map-label', help='Meaning of the mapped quantity; export default: Mapping value.')
  parser.add_argument('--map-unit', help='Unit of the mapped field, if known.')
  args = parser.parse_args(argv)
  if args.export_obj and args.do_slice:
    parser.error('-slice is interactive-only; do not combine it with --export-obj.')
  options = dict(map_min=args.map_min, map_max=args.map_max, colormap=args.colormap,
                  map_label=args.map_label, map_unit=args.map_unit)
  try:
    if args.export_obj:
      with redirect_stdout(sys.stderr):
        export_obj(args.ampcube, args.mapcube, args.isovalue,
                    filename=args.export_obj, **options)
    else:
      real3d(args.ampcube, args.mapcube, args.isovalue, slice=args.do_slice, **options)
  except Exception as exc:
    parser.exit(2, f'{parser.prog}: {type(exc).__name__}: {exc}\n')


if __name__ == '__main__':
  main()
