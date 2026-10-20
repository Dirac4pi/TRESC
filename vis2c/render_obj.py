#!/usr/bin/env python3
"""Render OBJ scenes with Blender 3.0+ in background mode.

Example:
blender --background --factory-startup --python-exit-code 1 \
  --python render_obj.py -- --input orbital.obj --output preview.png

The auto material mode rebuilds explicitly mapped scalar/orbital materials from a
matching .scene.json. Unmapped OBJ scenes keep their imported materials.
No geometry, contour value, scientific field, or source RGB/alpha is edited.
"""

import argparse
from contextlib import ExitStack
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import sys
import tempfile


# -----------------------------------------------------------------------------
# Arguments and provenance

def parse_args(argv=None):
  parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
  parser.add_argument('--input', required=True, type=Path)
  parser.add_argument('--output', required=True, type=Path,
                      help='PNG output; also writes OUTPUT.png.json. No overwrite.')
  parser.add_argument('--width', type=int, default=1200)
  parser.add_argument('--height', type=int, default=1000)
  parser.add_argument('--samples', type=int, default=64)
  parser.add_argument('--threads', type=int, default=4)
  parser.add_argument('--azimuth', type=float, default=45,
                      help='Degrees around source +Z, measured from +X.')
  parser.add_argument('--elevation', type=float, default=25,
                      help='Degrees above source XY; strictly between -89 and 89.')
  parser.add_argument('--margin', type=float, default=1.15)
  parser.add_argument('--transparent', action='store_true',
                      help='Transparent rather than white background.')
  parser.add_argument('--save-blend', type=Path)
  parser.add_argument('--scene-json', type=Path,
                      help='Default: INPUT basename with .scene.json extension.')
  parser.add_argument('--material-mode', choices=('auto', 'manifest', 'imported'),
                      default='auto', help='Rebuild verified mapped materials, or keep the importer result.')
  parser.add_argument('--source-color-space', choices=('srgb', 'linear'),
                      default='srgb', help='Color convention for manifest RGB only; imported colors are never reconverted.')
  parser.add_argument('--surface-shadows', choices=('transparent', 'physical'),
                      default='transparent', help='Manifest surfaces only: transparent shadows are a scientific-illustration preset.')
  parser.add_argument('--surface-roughness', type=float, default=0.36)
  parser.add_argument('--surface-specular', type=float, default=0.30)
  parser.add_argument('--world-strength', type=float, default=0.20,
                      help='Environment illumination; camera background remains white.')
  parser.add_argument('--key-energy', type=float, default=350)
  parser.add_argument('--fill-energy', type=float, default=70)
  parser.add_argument('--rim-energy', type=float, default=0,
                      help='Disabled by default for the scientific preset.')
  parser.add_argument('--no-denoise', action='store_true')
  parser.add_argument('--no-colorbar', action='store_true',
                      help='Omit the quantitative colorbar for mapped scalar surfaces.')
  parser.add_argument('--colorbar-ticks', type=int, default=5,
                      help='Mapped colorbar tick count, between 2 and 11.')
  if argv is None:
    argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else sys.argv[1:]
  args = parser.parse_args(argv)
  if min(args.width, args.height, args.samples, args.threads) < 1:
    parser.error('Dimensions, samples and threads must be positive.')
  numbers = (args.azimuth, args.elevation, args.margin, args.world_strength,
              args.key_energy, args.fill_energy, args.rim_energy,
              args.surface_roughness, args.surface_specular)
  if not all(math.isfinite(value) for value in numbers):
    parser.error('Numeric rendering settings must be finite.')
  if not -89 < args.elevation < 89 or args.margin <= 1:
    parser.error('Require -89 < elevation < 89 and margin > 1.')
  if min(args.world_strength, args.key_energy, args.fill_energy, args.rim_energy) < 0:
    parser.error('Light strengths must be nonnegative.')
  if not 0 <= args.surface_roughness <= 1 or not 0 <= args.surface_specular <= 1:
    parser.error('Surface roughness and specular must lie in [0, 1].')
  if not 2 <= args.colorbar_ticks <= 11:
    parser.error('--colorbar-ticks must lie between 2 and 11.')
  args.input = args.input.expanduser().resolve()
  args.output = args.output.expanduser().absolute()
  if not args.input.is_file() or args.input.suffix.lower() != '.obj':
    parser.error('Input must be an existing OBJ file.')
  if args.output.suffix.lower() != '.png':
    parser.error('Output must end with .png.')
  if args.scene_json:
    args.scene_json = args.scene_json.expanduser().resolve()
    if not args.scene_json.is_file():
      parser.error('--scene-json does not exist.')
  if args.save_blend:
    args.save_blend = args.save_blend.expanduser().absolute()
    if args.save_blend.suffix.lower() != '.blend':
      parser.error('--save-blend must end with .blend.')
  targets = [args.output, Path(str(args.output) + '.json')]
  if args.save_blend:
    targets.append(args.save_blend)
  for target in targets:
    if os.path.lexists(target):
      parser.error(f'Refusing existing output: {target}')
    if not target.parent.is_dir():
      parser.error(f'Output parent directory does not exist: {target.parent}')
  return args


def sha256(path):
  digest = hashlib.sha256()
  with Path(path).open('rb') as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b''):
      digest.update(block)
  return digest.hexdigest()


def material_dependencies(path):
  """Resolve the single-path mtllib entries produced by the vis2c exporters."""
  libraries = []
  with path.open(encoding='utf-8', errors='replace') as stream:
    for line in stream:
      fields = line.strip().split(maxsplit=1)
      if fields and fields[0] == 'mtllib' and len(fields) == 2:
        candidate = path.parent / fields[1].strip()
        if not candidate.is_file():
          raise FileNotFoundError(f'Missing referenced MTL: {candidate}')
        libraries.append(candidate.resolve())
  if not libraries:
    raise ValueError('OBJ has no mtllib reference; expected companion materials.')
  return sorted(set(libraries))


def _rgb_to_linear(rgb, convention):
  """Convert source display RGB once, without changing its intended sRGB palette."""
  rgb = tuple(float(value) for value in rgb)
  if len(rgb) != 3 or any(not math.isfinite(v) or not 0 <= v <= 1 for v in rgb):
    raise ValueError('Manifest RGB must contain three finite values in [0, 1].')
  if convention == 'linear':
    return rgb
  return tuple(v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb)


def load_manifest(args, libraries):
  """Use only explicit material-name mappings with matching OBJ/MTL hashes."""
  path = args.scene_json or args.input.with_suffix('.scene.json')
  if args.material_mode == 'imported':
    return None, {}, {'mode': 'imported', 'reason': 'Explicitly requested'}
  if not path.is_file():
    if args.material_mode == 'manifest':
      raise FileNotFoundError(f'Material manifest is required: {path}')
    return None, {}, {'mode': 'imported', 'reason': 'No scene manifest found'}
  manifest = json.loads(path.read_text(encoding='utf-8'))
  rows = manifest.get('actors', [])
  if not isinstance(rows, list):
    raise ValueError('Invalid manifest actor list.')
  mapped_rows = [row for row in rows if isinstance(row, dict) and row.get('mtl_material')]
  if not mapped_rows:
    if args.material_mode == 'manifest':
      raise ValueError('Manifest contains no explicit mtl_material mappings.')
    return manifest, {}, {
      'mode': 'imported', 'reason': 'Manifest has no verified material mapping',
      'manifest': str(path), 'manifest_sha256': sha256(path),
    }
  inventory = {row['name']: row for row in manifest.get('files', [])}
  for source in (args.input, *libraries):
    entry = inventory.get(source.name)
    if not entry or not entry.get('sha256') or sha256(source) != entry['sha256']:
      raise ValueError(f'Manifest does not match this OBJ/MTL bundle: {source.name}')
  styles = {}
  supported = {'orbital_positive', 'orbital_negative', 'isosurface', 'mapped_isosurface', 'molecular_geometry'}
  for row in mapped_rows:
    name = row['mtl_material']
    if name in styles or row.get('role') not in supported:
      raise ValueError(f'Duplicate or unsupported material mapping: {name}')
    rgb = tuple(float(v) for v in row['color'])
    linear = _rgb_to_linear(rgb, args.source_color_space)
    alpha = float(row['opacity'])
    if not math.isfinite(alpha) or not 0 <= alpha <= 1:
      raise ValueError(f'Invalid manifest opacity: {name}')
    if 'lut_rgba_255' in row:
      rgba = [float(v) for v in row['lut_rgba_255']]
      if (len(rgba) != 4 or any(not math.isfinite(v) or not 0 <= v <= 255 for v in rgba)
          or any(abs(rgb[i] - rgba[i] / 255) > 1e-7 for i in range(3))):
        raise ValueError(f'Inconsistent manifest LUT RGB: {name}')
      if 'actor_opacity' in row:
        expected = float(row['actor_opacity']) * rgba[3] / 255
        if not math.isfinite(expected) or abs(alpha - expected) > 1e-7:
          raise ValueError(f'Inconsistent manifest opacity: {name}')
    styles[name] = {
      'role': row['role'], 'source_rgb': list(rgb),
      'linear_rgb': list(linear), 'opacity': alpha,
    }
    if row['role'] == 'mapped_isosurface':
      mapping = row.get('color_mapping')
      if not isinstance(mapping, dict):
        raise ValueError('Mapped surfaces require quantitative color_mapping metadata.')
      limits, lut = mapping.get('range'), mapping.get('lut_rgba_255')
      if (not isinstance(limits, list) or len(limits) != 2
          or not all(math.isfinite(float(v)) for v in limits) or limits[0] >= limits[1]
          or not isinstance(lut, list) or len(lut) < 2
          or mapping.get('texture_width') != len(lut)):
        raise ValueError('Invalid LUT/range metadata.')
      if any(not isinstance(color, list) or len(color) != 4
              or any(not isinstance(v, (int, float)) or not math.isfinite(v) or not 0 <= v <= 255 for v in color)
              for color in lut):
        raise ValueError('Invalid LUT colors.')
      if mapping.get('texture_interpolation') != 'Linear' or mapping.get('texture_extension') != 'EXTEND':
        raise ValueError('Unsupported mapped texture sampling convention.')
      for key in ('color_texture', 'alpha_texture'):
        filename = row.get(key)
        if key == 'color_texture' and not filename:
          raise ValueError('Mapped surface has no color texture.')
        if not filename:
          continue
        if not isinstance(filename, str) or Path(filename).name != filename or Path(filename).suffix.lower() != '.png':
          raise ValueError('Texture reference must be a local PNG basename.')
        texture = args.input.parent / filename
        entry = inventory.get(filename)
        if (not texture.is_file() or texture.resolve().parent != args.input.parent
            or not entry or sha256(texture) != entry.get('sha256')):
          raise ValueError(f'Missing, external or mismatched texture: {filename}')
        styles[name][key] = str(texture.resolve())
        styles[name][key + '_sha256'] = entry['sha256']
      styles[name]['color_mapping'] = mapping
      styles[name]['backface_culling'] = bool(row.get('backface_culling', False))
    elif row.get('color_texture') or row.get('alpha_texture'):
      raise ValueError('A texture requires an explicit mapped_isosurface role.')

  return manifest, styles, {
    'mode': 'manifest', 'manifest': str(path), 'manifest_sha256': sha256(path),
    'source_color_space': args.source_color_space,
    'color_space_note': 'Explicit rendering convention, not inferred from screenshots or MTL',
  }


# -----------------------------------------------------------------------------
# Blender import and explicit scientific materials

def import_obj(bpy, path):
  if bpy.app.version >= (3, 3, 0):
    result = bpy.ops.wm.obj_import(filepath=str(path), forward_axis='Y', up_axis='Z')
  else:
    import addon_utils
    addon_utils.enable('io_scene_obj', default_set=False, persistent=False)
    result = bpy.ops.import_scene.obj(filepath=str(path), axis_forward='Y', axis_up='Z')
  if 'FINISHED' not in result:
    raise RuntimeError(f'OBJ import did not finish: {result}')
  meshes = [obj for obj in bpy.context.scene.objects if obj.type == 'MESH' and len(obj.data.vertices)]
  if not meshes:
    raise RuntimeError('No nonempty mesh was imported.')
  if any(not obj.material_slots or any(slot.material is None for slot in obj.material_slots) for obj in meshes):
    raise RuntimeError('Some imported meshes have no materials.')
  return meshes


def _principled(nodes, color, roughness, specular):
  node = nodes.new('ShaderNodeBsdfPrincipled')
  node.inputs['Base Color'].default_value = tuple(color) + (1.0,)
  node.inputs['Metallic'].default_value = 0
  node.inputs['Roughness'].default_value = roughness
  node.inputs['Alpha'].default_value = 1
  socket = node.inputs.get('Specular')
  if socket is None:
    socket = node.inputs.get('Specular IOR Level')
  if socket is None:
    raise RuntimeError('Unsupported Principled BSDF specular socket.')
  socket.default_value = specular
  return node


def _texture_node(bpy, nodes, links, filename, convention):
  """Load and pack a verified LUT texture, with explicit UV and color handling."""
  image = bpy.data.images.load(filename, check_existing=True)
  image.colorspace_settings.name = 'sRGB' if convention == 'srgb' else 'Non-Color'
  image.alpha_mode = 'STRAIGHT'
  image.pack()
  tex = nodes.new('ShaderNodeTexImage')
  tex.image, tex.interpolation, tex.extension = image, 'Linear', 'EXTEND'
  coords = nodes.new('ShaderNodeTexCoord')
  links.new(coords.outputs['UV'], tex.inputs['Vector'])
  return tex


def rebuild_materials(objects, styles, args):
  """Rebuild named materials; never infer roles from object order or color."""
  if not styles:
    return []
  materials = {slot.material.name: slot.material
                for obj in objects for slot in obj.material_slots if slot.material}
  unknown = set(materials) - set(styles)
  if unknown:
    raise ValueError('Imported material names lack manifest mappings: ' + ', '.join(sorted(unknown)))
  records = []
  for name, material in materials.items():
    style = styles[name]
    surface = style['role'].startswith('orbital_') or style['role'] in ('isosurface', 'mapped_isosurface')
    material.use_nodes = True
    nodes, links = material.node_tree.nodes, material.node_tree.links
    nodes.clear()
    output = nodes.new('ShaderNodeOutputMaterial')
    opaque = _principled(
      nodes, style['linear_rgb'],
      args.surface_roughness if surface else 0.40,
      args.surface_specular if surface else 0.25,
    )
    mapped = style['role'] == 'mapped_isosurface'
    alpha_socket = None
    if mapped:
      import bpy
      for obj in objects:
        if any(slot.material == material for slot in obj.material_slots) and obj.data.uv_layers.active is None:
          raise ValueError(f'Mapped mesh has no UV layer: {obj.name}')
      tex = _texture_node(bpy, nodes, links, style['color_texture'], args.source_color_space)
      links.new(tex.outputs['Color'], opaque.inputs['Base Color'])
      if style.get('alpha_texture'):
        alpha_tex = _texture_node(bpy, nodes, links, style['alpha_texture'], 'linear')
        gray = nodes.new('ShaderNodeRGBToBW')
        links.new(alpha_tex.outputs['Color'], gray.inputs['Color'])
        multiply = nodes.new('ShaderNodeMath')
        multiply.operation = 'MULTIPLY'
        multiply.inputs[1].default_value = style['opacity']
        links.new(gray.outputs[0], multiply.inputs[0])
        alpha_socket = multiply.outputs[0]

    transparent = nodes.new('ShaderNodeBsdfTransparent')
    transparent.inputs['Color'].default_value = (1, 1, 1, 1)
    # Alpha is applied exactly once. No refraction or volume absorption is added.
    mix = nodes.new('ShaderNodeMixShader')
    mix.inputs[0].default_value = style['opacity']
    if alpha_socket is not None:
      links.new(alpha_socket, mix.inputs[0])
    links.new(transparent.outputs[0], mix.inputs[1])
    links.new(opaque.outputs[0], mix.inputs[2])
    result = mix.outputs[0]
    if mapped and style.get('backface_culling'):
      geometry = nodes.new('ShaderNodeNewGeometry')
      backface = nodes.new('ShaderNodeMixShader')
      links.new(geometry.outputs['Backfacing'], backface.inputs[0])
      links.new(result, backface.inputs[1])
      links.new(transparent.outputs[0], backface.inputs[2])
      result = backface.outputs[0]
    shadow_bypass = surface and args.surface_shadows == 'transparent' 
    if shadow_bypass:
      path = nodes.new('ShaderNodeLightPath')
      shadow_mix = nodes.new('ShaderNodeMixShader')
      links.new(path.outputs['Is Shadow Ray'], shadow_mix.inputs[0])
      links.new(result, shadow_mix.inputs[1])
      links.new(transparent.outputs[0], shadow_mix.inputs[2])
      result = shadow_mix.outputs[0]
    links.new(result, output.inputs['Surface'])
    material.diffuse_color = tuple(style['linear_rgb']) + (style['opacity'],)
    records.append({
      'name': name, **style, 'shadow_bypass': shadow_bypass,
      'roughness': args.surface_roughness if surface else 0.40,
      'specular': args.surface_specular if surface else 0.25,
      'alpha_model': 'actor opacity times LUT alpha, then transparent/opaque mix' if mapped else
                      'mix white Transparent BSDF and opaque Principled; no refraction',
    })
  for role in ('orbital_positive', 'orbital_negative', 'isosurface', 'mapped_isosurface'):
    if any(style['role'] == role for style in styles.values()) and not any(r['role'] == role for r in records):
      raise ValueError(f'No imported material for expected role: {role}')
  return records


def make_world(bpy, scene, strength):
  """Keep a white camera background independent of environment illumination."""
  world = bpy.data.worlds.new('ScientificWorld')
  world.use_nodes = True
  nodes, links = world.node_tree.nodes, world.node_tree.links
  nodes.clear()
  output = nodes.new('ShaderNodeOutputWorld')
  environment = nodes.new('ShaderNodeBackground')
  environment.inputs['Color'].default_value = (1, 1, 1, 1)
  environment.inputs['Strength'].default_value = strength
  background = nodes.new('ShaderNodeBackground')
  background.inputs['Color'].default_value = (1, 1, 1, 1)
  background.inputs['Strength'].default_value = 1
  path = nodes.new('ShaderNodeLightPath')
  mix = nodes.new('ShaderNodeMixShader')
  links.new(path.outputs['Is Camera Ray'], mix.inputs[0])
  links.new(environment.outputs[0], mix.inputs[1])
  links.new(background.outputs[0], mix.inputs[2])
  links.new(mix.outputs[0], output.inputs['Surface'])
  scene.world = world


# -----------------------------------------------------------------------------
# Camera, lighting and rendering

def build_scene(bpy, objects, args):
  from mathutils import Matrix, Vector
  scene = bpy.context.scene
  bpy.context.view_layer.update()
  corners = [obj.matrix_world @ Vector(v) for obj in objects for v in obj.bound_box]
  if not all(math.isfinite(value) for corner in corners for value in corner):
    raise ValueError('Nonfinite model bounds.')
  low = Vector(tuple(min(v[i] for v in corners) for i in range(3)))
  high = Vector(tuple(max(v[i] for v in corners) for i in range(3)))
  center, span = (low + high) / 2, max(high - low)
  if span <= 0:
    raise ValueError('Degenerate model bounds.')
  transform = Matrix.Scale(4.0 / span, 4) @ Matrix.Translation(-center)
  for obj in objects:
    obj.matrix_world = transform @ obj.matrix_world
  bpy.context.view_layer.update()
  corners = [obj.matrix_world @ Vector(v) for obj in objects for v in obj.bound_box]
  az, el = math.radians(args.azimuth), math.radians(args.elevation)
  direction = Vector((math.cos(el)*math.cos(az), math.cos(el)*math.sin(az), math.sin(el)))
  data = bpy.data.cameras.new('RenderCamera')
  camera = bpy.data.objects.new('RenderCamera', data)
  scene.collection.objects.link(camera)
  camera.location = direction * 12
  camera.rotation_euler = (-camera.location).to_track_quat('-Z', 'Y').to_euler()
  rotation = camera.rotation_euler.to_matrix()
  right, up = rotation @ Vector((1, 0, 0)), rotation @ Vector((0, 1, 0))
  horizontal = max(abs(v.dot(right)) for v in corners)
  vertical = max(abs(v.dot(up)) for v in corners)
  data.type = 'ORTHO'
  data.sensor_fit = 'VERTICAL'
  data.ortho_scale = 2 * max(vertical, horizontal / (args.width / args.height)) * args.margin
  data.clip_start, data.clip_end = 0.01, 100
  scene.camera = camera

  def area(name, position, energy, size):
    if energy == 0:
      return
    light_data = bpy.data.lights.new(name, type='AREA')
    light_data.energy, light_data.shape, light_data.size = energy, 'DISK', size
    light = bpy.data.objects.new(name, light_data)
    scene.collection.objects.link(light)
    light.location = position
    light.rotation_euler = (-light.location).to_track_quat('-Z', 'Y').to_euler()

  area('Key', direction*6 + up*4 - right*3, args.key_energy, 4)
  area('Fill', direction*3 + right*5 + up, args.fill_energy, 5)
  area('Rim', -direction*4 + up*4, args.rim_energy, 4)
  make_world(bpy, scene, args.world_strength)
  scene.render.engine = 'CYCLES'
  scene.cycles.device = 'CPU'
  scene.cycles.samples = args.samples
  scene.cycles.seed = 0
  scene.cycles.use_denoising = not args.no_denoise
  scene.cycles.max_bounces = 12
  scene.cycles.transparent_max_bounces = 16
  scene.render.threads_mode, scene.render.threads = 'FIXED', args.threads
  scene.render.resolution_x, scene.render.resolution_y = args.width, args.height
  scene.render.resolution_percentage = 100
  scene.render.pixel_aspect_x = scene.render.pixel_aspect_y = 1
  scene.render.image_settings.file_format = 'PNG'
  scene.render.image_settings.color_mode = 'RGBA'
  scene.render.image_settings.color_depth = '8'
  scene.render.film_transparent = args.transparent
  scene.render.use_file_extension = True
  scene.view_settings.view_transform = 'Standard'
  try:
    scene.view_settings.look = 'None'
  except (TypeError, ValueError):
    print('WARNING: neutral look unavailable; recording the active look.', flush=True)
  scene.view_settings.exposure, scene.view_settings.gamma = 0, 1
  return scene, {
    'source_bounds': [list(low), list(high)],
    'source_to_display_matrix': [list(row) for row in transform],
    'camera_location': list(camera.location),
    'camera_euler_radians': list(camera.rotation_euler),
    'orthographic_scale': data.ortho_scale, 'sensor_fit': data.sensor_fit,
    'camera_source': 'Requested azimuth/elevation; not a recovered Mayavi camera',
  }


def _colorbar_values(low, high, count):
  values = [low + (high - low) * i / (count - 1) for i in range(count)]
  for precision in range(3, 16):
    labels = [f'{v:.{precision}g}' for v in values]
    if len(set(labels)) == count:
      return values, labels
  return values, [repr(v) for v in values]


def add_colorbar(bpy, scene, styles, args, layout):
  """Add a camera-facing quantitative LUT legend, without shading it by lights."""
  mapped = [style for style in styles.values() if style['role'] == 'mapped_isosurface']
  if not mapped or args.no_colorbar:
    return None
  if len(mapped) != 1:
    raise ValueError('This layout supports one quantitative colorbar per figure.')
  from mathutils import Vector
  style = mapped[0]
  mapping = style['color_mapping']
  camera = scene.camera
  rotation = camera.rotation_euler.to_matrix()
  right, up = rotation @ Vector((1, 0, 0)), rotation @ Vector((0, 1, 0))
  forward = rotation @ Vector((0, 0, -1))
  # Reserve the right 20% and translate the camera frame, not scientific geometry.
  camera.data.ortho_scale /= 0.8
  height = camera.data.ortho_scale
  width = height * args.width / args.height
  camera.location += right * (width * 0.1)
  bpy.context.view_layer.update()

  def emission_material(name, *, texture=None, rgb=(0, 0, 0)):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    nodes.clear()
    output = nodes.new('ShaderNodeOutputMaterial')
    emission = nodes.new('ShaderNodeEmission')
    emission.inputs['Color'].default_value = tuple(rgb) + (1,)
    emission.inputs['Strength'].default_value = 1
    if texture:
      tex = _texture_node(bpy, nodes, links, texture, args.source_color_space)
      links.new(tex.outputs['Color'], emission.inputs['Color'])
    transparent = nodes.new('ShaderNodeBsdfTransparent')
    transparent.inputs['Color'].default_value = (1, 1, 1, 1)
    path = nodes.new('ShaderNodeLightPath')
    mix = nodes.new('ShaderNodeMixShader')
    links.new(path.outputs['Is Camera Ray'], mix.inputs[0])
    links.new(transparent.outputs[0], mix.inputs[1])
    links.new(emission.outputs[0], mix.inputs[2])
    links.new(mix.outputs[0], output.inputs['Surface'])
    return mat

  color_mat = emission_material('MappingLegendPalette', texture=style['color_texture'])
  text_mat = emission_material('MappingLegendText')
  names = []

  def quad(name, x, y, w, h, mat, uv=None):
    origin = camera.location + forward * 1.2
    coords = [origin + right*px + up*py for px, py in
              ((x,y), (x+w,y), (x+w,y+h), (x,y+h))]
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata([tuple(v) for v in coords], [], [(0,1,2,3)])
    mesh.update()
    if uv:
      layer = mesh.uv_layers.new(name='UVMap')
      for i, value in enumerate(uv):
        layer.data[i].uv = value
    obj = bpy.data.objects.new(name, mesh)
    scene.collection.objects.link(obj)
    obj.data.materials.append(mat)
    names.append(name)
    return obj

  def text(name, body, x, y, size, max_width):
    curve = bpy.data.curves.new(name, type='FONT')
    curve.body, curve.size, curve.align_y = body, size, 'CENTER'
    obj = bpy.data.objects.new(name, curve)
    scene.collection.objects.link(obj)
    obj.location = camera.location + forward*1.1 + right*x + up*y
    obj.rotation_euler = camera.rotation_euler.copy()
    curve.materials.append(text_mat)
    bpy.context.view_layer.update()
    corners = [obj.matrix_world @ Vector(v) for v in obj.bound_box]
    extent = max(v.dot(right) for v in corners) - min(v.dot(right) for v in corners)
    if extent > max_width:
      obj.scale *= max_width / extent
    names.append(name)
    return obj

  n = mapping['texture_width']
  left, bottom, bar_width, bar_height = .31*width, -.25*height, .035*width, .50*height
  u0, u1 = .5/n, 1-.5/n
  quad('MappingLegendBar', left, bottom, bar_width, bar_height, color_mat,
        [(u0,.5), (u0,.5), (u1,.5), (u1,.5)])
  low, high = mapping['range']
  values, labels = _colorbar_values(low, high, args.colorbar_ticks)
  for i, label in enumerate(labels):
    y = bottom + bar_height * i / (len(labels)-1)
    quad(f'MappingTick{i}', left+bar_width, y-.0005*height, .008*width, .001*height, text_mat)
    text(f'MappingValue{i}', label, left+bar_width+.012*width, y, .025*height, .13*width)
  title = mapping.get('label') or 'Mapping value'
  if mapping.get('unit'):
    title += f' [{mapping["unit"]}]'
  text('MappingLegendTitle', title, .30*width, .30*height, .025*height, .19*width)
  layout['camera_location'] = list(camera.location)
  layout['orthographic_scale'] = camera.data.ortho_scale
  return {
    'range': list(mapping['range']), 'label': title, 'tick_values': values,
    'tick_labels': labels, 'objects': names,
    'appearance': 'Unlit LUT colors, without surface illumination or actor opacity',
    'reserved_width_fraction': .20,
  }


def inspect_png(bpy, path, width, height, transparent):
  """Decode the output and reject fully black/white/transparent images."""
  from array import array
  with path.open('rb') as stream:
    header = stream.read(24)
  if (len(header) != 24 or header[:8] != b'\x89PNG\r\n\x1a\n'
      or header[12:16] != b'IHDR' or struct.unpack('>II', header[16:24]) != (width, height)):
    raise RuntimeError('PNG header/dimensions validation failed.')
  image = bpy.data.images.load(str(path), check_existing=False)
  try:
    if tuple(image.size) != (width, height) or len(image.pixels) != width*height*4:
      raise RuntimeError('Unexpected decoded PNG shape/channels.')
    pixels = array('f', [0.0]) * len(image.pixels)
    image.pixels.foreach_get(pixels)
    rgb_min, rgb_max, total, alpha_min, alpha_max = float('inf'), -float('inf'), 0., 1., 0.
    visible = 0
    for i in range(0, len(pixels), 4):
      r, g, b, a = pixels[i:i+4]
      if not all(math.isfinite(v) for v in (r, g, b, a)):
        raise RuntimeError('Nonfinite decoded image pixels.')
      alpha_min, alpha_max = min(alpha_min, a), max(alpha_max, a)
      if a > 1e-6:
        visible += 1
        rgb_min, rgb_max = min(rgb_min, r, g, b), max(rgb_max, r, g, b)
        total += r + g + b
    if not visible:
      raise RuntimeError('The output is fully transparent.')
    if rgb_max <= 1e-7:
      raise RuntimeError('The output is fully black; inspect Blender/OCIO and scene diagnostics.')
    if rgb_min >= 1 - 1e-6:
      raise RuntimeError('The output is fully white; inspect camera and object visibility.')
    return {
      'decoded': True, 'visible_rgb_min': rgb_min, 'visible_rgb_max': rgb_max,
      'visible_rgb_mean': total / (3*visible), 'alpha_min': alpha_min,
      'alpha_max': alpha_max, 'visible_pixel_count': visible,
      'note': 'Pixel checks are not visual or scientific review; values are Blender-decoded pixels.',
    }
  finally:
    bpy.data.images.remove(image)


def publish(pairs):
  """Publish staged local files with no overwrite; roll back on failure."""
  written = []
  try:
    for source, destination in pairs:
      os.link(source, destination)
      written.append(destination)
  except Exception:
    for destination in reversed(written):
      destination.unlink()
    raise


def main():
  args = parse_args()
  libraries = material_dependencies(args.input)
  manifest, styles, material_info = load_manifest(args, libraries)
  try:
    import bpy
  except ImportError as exc:
    raise RuntimeError('Execute with Blender --background --python, not ordinary Python.') from exc
  print('Blender:', bpy.app.version_string, flush=True)
  print('Material strategy:', material_info, flush=True)
  bpy.ops.wm.read_factory_settings(use_empty=True)
  objects = import_obj(bpy, args.input)
  material_records = rebuild_materials(objects, styles, args)
  scene, layout = build_scene(bpy, objects, args)
  legend_record = add_colorbar(bpy, scene, styles, args, layout)
  if not styles:
    print('WARNING: keeping imported materials; mapped surface RGB/alpha and shadow controls were not applied.', flush=True)
  record = {
    'schema_version': 2, 'blender_version': bpy.app.version_string,
    'input_obj': str(args.input), 'obj_sha256': sha256(args.input),
    'mtl': [{'path': str(p), 'sha256': sha256(p)} for p in libraries],
    'render_script_sha256': sha256(__file__), 'output_png': str(args.output),
    'material_strategy': material_info, 'rebuilt_materials': material_records,
    'colorbar': legend_record,
    'texture_files': [
      {'path': style[key], 'sha256': style[key + '_sha256']}
      for style in styles.values() for key in ('color_texture', 'alpha_texture')
      if key in style
    ],
    'mesh_objects': [{'name': o.name, 'vertices': len(o.data.vertices),
                      'polygons': len(o.data.polygons)} for o in objects],
    'layout': layout, 'azimuth': args.azimuth, 'elevation': args.elevation,
    'width': args.width, 'height': args.height, 'samples': args.samples,
    'threads': args.threads, 'engine': 'CYCLES', 'device': 'CPU',
    'denoising': not args.no_denoise, 'transparent_background': args.transparent,
    'lighting': {'environment_strength': args.world_strength,
                  'camera_background_strength': 1, 'key_energy': args.key_energy,
                  'fill_energy': args.fill_energy, 'rim_energy': args.rim_energy},
    'color_management': {'view_transform': scene.view_settings.view_transform,
                          'look': scene.view_settings.look, 'exposure': 0, 'gamma': 1},
    'visual_review_performed': False,
    'notes': [
      'Scientific illustration preset, not a physical model of orbital material.',
      'Manifest fixed RGB or LUT textures define source colors; color convention is explicit.',
      'Transparent mapped surface shadows bypass shadow-ray attenuation, not camera-facing opacity.',
      'Front/back surface traversal may differ from VTK alpha compositing.',
      'No cube values, contour thresholds, source mesh vertices or normals are edited.',
      'Uniform display normalization is applied to the scene, not the source OBJ.',
      'Exact Mayavi view/lighting matching and automatic best-view selection are not implemented.',
    ],
  }
  with ExitStack() as stack:
    stage = Path(stack.enter_context(tempfile.TemporaryDirectory(
      prefix='.blender-render-', dir=args.output.parent,
    )))
    image = stage / args.output.name
    scene.render.filepath = str(image)
    print('Rendering with Cycles on CPU...', flush=True)
    result = bpy.ops.render.render(write_still=True)
    if 'FINISHED' not in result or not image.is_file():
      raise RuntimeError(f'Render did not finish: {result}')
    record['pixel_validation'] = inspect_png(
      bpy, image, args.width, args.height, args.transparent,
    )
    record['output_sha256'] = sha256(image)
    scene.render.filepath = str(args.output)
    pairs = [(image, args.output)]
    if args.save_blend:
      blend_stage = Path(stack.enter_context(tempfile.TemporaryDirectory(
        prefix='.blender-scene-', dir=args.save_blend.parent,
      )))
      blend = blend_stage / args.save_blend.name
      result = bpy.ops.wm.save_as_mainfile(filepath=str(blend), check_existing=False)
      if 'FINISHED' not in result or not blend.is_file():
        raise RuntimeError('Saving .blend failed.')
      record['blend_file'] = str(args.save_blend)
      pairs.append((blend, args.save_blend))
    meta_target = Path(str(args.output) + '.json')
    meta = stage / meta_target.name
    meta.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    pairs.append((meta, meta_target))
    if sha256(args.input) != record['obj_sha256'] or any(sha256(row['path']) != row['sha256'] for row in record['mtl']):
      raise RuntimeError('Input OBJ/MTL changed during rendering.')
    if any(sha256(row['path']) != row['sha256'] for row in record['texture_files']):
      raise RuntimeError('Input LUT texture changed during rendering.')
    publish(pairs)
  print('Saved PNG:', args.output, flush=True)
  print('Saved metadata:', meta_target, flush=True)


if __name__ == '__main__':
  main()
