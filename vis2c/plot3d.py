#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
plot 3D objects
author:Dirac4pi
env:vis2c
'''

from functools import lru_cache
from pathlib import Path
from os import environ
import numpy as np
from json import load as jsload
import qcelemental as qcel
qcel.CovalentRadii("ALVAREZ2008")
from mayavi import mlab
from scipy.ndimage import gaussian_filter
from tvtk.api import tvtk
from mendeleev.mendeleev import element
from typing import Optional, cast

dtol = 0.01 # tolerance of delaunay
maxpoints = 1e5 # max grid points

#------------------------------------------------------------------------------
@lru_cache(maxsize=118)
def _display_vdw_radius_pm(z: int) -> float:
  """Return H=120 pm or the Batsanov radius; reject missing data."""
  z = int(z)
  if z == 1:
    return 120.0
  if z in (2, 10, 18, 36, 54, 86, 118):
    raise ValueError(f"Noble-gas display radius is not configured: Z={z}.")
  el = element(z)
  radius = cast(Optional[float], el.vdw_radius_batsanov)
  if radius is None:
    raise ValueError(f"No Batsanov van der Waals radius for Z={z}.")
  return float(radius)

#-------------------------------------------------------------------------------
def sync_camera(fig1, fig2, sync_coe:float) -> None:
  """
  Synchronise camera parameters from fig1.scene to fig2.scene
  --
  """
  fig2.scene.camera.position = sync_coe*(fig1.scene.camera.position - \
  fig1.scene.camera.focal_point) + fig2.scene.camera.focal_point
  fig2.scene.camera.view_angle = fig1.scene.camera.view_angle
  fig2.scene.camera.view_up = fig1.scene.camera.view_up
  fig2.scene.render()

#------------------------------------------------------------------------------
def cmplx_orb_plot_cub(title: str, atoms, amp, pha, x, y, z, isovalue: float):
  '''
  Plot amplitude (contour) and phase (color) of scalar complex structured grid.
  --
  !!! For 2c MO, alpha and beta components should be plotted separately.\n
  title: Title of the figure\n
  atoms: Atom information dictionary\n
  amp: Amplitude of the function\n
  pha: Phase of the function\n
  x, y, z: Grid points\n
  isovalue: Isovalue of amplitude\n
  Return: Final isovalue used, figure object
  '''
  if isovalue <= 0.0:
    raise RuntimeError('isovalue should be positive')
  fig = mlab.figure(figure=f'{title}', size=(600, 600),
                    fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
  mlab.clf(fig)
  fig_try = mlab.figure(675)
  amp_smooth = gaussian_filter(amp, sigma=1)
  try:
    mlab.contour3d(x, y, z, amp_smooth, figure=fig_try, contours=[isovalue])
  except Exception as e:
    print(f'Error: {e}')
    if str(e).find('Contour instance must be') >= 0:
      # Extract safe isovalue from error message and scale by 0.5
      safe_val_str = str(e)[str(e).rfind('<=')+3 : str(e).rfind('<=')+14]
      isovalue = float(safe_val_str) * 0.5
      print('------\nIgnore previous error, adjusted isovalue to:', isovalue)
    else:
      raise RuntimeError()
  mlab.close(675)
  # Plot atoms
  # Force active figure to main figure before plotting atoms
  mlab.figure(fig)
  atom_list = atoms['index']
  atomx = atoms['x']
  atomy = atoms['y']
  atomz = atoms['z']
  atoms_plot(atom_list, atomx, atomy, atomz)
  src = mlab.pipeline.scalar_field(x, y, z, amp_smooth, figure=fig)
  # Explicitly name arrays to prevent VTK active scalar loss
  src.image_data.point_data.get_array(0).name = 'amplitude'
  src.image_data.point_data.add_array(pha.T.ravel())
  src.image_data.point_data.get_array(1).name = 'phase'
  src.update()
  src2 = mlab.pipeline.set_active_attribute(src, point_scalars='amplitude', \
                                            figure=fig)
  contourpha = mlab.pipeline.contour(src2, figure=fig)
  contourpha.filter.contours = [isovalue]
  contourpha2 = mlab.pipeline.set_active_attribute(contourpha, 
                                             point_scalars='phase', figure=fig)
  if 'Spin' in title or 'spin' in title:
    surf = mlab.pipeline.surface(contourpha2, colormap='RdBu', figure=fig,
                                 vmax=np.pi, vmin=0)
    surf.module_manager.scalar_lut_manager.data_range = [0, np.pi/2]
  else:
    surf = mlab.pipeline.surface(contourpha2, colormap='hsv', figure=fig,
                                 vmax=np.pi, vmin=-np.pi)
    surf.module_manager.scalar_lut_manager.data_range = [-np.pi, np.pi]
  # Increase opacity for richer colors, but keep backface_culling 
  # so the inner wall doesn't occlude the central atoms.
  surf.actor.property.opacity = 0.65
  surf.actor.property.backface_culling = True
  # Enhance diffuse lighting and reduce specular wash-out
  surf.actor.property.interpolation = 'phong'
  surf.actor.property.ambient  = 0.20
  surf.actor.property.diffuse  = 0.85
  surf.actor.property.specular = 0.20
  surf.actor.property.specular_power = 15
  cb = mlab.colorbar(object=surf, orientation='vertical', \
                     nb_labels=5)
  cb.label_text_property.color = (0, 0, 0)
  cb.title_text_property.color = (0, 0, 0)
  cb.scalar_bar.label_format = '%.3f'       # format of colorbar labels
  cb.scalar_bar.unconstrained_font_size = True 
  cb.label_text_property.font_size = 12     # size of colorbar labels
  cb.title_text_property.font_size = 12     # size of colorbar title
  cb.title_text_property.italic = False
  cb.title_text_property.bold = False
  cb.label_text_property.italic = False
  cb.label_text_property.bold = False
  mlab.view(figure=fig)
  fig.scene.show_axes = True
  return isovalue, fig

#-------------------------------------------------------------------------------
def cmplx_orb_plot_mog(title:str, atoms, amplitude, \
                       phase, x, y, z, isovalue:float)->float:
  '''
  plot amplitude(contour) and phase(color) of scalar complex unstructured grid
  --
  !!! for 2c MO, alpha and beta components should be plotted separately.\n
  title: title of figure\n
  atoms: atom information\n
  amplitude: amplitude of function\n
  phase: phase of function\n
  x: x grid points\n
  y: y grid points\n
  z: z grid points\n
  isovalue: isovalue of amplitude\n
  Return: final isovalue
  '''
  global dtol, maxpoints
  print('Unstructured data need to be Delaunay triangulated before extracting')
  print('isosurface, the time consumed increases significantly...')
  print('...please be patient...')
  if isovalue <= 0.0:
    raise RuntimeError('isovalue should be positive')
  n_points = len(x)
  if n_points > maxpoints:
    sample_size = int(maxpoints)
    indices = np.random.choice(n_points, sample_size, replace=False)
    x = x[indices]
    y = y[indices]
    z = z[indices]
    amplitude = amplitude[indices]
    phase = phase[indices]
    print(f"data points sample from {n_points} to {sample_size}")
  # -------------<plot amplitude with atoms>-------------
  fig1 = mlab.figure(figure=f'{title}(amplitude)', size=(400,400), \
    fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig1)
  # amplitude of MO not amplitude of linear combine of AO
  # dynamic adjustment of thresholds
  if isovalue > 0.5*amplitude.max():
    isovalue = 0.5 * (amplitude.max() + amplitude.min())
    print(f"Adjust isovalue to {isovalue} (mid-range)")
  elif isovalue < 2*amplitude.min():
    isovalue = 0.5 * (amplitude.max() + amplitude.min())
    print(f"Adjust isovalue to {isovalue} (mid-range)")
  points = mlab.pipeline.scalar_scatter(x, y, z, amplitude)
  delaunay = mlab.pipeline.delaunay3d(points)
  delaunay.filter.tolerance = dtol
  contour = mlab.pipeline.contour(delaunay)
  try:
    contour.filter.contours = [isovalue]
    surface = mlab.pipeline.surface(contour, figure=fig1, opacity=.7)
    surface.actor.mapper.scalar_visibility = True
    lut = surface.module_manager.scalar_lut_manager.lut
    lut.table = [(110, 110, 110, 180)] * 256
    surface.actor.property.interpolation = 'phong'
    surface.actor.property.ambient  = 0.15
    surface.actor.property.diffuse  = 0.75
    surface.actor.property.specular = 0.35
    surface.actor.property.specular_power = 20
  except Exception as e:
    raise RuntimeError(f'Error: {e}')
  atom_list = atoms['index']
  atomx = atoms['x']
  atomy = atoms['y']
  atomz = atoms['z']
  mlab.figure(fig1)
  atoms_plot(atom_list, atomx, atomy, atomz)
  # -------------<plot phase>-------------
  fig2 = mlab.figure(figure=f'{title}(phase), isovalue={isovalue}', 
    size=(400,400), fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig2)
  # create an unstructured point cloud source with amplitude as primary scalar.
  src = mlab.pipeline.scalar_scatter(x, y, z, amplitude)
  # add phase to point cloud
  src.mlab_source.dataset.point_data.add_array(phase.ravel())
  src.mlab_source.dataset.point_data.get_array(1).name = 'phase'
  src.update()
  # build pipeline: triangulation -> isosurface extraction
  delaunay = mlab.pipeline.delaunay3d(src)
  contour = mlab.pipeline.contour(delaunay)
  contour.filter.contours = [isovalue]
  # color mapping to phase
  contour_phase = mlab.pipeline.set_active_attribute(contour, \
                                                     point_scalars='phase')
  if title == 'spin':
    surf = mlab.pipeline.surface(contour_phase, colormap='RdBu', opacity=0.5)
    # spin phase is non-periodic and range in [0,pi/2]
    surf.module_manager.scalar_lut_manager.data_range = [0, np.pi/2]
  else:
    surf = mlab.pipeline.surface(contour_phase, colormap='hsv', opacity=0.5)
    surf.module_manager.scalar_lut_manager.data_range = [-np.pi, np.pi]
  mlab.colorbar(title='Phase', orientation='vertical')
  mlab.axes()
  # ---------<synchronise camera>---------
  sync_coe = (fig2.scene.camera.position - fig2.scene.camera.focal_point) \
           / (fig1.scene.camera.position - fig1.scene.camera.focal_point)
  fig1.scene.interactor.add_observer("InteractionEvent", lambda *args: \
    sync_camera(fig1, fig2, sync_coe))
  return isovalue

#-------------------------------------------------------------------------------
def real_orb_plot_cub(title:str, atoms, val, x, y, z, isovalue:float)->float:
  '''
  plot scalar real structured grid
  --
  title: title of figure\n
  atoms: atom information\n
  val: grid value of function\n
  x: x grid points\n
  y: y grid points\n
  z: z grid points\n
  isovalue: isovalue of val\n
  Return: final isovalue
  '''
  if isovalue <= 0.0:
    raise RuntimeError('isovalue should be positive')
  # -------------<plot orb>-------------
  fig1 = mlab.figure(figure=f'{title}', size=(400,400), \
    fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig1)
  # amplitude of MO not amplitude of linear combine of AO
  val_smoothed = gaussian_filter(val, sigma=1)
  if isovalue > 0.5*np.abs(val_smoothed).max():
    isovalue = 0.5 * (np.abs(val_smoothed).max() + np.abs(val_smoothed).min())
    print(f"Adjust isovalue to {isovalue} (mid-range)")
  elif isovalue < 2*np.abs(val_smoothed).min():
    isovalue = 0.5 * (np.abs(val_smoothed).max() + np.abs(val_smoothed).min())
    print(f"Adjust isovalue to {isovalue} (mid-range)")
  contourpos = mlab.contour3d(x, y, z, val_smoothed, \
                              contours=[isovalue], opacity=1.0)
  contourpos.actor.mapper.scalar_visibility = True
  lutpos = contourpos.module_manager.scalar_lut_manager.lut
  contourneg = mlab.contour3d(x, y, z, val_smoothed, \
                              contours=[-isovalue], opacity=1.0)
  contourneg.actor.mapper.scalar_visibility = True
  lutneg = contourneg.module_manager.scalar_lut_manager.lut
 
  #-- Scheme 1: Blue / Green,  Science Classic
  
  lutpos.table = [(0, 151, 141, 180)] * 256
  contourpos.actor.property.interpolation = 'phong'
  contourpos.actor.property.ambient = 0.30
  contourpos.actor.property.diffuse = 0.85
  contourpos.actor.property.specular = 0.35
  contourpos.actor.property.specular_power = 25
  contourpos.actor.property.specular_color = (1.0, 1.0, 1.0)

  lutneg.table = [(0, 91, 160, 180)] * 256
  contourneg.actor.property.interpolation = 'phong'
  contourneg.actor.property.ambient = 0.30
  contourneg.actor.property.diffuse = 0.85
  contourneg.actor.property.specular = 0.35
  contourneg.actor.property.specular_power = 25
  contourneg.actor.property.specular_color = (1.0, 1.0, 1.0)
  

  # -- Scheme 2: Red / Blue, Nat. Chem. & VMD Classic
  '''
  lutpos.table = [(214, 39, 40, 180)] * 256
  contourpos.actor.property.interpolation = 'phong'
  contourpos.actor.property.ambient  = 0.25
  contourpos.actor.property.diffuse  = 0.75
  contourpos.actor.property.specular = 0.35
  contourpos.actor.property.specular_power = 20
  
  lutneg.table = [(31, 119, 180, 180)] * 256
  contourneg.actor.property.interpolation = 'phong'
  contourneg.actor.property.ambient  = 0.3
  contourneg.actor.property.diffuse  = 0.8
  contourneg.actor.property.specular = 0.35
  contourneg.actor.property.specular_power = 20
  '''

  # -- Scheme 3: Blue / Amber, ACS & Avogadro Community
  '''
  lutpos.table = [(66, 134, 244, 180)] * 256
  contourpos.actor.property.interpolation = 'phong'
  contourpos.actor.property.ambient = 0.20
  contourpos.actor.property.diffuse = 0.80
  contourpos.actor.property.specular = 0.15
  contourpos.actor.property.specular_power = 35

  lutneg.table = [(226, 153, 34, 180)] * 256
  contourneg.actor.property.interpolation = 'phong'
  contourneg.actor.property.ambient = 0.20
  contourneg.actor.property.diffuse = 0.80
  contourneg.actor.property.specular = 0.15
  contourneg.actor.property.specular_power = 35
  '''

  # -- Scheme 4:  Custom Scheme
  '''
  lutpos.table = [(61, 159, 60, 180)] * 256
  contourpos.actor.property.interpolation = 'phong'
  contourpos.actor.property.ambient  = 0.25
  contourpos.actor.property.diffuse  = 0.75
  contourpos.actor.property.specular = 0.15
  contourpos.actor.property.specular_power = 25
  
  lutneg.table = [(54, 125, 176, 180)] * 256
  contourneg.actor.property.interpolation = 'phong'
  contourneg.actor.property.ambient  = 0.25
  contourneg.actor.property.diffuse  = 0.75
  contourneg.actor.property.specular = 0.15
  contourneg.actor.property.specular_power = 25
  '''

  # -------------<plot atoms>-------------
  atom_list = atoms['index']
  atomx = atoms['x']
  atomy = atoms['y']
  atomz = atoms['z']
  mlab.figure(fig1)
  atoms_plot(atom_list, atomx, atomy, atomz)
  fig1.scene.show_axes = True
  return isovalue

#-------------------------------------------------------------------------------
def real_func_plot_cub(title:str, atoms, amp, mapo, \
                       x, y, z, isovalue:float):
  '''
  plot scalar real data (non-orbital) via structured grid with colormap
  --
  title: title of figure\n
  atoms: atom information\n
  amp: grid value of amplitude function\n
  mapo: grid value of colormapping function\n
  x: x grid points\n
  y: y grid points\n
  z: z grid points\n
  isovalue: isovalue of val\n
  Return: final isovalue
  '''
  if isovalue <= 0.0:
    raise RuntimeError('isovalue should be positive')
  fig = mlab.figure(figure=f'{title}', size=(600, 600),
                    fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
  mlab.clf(fig)
  fig_try = mlab.figure(675)
  amp_smooth = gaussian_filter(amp, sigma=1)
  vmax=mapo.max()
  vmin=mapo.min()
  try:
    mlab.contour3d(x, y, z, amp_smooth, figure=fig_try, contours=[isovalue])
  except Exception as e:
    print(f'Error: {e}')
    if str(e).find('Contour instance must be') >= 0:
      # Extract safe isovalue from error message and scale by 0.5
      safe_val_str = str(e)[str(e).rfind('<=')+3 : str(e).rfind('<=')+14]
      isovalue = float(safe_val_str) * 0.5
      print('------\nIgnore previous error, adjusted isovalue to:', isovalue)
    else:
      raise RuntimeError()
  mlab.close(675)
  # Plot atoms
  # Force active figure to main figure before plotting atoms
  mlab.figure(fig)
  atom_list = atoms['index']
  atomx = atoms['x']
  atomy = atoms['y']
  atomz = atoms['z']
  atoms_plot(atom_list, atomx, atomy, atomz)
  src = mlab.pipeline.scalar_field(x, y, z, amp_smooth, figure=fig)
  # Explicitly name arrays to prevent VTK active scalar loss
  src.image_data.point_data.get_array(0).name = 'amplitude'
  src.image_data.point_data.add_array(mapo.T.ravel())
  src.image_data.point_data.get_array(1).name = 'phase'
  src.update()
  src2 = mlab.pipeline.set_active_attribute(src, point_scalars='amplitude', \
                                            figure=fig)
  contourpha = mlab.pipeline.contour(src2, figure=fig)
  contourpha.filter.contours = [isovalue]
  contourpha2 = mlab.pipeline.set_active_attribute(contourpha, 
                                             point_scalars='phase', figure=fig)
  # rainbow, jet, bwr colormap are usually suitable for non-periodic data,
  # while hsv is suitable for periodic data like phase
  surf = mlab.pipeline.surface(contourpha2, colormap='bwr', figure=fig,
                               vmax=vmax, vmin=vmin)
  surf.module_manager.scalar_lut_manager.data_range = [vmin, vmax]
  # Increase opacity for richer colors, but keep backface_culling
  # so the inner wall doesn't occlude the central atoms.
  surf.actor.property.opacity = 0.65
  surf.actor.property.backface_culling = True
  # Enhance diffuse lighting and reduce specular wash-out
  surf.actor.property.interpolation = 'phong'
  surf.actor.property.ambient  = 0.20
  surf.actor.property.diffuse  = 0.85
  surf.actor.property.specular = 0.20
  surf.actor.property.specular_power = 15
  cb = mlab.colorbar(object=surf, title='sign(lambda_2)*rho\n', \
                     orientation='vertical', nb_labels=5)
  cb.label_text_property.color = (0, 0, 0)
  cb.title_text_property.color = (0, 0, 0)
  cb.scalar_bar.label_format = '%.3f'       # format of colorbar labels
  cb.scalar_bar.unconstrained_font_size = True 
  cb.label_text_property.font_size = 12     # size of colorbar labels
  cb.title_text_property.font_size = 12     # size of colorbar title
  cb.title_text_property.italic = False
  cb.title_text_property.bold = False
  cb.label_text_property.italic = False
  cb.label_text_property.bold = False
  mlab.view(figure=fig)
  #fig.scene.show_axes = True
  return isovalue, fig

#-------------------------------------------------------------------------------
def slice_cub(title:str, atoms, x, y, z, dat, colormap:str,\
              vmax:float, vmin:float) -> None:
  '''
  2D slice based on 3D structured data
  --
  title: title of figure\n
  atoms: atom information\n
  x: x grid points\n
  y: y grid points\n
  z: z grid points\n
  dat: grid value of structured data\n
  colormap: name of the colormap of slice\n
  vmax: colorbar setting\n
  vmin: colorbar setting\n
  Return: None
  '''
  fig = mlab.figure(figure=f'{title}', size=(600, 600),
                    fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
  mlab.clf(fig)
  dat_smooth = gaussian_filter(dat, sigma=1)
  # Plot atoms
  # Force active figure to main figure before plotting atoms
  mlab.figure(fig)
  atom_list = atoms['index']
  atomx = atoms['x']
  atomy = atoms['y']
  atomz = atoms['z']
  atoms_plot(atom_list, atomx, atomy, atomz)
  src = mlab.pipeline.scalar_field(x, y, z, dat_smooth, figure=fig)
  # slice is initially centered in the z-direction.
  plane = mlab.pipeline.image_plane_widget(src, plane_orientation='z_axes', \
                                           slice_index=50, vmin=vmin, vmax=vmax)
  plane.module_manager.scalar_lut_manager.lut_mode = colormap 
  # mlab.outline(color=(0, 0, 0))
  cb = mlab.colorbar(object=plane, 
                   title=f"{title}\n", 
                   orientation='vertical',
                   nb_labels=5,
                   label_fmt="%.3f")
  cb.label_text_property.color = (0, 0, 0)
  cb.title_text_property.color = (0, 0, 0)
  cb.scalar_bar.label_format = '%.3f'       # format of colorbar labels
  cb.scalar_bar.unconstrained_font_size = True 
  cb.label_text_property.font_size = 12     # size of colorbar labels
  cb.title_text_property.font_size = 12     # size of colorbar title
  cb.title_text_property.italic = False
  cb.title_text_property.bold = False
  cb.label_text_property.italic = False
  cb.label_text_property.bold = False
  mlab.view(figure=fig)
  #fig.scene.show_axes = True
  return None

#-------------------------------------------------------------------------------
def real_orb_ampplot_cub(title:str, atoms, val, x, y, z, isovalue:float)->float:
  '''
  plot scalar real structured grid for orbital (amplitude only)
  --
  title: title of figure\n
  atoms: atom information\n
  val: grid value of function\n
  x: x grid points\n
  y: y grid points\n
  z: z grid points\n
  isovalue: isovalue of val\n
  Return: final isovalue
  '''
  if isovalue <= 0.0:
    raise RuntimeError('isovalue should be positive')
  # -------------<plot orb>-------------
  fig1 = mlab.figure(figure=f'{title}', size=(400,400), \
    fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig1)
  # amplitude of MO not amplitude of linear combine of AO
  val_amp = np.abs(val)
  val_smoothed = gaussian_filter(val_amp, sigma=1.0)
  # adjust isovalue if needed
  data = np.asarray(val_smoothed)
  if not np.all(np.isfinite(data)):
    raise ValueError("The smoothed field contains non-finite values.")
  if isovalue > 0:
    positive = data[data > 0]
    if positive.size == 0:
      raise ValueError("The field contains no positive values.")
    maximum = positive.max()
    if isovalue > 0.8 * maximum:
      smallest_positive = positive.min()
      isovalue = 0.5 * (maximum + smallest_positive)
      print(f"Adjusted positive isovalue to {isovalue:g}")
  elif isovalue < 0:
    negative = data[data < 0]
    if negative.size == 0:
      raise ValueError("The field contains no negative values.")
    minimum = negative.min()
    if isovalue < 0.8 * minimum:
      closest_to_zero_negative = negative.max()
      isovalue = 0.5 * (minimum + closest_to_zero_negative)
      print(f"Adjusted negative isovalue to {isovalue:g}")
  src = _structured_grid_source(
    x, y, z, val_smoothed, fig1, name="amplitude"
  )

  contour = mlab.contour3d(x, y, z, val_smoothed, figure=fig1, \
                           contours=[isovalue], opacity=1.0)
  contour.actor.mapper.scalar_visibility = True
  lut = contour.module_manager.scalar_lut_manager.lut
  lut.table = [(180, 180, 180, 180)] * 256
  contour.actor.property.interpolation = 'phong'
  contour.actor.property.ambient  = 0.15
  contour.actor.property.diffuse  = 0.75
  contour.actor.property.specular = 0.35
  contour.actor.property.specular_power = 20
  # -------------<plot atoms>-------------
  if title.startswith('real'): # real space function needs to plot atoms
    atom_list = atoms['index']
    atomx = atoms['x']
    atomy = atoms['y']
    atomz = atoms['z']
    mlab.figure(fig1)
    atoms_plot(atom_list, atomx, atomy, atomz)
  elif title.startswith('momentum'): # momentum space function don't need
    mlab.points3d(0.0, 0.0, 0.0, figure=fig1, \
                  scale_factor=0.3, resolution=20, \
                  color=(0.0,0.0,0.0), opacity=1.0, scale_mode='none')
  fig1.scene.show_axes = True
  return isovalue

#-------------------------------------------------------------------------------
def real_func_ampplot_cub(title:str, atoms, val, x, y, z, isovalue:float)->float:
  '''
  plot scalar real structured grid for real-space functions (amplitude only)
  amplitude can be positive or negative\n
  Note: This is capable of plotting a non-axis-aligned cube.
  --
  title: title of figure\n
  atoms: atom information\n
  val: grid value of function\n
  x: x grid points\n
  y: y grid points\n
  z: z grid points\n
  isovalue: isovalue of val\n
  Return: final isovalue
  '''
  # -------------<plot func>-------------
  fig1 = mlab.figure(figure=f'{title}', size=(400,400), \
    fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig1)
  atoms_plot(
    atoms["index"], atoms["x"], atoms["y"], atoms["z"]
  )
  val_smoothed = gaussian_filter(val, sigma=1.0)
  # adjust isovalue if needed
  data = np.asarray(val_smoothed)
  if not np.all(np.isfinite(data)):
    raise ValueError("The smoothed field contains non-finite values.")
  if isovalue > 0:
    positive = data[data > 0]
    if positive.size == 0:
      raise ValueError("The field contains no positive values.")
    maximum = positive.max()
    if isovalue > 0.8 * maximum:
      smallest_positive = positive.min()
      isovalue = 0.5 * (maximum + smallest_positive)
      print(f"Adjusted positive isovalue to {isovalue:g}")
  elif isovalue < 0:
    negative = data[data < 0]
    if negative.size == 0:
      raise ValueError("The field contains no negative values.")
    minimum = negative.min()
    if isovalue < 0.8 * minimum:
      closest_to_zero_negative = negative.max()
      isovalue = 0.5 * (minimum + closest_to_zero_negative)
      print(f"Adjusted negative isovalue to {isovalue:g}")
  src = _structured_grid_source(
    x, y, z, val_smoothed, fig1, name="amplitude"
  )

  contour_filter = mlab.pipeline.contour(src, figure=fig1)
  contour_filter.filter.contours = [isovalue]
  # Select one isosurface color.
  # surface_hex = "168A89"  # Deep teal
  # surface_hex = "267F9D"  # Blue teal
  # surface_hex = "7565A7"  # Blue violet
  surface_hex = "D89B24"  # Amber
  surface_rgb = tuple(channel / 255. for channel in hex2rgb(surface_hex))
  contour = mlab.pipeline.surface(
    contour_filter,
    figure=fig1,
    color=surface_rgb,
    opacity=1.0
  )
  contour.actor.mapper.scalar_visibility = False
  contour.actor.property.interpolation = 'phong'
  contour.actor.property.ambient = 0.25
  contour.actor.property.diffuse = 0.75
  contour.actor.property.specular = 0.20
  contour.actor.property.specular_power = 30
  contour.actor.property.specular_color = (1.0, 1.0, 1.0)
  fig1.scene.show_axes = True
  return isovalue

#-------------------------------------------------------------------------------
def hex2rgb(hexcolor:str) -> tuple[int,int,int]:
  '''
  color format conversion from Hex to RGB
  --
  Return: RGB color
  '''
  hexcolor = int(hexcolor, base=16) \
    if isinstance(hexcolor, str) else hexcolor
  rgb = ((hexcolor >> 16) & 0xff, (hexcolor >> 8) & 0xff, hexcolor & 0xff)
  return rgb

#-------------------------------------------------------------------------------
def atoms_plot(atom_list, atomx, atomy, atomz, opacity:float=1.0) -> None:
  """
  plot spherical model of atoms in current figure.
  --
  atom_list: atomic number of each atom.\n
  atomx: x coordinate of each atom.\n
  atomy: y coordinate of each atom.\n
  atomz: z coordinate of each atom.\n
  opacity: opacity of the ploted atoms\n
  Returns: None
  """
  # load standard data
  json_address = environ.get('TRESC') + '/vis2c/PubChemElements_all.json'
  with open(json_address,'r',encoding='utf-8') as f:
    std_dat = jsload(f)
  # get information of each atom
  atom_color=[]
  atom_symbol=[]
  atom_radius=[]
  atom_coradius=[]
  num_atom = 0
  for i in atom_list:
    num_atom += 1
    if i <= 0 or i > 118:
      raise RuntimeError('atomic number cannot be recognised.')
    for j in std_dat["Table"]["Row"]:
      if i == int(j['Cell'][0]):
        atom_symbol.append(j["Cell"][1])
        atom_color.append(j["Cell"][4])
        """
        Covalent radius is a suitable criterion but is not ideal for visualization.
        Bondi van der Waals radii for transition metals are too small.
        """
        atom_coradius.append(qcel.covalentradii.get(i)) # covalent (CSD)
        # atom_radius.append(float(j["Cell"][7])/150.)    # Bondi VdW
        atom_radius.append(_display_vdw_radius_pm(i) / 150.0) # Batsanov VdW
        break
  # plot atoms
  fig = mlab.gcf()
  brightness = {1: 0.88, 8: 0.90, 16: 0.85}
  for i in range(num_atom):
    base_rgb = np.array(hex2rgb(atom_color[i]), dtype=float) / 255.0
    factor = brightness.get(int(atom_list[i]), 1.0)
    RGB = tuple(base_rgb * factor)
    mlab.points3d(
      atomx[i], atomy[i], atomz[i], figure=fig,
      scale_factor=atom_radius[i], resolution=20,
      color=RGB, opacity=opacity, scale_mode='none'
    )
    for j in range(i): # atoms will not bond to themselves
      pointi = np.array([atomx[i], atomy[i], atomz[i]])
      pointj = np.array([atomx[j], atomy[j], atomz[j]])
      vec_ij = pointj - pointi
      distance = np.linalg.norm(vec_ij)
      if distance <= 1.15 * (atom_coradius[i]+atom_coradius[j]):
        dir_ij = vec_ij / distance # unit vector
        r_i = atom_radius[i]
        r_j = atom_radius[j]
        ratio = r_i / (r_i + r_j)
        split_point = pointi + dir_ij * (distance * ratio)
        color_i = tuple(np.array(hex2rgb(atom_color[i])) / 255.0)
        color_j = tuple(np.array(hex2rgb(atom_color[j])) / 255.0)
        linei = mlab.plot3d( 
          [pointi[0], split_point[0]],
          [pointi[1], split_point[1]],
          [pointi[2], split_point[2]],
          tube_radius=0.1,
          tube_sides=20,
          color=color_i
        )
        linei.actor.property.opacity = opacity
        linej = mlab.plot3d(
          [split_point[0], pointj[0]],
          [split_point[1], pointj[1]],
          [split_point[2], pointj[2]],
          tube_radius=0.1,
          tube_sides=20,
          color=color_j
        )
        linej.actor.property.opacity = opacity
  fig.scene.background = (1.0,1.0,1.0)

#------------------------------------------------------------------------------
def _structured_grid_source(x, y, z, values, figure,
                            name="amplitude", mapping=None,
                            mapping_name="mapping"):
  """Create a StructuredGrid with physical positions and optional mapping."""
  arrays = (np.asarray(x), np.asarray(y),
            np.asarray(z), np.asarray(values))
  shape = arrays[3].shape
  if len(shape) != 3 or any(a.shape != shape for a in arrays[:3]):
    raise ValueError("x, y, z and values must have identical 3D shapes.")
  if min(shape) < 2:
    raise ValueError("Each grid dimension needs at least 2 points.")
  if any(not np.all(np.isfinite(a)) for a in arrays):
    raise ValueError("Grid coordinates and values must be finite.")
  if mapping is not None:
    mapping = np.asarray(mapping)
    if mapping.shape != shape or not np.all(np.isfinite(mapping)):
      raise ValueError("mapping must match values and contain finite data.")
  points = np.empty(shape + (3,), dtype=np.float64)
  points[..., 0] = arrays[0]
  points[..., 1] = arrays[1]
  points[..., 2] = arrays[2]
  points = points.transpose(2, 1, 0, 3).copy().reshape(-1, 3)
  scalars = np.ascontiguousarray(
    arrays[3].transpose(2, 1, 0)
  ).ravel()
  grid = tvtk.StructuredGrid(dimensions=shape, points=points)
  grid.point_data.scalars = scalars
  grid.point_data.scalars.name = name
  if mapping is not None:
    mapping_values = np.ascontiguousarray(
      mapping.transpose(2, 1, 0)
    ).ravel()
    grid.point_data.add_array(mapping_values)
    grid.point_data.get_array(1).name = mapping_name
  return mlab.pipeline.add_dataset(grid, figure=figure)

#------------------------------------------------------------------------------
def plot_signed_cube_structured(title, atoms, values, x, y, z, isovalue):
  """Plot signed cube values on their actual grid coordinates, including rotations.

  Parameters:
    title: Mayavi window title.
    atoms: Structured atomic records with index, x, y and z fields in the
      same length unit as x, y and z.
    values: Signed 3D scalar field on the cube grid.
    x, y, z: Physical positions of each grid point, shaped like values.
    isovalue: Positive magnitude of the positive/negative contour levels.

  Returns:
    The created Mayavi figure. This function does not call mlab.show().
  """
  if not np.isfinite(isovalue) or isovalue <= 0:
    raise ValueError("isovalue must be positive and finite.")
  values = np.asarray(values)
  vmin, vmax = float(np.min(values)), float(np.max(values))
  if isovalue >= vmax or -isovalue <= vmin:
    raise ValueError(
      f"Both signed levels must lie within the data range [{vmin}, {vmax}]."
    )
  figure = mlab.figure(figure=title, size=(600, 600),
                       fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
  mlab.clf(figure)
  source = _structured_grid_source(x, y, z, values, figure)
  mlab.pipeline.iso_surface(source, contours=[isovalue],
                            color=(0.0, 0.59, 0.55), figure=figure)
  mlab.pipeline.iso_surface(source, contours=[-isovalue],
                            color=(0.0, 0.36, 0.63), figure=figure)
  mlab.figure(figure)
  atoms_plot(atoms['index'], atoms['x'], atoms['y'], atoms['z'])
  mlab.view(figure=figure)
  return figure


#===============================================================================
if __name__ == '__main__':
  print('plot3d.py is not designed to be run directly.')
  print('Please import the functions in this module to plot your data.')

