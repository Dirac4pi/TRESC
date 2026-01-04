'''
plot 3D objects
coding:UTF-8
env:vis2c
'''

from os import environ
import numpy as np
from json import load as jsload
import qcelemental as qcel
qcel.CovalentRadii("ALVAREZ2008")
from mayavi import mlab
from scipy.ndimage import zoom, gaussian_filter

dtol = 0.01 # tolerance of delaunay
maxpoints = 1e5 # max grid points

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

#-------------------------------------------------------------------------------
def cmplx_orb_plot_cub(\
    title:str, atoms, amp, pha, x, y, z, isovalue:float)->float:
  '''
  plot amplitude(contour) and phase(color) of scalar complex structured grid
  --
  !!! for 2c MO, alpha and beta components should be plotted separately.\n
  title: title of figure\n
  atoms: atom information\n
  amp: amplitude of function\n
  pha: phase of function\n
  x: x grid points\n
  y: y grid points\n
  z: z grid points\n
  isovalue: isovalue of amplitude\n
  Return: final isovalue
  '''
  if isovalue <= 0.0:
    raise RuntimeError('isovalue should be positive')
  # -------------<plot amplitude with atoms>-------------
  fig1 = mlab.figure(figure=f'{title}(amplitude)', size=(400,400), \
    fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig1)
  fig_try = mlab.figure(675)
  # amplitude of MO not amplitude of linear combine of AO
  amp_smooth = gaussian_filter(amp, sigma=1)
  try:
    mlab.contour3d(x, y, z, amp_smooth, figure=fig_try, contours=[isovalue])
  except Exception as e:
    print(f'Error: {e}')
    if str(e).find('Contour instance must be') >= 0:
      isovalue = float(str(e)[str(e).rfind('<=')+3 : \
        str(e).rfind('<=')+14]) * .5 
        # factor of 0.5 is usually appropriate
      contouramp = mlab.contour3d(x, y, z, amp_smooth, figure=fig1, \
        contours=[isovalue], opacity=.2, colormap='Greys')
      print('------')
      print('ignore previous error')
    else:
      raise RuntimeError()
  else:
    contouramp = mlab.contour3d(x, y, z, amp_smooth, figure=fig1, \
      contours=[isovalue], opacity=.7, colormap='Greys')
  finally:
    contouramp.actor.property.interpolation = 'phong'
    contouramp.actor.property.ambient  = 0.15
    contouramp.actor.property.diffuse  = 0.75
    contouramp.actor.property.specular = 0.35
    contouramp.actor.property.specular_power = 20
  mlab.close(675)
  atom_list = atoms['index']
  atomx = atoms['x']
  atomy = atoms['y']
  atomz = atoms['z']
  atoms_plot(atom_list, atomx, atomy, atomz)
  #-------------<plot phase>-------------
  fig2 = mlab.figure(figure=f'{title}(phase), isovalue={isovalue}', \
    size=(400,400), fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig2)
  # isosurface of amplitude
  # amplitude of MO are not linear combine of amplitude of AO
  src = mlab.pipeline.scalar_field(x,y,z,amp_smooth)
  # colour mapping of phase
  # use arctan2 to describe complex plane completly
  src.image_data.point_data.add_array(pha.T.ravel())
  src.image_data.point_data.get_array(1).name = 'phase'
  src.update()
  # contour of amplitude
  src2 = mlab.pipeline.set_active_attribute(src, point_scalars='scalar')
  contourpha = mlab.pipeline.contour(src2, figure=fig2)
  contourpha.filter.contours = [isovalue]
  # color map of phase based on contour of amplitude
  contourpha2 = mlab.pipeline.set_active_attribute(contourpha,
    point_scalars='phase', figure=fig2)
  mlab.pipeline.surface(contourpha2, colormap='hsv', opacity=.5, \
    vmax=np.pi, vmin=-np.pi, figure=fig2)
  mlab.colorbar(title='phase', orientation='vertical', nb_labels=5)
  mlab.view(figure=fig1)
  fig1.scene.show_axes = True
  mlab.view(figure=fig2)
  fig2.scene.show_axes = True
  # ---------<synchronise camera>---------
  sync_coe = (fig2.scene.camera.position - fig2.scene.camera.focal_point) \
    / (fig1.scene.camera.position - fig1.scene.camera.focal_point)
  fig1.scene.interactor.add_observer("InteractionEvent", lambda *args: \
    sync_camera(fig1, fig2, sync_coe))
  return isovalue

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
    surface = mlab.pipeline.surface(contour, figure=fig1, opacity=.7,\
                                    colormap='Greys')
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
    surf = mlab.pipeline.surface(contour_phase, colormap='jet', opacity=0.5)
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
  contourpos = mlab.contour3d(x, y, z, val_smoothed, contours=[isovalue], \
                              opacity=.7, colormap='winter')
  contourpos.actor.property.interpolation = 'phong'
  contourpos.actor.property.ambient  = 0.15
  contourpos.actor.property.diffuse  = 0.75
  contourpos.actor.property.specular = 0.35
  contourpos.actor.property.specular_power = 20
  contourneg = mlab.contour3d(x, y, z, val_smoothed, contours=[-isovalue], \
                              opacity=.7, colormap='autumn')
  contourneg.actor.property.interpolation = 'phong'
  contourneg.actor.property.ambient  = 0.15
  contourneg.actor.property.diffuse  = 0.75
  contourneg.actor.property.specular = 0.35
  contourneg.actor.property.specular_power = 20
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
def real_orb_ampplot_cub(title:str, atoms, val, x, y, z, isovalue:float)->float:
  '''
  plot scalar real structured grid (amplitude only)
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
  if isovalue > 0.5*val_smoothed.max():
    isovalue = 0.5 * (val_smoothed.max() + val_smoothed.min())
    print(f"Adjust isovalue to {isovalue} (mid-range)")
  elif isovalue < 2*val_smoothed.min():
    isovalue = 0.5 * (val_smoothed.max() + val_smoothed.min())
    print(f"Adjust isovalue to {isovalue} (mid-range)")
  contour = mlab.contour3d(x, y, z, val_smoothed, figure=fig1, \
                           contours=[isovalue], opacity=.7, colormap='Greys')
  contour.actor.property.interpolation = 'phong'
  contour.actor.property.ambient  = 0.15
  contour.actor.property.diffuse  = 0.75
  contour.actor.property.specular = 0.35
  contour.actor.property.specular_power = 20
  # -------------<plot atoms>-------------
  if (title.startswith('real')):
    atom_list = atoms['index']
    atomx = atoms['x']
    atomy = atoms['y']
    atomz = atoms['z']
    mlab.figure(fig1)
    atoms_plot(atom_list, atomx, atomy, atomz)
  elif (title.startswith('momentum')):
    mlab.points3d(0.0, 0.0, 0.0, figure=fig1, \
                  scale_factor=0.3, resolution=20, \
                  color=(0.0,0.0,0.0), opacity=1.0, scale_mode='none')
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
def atoms_plot(atom_list, atomx, atomy, atomz) -> None:
  '''
  plot spherical model of atoms in current figure.
  --
  atom_list: atomic number of each atom.\n
  atomx: x coordinate of each atom.\n
  atomy: y coordinate of each atom.\n
  atomz: z coordinate of each atom.\n
  Returns: None
  '''
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
        atom_coradius.append(qcel.covalentradii.get(i)) # covalent (CSD)
        atom_radius.append(float(j["Cell"][7])/150.) # Van der Waal
        break
  # plot atoms
  fig = mlab.gcf()
  for i in range(num_atom):
    RGB = hex2rgb(atom_color[i])
    R = RGB[0] / 255.
    G = RGB[1] / 255.
    B = RGB[2] / 255.
    mlab.points3d(atomx[i], atomy[i], atomz[i], figure=fig, \
      scale_factor=atom_radius[i], resolution=20, \
        color=(R,G,B), opacity=1., scale_mode='none')
    for j in range (num_atom):
      pointi = np.array([atomx[i], atomy[i], atomz[i]])
      pointj = np.array([atomx[j], atomy[j], atomz[j]])
      distance = np.linalg.norm(pointi - pointj)
      if distance <= 1.15 * (atom_coradius[i]+atom_coradius[j]):
        lx = np.linspace(pointi[0], pointj[0], 2)
        ly = np.linspace(pointi[1], pointj[1], 2)
        lz = np.linspace(pointi[2], pointj[2], 2)
        mlab.plot3d(lx, ly, lz, tube_radius=0.05, \
                    color=(0.3,0.3,0.3), tube_sides=20)
  fig.scene.background = (1.0,1.0,1.0)

