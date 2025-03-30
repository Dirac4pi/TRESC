'''
python module in vis2c
plot 3D objects
coding:UTF-8
env:2cvis
'''

from os import environ
import numpy as np
from json import load as jsload
import qcelemental as qcel
qcel.CovalentRadii("ALVAREZ2008")
from mayavi import mlab


def sync_camera(fig1, fig2, sync_coe:float) -> None:
  """
  Synchronise camera parameters from fig1.scene to fig2.scene.
  --
  """
  fig2.scene.camera.position = sync_coe*(fig1.scene.camera.position - \
  fig1.scene.camera.focal_point) + fig2.scene.camera.focal_point
  fig2.scene.camera.view_angle = fig1.scene.camera.view_angle
  fig2.scene.camera.view_up = fig1.scene.camera.view_up
  fig2.scene.render()

def cmplx_orb_plot_cub(title:str, atoms, real, img, x, y, z, isovalue:float)->float:


  '''
  plot amplitude(contour) and phase(color) of scalar complex structured grid.
  --
  !!! for 2c MO, alpha and beta components should be plotted separately.\n
  title: title of figure\n
  atoms: atom information\n
  real: real part of function\n
  img: imaginary part of function\n
  x: x grid points\n
  y: y grid points\n
  z: z grid points\n
  isovalue: isovalue of amplitude\n
  Return: final isovalue
  '''
  # -------------<plot amplitude with atoms>-------------
  fig1 = mlab.figure(figure=f'{title}(amplitude)', size=(400,400), \
    fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig1)
  fig_try = mlab.figure(675)
  # amplitude of MO not amplitude of linear combine of AO
  amplitude = np.sqrt(np.add(np.square(real), np.square(img)))
  try:
    mlab.contour3d(x, y, z, amplitude, figure=fig_try, contours=[isovalue])
  except Exception as e:
    print(f'Error: {e}')
    if str(e).find('Contour instance must be') >= 0:
      isovalue = float(str(e)[str(e).rfind('<=')+3 : \
        str(e).rfind('<=')+14]) * .5 
        # factor of 0.5 is usually appropriate
      contouramp = mlab.contour3d(x, y, z, amplitude, figure=fig1, \
        contours=[isovalue], opacity=.2, colormap='Greys')
      print('------')
      print('ignore previous error')
    else:
      raise RuntimeError()
  else:
    contouramp = mlab.contour3d(x, y, z, amplitude, figure=fig1, \
      contours=[isovalue], opacity=.2, colormap='Greys')
  finally:
    contouramp.actor.property.line_width = 1
    contouramp.actor.property.edge_visibility = True
  mlab.close(675)
  atom_list = atoms['index']
  atomx = atoms['x']
  atomy = atoms['y']
  atomz = atoms['z']
  atoms_plot(atom_list, atomx, atomy, atomz)
  
  # -------------<plot phase>-------------
  fig2 = mlab.figure(figure=f'{title}(phase), isovalue={isovalue}', \
    size=(400,400), fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig2)
  # isosurface of amplitude
  # amplitude of MO are not linear combine of amplitude of AO
  amplitude = np.sqrt(np.add(np.square(real), np.square(img)))
  src = mlab.pipeline.scalar_field(amplitude)
  # colour mapping of phase
  # use arctan2 to describe complex plane completly
  phase = np.arctan2(real,img)
  src.image_data.point_data.add_array(phase.T.ravel())
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

def cmplx_orb_plot_mog(title:str, atoms, real, img, x, y, z, isovalue:float)->float:
  '''
  plot amplitude(contour) and phase(color) of scalar complex unstructured grid.
  --
  !!! for 2c MO, alpha and beta components should be plotted separately.\n
  title: title of figure\n
  atoms: atom information\n
  real: real part of function\n
  img: imaginary part of function\n
  x: x grid points\n
  y: y grid points\n
  z: z grid points\n
  isovalue: isovalue of amplitude\n
  Return: final isovalue
  '''
  print('Unstructured data need to be Delaunay triangulated before extracting')
  print('isosurface, the time consumed increases significantly...')
  print('...please be patient...')
  n_points = len(x)
  if n_points > 1e5:  # downsampling triggered at over 100,000 points
    sample_size = int(1e5)
    indices = np.random.choice(n_points, sample_size, replace=False)
    x = x[indices]
    y = y[indices]
    z = z[indices]
    real = real[indices]
    img = img[indices]
    print(f"data points sample from {n_points} to {sample_size}")

  # -------------<plot amplitude with atoms>-------------
  fig1 = mlab.figure(figure=f'{title}(amplitude)', size=(400,400), \
    fgcolor=(0,0,0), bgcolor=(1,1,1))
  mlab.clf(fig1)
  # amplitude of MO not amplitude of linear combine of AO
  amplitude = np.sqrt(np.add(np.square(real), np.square(img)))
  # dynamic adjustment of thresholds
  if isovalue > amplitude.max():
    isovalue = 0.5 * amplitude.max()
    print(f"Adjust isovalue to {isovalue} (50% of max)")
  elif isovalue < amplitude.min():
    isovalue = 0.5 * (amplitude.max() + amplitude.min())
    print(f"Adjust isovalue to {isovalue} (mid-range)")


  points = mlab.pipeline.scalar_scatter(x, y, z, amplitude)
  delaunay = mlab.pipeline.delaunay3d(points)
  delaunay.filter.tolerance = 0.01
  contour = mlab.pipeline.contour(delaunay)
  try:
    contour.filter.contours = [isovalue]
    surface = mlab.pipeline.surface(contour, figure=fig1, opacity=0.2,\
                                    colormap='Greys')
    surface.actor.property.edge_visibility = True
    surface.actor.property.line_width = 0.1
  except Exception as e:
    raise RuntimeError(f'Error: {e}')
  else:
    surface = mlab.pipeline.surface(contour, figure=fig1, opacity=0.2,\
                                    colormap='Greys')
    surface.actor.property.edge_visibility = True
    surface.actor.property.line_width = 0.1
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
  phase = np.arctan2(real,img)
  src.mlab_source.dataset.point_data.add_array(phase.ravel())
  src.mlab_source.dataset.point_data.get_array(1).name = 'phase'
  src.update()

  # build pipeline: triangulation -> isosurface extraction
  delaunay = mlab.pipeline.delaunay3d(src)
  contour = mlab.pipeline.contour(delaunay)
  contour.filter.contours = [isovalue]

  # color mapping to phase
  contour_phase = mlab.pipeline.set_active_attribute(contour, point_scalars='phase')
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

def hex2rgb(hexcolor:str) -> tuple[int,int,int]:
  '''
  colour format conversion from Hex to RGB.
  --
  Return: RGB color
  '''
  hexcolor = int(hexcolor, base=16) \
    if isinstance(hexcolor, str) else hexcolor
  rgb = ((hexcolor >> 16) & 0xff, (hexcolor >> 8) & 0xff, hexcolor & 0xff)
  return rgb

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
  json_address = environ.get('TRESC') + '\\vis2c\PubChemElements_all.json'
  with open(json_address,'r',encoding='utf-8') as f:
    std_dat = jsload(f)
  # get information of each atom
  atom_color=[]
  atom_symbol=[]
  atom_radius=[]
  num_atom = 0
  for i in atom_list:
    num_atom += 1
    if i <= 0 or i > 118:
      raise RuntimeError('atomic number cannot be recognised.')
    for j in std_dat["Table"]["Row"]:
      if i == int(j['Cell'][0]):
        atom_symbol.append(j["Cell"][1])
        atom_color.append(j["Cell"][4])
        #atom_radius.append(qcel.covalentradii.get(i)) # covalent (CSD)
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
        color=(R, G, B), opacity=1., scale_mode='none')
    #mlab.text3d(atomx[i], atomy[i], atomz[i], \
    #    atom_symbol[i].strip(), line_width=.7, \
    #        color=(0,0,0), figure=fig, scale=(.3,.3,.3))

