#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
visualization of real-space function (e.g., electron density) and use another
function for color mapping (e.g., electrostatic potential, IGMH).
author:Dirac4pi
env:vis2c
'''

import viskit as vk
import plot3d as p3
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item
from os import path

doslice = False

#-------------------------------------------------------------------------------
def real3d(ampcube:str, mapcube:str, isovalue:float=0.05,\
           slice:bool=False)->None:
  '''
  visualization of real-space function with structured grid data
  --
  ampcube: address of amplitude.cub, contains amplitude data\n
  mapcube: address of colormap.cub, contains mapping data\n
  isovalue: isovalue of amplitude.\n
  slice: perform a slice.\n
  Returns: None
  '''
  if isovalue <= 0:
    raise RuntimeError(f'isovalue should be a positive float')
  if not ampcube.endswith('.cub') and not ampcube.endswith('.cube'):
    raise RuntimeError(f'ampcube seems not a cube file')
  if not mapcube.endswith('.cub') and not mapcube.endswith('.cube'):
    raise RuntimeError(f'mapcube seems not a cube file')
  ampbase = path.splitext(path.basename(ampcube))[0]
  mapbase = path.splitext(path.basename(mapcube))[0]
  # load data from cube files
  print(f'loading real space grid amplitude data from {ampcube}...')
  atoms, ampdat = vk.load_cube(ampcube)
  # grid data, ampdat and mapdat have the same grid
  x = ampdat['x']
  y = ampdat['y']
  z = ampdat['z']
  ar = ampdat['value']
  print(f'loading real space grid mapping data from {mapcube}...')
  atoms, mapdat = vk.load_cube(mapcube)
  mr = mapdat['value']
  # plot the real space function with colormap
  print('plotting...')
  if slice:
    p3.slice_cub(f"{ampbase}", atoms, x, y, z, ar,'viridis',ar.max(),ar.min())
    # p3.slice_cub(f"{mapbase}", atoms, x, y, z, mr, 'viridis', 0.05, -0.05)
  isovl, fig= p3.real_func_plot_cub(f'amplitude: {ampbase}, phase: {mapbase}', \
                                    atoms, ar, mr, x, y, z, isovalue)
  # UI regulation isovalue
  class IsoValueController(HasTraits):
    isov = Float(isovl, desc="isov", auto_set=False, enter_set=True)
    @observe('isov')
    def update_isov(self,event):
      new_isov = event.new
      mlab.close(f'amplitude: {ampbase}, phase: {mapbase}')
      isovl, fig= p3.real_func_plot_cub(\
    f'amplitude: {ampbase}, phase: {mapbase}', atoms, ar, mr, x, y, z, new_isov)
      mlab.draw()
      mlab.view()
    view = View(
      Item('isov', label="isovalue", show_label=True),
      width=300,
      height=200,
      resizable=True
    )
  controller = IsoValueController()
  controller.configure_traits()

#===============================================================================
if __name__ == "__main__":
  from sys import argv
  if len(argv) not in [3, 4]:
    print('Usage: real3d.py ampcube.cub mapcube.cub (-slice)')
    print('e.g. real3d.py amplitude.cub colormap.cub -slice')
  else:
    amp_cub = argv[1]
    map_cub = argv[2]
    if len(argv) == 4:
      if argv[3].lower() == '-slice':
        real3d(amp_cub, map_cub, isovalue=0.05, slice=True)
      else:
        exit(f'unkown arg: {argv[3]}')
    else:
      real3d(amp_cub, map_cub, isovalue=0.05, slice=False)

