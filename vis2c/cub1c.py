'''
visualization of real scalar orbital with structured grid data
coding:UTF-8
env:vis2c
'''

import viskit as vk
import plot3d as p3
from os import path
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item

def cub1c(molden:str, index:str, isovalue:float=0.05)->None:
  '''
  visualization of real scalar orbital with structured grid data
  --
  molden: title of .molden.input file.\n
  index: index of MO with spin.\n
  isovalue: isovalue of amplitude of MO.\n
  Returns: None
  '''
  if not molden.endswith('.molden.input') and not molden.endswith('.molden'):
    if not path.exists(molden.strip('.molden.input')):
      if not path.exists(molden.strip('.molden')):
        raise RuntimeError(f"can't find {molden}.molden.input \
                           or {molden}.molden")
      else:
        molden = molden + '.molden'
    else:
      molden = molden + '.molden.input'
  if not path.exists(molden):
    raise RuntimeError(f"can't find {molden}")
  #-----------------------------------------------------------------------------
  # generate cube file
  print('generating cub file...')
  vk.call_executable(['tshell.sh', '-cub1c', molden, index])
  #-----------------------------------------------------------------------------
  # load data from cube file
  print('loading grid data from cub file...')
  cub = index+'-real.cub'
  atoms, mat, nmo = vk.load_cube(cub, 1)
  if nmo == 0:
    raise RuntimeError('no orbital info in cube file '+cub)
  elif nmo != 1:
    raise RuntimeError('cube file '+cub+' should contain only 1 orbital.')
  # grid data
  x = mat['x']
  y = mat['y']
  z = mat['z']
  val = mat['value']
  #-----------------------------------------------------------------------------
  # plot the real scalar orbital
  print('plotting...')
  isovl = p3.real_orb_plot_cub(index, atoms, val, x, y, z, isovalue)
  # UI regulation isovalue
  class IsoValueController(HasTraits):
    isov = Float(isovl, desc="isov", auto_set=False, enter_set=True)
    @observe('isov')
    def update_isov(self,event):
      new_isov = event.new
      mlab.close(f'{index}')
      p3.real_orb_plot_cub(index, atoms, val, x, y, z, new_isov)
      mlab.draw()
      mlab.view()
    viewalpha = View(
      Item('isov', label="isovalue", show_label=True),
      width=300,
      height=200,
      resizable=True
    )
  controller = IsoValueController()
  controller.configure_traits()

if __name__ == "__main__":
  from sys import argv
  if len(argv) != 3:
    print('cub1c source.molden.input orb_index_with_spin (isovalue)')
    print('e.g. cub1c C6H6.molden.input 15a')
  else:
    cub1c(*argv[1:])