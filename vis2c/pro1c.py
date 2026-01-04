'''
visualization of real scalar orbital with structured
grid data in real space and momentum space, by the using
of dual space as an alternative to phase descriptions.
coding:UTF-8
env:vis2c
'''

import viskit as vk
import plot3d as p3
from os import path
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item

def pro1c(molden:str, index:str, isovalue:float=0.05)->None:
  '''
  visualization of real scalar orbital with structured
  grid data in real space and momentum space, by the using
  of dual space as an alternative to phase descriptions.
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
  # generate cube file (with real cube and momentum cube)
  print('generating cube file...')
  vk.call_executable(['tshell.sh', '-cub1c', molden, index])
  #-----------------------------------------------------------------------------
  # load data from real space cube file
  print('loading grid data from real space cub file...')
  cub = index+'-real.cub'
  atoms, mat, nmo = vk.load_cube(cub, 1)
  if nmo == 0:
    raise RuntimeError('no orbital info in cube file '+cub)
  elif nmo != 1:
    raise RuntimeError('cube file '+cub+' should contain only 1 orbital.')
  # grid data
  rx = mat['x']
  ry = mat['y']
  rz = mat['z']
  rval = mat['value']
  print('done')
  #-----------------------------------------------------------------------------
  # load data from momentum space cube file
  print('loading grid data from momentum space cub file...')
  cub = index+'-mmt.cub'
  atoms, mat, nmo = vk.load_cube(cub, 1)
  if nmo == 0:
    raise RuntimeError('no orbital info in cube file '+cub)
  elif nmo != 1:
    raise RuntimeError('cube file '+cub+' should contain only 1 orbital.')
  # grid data
  px = mat['x']
  py = mat['y']
  pz = mat['z']
  pval = mat['value']
  print('done')
  #-----------------------------------------------------------------------------
  # plot the real scalar orbital in real space and momentum space
  print('plotting in real space and momentum space ...')
  risovl = p3.real_orb_ampplot_cub('real space', atoms, rval, \
                                  rx, ry, rz, isovalue)
  pisovl = p3.real_orb_ampplot_cub('momentum space', atoms, pval, \
                                  px, py, pz, isovalue)
  # UI regulation isovalue
  class IsoValueController(HasTraits):
    risov = Float(risovl, desc="risov", auto_set=False, enter_set=True)
    pisov = Float(pisovl, desc="pisov", auto_set=False, enter_set=True)
    @observe('risov')
    def update_risov(self,event):
      new_risov = event.new
      mlab.close('real space')
      p3.real_orb_ampplot_cub('real space', atoms, rval, rx, ry, rz, new_risov)
      mlab.draw()
      mlab.view()
    @observe('pisov')
    def update_pisov(self,event):
      new_pisov = event.new
      mlab.close('momentum space')
      p3.real_orb_ampplot_cub('momentum space', atoms, pval, \
                              px, py, pz, new_pisov)
      mlab.draw()
      mlab.view()
    viewr = View(
      Item('risov', label="real space isovalue", show_label=True),
      width=400,
      height=300,
      resizable=True
    )
    viewp = View(
      Item('pisov', label="momentum space isovalue", show_label=True),
      width=400,
      height=300,
      resizable=True
    )
  controller = IsoValueController()
  controller.configure_traits()

if __name__ == "__main__":
  from sys import argv
  if len(argv) != 3:
    print('pro1c source.molden.input orb_index_with_spin (isovalue)')
    print('e.g. pro1c C6H6.molden.input 15a')
  else:
    pro1c(*argv[1:])