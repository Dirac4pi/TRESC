#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
visualization of real scalar orbital with structured
grid data in real space and momentum space, by the using
of dual space as an alternative to phase descriptions.
author:Dirac4pi
env:vis2c
'''

import viskit as vk
import plot3d as p3
from os import path
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item
import re

#-------------------------------------------------------------------------------
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
  check_index = lambda s: bool(re.fullmatch(r'[1-9]\d*[AB]', s, re.IGNORECASE))
  if check_index(index) is False:
    raise RuntimeError(f'index should be 24A, 25b, etc.')
  if float(isovalue) <= 0:
    raise RuntimeError(f'isovalue should be a positive float')
  # generate cube file (with real cube and momentum cube)
  print('calling TRESC:')
  vk.call_executable(['tshell.sh', '-pro1c', molden, index])
  # load data from real space cube file
  cub = index+'-real.cub'
  print(f'loading real space grid data from {cub}...')
  atoms, mat = vk.load_cube(cub)
  print(f'done')
  # grid data
  rx = mat['x']
  ry = mat['y']
  rz = mat['z']
  rval = mat['value']
  # load data from momentum space cube file
  cub = index+'-mmt.cub'
  print(f'loading momentum grid data from {cub}...')
  atoms, mat = vk.load_cube(cub)
  print('done')
  # grid data
  px = mat['x']
  py = mat['y']
  pz = mat['z']
  pval = mat['value']
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

#===============================================================================
if __name__ == "__main__":
  from sys import argv
  if len(argv) not in [3, 4]:
    print('Usage: pro1c.py source.molden.input orb_index_with_spin (isovalue)')
    print('e.g. pro1c.py C6H6.molden.input 15a')
  else:
    source_file = argv[1]
    orb_index = argv[2]
    if len(argv) == 4:
      isovalue = float(argv[3])
      pro1c(source_file, orb_index, isovalue)
    else:
      pro1c(source_file, orb_index)