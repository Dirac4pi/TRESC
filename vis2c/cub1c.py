#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
visualization of real scalar orbital with structured grid data
author: Dirac4pi
env:vis2c
'''

import viskit as vk
import plot3d as p3
import re
from os import path
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item

#-------------------------------------------------------------------------------
def cub1c(molden:str, index:str, isovalue:float=0.05, slice:bool=False)->None:
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
  check_index = lambda s: bool(re.fullmatch(r'[1-9]\d*[AB]', s, re.IGNORECASE))
  if check_index(index) is False:
    raise RuntimeError(f'index should be 24A, 25b, etc.')
  if float(isovalue) <= 0:
    raise RuntimeError(f'isovalue should be a positive float')
  # generate cube file
  print('calling TRESC:')
  vk.call_executable(['tshell.sh', '-cub1c', molden, index])
  # load data from cube file
  cub = index+'-real.cub'
  print(f'loading real space grid data from {cub}...')
  atoms, mat = vk.load_cube(cub)
  print('done')
  # grid data
  x = mat['x']
  y = mat['y']
  z = mat['z']
  val = mat['value']
  # plot the real scalar orbital
  print('plotting...')
  if slice:
    p3.slice_cub(f"Orbit Isosurfase Slice", atoms, x, y, z,val, \
                 'viridis', val.max(), val.min())
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

#===============================================================================
if __name__ == "__main__":
  from sys import argv
  if len(argv) not in [3, 4]:
    print('Usage: cub1c.py source.molden.input orb_index_with_spin (-slice)')
    print('e.g. cub1c.py C6H6.molden.input 15a -slice')
  else:
    source_file = argv[1]
    orb_index = argv[2]
    if len(argv) == 4:
      if argv[3].lower() == '-slice':
        cub1c(source_file, orb_index, isovalue=0.05, slice=True)
      else:
        exit(f'unkown arg: {argv[3]}')
    else:
      cub1c(source_file, orb_index, isovalue=0.05, slice=False)