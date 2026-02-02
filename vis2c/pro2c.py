'''
visualization of complex spinor orbital with structured
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
from numpy import sqrt, square

#-------------------------------------------------------------------------------
def pro2c(moldendir:str, index:int, isovalue:float=0.05)->None:
  '''
  visualization of complex spinor orbital with structured
  grid data in real space and momentum space, by the using
  of dual space as an alternative to phase descriptions.
  --
  moldendir: address of xxx.molden.d, contains molden file of 2c MOs\n
  index: index of MO.\n
  isovalue: isovalue of amplitude.\n
  Returns: None
  '''
  isovalue = float(isovalue)
  if not moldendir.endswith('.molden.d'):
    raise RuntimeError(f'moldendir should be xxx.molden.d')
  if not path.exists(moldendir):
    raise RuntimeError(f"can't find dir {moldendir}")
  # generate cube files (with real cube and momentum cube)
  print('calling TRESC:')
  vk.call_executable(['tshell.sh', '-pro2c', moldendir, str(index)])
  # load data from real space cube file
  alphareal = index+'-rar.cub'
  alphaimg  = index+'-rai.cub'
  betareal = index+'-rbr.cub'
  betaimg  = index+'-rbi.cub'
  print(f'loading real space grid alpha real part data from {alphareal}...')
  atoms, real_mat_alpha = vk.load_cube(alphareal)
  print('done')
  # grid data
  rx = real_mat_alpha['x']
  ry = real_mat_alpha['y']
  rz = real_mat_alpha['z']
  rar = real_mat_alpha['value']
  print(f'loading real space grid beta real part data from {betareal}...')
  atoms, real_mat_beta = vk.load_cube(betareal)
  print('done')
  rbr = real_mat_beta['value']
  print(f'loading real space grid alpha imaginary part data from {alphaimg}...')
  atoms, img_mat_alpha = vk.load_cube(alphaimg)
  print('done')
  rai = img_mat_alpha['value']
  print(f'loading real space grid beta imaginary part data from {betaimg}...')
  atoms, img_mat_beta = vk.load_cube(betaimg)
  print('done')
  rbi = img_mat_beta['value']
  # load data from momentum space cube file
  alphareal = index+'-mar.cub'
  alphaimg  = index+'-mai.cub'
  betareal = index+'-mbr.cub'
  betaimg  = index+'-mbi.cub'
  print(f'loading momentum grid alpha real part data from {alphareal}...')
  atoms, real_mat_alpha = vk.load_cube(alphareal)
  print('done')
  # grid data
  px = real_mat_alpha['x']
  py = real_mat_alpha['y']
  pz = real_mat_alpha['z']
  par = real_mat_alpha['value']
  print(f'loading momentum grid beta real part data from {betareal}...')
  atoms, real_mat_beta = vk.load_cube(betareal)
  print('done')
  pbr = real_mat_beta['value']
  print(f'loading momentum grid alpha imaginary part data from {alphaimg}...')
  atoms, img_mat_alpha = vk.load_cube(alphaimg)
  print('done')
  pai = img_mat_alpha['value']
  print(f'loading momentum grid beta imaginary part data from {betaimg}...')
  atoms, img_mat_beta = vk.load_cube(betaimg)
  print('done')
  pbi = img_mat_beta['value']
  # plot the complex spinor orbital in real space and momentum space
  print('plotting in real space and momentum space ...')
  rmod = sqrt(square(rar) + square(rai) + square(rbr) + square(rbi))
  pmod = sqrt(square(par) + square(pai) + square(pbr) + square(pbi))
  risovl = p3.real_orb_ampplot_cub('real space', atoms, rmod, \
                                   rx, ry, rz, isovalue)
  pisovl = p3.real_orb_ampplot_cub('momentum space', atoms, pmod, \
                                   px, py, pz, isovalue)
  # UI regulation isovalue
  class IsoValueController(HasTraits):
    risov = Float(risovl, desc="risov", auto_set=False, enter_set=True)
    pisov = Float(pisovl, desc="pisov", auto_set=False, enter_set=True)
    @observe('risov')
    def update_risov(self,event):
      new_risov = event.new
      mlab.close('real space')
      p3.real_orb_ampplot_cub('real space', atoms, rmod, rx, ry, rz, new_risov)
      mlab.draw()
      mlab.view()
    @observe('pisov')
    def update_pisov(self,event):
      new_pisov = event.new
      mlab.close('momentum space')
      p3.real_orb_ampplot_cub('momentum space', atoms, pmod, \
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
  if len(argv) != 3:
    print('pro2c source.molden.d orb_index (isovalue)')
    print('e.g. pro2c C6H6.molden.d 30')
  else:
    pro2c(*argv[1:])