'''
visualization of complex spinor orbital with structured grid data
coding:UTF-8
env:vis2c
'''

import viskit as vk
import plot3d as p3
from os import path
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item
from numpy import sqrt, square, arctan2

# True: alpha and beta visualization as phase, False: as components
spin_phase = True

def cub2c(moldendir:str, index:int, isovalue:float=0.05)->None:
  '''
  visualization of complex spinor orbital with structured grid data
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
  print(f'spin_phase = {spin_phase}')
  #-----------------------------------------------------------------------------
  # generate cube files
  print('generating cub files...')
  vk.call_executable(['tshell.sh', '-cub2c', moldendir, str(index)])
  #-----------------------------------------------------------------------------
  # load data from cube files
  print('loading grid data from cub files...')
  real = index+'-realreal.cub'
  img  = index+'-realimg.cub'
  atoms, real_mat_alpha, nmo = vk.load_cube(real, 1)
  if nmo == 0:
    raise RuntimeError('no orbital info in cube file '+real)
  elif nmo != 2:
    raise RuntimeError('cube file '+real+' should contain both \
      alpha and beta orbitals.')
  # grid data
  x = real_mat_alpha['x']
  y = real_mat_alpha['y']
  z = real_mat_alpha['z']
  ar = real_mat_alpha['value']
  atoms, real_mat_beta, nmo = vk.load_cube(real, 2)
  br = real_mat_beta['value']
  atoms, img_mat_alpha, nmo = vk.load_cube(img, 1)
  if nmo == 0:
    raise RuntimeError('no orbital info in cube file '+img)
  elif nmo != 2:
    raise RuntimeError('cube file '+img+' should contain both \
      alpha and beta orbitals.')
  ai = img_mat_alpha['value']
  atoms, img_mat_beta, nmo = vk.load_cube(img, 2)
  bi = img_mat_beta['value']
  #-----------------------------------------------------------------------------
  # plot the complex 2c orbital
  print('plotting...')
  if spin_phase:
    sa = sqrt(square(ar) + square(ai))
    sb = sqrt(square(br) + square(bi))
    mod = sqrt(square(ar) + square(ai) + square(br) + square(bi))
    sp = arctan2(sa, sb)
    mr = ar + br
    mi = ai + bi
    mp = arctan2(mr, mi)
    isovla = p3.cmplx_orb_plot_cub('spin', atoms, mod, sp, x, y, z, isovalue)
    isovlb = p3.cmplx_orb_plot_cub('orb', atoms, mod, mp, x, y, z, isovalue)
  else:
    moda = sqrt(square(ar) + square(ai))
    pha = arctan2(ar, ai)
    modb = sqrt(square(br) + square(bi))
    phb = arctan2(br, bi)
    isovla = p3.cmplx_orb_plot_cub('alpha', atoms, moda, pha, x, y, z, isovalue)
    isovlb = p3.cmplx_orb_plot_cub('beta', atoms, modb, phb, x, y, z, isovalue)
  # UI regulation isovalue
  class IsoValueController(HasTraits):
    isova = Float(isovla, desc="isova", auto_set=False, enter_set=True)
    isovb = Float(isovlb, desc="isovb", auto_set=False, enter_set=True)
    @observe('isova')
    def update_isova(self,event):
      old_isov = event.old
      new_isov = event.new
      if spin_phase:
        mlab.close('spin(amplitude)')
        mlab.close(f'spin(phase), isovalue={old_isov}')
        p3.cmplx_orb_plot_cub('spin', atoms, mod, sp, x, y, z, new_isov)
      else:
        mlab.close('alpha(amplitude)')
        mlab.close(f'alpha(phase), isovalue={old_isov}')
        p3.cmplx_orb_plot_cub('alpha', atoms, moda, pha, x, y, z, new_isov)
      mlab.draw()
      mlab.view()
    @observe('isovb')
    def update_isovb(self,event):
      old_isov = event.old
      new_isov = event.new
      if spin_phase:
        mlab.close('orb(amplitude)')
        mlab.close(f'orb(phase), isovalue={old_isov}')
        p3.cmplx_orb_plot_cub('orb', atoms, mod, mp, x, y, z, new_isov)
      else:
        mlab.close('beta(amplitude)')
        mlab.close(f'beta(phase), isovalue={old_isov}')
        p3.cmplx_orb_plot_cub('beta', atoms, modb, phb, x, y, z, new_isov)
      mlab.draw()
      mlab.view()
    viewalpha = View(
      Item('isova', label="isovalue(a)", show_label=True),
      width=300,
      height=200,
      resizable=True
    )
    viewbeta = View(
      Item('isovb', label="isovalue(b)", show_label=True),
      width=300,
      height=200,
      resizable=True
    )
  controller = IsoValueController()
  controller.configure_traits()

if __name__ == "__main__":
  from sys import argv
  if len(argv) != 3:
    print('cub2c source.molden.d orb_index (isovalue)')
    print('e.g. cub2c C6H6.molden.d 30')
  else:
    cub2c(*argv[1:])

