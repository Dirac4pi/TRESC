'''
visualizing 2-component complex orbital with unstructured data
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

def mog2c(moldendir:str, index:int, isovalue:float=0.05)->None:
  '''
  2-component complex MO visualization, unstructured grid net.
  --
  !!! 2 molden format file contain (alpha real & beta real),\n
  (alpha imagine & beta imagine) of a selected 2-component complex orbital.\n
  moldendir: address of xxx.molden.d, contains molden file of 2c MOs\n
  index: index of MO.\n
  isovalue: isovalue of amplitude of MO.\n
  Returns: None
  '''
  isovalue = float(isovalue)
  if not moldendir.endswith('.molden.d'):
    raise RuntimeError(f'moldendir should be xxx.molden.d')
  else:
    basename = moldendir.strip('.molden.d')
  if not path.exists(moldendir):
    raise RuntimeError(f"can't find dir {moldendir}")
  print(f'spin_phase = {spin_phase}')
  #-----------------------------------------------------------------------------
  # generate mog files
  print('generating mog files...')
  vk.call_executable(['tshell.sh', '-mog2c', moldendir, str(index)])
  #-----------------------------------------------------------------------------
  # load geometry
  print('loading grid data from mog files...')
  title = moldendir+'/'+basename+'-realpart1.molden.input'
  atoms = vk.load_geometry_molden(title)
  #-----------------------------------------------------------------------------
  # load mog
  x = vk.load_binary(basename+'.mogx')
  y = vk.load_binary(basename+'.mogy')
  if y.size != x.size:
    raise RuntimeError('x, y size not match')
  z = vk.load_binary(basename+'.mogz')
  if z.size != x.size:
    raise RuntimeError('x, z size not match')
  ar = vk.load_binary(str(index)+'.mogar')
  if ar.size != x.size:
    raise RuntimeError('x, ar size not match')
  ai = vk.load_binary(str(index)+'.mogai')
  if ai.size != x.size:
    raise RuntimeError('x, ai size not match')
  br = vk.load_binary(str(index)+'.mogbr')
  if br.size != x.size:
    raise RuntimeError('x, br size not match')
  bi = vk.load_binary(str(index)+'.mogbi')
  if bi.size != x.size:
    raise RuntimeError('x, bi size not match')
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
    isovla = p3.cmplx_orb_plot_mog('spin', atoms, mod, sp, x, y, z, isovalue)
    isovlb = p3.cmplx_orb_plot_mog('orb', atoms, mod, mp, x, y, z, isovalue)
  else:
    moda = sqrt(square(ar) + square(ai))
    pha = arctan2(ar, ai)
    modb = sqrt(square(br) + square(bi))
    phb = arctan2(br, bi)
    isovla = p3.cmplx_orb_plot_mog('alpha', atoms, moda, pha, x, y, z, isovalue)
    isovlb = p3.cmplx_orb_plot_mog('beta', atoms, modb, phb, x, y, z, isovalue)
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
        p3.cmplx_orb_plot_mog('spin', atoms, mod, sp, x, y, z, new_isov)
      else:
        mlab.close('alpha(amplitude)')
        mlab.close(f'alpha(phase), isovalue={old_isov}')
        p3.cmplx_orb_plot_mog('alpha', atoms, moda, pha, x, y, z, new_isov)
      mlab.draw()
      mlab.view()
    @observe('isovb')
    def update_isovb(self,event):
      old_isov = event.old
      new_isov = event.new
      if spin_phase:
        mlab.close('orb(amplitude)')
        mlab.close(f'orb(phase), isovalue={old_isov}')
        p3.cmplx_orb_plot_mog('orb', atoms, mod, mp, x, y, z, new_isov)
      else:
        mlab.close('beta(amplitude)')
        mlab.close(f'beta(phase), isovalue={old_isov}')
        p3.cmplx_orb_plot_mog('beta', atoms, modb, phb, x, y, z, new_isov)
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
    print('mog2c source.molden.d orb_index (isovalue)')
    print('e.g. mog2c C6H6.molden.d 30')
  else:
    mog2c(*argv[1:])

