'''
python module in vis2c
visualizing scalar real orbital
coding:UTF-8
env:vis2c
'''

import dataload as dl
import fileconv  as fc
import plot3D as p3
from os import path
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item

def mog1c(mofile:str, spin:str, index:int=1, isovalue:float=0.05)->None:
  '''
  scalar real MO visualization, irregular grid net.
  --
  mofile: title of .molden.input file.\n
  index: index of MO.\n
  isovalue: isovalue of amplitude of MO.\n
  Returns: None
  '''
  if not mofile.endswith('.molden.input'):
    mofile = mofile + '.molden.input'
  if not path.exists(mofile):
    raise RuntimeError(\
      f"can't find {mofile}, make sure it's .molden.input file")
  mofile = fc.mog_init(mofile)
  # print orb info
  alpha, beta = fc.load_orb(mofile)
  fc.print_orb(alpha, "Alpha")
  fc.print_orb(beta, "Beta")
  # generate .mo file
  print(f'generate {index}.mo')
  if mofile.endswith('.molden.input'):
    rMO = dl.load_mo_from_molden(mofile, spin, index)
    with open(str(index)+'.mo','w') as f:
      for i in range(len(rMO)):
        f.write(str(rMO[i])+'\n')
  else:
    raise RuntimeError('wavefunction file cannot be loaded in mog2c')
  
  print('call TRESC...')
  results = fc.call_executable(['TRESC.sh', '-1c', str(index)])
  print(f'stdout: {results}')

  atoms = dl.load_xyz('.xyz')

  x = dl.load_binary('.mogx')
  y = dl.load_binary('.mogy')
  if y.size != x.size:
    raise RuntimeError('x, y size not match')
  z = dl.load_binary('.mogz')
  if z.size != x.size:
    raise RuntimeError('x, z size not match')
  val = dl.load_binary(str(index)+'.mogv')
  if val.size != x.size:
    raise RuntimeError('x, val size not match')

  isovl = p3.real_orb_plot_mog('scalar', atoms, val, x, y, z, isovalue)
  # UI regulation isovalue
  class IsoValueController(HasTraits):
    isov = Float(isovl, desc="isov", auto_set=False, enter_set=True)
    @observe('isov')
    def update_isov(self,event):
      new_isov = event.new
      mlab.close('scalar')
      p3.real_orb_plot_mog('scalar', atoms, val, x, y, z, new_isov)
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
  import sys
  mog1c(*sys.argv[1:])