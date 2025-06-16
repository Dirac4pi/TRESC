'''
python module in vis2c
visualizing 2-component complex orbital
coding:UTF-8
env:vis2c
'''

import subprocess
import dataload as dl
import fileconv  as fc
import plot3D as p3
from os import path, rename
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item
from numpy import sqrt, square, arctan2

# True: alpha and beta visualization as phase, False: as components
spin_phase = True

def numbered_rename(src: str, dst: str) -> str:
  """automatically add serial number suffixes to avoid conflicts"""
  if not path.exists(src):
    raise FileNotFoundError(f"file not exists: {src}")
  
  base, ext = path.splitext(dst)
  counter = 1
  while path.exists(dst):
    dst = f"{base}_{counter}{ext}"
    counter += 1
  
  rename(src, dst)
  return dst

def cub2c(real:str, img:str, index:int=0, isovalue:float=0.05)->None:
  '''
  2-component complex MO visualization, uniform cube grid net.
  --
  !!! 2 Gaussian cube format file contain (alpha real & beta real),\n
  (alpha imagine & beta imagine) of a selected 2-component complex orbital.\n
  real: address of .cube file contains real part of MO.\n
  img: address of .cube file contains imaginary part of MO.\n
  index: index of MO.\n
  isovalue: isovalue of amplitude.\n
  Returns: None
  '''
  if not real.endswith('.molden.input'):
    real = real + '.molden.input'
  if not img.endswith('.molden.input'):
    img = img + '.molden.input'
  if not path.exists(real):
    raise RuntimeError(\
      f"can't find file {real}, may be it's not .molden.input file")
  if not path.exists(img):
    raise RuntimeError(\
      f"can't find file {img}, may be it's not .molden.input file")
  print(f'spin_phase = {spin_phase}')
  with open(real, 'r') as f:
    orbcount = 0
    while True:
      line = f.readline()
      line = line.lower()
      if line.startswith('spin=') and line.endswith('alpha\n'):
        orbcount += 1
      elif line.startswith('spin=') and line.endswith('beta\n'):
        break

  print('generate cube file')
  commandreal = f"""
  {real}
  200
  3
  {index},{index+orbcount}
  3
  2
  0
  q
  """
  commandimg = f"""
  {img}
  200
  3
  {index},{index+orbcount}
  3
  2
  0
  q
  """
  try:
    processreal = subprocess.run(
      'multiwfn',
      input=commandreal,
      text=True,
      capture_output=True,
      cwd="./",
      check=True
    )
    if "Error" in processreal.stderr:
      raise RuntimeError(f"Multiwfn error:\n{processreal.stderr}")
  except subprocess.CalledProcessError as e:
    raise RuntimeError(f"Multiwfn failed: {e.returncode}:\n{e.stderr}") from e
  real = numbered_rename('orbital.cub','real.cub')

  try:
    processimg = subprocess.run(
      'multiwfn',
      input=commandimg,
      text=True,
      capture_output=True,
      cwd="./",
      check=True
    )
    if "Error" in processimg.stderr:
      raise RuntimeError(f"Multiwfn error:\n{processimg.stderr}")
  except subprocess.CalledProcessError as e:
    raise RuntimeError(f"Multiwfn failed: {e.returncode}:\n{e.stderr}") from e
  img = numbered_rename('orbital.cub','img.cub')

  print('loading and plotting, please be patient...')
  atoms, real_mat_alpha, nmo = dl.load_cube(real, 1)
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
  atoms, real_mat_beta, nmo = dl.load_cube(real, 2)
  br = real_mat_beta['value']
  atoms, img_mat_alpha, nmo = dl.load_cube(img, 1)
  ai = img_mat_alpha['value']
  atoms, img_mat_beta, nmo = dl.load_cube(img, 2)
  bi = img_mat_beta['value']
  # plot the complex 2c orbital
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

def mog2c(real:str, img:str, index:int=1, isovalue:float=0.05)->None:
  '''
  2-component complex MO visualization, irregular grid net.
  --
  !!! 2 molden format file contain (alpha real & beta real),\n
  (alpha imagine & beta imagine) of a selected 2-component complex orbital.\n
  real: title of -real.molden.input file.\n
  img: title of -img.molden.input file.\n
  index: index of MO.\n
  isovalue: isovalue of amplitude of MO.\n
  Returns: None
  '''
  isovalue = float(isovalue)
  if not real.endswith('.molden.input'):
    real = real + '.molden.input'
  if not img.endswith('.molden.input'):
    img = img + '.molden.input'
  if not path.exists(real):
    raise RuntimeError(\
      f"can't find file {real}, may be it's not .molden.input file")
  if not path.exists(img):
    raise RuntimeError(\
      f"can't find file {img}, may be it's not .molden.input file")

  print(f'spin_phase = {spin_phase}')
  real = fc.mog_init(real)
  print(f'generate {index}.mo')
  if real.endswith('.molden.input'):
    realMOalpha = dl.load_mo_from_molden(real, 'alpha', index)
    realMObeta = dl.load_mo_from_molden(real, 'beta', index)
    imgMOalpha = dl.load_mo_from_molden(img, 'alpha', index)
    imgMObeta = dl.load_mo_from_molden(img, 'beta', index)
    with open(str(index)+'.mo','w') as f:
      for i in range(len(realMOalpha)):
        f.write('('+str(realMOalpha[i])+','+str(imgMOalpha[i])+')\n')
        f.write('('+str(realMObeta[i])+','+str(imgMObeta[i])+')\n')
  else:
    raise RuntimeError('wavefunction file cannot be loaded in mog2c')
  
  print('call TRESC...')
  fc.call_executable(['tshell.sh', '-2c', str(index)])

  atoms = dl.load_xyz('.xyz')

  x = dl.load_binary('.mogx')
  y = dl.load_binary('.mogy')
  if y.size != x.size:
    raise RuntimeError('x, y size not match')
  z = dl.load_binary('.mogz')
  if z.size != x.size:
    raise RuntimeError('x, z size not match')

  ar = dl.load_binary(str(index)+'.mogar')
  if ar.size != x.size:
    raise RuntimeError('x, ar size not match')
  ai = dl.load_binary(str(index)+'.mogai')
  if ai.size != x.size:
    raise RuntimeError('x, ai size not match')
  br = dl.load_binary(str(index)+'.mogbr')
  if br.size != x.size:
    raise RuntimeError('x, br size not match')
  bi = dl.load_binary(str(index)+'.mogbi')
  if bi.size != x.size:
    raise RuntimeError('x, bi size not match')

  # plot the complex 2c orbital
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
  import sys
  mog2c(*sys.argv[1:])

