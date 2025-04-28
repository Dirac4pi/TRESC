'''
python module in vis2c
visualizing 2-component complex orbital
coding:UTF-8
env:2cvis
'''

from .dataload import load_mo_from_molden
from .fileconv import mog_init, call_fortran
from .plot3D import cmplx_orb_plot_cub, cmplx_orb_plot_mog
from .dataload import load_binary, load_cube, load_xyz
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item
import os, subprocess, platform

def numbered_rename(src: str, dst: str) -> str:
  """automatically add serial number suffixes to avoid conflicts"""
  if not os.path.exists(src):
    raise FileNotFoundError(f"file not exists: {src}")
  
  base, ext = os.path.splitext(dst)
  counter = 1
  while os.path.exists(dst):
    dst = f"{base}_{counter}{ext}"
    counter += 1
  
  os.rename(src, dst)
  return dst

def cub2c(real:str, img:str, index:int=0, isovalue:float=0.1)->None:
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
  if not real.endswith('.molden'):
    real = real + '.molden'
  if not img.endswith('.molden'):
    img = img + '.molden'
  if not os.path.exists(real):
    raise RuntimeError(f"can't find file {real}, may be it's not .molden file")
  if not os.path.exists(img):
    raise RuntimeError(f"can't find file {img}, may be it's not .molden file")
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
  atoms, real_mat_alpha, nmo = load_cube(real, 1)
  if nmo == 0:
    raise RuntimeError('no orbital info in cube file '+real)
  elif nmo != 2:
    raise RuntimeError('cube file '+real+' should contain both \
      alpha and beta orbitals.')
  # grid data
  x = real_mat_alpha['x']
  y = real_mat_alpha['y']
  z = real_mat_alpha['z']
  real_val_alpha = real_mat_alpha['value']
  atoms, real_mat_beta, nmo = load_cube(real, 2)
  real_val_beta = real_mat_beta['value']
  atoms, img_mat_alpha, nmo = load_cube(img, 1)
  img_val_alpha = img_mat_alpha['value']
  atoms, img_mat_beta, nmo = load_cube(img, 2)
  img_val_beta = img_mat_beta['value']
  
  # plot the complex 2c orbital
  isovla = cmplx_orb_plot_cub('alpha', atoms, real_val_alpha, img_val_alpha, \
                              x, y, z, isovalue)
  isovlb = cmplx_orb_plot_cub('beta', atoms, real_val_beta, img_val_beta, \
                              x, y, z, isovalue)
  # UI regulation isovalue
  class IsoValueController(HasTraits):
    isova = Float(isovla, desc="isova", auto_set=False, enter_set=True)
    isovb = Float(isovlb, desc="isovb", auto_set=False, enter_set=True)
    @observe('isova')
    def update_isova(self,event):
      old_isov = event.old
      new_isov = event.new
      mlab.close('alpha(amplitude)')
      mlab.close(f'alpha(phase), isovalue={old_isov}')
      cmplx_orb_plot_cub('alpha', atoms, real_val_alpha, img_val_alpha, \
                        x, y, z, new_isov)
      mlab.draw()
      mlab.view()
    @observe('isovb')
    def update_isovb(self,event):
      old_isov = event.old
      new_isov = event.new
      mlab.close('beta(amplitude)')
      mlab.close(f'beta(phase), isovalue={old_isov}')
      cmplx_orb_plot_cub('beta', atoms, real_val_beta, img_val_beta, \
                        x, y, z, new_isov)
      mlab.draw()
      mlab.view()
    viewalpha = View(
      Item('isova', label="Alpha_isovalue", show_label=True),
      width=300,
      height=200,
      resizable=True
    )
    viewbeta = View(
      Item('isovb', label="Beta_isovalue", show_label=True),
      width=300,
      height=200,
      resizable=True
    )
  controller = IsoValueController()
  controller.configure_traits()

def mog2c(real:str, img:str, index:int=0, isovalue:float=0.1)->None:
  '''
  2-component complex MO visualization, irregular grid net.
  --
  !!! 2 molden format file contain (alpha real & beta real),\n
  (alpha imagine & beta imagine) of a selected 2-component complex orbital.\n
  real: title of -real.molden file.\n
  img: title of -img.molden file.\n
  index: index of MO.\n
  isovalue: isovalue of amplitude of MO.\n
  Returns: None
  '''
  if not real.endswith('.molden'):
    real = real + '.molden'
  if not img.endswith('.molden'):
    img = img + '.molden'
  if not os.path.exists(real):
    raise RuntimeError(f"can't find file {real}, may be it's not .molden file")
  if not os.path.exists(img):
    raise RuntimeError(f"can't find file {img}, may be it's not .molden file")

  real = mog_init(real)
  print(f'generate {index}.mo')
  if real.endswith('.molden'):
    realMOalpha = load_mo_from_molden(real, 'alpha', index)
    realMObeta = load_mo_from_molden(real, 'beta', index)
    imgMOalpha = load_mo_from_molden(img, 'alpha', index)
    imgMObeta = load_mo_from_molden(img, 'beta', index)
    with open(str(index)+'.mo','w') as f:
      for i in range(len(realMOalpha)):
        f.write('('+str(realMOalpha[i])+','+str(imgMOalpha[i])+')\n')
        f.write('('+str(realMObeta[i])+','+str(imgMObeta[i])+')\n')
  else:
    raise RuntimeError('wavefunction file cannot be loaded in mog2c')
  
  print('call TRESC')
  if platform.system() == "Linux":
    results = call_fortran(['tresc', '-2c', str(index)])
  elif platform.system() == "Windows":
    results = call_fortran(['thomas', '-2c', str(index)])
  print(f'stdout: {results}')

  atoms = load_xyz('.xyz')

  x = load_binary('.mogx')
  size = x.size
  y = load_binary('.mogy')
  if y.size != size:
    raise RuntimeError('x, y size not match')
  z = load_binary('.mogz')
  if z.size != size:
    raise RuntimeError('x, z size not match')

  ar = load_binary(str(index)+'.mogar')
  if ar.size != size:
    raise RuntimeError('x, ar size not match')
  ai = load_binary(str(index)+'.mogai')
  if ai.size != size:
    raise RuntimeError('x, ai size not match')
  br = load_binary(str(index)+'.mogbr')
  if br.size != size:
    raise RuntimeError('x, br size not match')
  bi = load_binary(str(index)+'.mogbi')
  if bi.size != size:
    raise RuntimeError('x, bi size not match')

  # plot the complex 2c orbital
  isovla = cmplx_orb_plot_mog('alpha', atoms, ar, ai, x, y, z, isovalue)
  isovlb = cmplx_orb_plot_mog('beta', atoms, br, bi, x, y, z, isovalue)
  # UI regulation isovalue
  class IsoValueController(HasTraits):
    isova = Float(isovla, desc="isova", auto_set=False, enter_set=True)
    isovb = Float(isovlb, desc="isovb", auto_set=False, enter_set=True)
    @observe('isova')
    def update_isova(self,event):
      old_isov = event.old
      new_isov = event.new
      mlab.close('alpha(amplitude)')
      mlab.close(f'alpha(phase), isovalue={old_isov}')
      cmplx_orb_plot_mog('alpha', atoms, ar, ai, x, y, z, new_isov)
      mlab.draw()
      mlab.view()
    @observe('isovb')
    def update_isovb(self,event):
      old_isov = event.old
      new_isov = event.new
      mlab.close('beta(amplitude)')
      mlab.close(f'beta(phase), isovalue={old_isov}')
      cmplx_orb_plot_mog('beta', atoms, br, bi, x, y, z, new_isov)
      mlab.draw()
      mlab.view()
    viewalpha = View(
      Item('isova', label="Alpha_isovalue", show_label=True),
      width=300,
      height=200,
      resizable=True
    )
    viewbeta = View(
      Item('isovb', label="Beta_isovalue", show_label=True),
      width=300,
      height=200,
      resizable=True
    )
  controller = IsoValueController()
  controller.configure_traits()

