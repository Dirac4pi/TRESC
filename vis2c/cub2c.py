#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
visualization of complex spinor orbital with structured grid data
author:Dirac4pi
env:vis2c
'''

import viskit as vk
import plot3d as p3
from os import path
from mayavi import mlab
from traits.api import HasTraits, Float, observe, Any
from traitsui.api import View, Item
from numpy import sqrt, square, arctan2

# True: alpha and beta visualization as phase, False: as components
spin_phase = True

#-------------------------------------------------------------------------------
def cub2c(moldendir:str, index:int, isovalue:float=0.05,\
          slice:bool=False)->None:
  '''
  visualization of complex spinor orbital with structured grid data
  --
  moldendir: address of xxx.molden.d, contains molden file of 2c MOs\n
  index: index of MO.\n
  isovalue: isovalue of amplitude.\n
  slice: perform a slice.\n
  Returns: None
  '''
  global spin_phase
  if not moldendir.endswith('.molden.d') and \
     not moldendir.endswith('.molden.dir') and \
     not moldendir.endswith('.molden.d/') and \
     not moldendir.endswith('.molden.dir/'):
    raise RuntimeError(f'moldendir should be xxx.molden.d')
  if not path.exists(moldendir):
    raise RuntimeError(f"can't find dir {moldendir}")
  if int(index) <= 0:
    raise RuntimeError(f'index should be a positive integer')
  if float(isovalue) <= 0:
    raise RuntimeError(f'isovalue should be a positive float')
  print(f'spin_phase = {spin_phase}')
  # generate cube files
  print('calling TRESC:')
  vk.call_executable(['tshell.sh', '-cub2c', moldendir, str(index)])
  # load data from cube files
  alphareal = index+'-rar.cub'
  alphaimg  = index+'-rai.cub'
  betareal = index+'-rbr.cub'
  betaimg  = index+'-rbi.cub'
  print(f'loading real space grid alpha real part data from {alphareal}...')
  atoms, real_mat_alpha = vk.load_cube(alphareal)
  print('done')
  # grid data
  x = real_mat_alpha['x']
  y = real_mat_alpha['y']
  z = real_mat_alpha['z']
  ar = real_mat_alpha['value']
  print(f'loading real space grid beta real part data from {betareal}...')
  atoms, real_mat_beta = vk.load_cube(betareal)
  print('done')
  br = real_mat_beta['value']
  print(f'loading real space grid alpha imaginary part data from {alphaimg}...')
  atoms, img_mat_alpha = vk.load_cube(alphaimg)
  print('done')
  ai = img_mat_alpha['value']
  print(f'loading real space grid beta imaginary part data from {betaimg}...')
  atoms, img_mat_beta = vk.load_cube(betaimg)
  print('done')
  bi = img_mat_beta['value']
  # plot the complex 2c orbital
  print('plotting...')
  sa = sqrt(square(ar) + square(ai))
  sb = sqrt(square(br) + square(bi))
  mod = sqrt(square(ar) + square(ai) + square(br) + square(bi))
  sp = arctan2(sa, sb)
  mr = ar + br
  mi = ai + bi
  mp = arctan2(mr, mi)
  if slice:
    p3.slice_cub(f"Spin-phase Slice", atoms, x, y, z, sp, \
                 'viridis', sp.max(), sp.min())
    p3.slice_cub(f"Orbit Isosurfase Slice", atoms, x, y, z, mod, \
                 'viridis', mod.max(), mod.min())
  if spin_phase:
    isovla, figa = p3.cmplx_orb_plot_cub(\
    'Orbit Isosurfase with Spin-phase Mapping', atoms, mod, sp, x, y, z, isovalue)
    isovlb, figb = p3.cmplx_orb_plot_cub(\
    'Orbit Isosurfase with Orbital-phase Mapping', atoms, mod, mp, x, y, z, isovalue)
  else:
    moda = sqrt(square(ar) + square(ai))
    pha = arctan2(ar, ai)
    modb = sqrt(square(br) + square(bi))
    phb = arctan2(br, bi)
    isovla, figa = p3.cmplx_orb_plot_cub('alpha', atoms, moda, pha, \
                                         x, y, z, isovalue)
    isovlb, figb = p3.cmplx_orb_plot_cub('beta', atoms, modb, phb, \
                                         x, y, z, isovalue)
  # synchronise camera
  initial_sync_coe=(figb.scene.camera.position - figb.scene.camera.focal_point)\
                  / (figa.scene.camera.position - figa.scene.camera.focal_point)
  # UI regulation isovalue
  class IsoValueController(HasTraits):
    isova = Float(isovla, desc="isova", auto_set=False, enter_set=True)
    isovb = Float(isovlb, desc="isovb", auto_set=False, enter_set=True)
    fig_a = Any(figa)
    fig_b = Any(figb)
    sync_coe = Any(initial_sync_coe)
    def __init__(self, **traits):
      super(IsoValueController, self).__init__(**traits)
      self._bind_sync()
    def _bind_sync(self):
      self.fig_a.scene.interactor.add_observer("InteractionEvent", \
                                               self.sync_action)
    def sync_action(self, *args):
      if self.fig_a and self.fig_b:
        p3.sync_camera(self.fig_a, self.fig_b, self.sync_coe)
    @observe('isova')
    def update_isova(self, event):
      new_isov = event.new
      title = 'Orbit Isosurfase with Spin-phase Mapping' if spin_phase else 'alpha'
      cam_pos = self.fig_a.scene.camera.position
      cam_foc = self.fig_a.scene.camera.focal_point
      cam_ang = self.fig_a.scene.camera.view_angle
      cam_up  = self.fig_a.scene.camera.view_up
      mlab.close(f'{title}')
      if spin_phase:
        _, self.fig_a = p3.cmplx_orb_plot_cub(f'{title}', atoms, mod, sp,\
                                              x, y, z, new_isov)
      else:
        _, self.fig_a = p3.cmplx_orb_plot_cub(f'{title}', atoms, moda, pha,\
                                              x, y, z, new_isov)
      self.fig_a.scene.camera.position = cam_pos
      self.fig_a.scene.camera.focal_point = cam_foc
      self.fig_a.scene.camera.view_angle = cam_ang
      self.fig_a.scene.camera.view_up = cam_up
      self._bind_sync()
      mlab.draw()
    @observe('isovb')
    def update_isovb(self, event):
      new_isov = event.new
      title = 'Orbit Isosurfase with Orbital-phase Mapping' if spin_phase else 'beta'
      cam_pos = self.fig_b.scene.camera.position
      cam_foc = self.fig_b.scene.camera.focal_point
      cam_ang = self.fig_b.scene.camera.view_angle
      cam_up  = self.fig_b.scene.camera.view_up
      mlab.close(f'{title}')
      if spin_phase:
        _, self.fig_b = p3.cmplx_orb_plot_cub(f'{title}', atoms, mod, mp,\
                                              x, y, z, new_isov)
      else:
        _, self.fig_b = p3.cmplx_orb_plot_cub(f'{title}', atoms, modb, phb,\
                                              x, y, z, new_isov)
      self.fig_b.scene.camera.position = cam_pos
      self.fig_b.scene.camera.focal_point = cam_foc
      self.fig_b.scene.camera.view_angle = cam_ang
      self.fig_b.scene.camera.view_up = cam_up
      # no need to call self._bind_sync(), the listener is attached to fig_a
      # fig_a remains unchanged, while self.sync_action will automatically
      # access the newly created self.fig_b.
      mlab.draw()
    viewalpha = View(
      Item('isova', label="isovalue(a)", show_label=True),
      width=300, height=200, resizable=True
    )
    viewbeta = View(
      Item('isovb', label="isovalue(b)", show_label=True),
      width=300, height=200, resizable=True
    )
  controller = IsoValueController()
  controller.configure_traits()

#===============================================================================
if __name__ == "__main__":
  from sys import argv
  if len(argv) not in [3, 4]:
    print('Usage: cub2c.py source.molden.d orb_index (-slice)')
    print('e.g. cub2c.py C6H6.molden.d -slice')
  else:
    source_file = argv[1]
    orb_index = argv[2]
    if len(argv) == 4:
      if argv[3].lower() == '-slice':
        cub2c(source_file, orb_index, isovalue=0.05, slice=True)
      else:
        exit(f'unkown arg: {argv[3]}')
    else:
      cub2c(source_file, orb_index, isovalue=0.05, slice=False)

