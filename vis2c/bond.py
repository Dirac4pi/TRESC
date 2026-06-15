#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
Visualizing Bond Populations in Molecules.
author:007
env:vis2c
'''

import numpy as np
import viskit as vk
import plot3d as p3
import sys
from os import path, environ
from json import load as jsload
import qcelemental as qcel
qcel.CovalentRadii("ALVAREZ2008")
from mayavi import mlab

#-------------------------------------------------------------------------------
def bond_plot(atom_list, atomx, atomy, atomz, bondmat) -> None:
  """
  Visualizing Bond Populations in Molecules.
  --
  atom_list: atomic number of each atom.\n
  atomx: x coordinate of each atom.\n
  atomy: y coordinate of each atom.\n
  atomz: z coordinate of each atom.\n
  bondmat: File path for storing the bond population matrix\n
  Returns: None
  """
  # load standard data
  json_address = environ.get('TRESC') + '/vis2c/PubChemElements_all.json'
  with open(json_address,'r',encoding='utf-8') as f:
    std_dat = jsload(f)
  # get information of each atom
  atom_color=[]
  atom_symbol=[]
  atom_radius=[]
  atom_coradius=[]
  num_atom = 0
  for i in atom_list:
    num_atom += 1
    if i <= 0 or i > 118:
      raise RuntimeError('atomic number cannot be recognised.')
    for j in std_dat["Table"]["Row"]:
      if i == int(j['Cell'][0]):
        atom_symbol.append(j["Cell"][1])
        atom_color.append(j["Cell"][4])
        # Angstrom in .xyz and Bohr in .cub, they're different
        atom_coradius.append(qcel.covalentradii.get(i) * \
                             qcel.constants.bohr2angstroms) # covalent (CSD)
        atom_radius.append(float(j["Cell"][7])/225.) # Van der Waal
        break
  # prepare a set of drawable bonds.
  bondset = []
  for i in range(num_atom):
    for j in range(i): # atoms will not bond to themselves
      pointi = np.array([atomx[i], atomy[i], atomz[i]])
      pointj = np.array([atomx[j], atomy[j], atomz[j]])
      vec_ij = pointj - pointi
      distance = np.linalg.norm(vec_ij)
      if distance <= 1.15 * (atom_coradius[i]+atom_coradius[j]):
        bondset.append(bondmat[i][j])
  # plot atoms
  fig = mlab.gcf()
  for i in range(num_atom):
    RGB = tuple(np.array(p3.hex2rgb(atom_color[i])) / 255.0)
    mlab.points3d(atomx[i], atomy[i], atomz[i], figure=fig, \
                  scale_factor=atom_radius[i], resolution=20, \
                  color=RGB, opacity=1., scale_mode='none')
    for j in range(i): # atoms will not bond to themselves
      pointi = np.array([atomx[i], atomy[i], atomz[i]])
      pointj = np.array([atomx[j], atomy[j], atomz[j]])
      vec_ij = pointj - pointi
      distance = np.linalg.norm(vec_ij)
      if distance <= 1.15 * (atom_coradius[i]+atom_coradius[j]):
        current_bo = bondmat[i][j] 
        x_pts = [atomx[i], atomx[j]]
        y_pts = [atomy[i], atomy[j]]
        z_pts = [atomz[i], atomz[j]]
        s_pts = [current_bo, current_bo] 
        bond_plot = mlab.plot3d(x_pts, y_pts, z_pts, s_pts,
                                tube_radius=0.1,
                                colormap='Greys',
                                vmin=min(bondset),
                                vmax=max(bondset),
                                figure=fig)
  cb = mlab.colorbar(object=bond_plot, title="Bond Order\n",\
                orientation="vertical", nb_labels=5)
  cb.label_text_property.color = (0, 0, 0)
  cb.title_text_property.color = (0, 0, 0)
  cb.scalar_bar.label_format = '%.3f'       # format of colorbar labels
  cb.scalar_bar.unconstrained_font_size = True 
  cb.label_text_property.font_size = 12     # size of colorbar labels
  cb.title_text_property.font_size = 12     # size of colorbar title
  cb.title_text_property.italic = False
  cb.title_text_property.bold = False
  cb.label_text_property.italic = False
  cb.label_text_property.bold = False
  fig.scene.background = (1.0,1.0,1.0)
  mlab.show()



#===============================================================================
if __name__ == "__main__":
  if len(sys.argv) != 3:
    print(f"usage: bond.py xxx.xyz xxx.[txt/dat/mat/bnd]")
    sys.exit(1)
  if sys.argv[1].endswith('.xyz'):
    xyz_path = sys.argv[1]
  else:
    print(f"unsupported file format: {sys.argv[1]}")
    sys.exit(1)
  if sys.argv[2].endswith('.txt') or sys.argv[2].endswith('.dat'):
    bond_path = sys.argv[2]
  elif sys.argv[2].endswith('.mat') or sys.argv[2].endswith('.bnd'):
    bond_path = sys.argv[2]
  else:
    print(f"unsupported file format: {sys.argv[2]}")
    sys.exit(1)
  if not path.isfile(xyz_path):
    print(f"can not find {xyz_path}")
    sys.exit(1)
  if not path.isfile(bond_path):
    print(f"can not find {bond_path}")
    sys.exit(1)
  print(f'loading geometry from {xyz_path} ...')
  atoms, x, y, z = vk.load_geometry_xyz(xyz_path)
  print(f'loading bond population matrix from {bond_path} ...')
  bondmat = vk.load_matrix(bond_path)
  if len(bondmat) != len(bondmat[0]):
    print(f'Bond matrix is ​​not a square matrix')
    exit(1)
  if len(atoms) != len(bondmat):
    print(f'Molecule in {xyz_path} appears to be different from {bond_path}')
    exit(1)
  print('plotting ...')
  bond_plot(atoms, x, y, z, bondmat)


