#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
Visualizing Atomic Populations in Molecules.
author:Dirac4pi
env:vis2c
'''

from typing import Tuple, List
from os import environ, path
import numpy as np
from json import load as jsload
import qcelemental as qcel
qcel.CovalentRadii("ALVAREZ2008")
from mayavi import mlab
import viskit as vk
import sys

#-------------------------------------------------------------------------------
def datainline(line: str) -> Tuple[int, int, List[float], List[float]]:
  """
  Extract atomic population information from a string.
  --
  Returns:
    tuple: (population_count, atomic_number, [x, y, z], [populations])
    Returns (-1, -1, [], []) if the format is invalid.
  """
  parts = line.split()
  if len(parts) < 4:
    return -1, -1, [], []
  first_item = parts[0]
  try:
    int(first_item)
  except ValueError:
    try:
      float(first_item)
      return -1, -1, [], []
    except ValueError:
      index = vk.get_element_info(first_item)
  else:
    index = int(first_item)
  tail_items = parts[1:]
  for item in tail_items:
    try:
      float(item)
    except ValueError:
      return -1, -1, [], []
  coords = [float(x) for x in parts[1:4]]
  populations = [float(x) for x in parts[4:]]
  pop_count = len(populations)
  return pop_count, index, coords, populations

#-------------------------------------------------------------------------------
def pop_plot(index: List[int], coords: List[List[float]], \
             populations: List[float]) -> None:
  '''
  visualize atomic populations in molecules.
  --
  index: atomic numbers of each atom.\n
  coords: coordinates of each atom.\n
  populations: population values for each atom.\n
  Returns: None
  '''
  # load standard data
  json_address = environ.get('TRESC') + '/vis2c/PubChemElements_all.json'
  with open(json_address,'r',encoding='utf-8') as f:
    std_dat = jsload(f)
  # get information of each atom
  atom_symbol=[]
  atom_radius=[]
  atom_coradius=[]
  num_atom = 0
  for i in index:
    num_atom += 1
    if i <= 0 or i > 118:
      raise RuntimeError('atomic number cannot be recognised.')
    for j in std_dat["Table"]["Row"]:
      if i == int(j['Cell'][0]):
        atom_symbol.append(j["Cell"][1])
        # Angstrom in .xyz and Bohr in .cub, they're different
        atom_coradius.append(qcel.covalentradii.get(i) * \
                             qcel.constants.bohr2angstroms) # covalent (CSD)
        atom_radius.append(float(j["Cell"][7])/225.) # Van der Waal
        break
  # plot atomic populations
  fig = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
  max_abs = max(abs(min(populations)), abs(max(populations)))*1.3
  for i in range(num_atom):
    # colormap can choose in 'coolwarm', 'bwr', 'seismic'
    atom_plot = mlab.points3d(coords[i][0], coords[i][1], coords[i][2], \
                              populations[i], \
                  figure=fig, scale_factor=atom_radius[i], resolution=20, \
                  colormap='bwr', opacity=1., scale_mode='none', \
                  vmin=-max_abs, vmax=max_abs)
    # invert the colormap to make positive values red and negative values blue
    lut = atom_plot.module_manager.scalar_lut_manager.lut.table.to_array()
    ilut = lut[::-1]
    atom_plot.module_manager.scalar_lut_manager.lut.table = ilut
    for j in range(i): # atoms will not bond to themselves
      pointi = np.array([coords[i][0], coords[i][1], coords[i][2]])
      pointj = np.array([coords[j][0], coords[j][1], coords[j][2]])
      vec_ij = pointj - pointi
      distance = np.linalg.norm(vec_ij)
      if distance <= 1.15 * (atom_coradius[i]+atom_coradius[j]):
        mlab.plot3d(
          [pointi[0], pointj[0]],
          [pointi[1], pointj[1]],
          [pointi[2], pointj[2]],
          tube_radius=0.05,
          tube_sides=20,
          color=(0.1,0.1,0.1)
        )
  fig.scene.background = (1.0,1.0,1.0)
  # plot the colorbar
  cb = mlab.colorbar(orientation='vertical')
  cb.label_text_property.color = (0, 0, 0)
  cb.title_text_property.color = (0, 0, 0)
  cb.scalar_bar.label_format = '%.3f'       # format of colorbar labels
  cb.scalar_bar.unconstrained_font_size = True 
  cb.label_text_property.font_size = 12     # size of colorbar labels
  cb.title_text_property.font_size = 12     # size of colorbar title
  cb.title_text_property.italic = False
  cb.title_text_property.bold = False
  cb.label_text_property.italic = False
  cb.label_text_property.bold = True
  mlab.show()

#-------------------------------------------------------------------------------
def pdb2xyzq(pdb_path: str) -> str:
  '''
  Convert .pdb file to .xyzq file with atomic population information.
  --
  file_path: path of the .pdb file.\n
  Returns: file name of the created .xyzq file.
  '''
  base_name = path.splitext(pdb_path)[0]
  xyzq_path = f"{base_name}.xyzq"
  with open(pdb_path, 'r') as pdb_in, open(xyzq_path, 'w') as xyzq_out:
    for line_num, line in enumerate(pdb_in, start=1):
      if line.startswith(("ATOM", "HETATM")):
        clean_line = line.rstrip('\n\r')
        if len(clean_line) < 78:
          error_msg = (f"Format Error at line {line_num}: "
                       f"ATOM/HETATM record is shorter than 78 columns. "
                       f"Actual length: {len(clean_line)}")
          raise ValueError(error_msg)
        string_a = clean_line[30:72]
        string_b = clean_line[76:78].strip()
        string_c = f"{string_b}    {string_a}\n"
        xyzq_out.write(string_c)
  return xyzq_path

#===============================================================================
if __name__ == "__main__":
  if len(sys.argv) != 2:
    print(f"usage: atomic.py xxx.[xyzq/chg/txt/dat/pdb]")
    sys.exit(1)
  if sys.argv[1].endswith('.xyzq') or sys.argv[1].endswith('.txt'):
    file_path = sys.argv[1]
  elif sys.argv[1].endswith('.chg') or sys.argv[1].endswith('.dat'):
    file_path = sys.argv[1]
  elif sys.argv[1].endswith('.xyz'):
    file_path = sys.argv[1]
  elif sys.argv[1].endswith('.pdb'):
    print(f"this seems to be a .pdb file, will convert it to .xyzq")
    file_path = pdb2xyzq(sys.argv[1])
  else:
    print(f"unsupported file format: {sys.argv[1]}")
    sys.exit(1)
  if not path.isfile(file_path):
    print(f"can not find {file_path}")
    sys.exit(1)
  print(\
'We are plotting atomic populations, which column would you like to process?')
  column_index = input("Enter the column number: ")
  try:
    column_index = int(column_index)
  except ValueError:
    print("Invalid column number. Please enter a valid integer.")
    sys.exit(1)
  if column_index < 1:
    print("Column number must be positive")
    sys.exit(1)
  start = False
  pop_init = 0
  atom_index = []
  atom_coords = []
  atom_pop = []
  with open(file_path, 'r', encoding='utf-8') as file:
    for line_number, line in enumerate(file, start=1):
      pop_count, index, coords, populations = datainline(line)
      if not start:
        if pop_count == -1:
          continue
        else:
          start = True
          if column_index > pop_count:
            print(f"Column index {column_index} exceeds the "+\
                  f"number of population columns ({pop_count})")
            sys.exit(1)
          pop_init = pop_count
      else:
        if pop_count == -1:
          break
        if pop_count != pop_init:
          print(f"Warning: a data error may exist at line {line_number}.")
          break
      atom_index.append(index)
      atom_coords.append(coords)
      atom_pop.append(populations[column_index-1])
  if len(atom_index) == 0:
    print("Error: no valid data found.")
    sys.exit(1)
  pop_plot(atom_index, atom_coords, atom_pop)

