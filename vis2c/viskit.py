#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
vis2c toolkit
author:Dirac4pi
env:vis2c
'''

from os import path
import numpy as np
import subprocess
from typing import Tuple, List, Union

ELEMENTS = {
  'H' : 1, 'He': 2, 'Li': 3, 'Be': 4, 'B' : 5, 'C' : 6, 'N' : 7, 'O' : 8,
  'F' : 9, 'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P' :15, 'S' :16,
  'Cl':17, 'Ar':18, 'K' :19, 'Ca':20, 'Sc':21, 'Ti':22, 'V' :23, 'Cr':24,
  'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32,
  'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y' :39, 'Zr':40,
  'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48,
  'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I' :53, 'Xe':54, 'Cs':55, 'Ba':56,
  'La':57, 'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64,
  'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71, 'Hf':72,
  'Ta':73, 'W' :74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80,
  'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 'Fr':87, 'Ra':88,
  'Ac':89, 'Th':90, 'Pa':91, 'U' :92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96,
  'Bk':97, 'Cf':98, 'Es':99
}
ATOMIC_NUMBERS = {v: k for k, v in ELEMENTS.items()}

#-------------------------------------------------------------------------------
def get_element_info(query: Union[int, str]) -> Union[int, str, None]:
  """
  Bi-directional lookup for chemical elements.
  --
  - If query is an int (e.g., 8), returns the symbol (e.g., 'O').
  -If query is a str (e.g., 'O'), returns the atomic number (e.g., 8).
  - Returns None if the element is not found.
  """
  if isinstance(query, int):
    return ATOMIC_NUMBERS.get(query)
  if isinstance(query, str):
    query = query.strip()
    if query.isdigit():
      return ATOMIC_NUMBERS.get(int(query))
    symbol = query.capitalize()
    return ELEMENTS.get(symbol)
  return None

#-------------------------------------------------------------------------------
def isfloat(string:str) -> bool:
  '''
  determine whether the sting is a float or not
  --
  string: input string\n
  return: true/false
  '''
  try:
    float(string)
  except:
    return False
  else:
    return True

#-------------------------------------------------------------------------------
def iscood(line:str) -> bool:
  '''
  determine whether the sting is a coordinate or not
  --
  line: input string\n
  return: true/false
  '''
  s = line.split()
  if len(s) != 4:
    return False
  elif s[0][0] == '!':
    return False
  elif not s[0][0].isalpha():
    return False
  elif isfloat(s[1]) and isfloat(s[2]) and isfloat(s[3]):
    return True
  else:
    return False

#-------------------------------------------------------------------------------
def load_geometry_molden(title:str):
  '''
  load geometry from MOLDEN files
  '''
  try:
    open(title,'r')
  except FileNotFoundError:
    raise RuntimeError(title+' not found')
  with open(title,'r') as f:
    ncenter = 0
    for line in f:
      ncenter += 1
      if '[Atoms]' in line:
        ncenter = 0
      if '[GTO]' in line:
        break
    ncenter -= 1
  with open(title,'r') as f:
    atoms = np.zeros(ncenter, dtype=[('index','i4'), \
    ('charge','f8'), ('x','f8'), ('y','f8'), ('z','f8')])
    start = False
    for line in f:
      if start:
        parts = line.strip().split()
        atoms[i]['index'] = parts[2]
        atoms[i]['x'],atoms[i]['y'],atoms[i]['z'] = \
          np.array(parts[3:6],dtype=float)
        i += 1
        if i == ncenter:
          return atoms
      if '[Atoms]' in line:
        start = True
        i = 0

#-------------------------------------------------------------------------------
def load_geometry_xyz(file_path: str):
  '''
  load geometry from XYZ files
  '''
  parsed = []
  with open(file_path, 'r', encoding='utf-8') as f:
    lines = f.readlines()[2:]
  for line in lines:
    p = line.split()
    if len(p) >= 4:
      z = int(p[0]) if p[0].isdigit() else ELEMENTS.get(p[0].capitalize(), 0)
      parsed.append((z, float(p[1]), float(p[2]), float(p[3])))
  if not parsed:
    return [], [], [], []
  z_out, x_out, y_out, z_out_coord = map(list, zip(*parsed))
  return z_out, x_out, y_out, z_out_coord

#-------------------------------------------------------------------------------
def load_cube(cube_address:str):
  '''
  load data from Gaussian cube format file.
  --
  cube_address: address of .cube file\n
  return: atomic info, cube info
  '''
  try:
    f = open(cube_address, 'r')
  except FileNotFoundError:
    if cube_address.endswith('.cub'):
      cube_address = cube_address + 'e'
    elif cube_address.endswith('.cube'):
      cube_address = cube_address.removesuffix('e')
  else:
    f.close()
  with open(cube_address, 'r') as f:
    lines = f.readlines()
  # read natom and coordination
  iline = 0
  while True:
    try:
      line = lines[iline].split()
      ncenter = int(line[0])
      orgx, orgy, orgz = map(float, line[1:4])
    except ValueError:
      iline += 1
      continue
    else:
      break
  # read grid information
  iline += 1
  line = lines[iline].split()
  n1 = int(line[0])
  v1x, v1y, v1z = map(float, line[1:4])
  iline += 1
  line = lines[iline].split()
  n2 = int(line[0])
  v2x, v2y, v2z = map(float, line[1:4])
  iline += 1
  line = lines[iline].split()
  n3 = int(line[0])
  v3x, v3y, v3z = map(float, line[1:4])
  # allocate memory
  atoms = np.zeros(ncenter, dtype=[('index','i4'), \
                                  ('charge','f8'), \
                                  ('x','f8'), \
                                  ('y','f8'), \
                                  ('z','f8')])
  cubmat = np.zeros((n1,n2,n3), dtype=[('value','f8'), \
                    ('x','f8'), ('y','f8'), ('z','f8')])
  # read atomic information
  for i in range(ncenter):
    iline += 1
    line = lines[iline].split()
    atoms[i]['index'] = int(line[0])
    atoms[i]['charge'], atoms[i]['x'], \
      atoms[i]['y'], atoms[i]['z'] = map(float, line[1:5])
  # read grid data
  for i in range(n1):
    for j in range(n2):
      ilinestart = iline
      iline += 1
      line = lines[iline].split()
      for k in range(n3):
        if k//6 >= iline - ilinestart:
          iline += 1
          line = lines[iline].split()
        cubmat[i, j, k]['value'] = float(line[k%6])
        cubmat[i, j, k]['x'] = orgx + (i)*v1x + (j)*v2x + (k)*v3x
        cubmat[i, j, k]['y'] = orgy + (i)*v1y + (j)*v2y + (k)*v3y
        cubmat[i, j, k]['z'] = orgz + (i)*v1z + (j)*v2z + (k)*v3z
  return atoms, cubmat

#-------------------------------------------------------------------------------
def load_binary(title:str):
  """
  load data from binary file
  --
  default as 8 byte float (double precision in Fortran)
  --
  title: title of binary file\n
  return: Contents of the binary file\n
  """
  try:
    f = open(title, 'rb')
  except FileNotFoundError:
    raise RuntimeError("can't find binary file "+title)
  else:
    f.close()
  dat = np.fromfile(title, dtype=np.float64, offset=8)
  return dat

#-------------------------------------------------------------------------------
def call_executable(command: list):
  """
  call executable command-line programs and handling errors
  --
  command: command list(such as ["./program", "arg1", "arg2"])
  return: standard output
  """
  if not command:
    raise RuntimeError("No command reserved")
  try:
    result = subprocess.run(command, check=True, text=True)
  except subprocess.CalledProcessError as e:
    raise RuntimeError(f"Process failed with return code {e.returncode}") from e
  except FileNotFoundError as e:
    raise RuntimeError(f"Can't find {command[0]}") from e
  except Exception as e:
    raise RuntimeError(f"Error: {str(e)}") from e

#-------------------------------------------------------------------------------
def load_matrix(file_path):
  """
  Parses a block-formatted float matrix from a given file.
  --
  file_path (str): The path to the text file.\n
  Returns: (np.ndarray) An n*n numpy array containing the matrix data.
  """
  matrix_dict = {}
  current_cols = []
  max_dim = 0
  if not path.exists(file_path):
    print(f'{file_path} does not exist')
    exit(1)
  with open(file_path, 'r') as f:
    for line in f:
      tokens = line.strip().split()
      if not tokens:
        continue
      if all(t.isdigit() for t in tokens):
        current_cols = [int(t) for t in tokens]
        continue
      if tokens[0].isdigit() and len(tokens) > 1:
        try:
          row_idx = int(tokens[0])
          values = [float(t) for t in tokens[1:]]
          for i, val in enumerate(values):
            if i < len(current_cols):
              col_idx = current_cols[i]
              matrix_dict[(row_idx, col_idx)] = val
              max_dim = max(max_dim, row_idx, col_idx)
        except ValueError:
          continue
  matrix = np.zeros((max_dim, max_dim), dtype=float)
  # subtracting 1 to convert 1-based file indices to 0-based Python indices
  for (r, c), val in matrix_dict.items():
    matrix[r - 1, c - 1] = val
  return matrix


#===============================================================================
if __name__ == '__main__':
  print('viskit.py is not designed to be run directly.')
  print('Please import the functions in this module to plot your data.')
