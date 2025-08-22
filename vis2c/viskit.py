'''
vis2c toolkit
coding:UTF-8
env:vis2c
'''

import numpy as np
import subprocess

elements = {
    1: {"symbol": "H", "name": "Hydrogen", "mass": 1.008},
    2: {"symbol": "He", "name": "Helium", "mass": 4.0026},
    3: {"symbol": "Li", "name": "Lithium", "mass": 6.94},
    4: {"symbol": "Be", "name": "Beryllium", "mass": 9.0122},
    5: {"symbol": "B", "name": "Boron", "mass": 10.81},
    6: {"symbol": "C", "name": "Carbon", "mass": 12.011},
    7: {"symbol": "N", "name": "Nitrogen", "mass": 14.007},
    8: {"symbol": "O", "name": "Oxygen", "mass": 15.999},
    9: {"symbol": "F", "name": "Fluorine", "mass": 18.998},
    10: {"symbol": "Ne", "name": "Neon", "mass": 20.180},
    11: {"symbol": "Na", "name": "Sodium", "mass": 22.990},
    12: {"symbol": "Mg", "name": "Magnesium", "mass": 24.305},
    13: {"symbol": "Al", "name": "Aluminium", "mass": 26.982},
    14: {"symbol": "Si", "name": "Silicon", "mass": 28.085},
    15: {"symbol": "P", "name": "Phosphorus", "mass": 30.974},
    16: {"symbol": "S", "name": "Sulfur", "mass": 32.06},
    17: {"symbol": "Cl", "name": "Chlorine", "mass": 35.45},
    18: {"symbol": "Ar", "name": "Argon", "mass": 39.948},
    19: {"symbol": "K", "name": "Potassium", "mass": 39.098},
    20: {"symbol": "Ca", "name": "Calcium", "mass": 40.078},
    21: {"symbol": "Sc", "name": "Scandium", "mass": 44.956},
    22: {"symbol": "Ti", "name": "Titanium", "mass": 47.867},
    23: {"symbol": "V", "name": "Vanadium", "mass": 50.942},
    24: {"symbol": "Cr", "name": "Chromium", "mass": 51.996},
    25: {"symbol": "Mn", "name": "Manganese", "mass": 54.938},
    26: {"symbol": "Fe", "name": "Iron", "mass": 55.845},
    27: {"symbol": "Co", "name": "Cobalt", "mass": 58.933},
    28: {"symbol": "Ni", "name": "Nickel", "mass": 58.693},
    29: {"symbol": "Cu", "name": "Copper", "mass": 63.546},
    30: {"symbol": "Zn", "name": "Zinc", "mass": 65.38},
    31: {"symbol": "Ga", "name": "Gallium", "mass": 69.723},
    32: {"symbol": "Ge", "name": "Germanium", "mass": 72.630},
    33: {"symbol": "As", "name": "Arsenic", "mass": 74.922},
    34: {"symbol": "Se", "name": "Selenium", "mass": 78.971},
    35: {"symbol": "Br", "name": "Bromine", "mass": 79.904},
    36: {"symbol": "Kr", "name": "Krypton", "mass": 83.798},
    37: {"symbol": "Rb", "name": "Rubidium", "mass": 85.468},
    38: {"symbol": "Sr", "name": "Strontium", "mass": 87.62},
    39: {"symbol": "Y", "name": "Yttrium", "mass": 88.906},
    40: {"symbol": "Zr", "name": "Zirconium", "mass": 91.224},
    41: {"symbol": "Nb", "name": "Niobium", "mass": 92.906},
    42: {"symbol": "Mo", "name": "Molybdenum", "mass": 95.95},
    43: {"symbol": "Tc", "name": "Technetium", "mass": 98.0},
    44: {"symbol": "Ru", "name": "Ruthenium", "mass": 101.07},
    45: {"symbol": "Rh", "name": "Rhodium", "mass": 102.91},
    46: {"symbol": "Pd", "name": "Palladium", "mass": 106.42},
    47: {"symbol": "Ag", "name": "Silver", "mass": 107.87},
    48: {"symbol": "Cd", "name": "Cadmium", "mass": 112.41},
    49: {"symbol": "In", "name": "Indium", "mass": 114.82},
    50: {"symbol": "Sn", "name": "Tin", "mass": 118.71},
}

ang2bohr = 1.8897259886

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
def load_cube(cube_address:str, imo:int=0):
  '''
  load data from Gaussian cube format file.
  --
  cube_address: address of .cube file\n
  imo: serial number of the orbital to be loaded\n
  return: atomic info, cube info, number of orbitals
  '''
  try:
    f = open(cube_address, 'r')
  except FileNotFoundError:
    if cube_address.endswith('.cub'):
      cube_address = cube_address + 'e'
    elif cube_address.endswith('.cube'):
      cube_address = cube_address.rstrip('e')
  else:
    f.close()
  with open(cube_address, 'r') as f:
    # read natom and coordination
    while True:
      try:
        ncenter, orgx, orgy, orgz = map(float, f.readline().split())
      except ValueError:
        continue
      else:
        break
    # read grid information
    n1, v1x, v1y, v1z = map(float, f.readline().split())
    n2, v2x, v2y, v2z = map(float, f.readline().split())
    n3, v3x, v3y, v3z = map(float, f.readline().split())
    # determine if it contains MOs
    nmo = 0
    if ncenter < 0:
      nmo = 1
      ncenter = abs(ncenter)  # reduce to positive
    # allocate memory
    atoms = np.zeros(int(ncenter), dtype=[('index','i4'), \
      ('charge','f8'), ('x','f8'), ('y','f8'), ('z','f8')])
    cubmat = np.zeros((int(n1),int(n2),int(n3)), dtype=[('value','f8'), \
      ('x','f8'), ('y','f8'), ('z','f8')])
    # read atomic information
    for i in range(int(ncenter)):
      atoms[i]['index'], atoms[i]['charge'], atoms[i]['x'], \
        atoms[i]['y'], atoms[i]['z'] = map(float, f.readline().split())
    if nmo == 1:  # read the number of MOs if exist
      nmo = int(f.readline().split()[0])
      if nmo < imo:
        raise RuntimeError('read_cube: nmo < imo in '+cube_address)
      elif imo <= 0:
        raise RuntimeError('read_cube: nmo > 0 but imo <= 0 in '+cube_address)
      # read grid data
      for i in range(int(n1)):
        for j in range(int(n2)):
          iline = 1
          line = f.readline().split()
          for k in range(imo, int(n3*nmo)+1, nmo):
            if k//6 >= iline:
              line = f.readline().split()
              iline += 1
            if imo == nmo:
              cubmat[i,j,k//nmo-1]['value'] = float(line[k%6-1])
              #line[-1] is same as line[5]
            else:
              cubmat[i,j,k//nmo]['value'] = float(line[k%6-1])
    else:   # read real-space function
      if imo != 0:
        Warning('no orbital in cube file but imo != 0, file is '+cube_address)
      # read grid data
      for i in range(int(n1)):
        for j in range(int(n2)):
          iline = 1
          line = f.readline().split()
          for k in range(1, int(n3)+1):
            if k // 6 >= iline:
              line = f.readline().split()
              iline += 1
            cubmat[i,j,k-1]['value'] = float(line[k%6-1])
            #line[-1] is same as line[5]
    # assign coordinate information to each grid point
    for i in range(int(n1)):
      for j in range(int(n2)):
        for k in range(int(n3)):
          cubmat[i, j, k]['x'] = orgx + (i)*v1x + (j)*v2x + (k)*v3x
          cubmat[i, j, k]['y'] = orgy + (i)*v1y + (j)*v2y + (k)*v3y
          cubmat[i, j, k]['z'] = orgz + (i)*v1z + (j)*v2z + (k)*v3z
  return atoms, cubmat, nmo

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
