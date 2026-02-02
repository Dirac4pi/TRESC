'''
vis2c toolkit
coding:UTF-8
env:vis2c
'''

import numpy as np
import subprocess

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
      cube_address = cube_address.rstrip('e')
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
