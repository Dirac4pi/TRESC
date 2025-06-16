'''
python module in vis2c
sbatch convert .gjf to .xyz
coding:UTF-8
env:base
'''

import os, subprocess
import sys
from cclib.io import ccread
from cclib.parser.utils import convertor

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

def gjf2xyz(doc:str='') -> None:
  '''
  convert gjf to xyz
  --
  doc: add comment at the end of all xyz file generated.
  '''
  wd = os.getcwd()
  gjfs = []
  if doc != '':
    with open(doc, 'r') as f:
      keywords = f.read()
  else:
    keywords = ''

  for file in os.listdir(wd):
    if file.endswith('.gjf'):
      gjfs.append(file)
  for gjf in gjfs:
    with open(gjf, 'r', encoding='utf-8') as rgjf:
      xyz = gjf.rstrip('.gjf') + '.xyz'
      wxyz = open(xyz, 'w', encoding='utf-8')
      cood = []
      natom = 0
      for line in rgjf:
        if line.startswith('!'):
          continue
        elif iscood(line):
          i = line.find('(')
          if i != -1:
            j = line.find(')')
            line = line[0:i]+line[j+1:-1]
          cood.append(line)
          natom += 1
      wxyz.write(str(natom)+'\n')
      wxyz.write('\n')
      wxyz.writelines(cood)
      wxyz.write('\n')
      wxyz.write(keywords)
      wxyz.close()

def mog_init(title:str) -> str:
  """
  generate geometry(.xyz) and basis set(.gbs) file
  --
  title: file name of .chk, .molden.input or .gbw file\n
  number: orbital number\n
  return: modified title
  """
  # ----------<formatting>----------
  if title.endswith('.molden.input'):
    print('This is a molden file')
  elif title.endswith('.chk'):
    print('This is a Gaussian checkpoint file')
    print('formatting...')
    try:
      subprocess.run(
        ["formchk", title, title[0:-4]+'.fchk'],
        capture_output=True,
        text=True,
        check=True
      )
    except subprocess.CalledProcessError as e:
      raise RuntimeError(f"formchk failed: {e.stderr}")
    except Exception as e:
      raise RuntimeError(f"error: {str(e)}")
    title = title[0:-4]
    print('chk file transfer to fchk file')

    command = f"""
    {title+'.fchk'}
    100
    2
    6
    ./{title+'.molden.input'}
    0
    q
    """
    title = title + '.molden.input'
    try:
      process = subprocess.run(
        'multiwfn',
        input=command,
        text=True,
        capture_output=True,
        cwd="./",
        check=True
      )
      if "Error" in process.stderr:
        raise RuntimeError(f"Multiwfn error:\n{process.stderr}")
    except subprocess.CalledProcessError as e:
      raise RuntimeError(f"Multiwfn failed: {e.returncode}:\n{e.stderr}") from e

  elif title.endswith('.fch') or title.endswith('.fchk'):
    print('This is a formatted Gaussian checkpoint file')
    command = f"""
    {title}
    100
    2
    6
    ./{title[0:-5] + '.molden.input'}
    0
    q
    """
    title = title[0:-5] + '.molden.input'
    try:
      process = subprocess.run(
        'multiwfn',
        input=command,
        text=True,
        capture_output=True,
        cwd="./",
        check=True
      )
      if "Error" in process.stderr:
        raise RuntimeError(f"Multiwfn error:\n{process.stderr}")
    except subprocess.CalledProcessError as e:
      raise RuntimeError(f"Multiwfn failed: {e.returncode}:\n{e.stderr}") from e
  elif title.endswith('.gbw'):
    print('This is an ORCA wavefunction file')
    print('formatting...')
    try:
      base_name = os.path.splitext(title)[0]
      output_molden_input = f"{base_name}.molden.input"
      output_molden = f"{base_name}.molden.input"
      subprocess.run(
        ['orca_2mkl', base_name, "-molden"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
      )
    except subprocess.CalledProcessError as e:
      raise RuntimeError(f"orca_2mkl failed: {e.stderr}")
    except Exception as e:
      raise RuntimeError(f"error: {str(e)}")
    os.rename(output_molden_input, output_molden)
    title = output_molden
    print('gbw file transfer to molden file')
  else:
    raise RuntimeError("mog_init error, can't recognize file "+title)
  if not os.path.exists(title):
    raise RuntimeError("can't find file "+title)
  
  # ----------<generate .xyz>----------
  print('generate .xyz')
  try:
    subprocess.run(
      ["obabel", "-imolden", title, "-O", '.xyz'],
      check=True,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True
    )
  except subprocess.CalledProcessError as e:
    raise RuntimeError(f"obabel failed: {e.stderr}")
  except Exception as e:
    raise RuntimeError(f"error: {str(e)}")

  # ----------<generate .gbs>----------
  print('generate .gbs')
  with open(title,'r') as molden:
    moldenlines = molden.readlines()
  in_Atoms = False
  in_GTO = False
  atomsymbol = set()
  atomindex = set()
  atomlist = [None] * 200
  with open('.gbs','w') as gbs:
    gbs.write('! generated by mog_init\n')
    gbs.write('\n')
    gbs.write('\n')
    for line in moldenlines:
      if line.strip() == '[GTO]':
        in_GTO = True
        in_Atoms = False
        GTOindex = 1
      elif line.startswith('[Atoms]'):
        in_Atoms = True
      elif in_Atoms and not line.startswith('['):
        if not line.split()[0] in atomsymbol:
          atomsymbol.add(line.split()[0])
          atomindex.add(line.split()[1])
          atomlist[int(line.split()[1])] = line.split()[0]
      elif in_GTO and line.strip() == '':
        if str(GTOindex) in atomindex:
          gbs.write('****\n')
        GTOindex += 1
      elif in_GTO and not line.startswith('['):
        if str(GTOindex) in atomindex and line.split()[1] != '0':
          gbs.write(line.strip().upper()+'\n')
        if str(GTOindex) in atomindex and line.split()[1] == '0':
          gbs.write('-'+atomlist[GTOindex]+'  0\n')
      elif in_GTO and line.startswith('['):
        break
  print('mog_grid initialized')
  return title

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

def convert_to_orca(input_file, output_file):
  try:
    data = ccread(input_file)
  except Exception as e:
    raise ValueError(f"Failed to parse {input_file}: {str(e)}")
  if not hasattr(data, 'atomcoords'):
    raise RuntimeError("No coordinate information found in the input file")
  orca_template = """\
! wB97M-V def2-TZVP def2/J RIJCOSX strongSCF noautostart miniprint nopop
%maxcore 4000
%pal nprocs 32 end
* xyz {charge} {mult}
{coordinates}
*
"""
  charge = getattr(data, 'charge', 0)
  mult = getattr(data, 'mult', 1)

  coordinates = "\n".join(
    f"{data.atomnos[i]} {x:.6f} {y:.6f} {z:.6f}"
    for i, (x, y, z) in enumerate(data.atomcoords[-1])
  )
  with open(output_file, 'w') as f:
    f.write(orca_template.format(
      charge=charge,
      mult=mult,
      coordinates=coordinates
    ))
  print(f"Successfully generated ORCA input: {output_file}")

def load_orb(filename):
  """Fast parser for Molden files to extract orbital energies and occupations.
  """
  alpha, beta = [], []
  current_spin = None
  energy = None
  occ = None
  index = 0
  with open(filename, 'r') as f:
    for line in f:
      line = line.strip()
      # Detect MO section
      if line == '[MO]':
        index += 1
        current_spin = None
        energy = None
        occ = None
      # Detect spin state
      elif 'Spin=' in line:
        if 'Alpha' in line:
          current_spin = 'alpha'
        elif 'Beta' in line:
          current_spin = 'beta'
      # Parse energy
      elif line.startswith('Ene='):
        energy = float(line.split('=')[1])
      # Parse occupation
      elif line.startswith('Occup='):
        occ = float(line.split('=')[1])
        # Only store if we have both energy and occupation
        if energy is not None and current_spin is not None:
          if current_spin == 'alpha':
            alpha.append((index, energy, occ))
          elif current_spin == 'beta':
            beta.append((index, energy, occ))
  return alpha, beta

def print_orb(orbitals, spin, lumo_cutoff=10):
  """Print orbital information up to LUMO+10 with correct HOMO marking"""
  if not orbitals:
    print(f"No {spin} orbitals found")
    return
  # Find HOMO (last occupied orbital with occupation > 0.1)
  homo_index = 0
  for i, (idx, en, occ) in enumerate(orbitals):
    if occ < 0.1:  # Consider orbitals with < 0.1 occupation as unoccupied
      homo_index = i-1
      break
  else:
    # All orbitals are occupied (unlikely but possible)
    homo_index = len(orbitals) - 1
  print(homo_index)
  # Determine print range (HOMO-10 to LUMO+10)
  start = max(0, homo_index - 10)
  end = min(len(orbitals), homo_index + lumo_cutoff + 1)
  print(f"\n{spin} Orbitals (showing HOMO-10 to LUMO+{lumo_cutoff}):")
  print("-" * 50)
  print(f"{'Index':<6}{'Energy(Ha)':<15}{'Occupation':<10}")
  print("-" * 50)
  ii = start
  for idx, en, occ in orbitals[start:end]:
    print(f"{ii:<6}{en:<15.6f}{occ:<10.3f}")
    ii += 1

if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: python convert_to_orca.py input.[gjf|com|xyz] output.inp")
    sys.exit(1)
  input_file = sys.argv[1]
  output_file = sys.argv[2]
  if input_file.endswith('.gjf'):
    gjf2xyz()
    input_file = input_file.rstrip('.gjf') + '.xyz'
  try:
    convert_to_orca(input_file, output_file)
  except Exception as e:
    print(f"Error: {str(e)}")
    sys.exit(1)
