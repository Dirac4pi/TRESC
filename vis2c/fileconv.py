'''
python module in vis2c
sbatch convert .gjf to .xyz
coding:UTF-8
env:base
'''

import os, subprocess


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
          cood.append(line+'\n')
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
  title: file name of .chk, .molden or .gbw file\n
  number: orbital number\n
  return: modified title
  """
  # ----------<formatting>----------
  if title.endswith('.molden'):
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
    ./{title+'.molden'}
    0
    q
    """
    title = title + '.molden'
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
    ./{title[0:-5] + '.molden'}
    0
    q
    """
    title = title[0:-5] + '.molden'
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
      output_molden = f"{base_name}.molden"
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
      ["obabel", title, "-O", '.xyz'],
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

def call_fortran(command: list):
  """
  call Fortran command-line programs and handling errors
  --
  command: command list(such as ["./program", "arg1", "arg2"])
  return: standard output
  """
  if list == []:
    raise RuntimeError('None command reserved')
  try:
    result = subprocess.run(
      command,
      check=True,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True
    )
    return result.stdout
  
  except subprocess.CalledProcessError as e:
    # Fortran: stop(1) returns error
    error_message = f"""
    Fortran error stop({e.returncode}):
    ----------------------------
    (stderr):
    {e.stderr}
    ----------------------------
    (stdout):
    {e.stdout}
    """
    raise RuntimeError(error_message) from e

  except FileNotFoundError as e:
    raise RuntimeError(f"can't found {command[0]}") from e

  except Exception as e:
    raise RuntimeError(f"error: {str(e)}") from e
