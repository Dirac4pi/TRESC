'''
plot 2D objects
coding:UTF-8
env:vis2c
'''

import numpy as np
import matplotlib.pyplot as plt


def fill(digit:int, num:int=2, element:str='0') -> str:
  '''
  fill in the blank of a digit
  --
  digit: digit.\n
  num: number of digit.\n
  element: fill element.\n
  Returns: filled digit
  '''
  if len(element) != 1:
    raise RuntimeError('fill too long or too short')
  if len(str(digit)) > num:
    raise RuntimeError('digit too long')
  filled = str(digit)
  i = 1
  while i <= num-len(str(digit)):
    filled = element + filled
    i += 1
  return filled

def batch_plot(start:str, flag:str) -> None:
  '''
  Plot data from batch calculations
  --
  start: start of the objective\n
  flag: flag of the objective\n
  '''
  initi = 1
  initj = 1
  fini = 20
  finj = 20
  i = initi
  while i <= fini:
    dat = []
    j = initj
    while j <= finj:
      try:
        filename = fill(i) + fill(j) + '.esc'
        file = open(filename,'r')
      except FileNotFoundError:
        raise RuntimeError('file ' + filename + ' not found')
      while True:
        record = file.readline()
        if record.startswith(start):
          pos = record.find(flag)
          tgt = record[pos+len(flag):-1].strip()
          dat.append(float(tgt))
          break
        elif record == '':
          raise RuntimeError(filename+" can't find "+start)
      file.close()
      j += 1
    count = np.linspace(1, (finj-initj+1), (finj-initj+1))
    plt.plot(count, dat, linewidth=2, label='serial='+str(i))
    i += 1
  plt.legend()
  plt.xlabel('xlabel')
  plt.ylabel(start.strip())
  plt.grid()
  plt.show()


def scf_plot(filename:str) -> None:
  '''
  visualisation of iteration for ongoing or completed SCF computation
  --
  filename: file name contain .esc
  '''
  eng = []
  rmsdp = []
  damp = []
  loop = 0
  flag = 0
  converged = False
  try:
    file = open(filename,'r')
  except FileNotFoundError:
    raise RuntimeError('file ' + filename + ' does nor exist')
  for line in file:
    if line.startswith('  SCF energy (A.U.)'):
      i = line.find('-')
      j = line.find(';')
      eng.append(float(line[i+1:j-1]))
      loop += 1
    elif line.startswith('  -- RMSDP'):
      rmsdp.append(float(line.split()[2]))
    elif line.startswith('  DIIS information'):
      flag = 1
    elif flag == 1:
      if line.startswith('  -- no DIIS acceleration'):
        flag = 2
      else:
        flag = 0
        damp.append(0.)
    elif flag == 2:
      if line.startswith('  -- undamped'):
        damp.append(0.)
      elif line.startswith('  -- fallback'):
        damp.append(damp[-1] + (1.-damp[-1])/2.)
      else:
        damp.append(float(line.split()[2]))
      flag = 0
    elif line.find('SCF succeed!') != -1:
      converged = True
      
  file.close()
  count1 = np.linspace(1, loop, loop)
  if converged:
    count2 = np.linspace(2, loop-1, loop-2)
  else:
    count2 = np.linspace(2, loop, loop-1)
  if converged:
    countd = np.linspace(1, loop-1, loop-1)
  fig, axs = plt.subplots(2, 1, figsize=(8, 6))
  ax1 = axs[0]
  ax1.plot(count1, eng, linewidth=2, color='darkblue', \
  marker='x', markersize=10)
  ax1.set_xlabel('iteration', fontsize=15)
  ax1.set_ylabel('SCF energy', color='darkblue', fontsize=15)
  ax1.tick_params(axis='y', labelcolor='darkblue')
  ax2 = ax1.twinx()
  plt.grid()
  ax2.plot(countd, damp, linewidth=2, color='darkorange',\
   marker='x', markersize=10)
  ax2.set_ylabel('damp', color='darkorange', fontsize=15)
  ax2.tick_params(axis='y', labelcolor='darkorange')
  if loop >= 2:
    ax3 = axs[1]
    ax3.plot(count2, rmsdp, linewidth=2, color='darkblue', \
    marker='x', markersize=10)
    ax3.set_xlabel('iteration', fontsize=15)
    ax3.set_ylabel('RMSDP', color='darkblue', fontsize=15)
    ax3.tick_params(axis='y', labelcolor='darkblue')
    ax4 = ax3.twinx()
    ax4.plot(countd, damp, linewidth=2, color='darkorange',\
    marker='x', markersize=10)
    ax4.set_ylabel('damp', color='darkorange', fontsize=15)
    ax4.tick_params(axis='y', labelcolor='darkorange')
    plt.grid()
  fig.tight_layout()
  plt.show()

