'''
python module in 2cvis
Transform all Gaussian basis set file(s) in specified dir to Gaussian
system libarary basis set file(s).
coding:UTF-8
env:vis2c
'''

import os

def gbsmod(dir:str='./') -> None:
  '''
  modify all .gbs file in working dir
  --
  dir: specified dir to be modified
  '''
  for file in os.listdir(dir):
    if file.endswith('.gbs'):
      print(dir + file)
      content = ''
      with open(dir+file, 'r', encoding="utf-8") as f:
        for line in f:
          if line.endswith(' 0\n') and not line.startswith('-'):
            line = '-' + line
          content += line
      with open(dir+file, 'w', encoding="utf-8") as f:
        f.write(content)
      new_file = file.replace(file, file.split('.',1)[0]+'.gbs')
      os.rename(dir+file, dir+new_file)

if __name__ == '__main__':
  gbsmod()