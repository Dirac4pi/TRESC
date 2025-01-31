'''
Transform all Gaussian basis set file(s) in specified folder to Gaussian
system libarary basis set file(s).
env:base
'''

import os

path = './'
for file in os.listdir(path):
    if file.endswith('.gbs'):
        print(path + file)
        content = ''
        with open(path+file, 'r', encoding="utf-8") as f:
            for line in f:
                if line.endswith(' 0\n') and not line.startswith('-'):
                    line = '-' + line
                content += line
        with open(path+file, 'w', encoding="utf-8") as f:
            f.write(content)
        new_file = file.replace(file, file.split('.',1)[0]+'.gbs')
        os.rename(path+file, path+new_file)
