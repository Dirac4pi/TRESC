'''
script to batch convert .gjf to .xyz
env:base
'''

import os

def isfloat(string):
    '''
    Parameters
    ----------
    string : string

    Returns
    -------
    bool

    '''
    try:
        float(string)
    except:
        return False
    else:
        return True

def iscood(line):
    '''
    Parameters
    ----------
    line : string

    Returns
    -------
    bool

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

def gjf2xyz(keywords=''):
    '''
    Parameters
    ----------
    keywords : string

    Returns
    -------
    None.

    '''
    wd = os.getcwd()
    gjfs = []
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
            wxyz.writelines(keywords)
            wxyz.close()
            

def run():
    keywords = []
    keywords.append('%Atoms\n')
    keywords.append('spin=3\n')
    keywords.append('basis=dkh-def2-svp\n')
    keywords.append('endAtoms\n')
    
    keywords.append('%Hamiltonian\n')
    keywords.append('dkh2\n')
    keywords.append('threads=8\n')
    keywords.append('finitenuc\n')
    keywords.append('endHamiltonian\n')
    
    keywords.append('%SCF\n')
    keywords.append('keepspin\n')
    keywords.append('emd4\n')
    keywords.append('endSCF\n')
    
    gjf2xyz(keywords)

if __name__ == '__main__':
    run()
