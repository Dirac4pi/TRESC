'''
script to plot batch calculation
'''

import numpy as np
import matplotlib.pyplot as plt


def fill(digit, num=2, element='0'):
    '''
    Parameters
    ----------
    digit : integer
        digit.
    num : integer
        number of digit. The default is 2.
    element : string
        fill element. The default is '0'.

    Returns
    -------
    string.

    '''
    if len(element) != 1:
        exit('fill too long or too short')
    if len(str(digit)) > num:
        exit('digit too long')
    filled = str(digit)
    i = 1
    while i <= num-len(str(digit)):
        filled = element + filled
        i += 1
    return filled

def batch_plot(start, flag):
    '''
    Parameters
    ----------
    start : string
        start of the objective
    flag : string
        flag of the objective

    Returns
    -------
    None.

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
                exit('file ' + filename + ' not found')
            while True:
                record = file.readline()
                if record.startswith(start):
                    pos = record.find(flag)
                    tgt = record[pos+len(flag):-1].strip()
                    dat.append(float(tgt))
                    break
                elif record == '':
                    exit(filename+" can't find "+start)
            file.close()
            j += 1
        count = np.linspace(1, (finj-initj+1), (finj-initj+1))
        plt.plot(count, dat, linewidth=2, label='数据系列标签='+str(i))
        i += 1
    plt.legend()
    plt.xlabel('横坐标标签')
    plt.ylabel(start.strip())
    plt.grid()
    plt.show()
    
def run():
    batch_plot('   total electronic energy / Eh', '...')

if __name__ == '__main__':
    run()