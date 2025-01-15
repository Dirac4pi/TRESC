'''
script for SCF monitoring
'''

import numpy as np
import matplotlib.pyplot as plt


def scf_plot(filename):
    '''
    Parameters
    ----------
    filename : string

    Returns
    -------
    None.

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
        exit('file ' + filename + ' does nor exist')
    for line in file:
        if line.startswith('   SCF energy (A.U.)'):
            i = line.find('-')
            j = line.find(';')
            eng.append(float(line[i+1:j-1]))
            loop += 1
        elif line.startswith('   --- RMSDP'):
            rmsdp.append(float(line.split()[2]))
        elif line.startswith('   DIIS information'):
            flag = 1
        elif flag == 1:
            if line.startswith('   --- no DIIS acceleration'):
                flag = 2
            else:
                flag = 0
                damp.append(0.)
        elif flag == 2:
            if line.startswith('   --- undamped'):
                damp.append(0.)
            elif line.startswith('   --- fallback'):
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
    ax1.plot(count1, eng, linewidth=2, color='darkblue', marker='x', markersize=10)
    ax1.set_xlabel('iteration', fontsize=15)
    ax1.set_ylabel('SCF energy', color='darkblue', fontsize=15)
    ax1.tick_params(axis='y', labelcolor='darkblue')
    ax2 = ax1.twinx()
    plt.grid()
    ax2.plot(countd, damp, linewidth=2, color='darkorange', marker='x', markersize=10)
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

def run():
    scf_plot('FeO.esc')
        

if __name__ == '__main__':
    run()
