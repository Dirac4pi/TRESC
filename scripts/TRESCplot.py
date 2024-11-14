# This script is to plot the molecular energy curve of each SCF loop based 
# on defferent calculation conditions.
import numpy as np
import sys
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
dattyp = 2 # 1 = variance, 2 = inaccuracy
initi = 1
initj = 4
fini = 10
finj = 30
i = initi
while i <= fini:
    dat1 = [] # variance
    dat2 = [] # inaccuracy
    j = initj
    while j <= finj:
        arr = []
        try:
            if len(str(i)) == 1 and len(str(j)) == 1:
                filename = '0' + str(i) + '0' + str(j) + '.tot'
            elif len(str(i)) == 1 and len(str(j)) != 1:
                filename = '0' + str(i) + str(j) + '.tot'
            elif len(str(i)) != 1 and len(str(j)) == 1:
                filename = str(i) + '0' + str(j) + '.tot'
            else:
                filename = str(i) + str(j) + '.tot'
            file = open(filename,'r')
        except FileNotFoundError:
            print('check the limit of loop')
            sys.exit()
        while True:
            record = file.readline()
            if record.startswith('   SCF loop  96'): # larger = more loops
                break
        while len(record) != 0:
            record = file.readline()
            if record.startswith('   molecular energy (A.U.)'):
                splrec = record.split()
                arr.append(float(splrec[4].rstrip(';')))
        dat1.append(np.var(arr))
        err = 0.0
        for k in range(1,5):
            err += abs(arr[k] + 1.133024)
        err = err / 5.0
        dat2.append(err)
        j += 1
    count = np.linspace(1, (finj-initj+1), (finj-initj+1))
    count_smooth = np.linspace(1, (finj-initj+1), 500) # larger = smoother
    dat1_smooth = make_interp_spline(count, dat1)(count_smooth)
    dat2_smooth = make_interp_spline(count, dat2)(count_smooth)
    if dattyp == 1:
        plt.plot(count_smooth, dat1_smooth, linewidth=2, label='nudge='+str(i))
    elif dattyp == 2:
        plt.plot(count_smooth, dat2_smooth, linewidth=2, label='nudge='+str(i))
    i += 1
plt.legend()
plt.xlabel('subsp')
if dattyp == 1:
    plt.ylabel('variance')
elif dattyp == 2:
    plt.ylabel('inaccuracy')
plt.grid() 
plt.show()