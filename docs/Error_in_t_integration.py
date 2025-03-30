'''
investigate the error on t-expansion series of Coulomb integral
env:base
'''

import math
import copy

choice=2 # 1: direct integration; 2: Taylor expansion of integration at X=0
X=2.3353 # The smaller the value, the greater the error
it=30 # The number of iterations is positively related to angular momentum. 
# The more iterations, the greater the error.
d=20 # Taylor expansion series at X = 0
a=[0]*it
NX_mic=[0]*d
NX_mic_pre=[0]*d
NX_mic_pre_pre=[0]*d
if choice==1:
# direct integration, error caused by numerical errors in calculating limits
  a[0]=math.sqrt(math.pi/X) * math.erf(math.sqrt(X)) / 2.0
  print(1, a[0])
  a[1]=(1.0 - math.exp(-X)) / (2.0 * X)
  print(2, a[1])
  i=2
  while i<=(it-1):
    a[i]=(-math.exp(-X) + (i-1) * a[i-2]) / (2.0 * X)
    print(i+1, a[i])
    i+=1
else:
# Taylor expansion of integration at X=0, error caused by offset of X from 0
  i=0
  while i<=(it-1):
    if i==0:
      j=0
      while j<d:
        NX_mic[j]=(-1)**(j+1)*(1/((2*(j+1)+1)*math.factorial(j+1)))
        j+=1
      a[i]= 1.0 - X/3.0 + X*X/10.0 - X**3/42.0 + X**4/216.0 -\
        X**5/1320.0 + X**6/9360.0
    elif i==1:
      NX_mic_pre=copy.copy(NX_mic)
      j=0
      while j<d:
        NX_mic[j]=(-1)**(j+1)*(1/(math.factorial(j+2)))
        j+=1
      a[i]=(1.0 - X/2.0 + X*X/6.0 - X**3/24.0 + X**4/120.0 -\
          X**5/720.0 + X**6/5040.0) / (2.0)
    else:
      NX_mic_pre_pre=copy.copy(NX_mic_pre)
      NX_mic_pre=copy.copy(NX_mic)
      j=0
      while j<(d-1):
        NX_mic[j]=((-1)**(j+1)*(1/(2*math.factorial(j+2)))+\
               (NX_mic_pre_pre[j+1])/2)*(i+1)
        j+=1
      NX_mic[d-1]=((-1)**(j+1)*(1/(2*math.factorial(j+2))))*(i+1)
      j=0
      while j<d:
        a[i]+=(-1)**j*X**j*(1/math.factorial(j+1))/2
        j+=1
      j=0
      while j<d:
        a[i]+=X**j*NX_mic_pre_pre[j]/2
        j+=1
    print(i+1,'  ', a[i],'  ', abs(a[i]-1/(i+1))/(1/(i+1)))
    i+=1
    
"""
To obtain sufficient accuracy within 50 rounds of iterations, 
the truncation X is set to 4 and the number of 
Taylor expansion series is set to 30!

To obtain sufficient accuracy within 40 rounds of iterations, 
the truncation X is set to 3 and the number of 
Taylor expansion series is set to 25!

To obtain sufficient accuracy within 30 rounds of iterations, 
the truncation X is set to 1 and the number of 
Taylor expansion series is set to 20!

To obtain sufficient accuracy within 25 rounds of iterations, 
the truncation X is set to 1 and the number of 
Taylor expansion series is set to 15!

To obtain sufficient accuracy within 15 rounds of iterations, 
the truncation X is set to 0.1 and the number of 
Taylor expansion series is set to 8!

To obtain sufficient accuracy within 10 rounds of iterations, 
the truncation X is set to 0.01 and the number of 
Taylor expansion series is set to 5!

To obtain sufficient accuracy within 5 rounds of iterations, 
the truncation X is set to 0.001 and the number of 
Taylor expansion series is set to 2!
"""
