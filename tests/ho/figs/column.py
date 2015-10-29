import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py as h5py
import sys as sys

if(len(sys.argv) > 1):
  nx = int(sys.argv[1])
  nr = int(sys.argv[2])
else:
  nx = 150
  nr  = 100

cin = 1.9941e-07 #M
kp  = 2.5e5
length = 15.0
radius = 0.028

xx = np.linspace(0, length, nx)

f0=h5py.File('bv01_%d_%d_restart.h5' % (nx, nr), 'r')
c0 = f0['restart'][0:nx]
c0 = c0 / cin
print 'c0'
f0.close()

f1=h5py.File('bv03_%d_%d_restart.h5' % (nx, nr), 'r')
c1 = f1['restart'][0:nx]
c1 = c1 / cin
print 'c1'
f1.close()

f2=h5py.File('bv10_%d_%d_restart.h5' % (nx, nr), 'r')
c2 = f2['restart'][0:nx]
c2 = c2 / cin
print 'c2'
f2.close()

plt.plot(xx, c0, 'b-', xx, c1, 'r-', xx, c2, 'g-')
fig = plt.gcf()
fig.set_size_inches(6, 4.5)
plt.savefig('column_%d_%d.png'%(nx, nr))
plt.show()
