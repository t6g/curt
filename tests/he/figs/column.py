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

f0=h5py.File('bv02_%d_%d_restart.h5' % (nx, nr), 'r')
c0 = f0['restart'][0:nx]
c0 = c0 / cin
print 'c0'
f0.close()

f1=h5py.File('bv04_%d_%d_restart.h5' % (nx, nr), 'r')
c1 = f1['restart'][0:nx]
c1 = c1 / cin
print 'c1'
f1.close()

f2=h5py.File('bv10_%d_%d_restart.h5' % (nx, nr), 'r')
c2 = f2['restart'][0:nx]
c2 = c2 / cin
print 'c2'
f2.close()

f3=h5py.File('bv20_%d_%d_restart.h5' % (nx, nr), 'r')
c3 = f3['restart'][0:nx]
c3 = c3 / cin
print 'c3'
f3.close()

f4=h5py.File('bv30_%d_%d_restart.h5' % (nx, nr), 'r')
c4 = f4['restart'][0:nx]
c4 = c4 / cin
print 'c4'
f4.close()

f5=h5py.File('bv40_%d_%d_restart.h5' % (nx, nr), 'r')
c5 = f5['restart'][0:nx]
c5 = c5 / cin
print 'c5'
f5.close()

f6=h5py.File('bv50_%d_%d_restart.h5' % (nx, nr), 'r')
c6 = f6['restart'][0:nx]
c6 = c6 / cin
print 'c6'
f6.close()

plt.plot(xx, c0, 'b-', xx, c1, 'r-', xx, c2, 'g-', xx, c3, 'm-', xx, c4, 'c-', xx, c5, 'y-', xx, c6, 'k-')
fig = plt.gcf()
fig.set_size_inches(6, 4.5)
plt.savefig('column_%d_%d.png'%(nx, nr))
plt.show()
