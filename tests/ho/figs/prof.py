import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py as h5py
import sys as sys

if(len(sys.argv) > 1):
  nx = int(sys.argv[1])
  nr = int(sys.argv[2])
else:
  nx = 200
  nr  = 200

cin = 1.9941e-07 #M
kp  = 2.5e5
length = 20.0
radius = 0.028

z, r = np.meshgrid(np.linspace(0, length, nx), np.linspace(0, radius*10.0, nr + 1))

f0=h5py.File('bv01_%d_%d_restart.h5'%(nx, nr), 'r')
c0 = f0['restart'][:]
c0 = c0 / cin
#c0[nx:nx*(nr+1), :] = c0[nx:nx*(nr+1), :] * kp
cmax = max(1.0, c0.max()) 
cmin = 0.0

ax0 = plt.subplot(1, 3, 1)
plt.pcolormesh(r, z, np.asarray(c0.reshape(nr+1, nx)), vmin=cmin, vmax=cmax)
plt.xlim([0, 0.28])
plt.ylim([20, 0])
plt.ylabel('x (cm)')
plt.xlabel('r (mm)')
plt.setp(ax0.get_xticklabels(), visible=False)
plt.xticks([0, .1, .2])
plt.title('(a) 0.1 bv/min')
f0.close()

f1=h5py.File('bv03_%d_%d_restart.h5'%(nx, nr), 'r')
c1 = f1['restart'][:]
c1 = c1 / cin
#c1[nx:nx*(nr+1), :] = c1[nx:nx*(nr+1), :] * kp
ax1 = plt.subplot(1, 3, 2)
plt.pcolormesh(r, z, np.asarray(c1.reshape(nr+1, nx)), vmin=cmin, vmax=cmax)
plt.xlim([0, 0.28])
plt.ylim([20, 0])
plt.setp(ax1.get_yticklabels(), visible=False)
plt.xticks([0, .1, .2])
plt.title('(b) 0.3 bv/min')
f1.close()
plt.xlabel('r (mm)')

f2=h5py.File('bv10_%d_%d_restart.h5'%(nx, nr), 'r')
c2 = f2['restart'][:]
c2 = c2 / cin
#c2[nx:nx*(nr+1), :] = c2[nx:nx*(nr+1), :] * kp
ax2 = plt.subplot(1, 3, 3)
ct = plt.pcolormesh(r, z, np.asarray(c2.reshape(nr+1, nx)), vmin=cmin, vmax=cmax)
plt.xlim([0, 0.28])
plt.ylim([20, 0])
plt.setp(ax2.get_yticklabels(), visible=False)
plt.xticks([0, .1, .2])
plt.title('(c) 1 bv/min')
f2.close()
plt.xlabel('r (mm)')

fig = plt.gcf()
fig.subplots_adjust(left=0.1, right=0.88, wspace=0.05, hspace=0.15)
cbar_ax = fig.add_axes([0.90, 0.1, 0.03, 0.8])
fig.colorbar(ct, cax=cbar_ax) #, ticks=[0, .2, .4, .6, .8, 1.0])

fig.set_size_inches(9, 3)
plt.savefig('prof.png')
#plt.savefig('prof.pdf')
plt.show()

