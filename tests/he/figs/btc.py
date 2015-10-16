import numpy as np
import math  as math
import matplotlib.pyplot as plt
import h5py as h5py
import sys as sys

def readh5file(nx, nr):
  cc  = []
  tt  = []
  t0  = 0.0

  f0  = h5py.File('bv02_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File('bv04_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File('bv10_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File('bv20_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File('bv30_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File('bv40_%d_%d.h5'%(nx, nr), 'r')
 
  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File('bv50_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  res = {}
  res['t'] = np.asarray(tt)
  res['c'] = np.asarray(cc)

  return res

if(len(sys.argv) > 1):
  nx = int(sys.argv[1])
  nr = int(sys.argv[2])
else:
  nx = 150
  nr  = 100

cin = 1.9941e-07;

cc = readh5file(nx, nr)
plt.plot(cc['t']/60.0, cc['c']/cin, 'b-')
plt.xlabel('Time (min)')
plt.ylabel('c/c$_0$')

fig = plt.gcf()
fig.set_size_inches(6, 4.5)
plt.savefig('btc_%d_%d.png'%(nx, nr))
plt.show()
