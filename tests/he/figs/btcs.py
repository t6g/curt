import numpy as np
import math  as math
import matplotlib.pyplot as plt
import h5py as h5py
import sys as sys

def readh5file(directory, nx, nr):
  cc  = []
  tt  = []
  t0  = 0.0

  f0  = h5py.File(directory + 'bv02_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File(directory + 'bv04_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File(directory + 'bv10_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File(directory + 'bv20_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File(directory + 'bv30_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File(directory + 'bv40_%d_%d.h5'%(nx, nr), 'r')
 
  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File(directory + 'bv50_%d_%d.h5'%(nx, nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  res = {}
  res['t'] = np.asarray(tt)
  res['c'] = np.asarray(cc)

  return res

cin = 1.9941e-07;
c1 = readh5file('../150x100/', 150, 100)
c2 = readh5file('../300x1000/', 300, 1000)
c3 = readh5file('../600x2000/', 600, 2000)


plt.semilogy(c1['t']/60.0, c1['c']/cin, 'b-', c2['t']/60.0, c2['c']/cin, 'r-', c3['t']/60.0, c3['c']/cin, 'g-')
plt.xlabel('Time (min)')
plt.ylabel('c/c$_0$')

fig = plt.gcf()
fig.set_size_inches(6, 4.5)
plt.savefig('btc.png')
plt.show()
