import numpy as np
import math  as math
import matplotlib.pyplot as plt
import h5py as h5py

def readh5file(directory, nx, nr):
  cc  = []
  tt  = []
  t0  = 0.0

  f0  = h5py.File(directory + '/bv01_%d_%d.h5'%(nx,nr), 'r')
  print directory + '/bv01_%d_%d.h5'%(nx,nr)

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File(directory + '/bv03_%d_%d.h5'%(nx,nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t0 = max(tt)
  f0 = h5py.File(directory + '/bv10_%d_%d.h5'%(nx,nr), 'r')

  for i in f0.keys():
    cc.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  t  = sorted(tt)
  it = sorted(range(len(tt)), key=lambda k: tt[k])

  co = np.zeros_like(cc)
  for i in range(len(tt)):
    co[i] = cc[it[i]]

  res = {}
  res['t'] = np.asarray(t)
  res['c'] = np.asarray(co)

  return res

r2 = readh5file('../200x500', 200, 500)
r3 = readh5file('../500x500', 500, 500)
r4 = readh5file('../500x1000', 500, 1000)
r1 = readh5file('../500x10000', 500, 10000)

cin = 1.9941e-07;
lx = 0.08
ly = 0.88
s2d = 1.0/86400.0

#ax1 = plt.subplot(2, 2, 1)
plt.semilogy(r1['t']*s2d, r1['c']/cin, 'b-', r2['t']*s2d, r2['c']/cin, 'r-', r3['t']*s2d, r3['c']/cin, 'g-', r4['t']*s2d, r4['c']/cin, 'm-')
lgd = plt.legend(('nx = 500, nr = 10000$', 'nx = 200, nr = 500', 'nx = 500, nr = 500', 'nx = 500, nr = 1000'),loc=4,numpoints=1)
lgd.draw_frame(False)
txt = lgd.get_texts()
plt.setp(txt, fontsize='small') 
plt.plot([1, 1], [1e-3, 1.0], 'k:')
plt.plot([2, 2], [1e-3, 1.0], 'k:')
plt.plot([3, 3], [1e-3, 1.0], 'k:')
plt.plot([0, 5], [0.01, 0.01], 'r--')
plt.plot([0, 5], [0.05, 0.05], 'y--')
plt.ylabel('c/c$_0$')
plt.xlim([0, 3])
#plt.ylim([1e-30, 1])
#plt.text(lx, ly, '(a)', transform=ax1.transAxes)
#plt.setp(ax1.get_xticklabels(), visible=False)

fig = plt.gcf()
#fig.subplots_adjust(wspace=0.08, hspace=0.08)
fig.set_size_inches(10, 7.5)
plt.savefig('cbtc.png')
#plt.savefig('hist.pdf')
plt.show()

