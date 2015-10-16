import numpy as np
import math  as math
import matplotlib.pyplot as plt
import h5py as h5py

def readh5file(nx, nr):
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

if(len(sys.argv) > 1):
  nx = int(sys.argv[1])
  nr = int(sys.argv[2])
else:
  nx = 150
  nr  = 100

cin = 1.9941e-07;
t = np.arange(0.3, 21.3, 0.3)


lx = 0.08
ly = 0.88

ax1 = plt.subplot(2, 2, 1)
plt.plot(t, de0['c']/cin, 'b-', t, base['c']/cin, 'r-', t, de1['c']/cin, 'g-')
lgd = plt.legend(('D$_e$ = 1 $\\times$ 10$^{-6}$', 'D$_e$ = 2 $\\times$ 10$^{-6}$', 'D$_e$ = 5 $\\times$ 10$^{-6}$'),loc=6,numpoints=1)
lgd.draw_frame(False)
txt = lgd.get_texts()
plt.setp(txt, fontsize='small') 
plt.ylabel('c/c$_0$')
plt.xlim([0, 20])
plt.ylim([0, 0.15])
plt.yticks([0, 0.05, 0.1, 0.15])
plt.text(lx, ly, '(a)', transform=ax1.transAxes)
plt.setp(ax1.get_xticklabels(), visible=False)


ax2 = plt.subplot(2, 2, 2)
plt.plot(t, rd0['c']/cin, 'b-', t, base['c']/cin, 'r-', t, rd1['c']/cin, 'g-')
lgd = plt.legend(('r = 0.20 mm', 'r = 0.28 mm', 'r = 0.36 mm'),loc=6,numpoints=1)
lgd.draw_frame(False)
txt = lgd.get_texts()
plt.setp(txt, fontsize='small') 
plt.xlim([0, 20])
plt.ylim([0, 0.15])
plt.yticks([0, 0.05, 0.1, 0.15])
plt.text(lx, ly, '(b)', transform=ax2.transAxes)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

ax3 = plt.subplot(2, 2, 3)
plt.plot(t, kp0['c']/cin, 'b-', t, base['c']/cin, 'r-', t, kp1['c']/cin, 'g-')
lgd = plt.legend(('k$_p$ = 2.5 $\\times$ 10$^{5}$', 'k$_p$ = 5.0 $\\times$ 10$^{5}$', 'k$_p$ = 5.0 $\\times$ 10$^{6}$'),loc=6,numpoints=1)
lgd.draw_frame(False)
txt = lgd.get_texts()
plt.setp(txt, fontsize='small') 
plt.xlim([0, 20])
plt.ylim([0, 0.15])
plt.yticks([0, 0.05, 0.1, 0.15])
plt.xlabel('Bed volume')
plt.ylabel('c/c$_0$')
plt.text(lx, ly, '(c)', transform=ax3.transAxes)

ax4 = plt.subplot(2, 2, 4)
plt.plot(t, dp0['c']/cin, 'b-', t, base['c']/cin, 'r-', t, dp1['c']/cin, 'g-')
lgd = plt.legend(('$\lambda$ = 0.25 mm', '$\lambda$ = 1.00 mm', '$\lambda$ = 5.00 mm'),loc=6,numpoints=1)
lgd.draw_frame(False)
txt = lgd.get_texts()
plt.setp(txt, fontsize='small') 
plt.xlim([0, 20])
plt.ylim([0, 0.15])
plt.yticks([0, 0.05, 0.1, 0.15])
plt.xlabel('Bed volume')
plt.text(lx, ly, '(d)', transform=ax4.transAxes)
plt.setp(ax4.get_yticklabels(), visible=False)

fig = plt.gcf()
fig.subplots_adjust(wspace=0.08, hspace=0.08)
fig.set_size_inches(8, 6)
plt.savefig('hist2.png')
#plt.savefig('hist.pdf')
#plt.show()

