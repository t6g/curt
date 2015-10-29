import numpy as np
import math  as math
import matplotlib.pyplot as plt
import h5py as h5py

def readh5file(h5file, nx, nr):
  length = 5.0
  cc  = []
  tt  = []
  t0  = 0.0

  f0  = h5py.File(h5file, 'r')

  j = 0 
  for i in f0.keys():
    if j == 0:
      cc = f0[i][0:nx, 0]
      j = j + 1
    else:
      cc = np.column_stack((cc, f0[i][0:nx, 0]))
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  res = {}
  res['t'] = np.asarray(tt)
  res['c'] = np.asarray(cc)
  res['x'] = np.linspace(0, length, nx)

  return res

"""
Ogata-Banks solution
\frac{c}{c_0} = \frac{1}{2}\left[erfc\left(\frac{x-vt}{2\sqrt{Dt}}\right) - \exp\left(\frac{xv}{D}\right)erfc\left(\frac{x+vt}{2\sqrt{Dt}}\right)\right]

"""

def Ogata_Banks(x, t, v, D):
  c = 0.5 * ( math.erfc((x - v * t) / 2.0 / math.sqrt(D * t)) - \
              math.exp(x * v / D) * math.erfc((x + v * t)/2.0 / math.sqrt(D * t)))
  return c

r1 = readh5file('cde_20_5.h5', 10, 5)
r2 = readh5file('cde_100_5.h5', 100, 5)
r3 = readh5file('cde_1000_5.h5', 1000, 5)

cin = 1.9941e-07;

c1 = [] 
c2 = [] 
c3 = [] 
c4 = [] 

#t = 1, 2, 3, 4 s
for x in r3['x']:
  ct = Ogata_Banks(x, 1.0, 1.0, 0.1)
  c1.append(ct)
  ct = Ogata_Banks(x, 2.0, 1.0, 0.1)
  c2.append(ct)
  ct = Ogata_Banks(x, 3.0, 1.0, 0.1)
  c3.append(ct)
  ct = Ogata_Banks(x, 4.0, 1.0, 0.1)
  c4.append(ct)

lx = 0.08
ly = 0.90

ax1 = plt.subplot(2, 2, 1)
i = 1
j = 4
plt.plot(r1['x'], r1['c'][:, i]/cin, 'b-', r2['x'], r2['c'][:, j]/cin, 'r-', r3['x'], r3['c'][:, j]/cin, 'g-', r3['x'], c1, 'm-')
plt.ylabel('c/c$_0$')
plt.xlim([0, 5])
plt.ylim([0, 1.2])
plt.text(lx, ly, '(a) t = 1', transform=ax1.transAxes)


ax2 = plt.subplot(2, 2, 2)
i = 3
j = 9
plt.plot(r1['x'], r1['c'][:, i]/cin, 'b-', r2['x'], r2['c'][:, j]/cin, 'r-', r3['x'], r3['c'][:, j]/cin, 'g-', r3['x'], c2, 'm-')
plt.xlim([0, 5])
plt.ylim([0, 1.2])
plt.setp(ax2.get_yticklabels(), visible=False)
plt.text(lx, ly, '(b) t = 2', transform=ax2.transAxes)

ax3 = plt.subplot(2, 2, 3)
i = 5
j = 14
plt.plot(r1['x'], r1['c'][:, i]/cin, 'b-', r2['x'], r2['c'][:, j]/cin, 'r-', r3['x'], r3['c'][:, j]/cin, 'g-', r3['x'], c3, 'm-')
plt.xlim([0, 5])
plt.ylim([0, 1.2])
plt.setp(ax2.get_yticklabels(), visible=False)
plt.text(lx, ly, '(c) t = 3', transform=ax3.transAxes)
plt.xlabel('Time (s)')
plt.ylabel('c/c$_0$')

ax4 = plt.subplot(2, 2, 4)
i = 7
j = 19
plt.plot(r1['x'], r1['c'][:, i]/cin, 'b-', r2['x'], r2['c'][:, j]/cin, 'r-', r3['x'], r3['c'][:, j]/cin, 'g-', r3['x'], c4, 'm-')
plt.xlim([0, 5])
plt.ylim([0, 1.2])
plt.setp(ax2.get_yticklabels(), visible=False)
plt.text(lx, ly, '(c) t = 3', transform=ax3.transAxes)
plt.xlabel('Time (s)')
plt.setp(ax4.get_yticklabels(), visible=False)
lgd = plt.legend(('nx = 20', 'nx = 100', 'nx = 1000', 'Analytical'),loc=3)
lgd.draw_frame(False)
#txt = lgd.get_texts()
#plt.setp(txt, fontsize='small') 

fig = plt.gcf()
fig.subplots_adjust(left=0.08, right=0.95, wspace=0.05, hspace=0.08)
fig.set_size_inches(12, 10)
plt.savefig('prof.pdf')
plt.show()
