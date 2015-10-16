import numpy as np
import math  as math
import matplotlib.pyplot as plt
import h5py as h5py

def readh5file(h5file, nx, nr):
  c1  = []
  c2  = []
  c3  = []
  tt  = []
  t0  = 0.0

  f0  = h5py.File(h5file, 'r')
 
  for i in f0.keys():
    c1.append(f0[i][:, 0][nx/2-1])
    c2.append(f0[i][:, 0][nx-1])
    t = i.split()
    tt.append(t0 + float(t[2]))
  f0.close()

  res = {}
  res['t'] = np.asarray(tt)
  res['c1'] = np.asarray(c1)
  res['c2'] = np.asarray(c2)

  return res

"""
Ogata-Banks solution
\frac{c}{c_0} = \frac{1}{2}\left[erfc\left(\frac{x-vt}{2\sqrt{Dt}}\right) - \exp\left(\frac{xv}{D}\right)erfc\left(\frac{x+vt}{2\sqrt{Dt}}\right)\right]

"""

def Ogata_Banks(x, t, v, D):
  c = 0.5 * ( math.erfc((x - v * t) / 2.0 / math.sqrt(D * t)) - \
              math.exp(x * v / D) * math.erfc((x + v * t)/2.0 / math.sqrt(D * t)))
  return c

r1 = readh5file('cde_10_5.h5', 10, 5)
r2 = readh5file('cde_100_5.h5', 100, 5)
r3 = readh5file('cde_1000_5.h5', 1000, 5)

cin = 1.9941e-07;

c1 = [] 
c2 = [] 

for t in r3['t']:
  ct = Ogata_Banks(2.5, t, 1.0, 0.1)
  c1.append(ct)
  ct = Ogata_Banks(5.0, t, 1.0, 0.1)
  c2.append(ct)

lx = 0.08
ly = 0.90

ax1 = plt.subplot(1, 2, 1)
plt.plot(r1['t'], r1['c1']/cin, 'b-', r2['t'], r2['c1']/cin, 'r-', r3['t'], r3['c1']/cin, 'g-', r3['t'], c1, 'm-')
lgd = plt.legend(('nx = 20', 'nx = 100', 'nx = 1000', 'Analytical'),loc=4)
lgd.draw_frame(False)
#txt = lgd.get_texts()
#plt.setp(txt, fontsize='small') 
plt.ylabel('c/c$_0$')
plt.xlim([0, 10])
plt.ylim([0, 1.2])
plt.xlabel('Time (s)')
plt.text(lx, ly, '(a) X = 0.5L', transform=ax1.transAxes)


ax2 = plt.subplot(1, 2, 2)
plt.plot(r1['t'], r1['c2']/cin, 'b-', r2['t'], r2['c2']/cin, 'r-', r3['t'], r3['c2']/cin, 'g-', r3['t'], c2, 'm-')
plt.xlim([0, 10])
plt.ylim([0, 1.2])
plt.setp(ax2.get_yticklabels(), visible=False)
plt.xlabel('Time (s)')
plt.xlabel('Time (s)')
plt.text(lx, ly, '(b) X = 1.0L', transform=ax2.transAxes)

fig = plt.gcf()
fig.subplots_adjust(left=0.08, right=0.95, wspace=0.05)
fig.set_size_inches(10, 4.8)
plt.savefig('comp.png')
plt.show()

