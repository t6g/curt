import numpy as np
import math  as math
import scipy.stats as stats

"""
solve the diffusion equation in a sphere with constant c at the surface
  Crank 1975 The Mathematics of Diffusion Eq. 6.18
"""

def sdl_analytical_rt(r, t, a, D, tol = 1.0e-20, nmax = 20):
  pi = 4.0 * math.atan(1.0)
  c  = 1.0
  mp = 2.0 * a / r /pi
  n  = 1
  
  dc = 1.0e10

  # nmax is used to avoid pre-mature exit from loop
  while abs(mp*dc) > tol or n < nmax:
    fn = float(n) 
    dc = (-1)**n / fn * math.sin(fn * pi * r / a) * math.exp(-D * fn * fn * pi * pi * t / a / a)
    #print 'iteraion = ', n, 'concentration = ', c,  'update=', mp*dc
    c = c + mp * dc
    n = n + 1

  return c

def sdl_analytical(array_r, array_t, a, D, tol = 1.0e-20, nmax = 20):
  c = np.zeros((len(array_r), len(array_t)))
  i = 0
  for r in array_r:
    j = 0
    for t in array_t:
      c[i, j] = sdl_analytical_rt(array_r[i], array_t[j], a, D, tol, nmax)
      j = j + 1
    i = i + 1
  return c

def sdl_numerical(radius_sphere, diffusion_coefficient, kp, simulation_duration, nr, nt):
  pi = 4.0 * math.atan(1.0)

  dr = radius_sphere / float(nr)
  dt = simulation_duration / float(nt)

  volu = np.zeros(nr)
  area = np.zeros(nr)
  radi = np.zeros(nr)
  dist = np.zeros(nr)

  r = 0.0
  alpha = 0.0

  for i in range(nr):
    vol0    = 4.0 * pi * r * r * r / 3.0
    r = r + dr
    radi[i] = r
    dist[i] = dr
    volu[i] = 4.0 * pi * r * r * r / 3.0 - vol0
    area[i] = 4.0 * pi * r * r

  dist[nr-1] = dist[nr-1]/2.0 

  for i in range(nr):
    volu[i] = volu[i] * (1.0 + kp)

  #matrix
  Kmatrix = np.zeros([nr, nr])

  #inner cell 0 
  Kmatrix[0, 0] = -diffusion_coefficient * area[0] / dist[0] / volu[0] * dt
  Kmatrix[0, 1] = diffusion_coefficient * area[0] / dist[0] / volu[0] * dt
  
  for i in range(nr-1):
    Kmatrix[i+1, i] = diffusion_coefficient * area[i] / dist[i] / volu[i+1] * dt
    Kmatrix[i+1, i+1] = -diffusion_coefficient * (area[i]/dist[i] + area[i+1] / dist[i+1]) / volu[i+1] * dt
    if i == nr - 2:
      alpha = diffusion_coefficient * area[i+1] / dist[i+1] / volu[i+1] * dt
    if i < nr - 2:
      Kmatrix[i+1, i+2] = diffusion_coefficient * area[i+1] / dist[i+1] / volu[i+1] * dt

  #outer cell
  #Kmatrix[nr, nr-1] = diffusion_coefficient * area[nr-1] / dist[nr-1] / volu[nr]
  #Kmatrix[nr, nr]   = -diffusion_coefficient * area[nr-1] / dist[nr-1] / volu[nr]

  invk = np.linalg.inv(np.matrix(np.identity(nr), copy=False) - Kmatrix)
  print np.matrix(np.identity(nr), copy=False) - Kmatrix

  # initial
  tt = np.asarray(np.arange(0.0, simulation_duration, dt))
  cc = np.zeros((nr, len(tt)))
  c0 = np.zeros((nr, 1))
  c1 = np.zeros_like(c0)
  

  for i in range(len(tt)-1):
    c0[nr-1, 0] = c0[nr-1, 0] + alpha
    c1 = invk * c0
    cc[:, i+1] = np.asarray(c1)[:, 0]
    c0 = c1 

  results = {}
  results['c'] = cc
  results['t'] = tt 
  results['r'] = radi
  return results


