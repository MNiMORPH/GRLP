import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse import spdiags, identity
from scipy.sparse.linalg import spsolve, isolve

dx = 1000.
nx = 75
x = dx * np.arange(nx) + 5*dx
q_R = 0.001E-5/3600. #700 / 1000. / 3.15E7 # 700 mm/yr
P_xA = 7/4. # 1/Hack exponent
P_AQ = 1.#0.8
P_xB = 0.8
z_init = np.arange(nx*10,0.,-10) # copy/paste

#S0 = (z[1] - z[0])/dx * Qs_mult**(6/7.)
#Qs0 = 0.041 * Q[0] * np.abs(S0)**(7/6.)

def calc(nt, U, Qs_mult=1, x=None, z=None, bcr=0, dt=1E9, S0=-1E-2, P_xA=7/4., P_AQ=None, P_xB=None):
  dx = 1000.
  nx = 75
  if x is None:
      x = dx * np.arange(nx) + 5*dx
  #R = 700 / 1000. / 3.15E7 # 700 mm/yr
  #P_xA = 7/4.
  #P_AQ = .8#0.8
  #k_xB = 1.
  #P_xB = 0.3#0.5
  k_xB = 10./np.max(x**P_xB)
  #Q = (x + 100*dx)**1.5 * R #/ 1E5
  #Ql = (x[0] + 100*dx - dx)**1.5 * R
  #Qr = (x[-1] + 100*dx + dx)**1.5 * R
  coeff = q_R * 10E6**(1-P_AQ) # n/3600 cm/hr, storm size = 10 km2 * 103
  Q = (x)**(P_xA * P_AQ) * coeff
  Ql = (x[0] - dx)**(P_xA * P_AQ) * coeff #R*x[0]**P_AQ
  Qr = (x[-1] + dx)**(P_xA * P_AQ) * coeff #R*x[-1]**P_AQ
  #Ql = (x[0] - dx)**(P_xA * P_AQ) * 2E-5 * 2.#R*x[0]**P_AQ
  #Qr = (x[-1] + dx)**(P_xA * P_AQ) * 2E-5 #R*x[-1]**P_AQ
  Q_ext = np.hstack((Ql, Q, Qr))
  Q_ext_0 = Q_ext.copy()
  #D = .12 * np.ones(nx)
  #D = np.linspace(0.2, 0.04, nx)
  if z is None:
      z = np.arange(nx*10,0.,-10)
  z0 = z.copy()
  #Q = 0.1 * (x/1000.)**1.5
  #dQdx = 0.1E-3 * 1.5 * (x/1000.)**0.5
  dQ = (Q_ext[2:] - Q_ext[:-2]) # What about 2 before dx??? down below?
  #dQdx = (Q_ext[2:] - Q_ext[:-2]) / (2*dx)
  
  B = k_xB*x**P_xB
  Bl = k_xB*(x[0] - dx)**P_xB
  Br = k_xB*(x[-1] + dx)**P_xB
  B_ext = np.hstack((Bl, B, Br))

  dB = (B_ext[2:] - B_ext[:-2]) # What about 2 before dx??? down below?
  #dBdx = (B_ext[2:] - B_ext[:-2]) / (2*dx)

  
  lambda_p = 0.35
  rho_s = 2650.
  rho = 1000.
  g = 9.805
  epsilon = 0.2
  tau_star_c = 0.0495
  phi = 3.97

  #dt = 1E9 # second

  #k0 = (2.88/(1-lambda_p)) * ((rho_s - rho)/rho)**(13/6.) * g**.5 * (1 + epsilon)**2.5 * epsilon**1.5 * tau_star_c**(19/6.)

  #k0 = 0.67 / ((rho_s - rho)/rho)**(7/6.) * epsilon**1.5/(1 + epsilon)**(5/3.) * (1/tau_star_c**(5/3.))
  kQs = 0.17 * phi * epsilon**(3/2.) / ( ((rho_s - rho)/rho)**(7/6.) * (1 + epsilon)**(5/3.) * tau_star_c**(1/6.) )
  
  C0 = kQs/(1-lambda_p) * dt / dx**2

  # BC for implicit
  # explicit has padded array
  # this is shoddy but I just want to see if something works!
  #bcl = np.max(z) + np.mean(np.diff(z))
  #bcl0 = (z[1] - z[0]) * 2 * ( (7/3.) * Q[0]/dx - (Q_ext[2] - Q_ext[0])/(2*dx))
  # Ghost nodes on left -- assuming same as starting slope
  #S0 = (z[1] - z[0])/dx * Qs_mult**(6/7.)
  #print S0
  S0 *= Qs_mult
  # CHECK!
  #bcl0 = -2 * dx * S0 * C0 * ( (14/3.) * Q[0] - dQ[0]/dx)# + Q[0]*dB[0]/B[0])# * 200
  #bcl0 = -2 * dx * S0 * ( (7/3.) * Q[0]/dx - dQdx[0])# * 200
  # bcl0 = -2 * dx * S0 * dB[0] + z_init[1]#* ( (7/3.) * Q[0]/dx - dQdx[0])# * 200 # wrong -- needs extra /(2 dx)
  #bcl0 = -4 * dx**2 * S0 / dB[0] + z_init[1]#* ( (7/3.) * Q[0]/dx - dQdx[0])# * 200
  #print -2 * dx * S0
  #print S0
  #bcl = bcl0.copy()
  #bcr = 0
  #bcarray0 = bcarray = np.hstack((bcl, np.zeros(len(z)-2), bcr ))

  # Before main loop: set up - just has to be in the ballpark
  C1 = C0 * np.sign(np.diff(z)[0]) * (np.abs(np.diff(z)[0]/dx))**(1/6.) * Q / B# * sign_S
  left = -C1 * ( (7/6.) - dQ/Q/4. + dB/B/4.)

  bcl = z[0] + 2*dx*S0*left[0]
  bcr = 0
  z_ext = np.hstack(( bcl, z, bcr ))

  niter = 3

  #plt.plot(x, z, 'k-')
  #z[0] += 100
  for ti in range(int(nt)):
    zi = z_ext.copy()
    #bcr += U * dt / 3.15E7
    for i in range(niter):
      # 0 discharge if river reverses
      Q_ext = Q_ext_0.copy()
      #Q_ext[np.hstack((False, np.diff(z) > 0))] = 0
      #dQdx = (Q_ext[2:] - Q_ext[:-2]) / (2*dx) # COMMENT?
      # Sign check, pull sign out here
      # Just assuming negative for now
      # Should have extra - for implicit solution
      dzdt_0 = np.abs( (zi[2:] - zi[:-2]) / (2*dx) )**(1/6.)
      sign_S = - np.sign(zi[2:] - zi[:-2])
      #B = 1#x**2#1.#x**0.5
      #B = Q**.5 * 20
      C1 = C0 * dzdt_0 / B# * sign_S
      left = -C1 * ( (7/6.) - dQ/Q/4. + dB/B/4.)
      center = C1 * 2 * ( (7/6.) ) + 1
      right = -C1 * ( (7/6.) + dQ/Q/4. - dB/B/4. )
      #right[0] = -2 * C1[0] * (14/3. * Q[0] + dQ[0] - Q[0]*dB[0]/B[0])
      #right[0] *= 2
      #right[0] = -2 * C1[0] * 14/3 * Q[0]
      # RHS b.c. from ghost node approach
      right[0] = -2 * C1[0] * 7/6. * Q[0] # Should be the same as below
      #right[0] = left[0] + right[0]
      # Boundary condition update: RHS from ghost-node approach
      # Could also approach this with an array - subtract from RHS
      bcl = z[0] + 2*dx*S0*left[0] # left[0] already negative, - to +
      left = np.roll(left, -1)
      right = np.roll(right, 1)
      # NEUMANN
      diagonals = np.vstack((left, center, right))
      offsets = np.array([-1, 0, 1])
      A = spdiags(diagonals, offsets, len(z), len(z), format='csr')
      RHS = np.hstack((bcl, z[1:-1], bcr))
      zi[1:-1] = spsolve(A, RHS)
      #zi[1:-1] = spsolve(A, z + bcarray)
      #zi[1:-1] = spsolve(A, np.hstack((bcl, z[1:-1] + bcarray[1:-1], bcr)))
    z = zi[1:-1]
    z_ext = zi.copy()
  return x, z
  
def analytical0upliftDirichlet(xext, x0, x1, z0, z1):
  # Analytical: no uplift
  z_Uequals0 = (z1 - z0) * \
               (xext**(-1/14.) - x0**(-1/14.)) / \
               (x1**(-1/14.) - x0**(-1/14.)) \
               + z0
  return z_Uequals0

def running_mean(x, y, dx):
    _x = x[0]
    xout = []
    yout = []
    while _x < x[-1]:
        _y = np.nanmean(y[(x >= _x) * (x < _x+dx)])
        _x = np.nanmean(x[(x >= _x) * (x < _x+dx)])
        xout.append(_x)
        yout.append(_y)
        _x += dx
    return xout, yout
        
plt.ion()
plt.figure(figsize=(12,6))
dt = 1E12
x0, z0 = calc(500, 0.0, z=z_init, dt=dt, P_xA=P_xA, P_AQ=P_AQ, P_xB=P_xB)
# Uplift pulse
plt.plot(x0/1000., z0, '0.6', linewidth=6)
plt.xlabel('Downstream distance [km]', fontsize=26)
plt.ylabel('Elevation [m]', fontsize=26)
plt.tick_params(axis='both', which='major', labelsize=16)
#plt.ylim(0, 8000)
plt.tight_layout()
plt.pause(0.5)
x = x0
z = z0
"""
U = 0#-1E-2/3.15E7 #-3.1746031746031744e-08
for i in range(20):
    z[-1] += U * dt
    x, z = calc(50, 0.0, z=z, Qs_mult=1, dt=dt)
    plt.cla()
    plt.plot(x0/1000., z0, '0.6', linewidth=6)
    plt.plot(x/1000., z, 'k-', linewidth=2)
    plt.xlabel('Downstream distance [km]', fontsize=26)
    plt.ylabel('Elevation [m]', fontsize=26)
    plt.tick_params(axis='both', which='major', labelsize=16)
    #plt.ylim(0, 8000)
    plt.text(x=50, y=80, s='dt = 10$^'+str(int(np.log10(dt)))+'$ s', fontsize=26)
    plt.tight_layout()
    plt.pause(0.01)
"""

# Analytical
P_xQ = P_xA * P_AQ
e = 1 + 6*(P_xB - P_xQ)/7.
zanalytical = (z[-1]-z[0]) * (x**e - x[0]**e)/(x[-1]**e - x[0]**e) + z[0]
plt.plot(x/1000., zanalytical, 'k', linewidth=2)

from scipy.stats import linregress

for _x, _z in [[x0, z0]]:
  S = np.abs( (_z[2:] - _z[:-2]) / (2*dx) )
  Ascaled = 100* _x[1:-1]**P_xA

  logS = np.log10(S)
  logA = np.log10(Ascaled)

  out = linregress(logA, logS)
  print "Concavity = ", -out.slope
  print "R2 = ", out.rvalue**2.
