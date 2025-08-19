import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse import spdiags, identity
from scipy.sparse.linalg import spsolve, isolve

dx = 1000.
nx = 75
x = dx * np.arange(nx) + 5*dx
R = 700 / 1000. / 3.15E7 # 700 mm/yr
P_xA = 7/4. # 1/Hack exponent
#P_AQ = .8#0.8
P_AQ = 0.7
z_init = np.arange(nx*10,0.,-10) # copy/paste

def calc(nt, U, Qs_mult=1):
  dx = 1000.
  nx = 75
  x = dx * np.arange(nx) + 5*dx
  #R = 700 / 1000. / 3.15E7 # 700 mm/yr
  P_xA = 7/4.
  P_AQ = .7#0.8
  #Q = (x + 100*dx)**1.5 * R #/ 1E5
  #Ql = (x[0] + 100*dx - dx)**1.5 * R
  #Qr = (x[-1] + 100*dx + dx)**1.5 * R
  coeff = 0.041 * 0.01/3600. * 10E6**(1-P_AQ) # n/3600 cm/hr, storm size = 10 km2
  Q = (x)**(P_xA * P_AQ) * coeff
  Ql = (x[0] - dx)**(P_xA * P_AQ) * coeff #R*x[0]**P_AQ
  Qr = (x[-1] + dx)**(P_xA * P_AQ) * coeff #R*x[-1]**P_AQ
  #Ql = (x[0] - dx)**(P_xA * P_AQ) * 2E-5 * 2.#R*x[0]**P_AQ
  #Qr = (x[-1] + dx)**(P_xA * P_AQ) * 2E-5 #R*x[-1]**P_AQ
  Q_ext = np.hstack((Ql, Q, Qr))
  #D = .12 * np.ones(nx)
  #D = np.linspace(0.2, 0.04, nx)
  z = np.arange(nx*10,0.,-10)
  z0 = z.copy()
  #Q = 0.1 * (x/1000.)**1.5
  #dQdx = 0.1E-3 * 1.5 * (x/1000.)**0.5
  dQdx = (Q_ext[2:] - Q_ext[:-2]) / (2*dx)

  lambda_p = 0.35
  rho_s = 2650.
  rho = 1000.
  g = 9.805
  epsilon = 0.2
  tau_star_c = 0.0495
  phi = 3.97

  dt = 1E10 # second

  #k0 = (2.88/(1-lambda_p)) * ((rho_s - rho)/rho)**(13/6.) * g**.5 * (1 + epsilon)**2.5 * epsilon**1.5 * tau_star_c**(19/6.)

  #k0 = 0.67 / ((rho_s - rho)/rho)**(7/6.) * epsilon**1.5/(1 + epsilon)**(5/3.) * (1/tau_star_c**(5/3.))
  k0 = 0.17 * phi * epsilon**(3/2.) / ( ((rho_s - rho)/rho)**(7/6.) * (1 + epsilon)**(5/3.) * tau_star_c**(1/6.) )
  
  C0 = k0/(1-lambda_p) * dt / (2*dx)

  # BC for implicit
  # explicit has padded array
  # this is shoddy but I just want to see if something works!
  #bcl = np.max(z) + np.mean(np.diff(z))
  #bcl0 = (z[1] - z[0]) * 2 * ( (7/3.) * Q[0]/dx - (Q_ext[2] - Q_ext[0])/(2*dx))
  # Ghost nodes on left -- assuming same as starting slope
  S0 = (z[1] - z[0])/dx * Qs_mult**(6/7.)
  Qs0 = 0.041 * Q[0] * np.abs(S0)**(7/6.)
  bcl0 = -2 * dx * S0 * ( (7/3.) * Q[0]/dx - dQdx[0])# * 200
  bcl = bcl0.copy()
  bcr = 0
  bcarray0 = np.hstack((bcl, np.zeros(len(z)-2), bcr ))

  z_ext = np.hstack(( bcl, z, bcr ))

  niter = 3

  #plt.plot(x, z, 'k-')
  #z[0] += 100
  for ti in range(int(nt)):
    zi = z_ext.copy()
    bcr += U * dt / 3.15E7
    if ti == 100:
      pass
      #bcr += 10
      #S0 *= 10
      #bcl = -2 * dx * S0 * ( (7/3.) * Q[0]/dx - dQdx[0])# * 200
      #bcarray0 = np.hstack(( bcl, np.zeros(len(z)-2), bcr ))
    #if ti >= 400:
      #bcr += .2
    bcarray0 = np.hstack(( bcl, np.zeros(len(z)-2), bcr ))
    for i in range(niter):
      # Sign check, pull sign out here
      # Just assuming negative for now
      # Should have extra - for implicit solution
      dzdt_0 = np.abs( (zi[2:] - zi[:-2]) / (2*dx) )**(1/6.)
      sign_S = - np.sign(zi[2:] - zi[:-2])
      #B = 1#x**2#1.#x**0.5
      B = 10.#Q*2
      C1 = C0 * dzdt_0 / B# * sign_S # USED TO BE *B; MISTAKE?
      #bcl = bcl0 * C1[0]
      #bcarray0 = np.hstack(( bcl, np.zeros(len(z)-2), bcr ))
      left = -C1 * ( (7/3.) * Q/dx - dQdx )
      center = C1 * 2 * ( (7/3.) * Q/dx ) + 1
      right = -C1 * ( (7/3.) * Q/dx + dQdx )
      right[0] = -2 * C1[0] * 7/3. * Q[0]/dx
      #right[0] *= 2
      left = np.roll(left, -1)
      right = np.roll(right, 1)
      # NEUMANN
      diagonals = np.vstack((left, center, right))
      offsets = np.array([-1, 0, 1])
      A = spdiags(diagonals, offsets, len(z), len(z), format='csr')
      #spsolve(A, zi[1:-1] + C1 * (-7/6. * Q[1:-1]/(dx**2) - dQdx[1:-1]/(2*dx)) * bcarray)
      #zi[1:-1] = spsolve(A, zi[1:-1] + C1*bcarray*-1)
      # FIX Q, DQDX
      bcarray = bcarray0.copy()
      bcarray[0] *= C1[0]# * ( (7/3.) * Q[0]/dx - dQdx[0] )
      bcarray[-1] *= C1[-1] * ( (7/3.) * Q[-1]/dx + dQdx[-1] )
      zi[1:-1] = spsolve(A, z + bcarray)
      #spsolve(A, z + C1 * (-7/6. * Q/(dx**2) - dQdx/(2*dx)) * bcarray)
      #zi[1:-1] = spsolve(A, z + C1 * (-7/6. * Q/(dx**2) - dQdx/(2*dx)) * bcarray)
    z = zi[1:-1]
    print z[0]
  print ""
  print ""
  #print left
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
        
  
#xS, zS = calc(200, 0.001)
x0, z0 = calc(200, 0.0)
#xU, zU = calc(200, -0.001)

x = dist = x0.copy()
z = z0.copy()

"""
Ascaled = 100* x**(7/4.)

S = -np.diff(z) / np.diff(dist)
AS = (Ascaled[1:] + Ascaled[:-1])/2.
distS = (dist[1:] + dist[:-1])/2.
zS = (z[1:] + z[:-1])/2.

xm, Sm = running_mean(distS, S, 300)
xm, ASm = running_mean(distS, AS, 300)

plt.figure()
for power in np.linspace(0.2,0.25, 17):
    chi = np.cumsum((1E6/AS[::-1])**power * -np.diff(dist[::-1]))
    plt.plot(chi,zS[::-1])
"""


#zA = analytical0uplift(x, x[10], x[-1], z[10], z[-1])

#x,z = xU,zU

#zmin = np.min(np.hstack((zS, z0, zU)))

plt.ion()
#plt.plot(x/1000.,  zA + z0[0] - zA[0], '.5', linewidth=8)
plt.figure()
#plt.plot(xU/1000., zU + z0[0] - zU[0], 'r-', linewidth=4)
#plt.plot(xS/1000., zS+ z0[0] - zS[0], 'b-', linewidth=4)
plt.plot(x0/1000., z0, 'k-', linewidth=4)
#plt.plot(xU/1000., zU, 'r-', linewidth=4)
#plt.plot(xS/1000., zS, 'b-', linewidth=4)
#plt.plot(x0/1000., z0, 'k-', linewidth=4)
plt.xlabel('Downstream distance [km]', fontsize=24)
plt.ylabel('Elevation [m]', fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.show()

from scipy.stats import linregress

for x, z in [[xS, zS], [x0, z0], [xU, zU]]:
  S = np.abs( (z[2:] - z[:-2]) / (2*dx) )
  Ascaled = 100* x[1:-1]**P_xA

  logS = np.log10(S)
  logA = np.log10(Ascaled)

  out = linregress(logA, logS)
  print "Concavity = ", -out.slope
