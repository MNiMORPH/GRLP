import numpy as np
from matplotlib import pyplot as plt
#plt.ion()
from scipy.stats import linregress
#plt.ion()

z0 = 1
z1 = 0
x0 = 3000.
dx = 1000.
x = np.arange(0,1E5+1,dx) + x0
x1 = np.max(x)

P_xB = .8
P_xA = 7/4.
P_AQ = 1.
P_xQ = P_xA * P_AQ

e = 1 + 6*(P_xB - P_xQ)/7.

z = (z1-z0) * (x**e - x0**e)/(np.max(x)**e - x0**e) + z0


S = np.abs( (z[2:] - z[:-2]) / (2*dx) )
Ascaled = 100* x[1:-1]**P_xA

logS = np.log10(S)
logA = np.log10(Ascaled)

out = linregress(logA, logS)
print( "Concavity = ", -out.slope )
print( "R2 = ", out.rvalue**2. )

plt.plot(x,z); plt.show()
