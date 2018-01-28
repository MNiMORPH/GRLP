from matplotlib import pyplot as plt
import numpy as np
plt.ion()

P_AQ = np.linspace(0,1,101)

def calc_PxB(P_AQ, theta):
    P_xB = 7/4. * P_AQ - 49/24. * theta
    return P_xB
    
P_xB_05 = calc_PxB(P_AQ, 0.5)
P_xB_04 = calc_PxB(P_AQ, 0.4)
P_xB_00 = calc_PxB(P_AQ, 0.0)
#P_xB_07 = calc_PxB(P_AQ, 0.7)

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.fill_between(P_AQ[P_AQ <= 0.5], P_xB_04[P_AQ <= 0.5], P_xB_5[P_AQ <= 0.5], color='0.8')
ax.fill_between(P_AQ, P_xB_00, color='0.85')
#ax.fill_between(P_AQ[P_AQ >= 0.5], P_xB_05[P_AQ >= 0.5], P_xB_04[P_AQ >= 0.5], color='0.5')
ax.fill_between(P_AQ, P_xB_05, P_xB_04, color='0.5')
#for theta in [0.0, 0.1, 0.2, 0.3, 0.6, 0.7, 0.8, 0.9, 1]:
#    plt.plot(P_AQ, calc_PxB(P_AQ, theta), 'k-', linewidth=1)
plt.plot(P_AQ, P_xB_04, 'k-', linewidth=2)
plt.plot(P_AQ, P_xB_05, 'k-', linewidth=2)
#plt.plot(P_AQ, P_xB_07, 'k--', linewidth=2)
#ax.fill_between(P_AQ, , color='0.4')
ax.set_ylim(0,1)
ax.set_xlabel(r'Drainage area to discharge exponent, $P_{AQ}$', fontsize=20)
ax.set_ylabel('Distance to valley width\n'+r'exponent, $P_{xB}$', fontsize=20)
plt.tight_layout()
plt.show()
