import numpy as np
from matplotlib import pyplot as plt

import grlp
#reload(grlp)

S0 = 1E-2
P_xB = 0.2
#z1 = 1000
z1 = 0

Qamp = 0.5
dt = 3.15E7 * 1E3
Qperiod = 100000.*3.15E7
#nt = 2 * int(np.ceil(Qperiod/dt))
nt = int(np.ceil(Qperiod*5/dt))+1
Bmax = 100.

lp = grlp.LongProfile()
self = lp

self.bcr = z1
self.intermittency = 0.3

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=1000, nx=360, x0=4E4)
lp.set_z(S0=-S0, z1=z1)
lp.set_A(k_xA=1.)
lp.set_Q(q_R=0.01, A_R=1E4)
lp.set_B(k_xB=Bmax/np.max(lp.x**P_xB), P_xB=P_xB)
#lp.set_B(k_xB=100., P_xB=0.)
lp.set_uplift_rate(0)
lp.set_niter(3)
lp.set_z_bl(z1)
Qs0 = lp.k_Qs * lp.Q[0] * (S0)**(7/6.)
lp.set_Qs_input_upstream(Qs0)
#lp.set_bcr_Dirichlet()
#lp.set_uplift_rate(0.01/3.15E7)



x0 = lp.x.copy()
Q0 = lp.Q.copy()

"""
lp.Q = Q0*.5
lp.set_Qs_input_upstream(Qs0)
lp.evolve_threshold_width_river(150, 1E11)
z_min_eq = lp.z.copy()

lp.Q = Q0 * 1.5
lp.set_Qs_input_upstream(Qs0)
lp.evolve_threshold_width_river(150, 1E11)
z_max_eq = lp.z.copy()
"""

print( "spin-up" )
lp.Q = Q0
lp.set_Qs_input_upstream(Qs0)
lp.evolve_threshold_width_river(1000, 3.15E7 * 1000.)
z0 = lp.z.copy()

lp.analytical_threshold_width(P_xB=P_xB)
lp.compute_Q_s()

Qs0_mean = np.mean(lp.Q_s)

ts = np.arange(nt)
t = dt * ts# - 3 * Qperiod
Qmult = np.sin(t / Qperiod * 2 * np.pi) * Qamp + 1

zmax = z0.copy()
zmin = z0.copy()

zall = []
Q_s_out_list = []
S_out_list = []
mean_curvature = []
mean_driver = []
Qs_all = []

print( "time-series" )
for i in range(nt):
    lp.Q = Q0 * Qmult[i]
    lp.set_Qs_input_upstream(Qs0)
    lp.compute_Q_s() # slope too
    lp.evolve_threshold_width_river(5, dt/5.)
    zmax = np.maximum(zmax, lp.z)
    zmin = np.minimum(zmin, lp.z)   
    zall.append(lp.z.copy())
    Q_s_out_list.append(lp.Q_s[-1])
    S_out_list.append(lp.S[-1])
    mean_curvature.append(np.mean(np.diff(lp.S)))
    #mean_driver.append( (7/6.*np.mean( 1/((lp.S[:-1]+lp.S[1:])/2.)  * np.diff(lp.S)/lp.dx)
    #                     + np.mean(1/((lp.Q[:-1]+lp.Q[1:])/2.) * np.diff(lp.Q)/lp.dx)
    #                     - np.mean( 1/((lp.B[:-1]+lp.B[1:])/2.) * np.diff(lp.B)/lp.dx))
    #                    * np.mean(lp.Q/lp.B * lp.S * np.abs(lp.S)**(7/6.)) )
    Sint = np.diff( (lp.z_ext[:-1]+lp.z_ext[1:])/2. )
    mean_driver.append( np.mean( (7/6.*( 1/((Sint[:-1]+Sint[1:])/2.)  * np.diff(Sint)/lp.dx)
                         + (1/((lp.Q[:-1]+lp.Q[1:])/2.) * np.diff(lp.Q)/lp.dx)
                         - ( 1/((lp.B[:-1]+lp.B[1:])/2.) * np.diff(lp.B)/lp.dx))
                        * (lp.Q[:-1]+lp.Q[1:])/2./((lp.B[:-1]+lp.B[1:])/2.) * ((Sint[:-1]+Sint[1:])/2.) * np.abs(((Sint[:-1]+Sint[1:])/2.))**(7/6.)) )
    Qs_all.append(lp.Q_s)

"""
plt.plot(Q_s_out_list_20k, 'r')
plt.plot(Q_s_out_list_40k, 'b')
plt.plot(Q_s_out_list_100k, 'k')

S_out_list_20k = np.array(S_out_list)

plt.plot(S_out_list_20k, 'r')
plt.plot(S_out_list_40k, 'b')
plt.plot(S_out_list_100k, 'k')
plt.plot(S_out_list_400k, 'g')
"""



print( "plotting" )

"""
plt.ion()
fig = plt.figure(figsize=(12,9))
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
ax1.set_xlabel('Downstream distance [km]', fontsize=26)
ax1.set_ylabel('Elevation [m]', fontsize=26)
ax2.set_xlabel('Time [kyr]', fontsize=26)
ax2.set_ylabel('Discharge multiplier [-]', fontsize=26)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax1.set_ylim((0, 500))
ax2.set_xlim((0, 40))
plt.tight_layout()

for i in np.hstack((np.zeros(5), range(nt))).astype(int):
    i = int(i)
    if i%1 == 0:
        #fig = plt.figure(figsize=(12,9))
        #ax1 = plt.subplot(211)
        #ax2 = plt.subplot(212)
        ax1.cla()
        ax2.cla()
        #fig.clf()
        #ax1.set_xlabel('Downstream distance [km]', fontsize=26)
        #ax1.set_ylabel('Elevation [m]', fontsize=26)
        #ax2.set_xlabel('Time [kyr]', fontsize=26)
        #ax2.set_ylabel('Discharge multiplier [-]', fontsize=26)
        #ax1.tick_params(axis='both', which='major', labelsize=16)
        #ax2.tick_params(axis='both', which='major', labelsize=16)
        ax1.plot(lp.x/1000., z0, '0.7', linewidth=6)
        ax1.plot(lp.x/1000., zall[i], '0.3', linewidth=2)
        ax2.plot(t/3.15E10, Qmult, 'k-', linewidth=3)
        ax2.plot(t[i]/3.15E10, Qmult[i], 'ko')
        #ax1.set_ylim((0, 500))
        ax1.set_ylim((0, 800))
        #ax1.set_ylim((500, 1100))
        ax1.set_xlim((40, 200))
        #ax2.set_xlim((0, 40))
        print i
        fig.tight_layout()
        plt.pause(0.01)
        #plt.savefig('LP_'+'%06d' %(t[i]/3.15E7)+'.png')
        #plt.close()
#ax1.plot(lp.x/1000., zmax, '--', color='.0')
#ax1.plot(lp.x/1000., zmin, '--', color='.0')
#ax1.plot(lp.x/1000., z_max_eq, '--', color='.5')
#ax1.plot(lp.x/1000., z_min_eq, '--', color='.5')

# long profile with envelope of total change
"""

plt.figure()
plt.plot(lp.x, np.max(zall[-20:], axis=0) - np.min(zall[-20:], axis=0), 'k-', linewidth=4)
plt.xlabel('Downstream distance [km]', fontsize=16)
plt.ylabel('Maximum elevation difference\nAcross one cycle [m]', fontsize=16)
plt.title( str(int(np.round(Qperiod/(1000*3.15E7)))) + ' kyr cycle', fontsize=16, fontweight='bold')
plt.tight_layout()

"""
plt.figure()
plt.plot(lp.x, dz20, '.7', linewidth=4, label='20 kyr cycle')
plt.plot(lp.x, dz40, '.4', linewidth=4, label='40 kyr cycle')
plt.plot(lp.x, dz100, '.1', linewidth=4, label='100 kyr cycle')
plt.xlabel('Downstream distance [km]', fontsize=16)
plt.ylabel('Maximum elevation difference\nAcross one cycle [m]', fontsize=16)
#plt.title( str(int(np.round(Qperiod/(1000*3.15E7)))) + ' kyr cycle', fontsize=16, fontweight='bold')
plt.legend()
plt.tight_layout()
"""


"""
# Glacial inputs
lp.Q += 1E3
lp.set_Qs_input_upstream(Qs0*2.)
for i in range(10):
    lp.evolve_threshold_width_river(1, 1E9)
    plt.plot(lp.x/1000., lp.z, '0.3', linewidth=2)
    plt.pause(.1)
# Now only excess water and a bit of sediment
#lp.Q += 5E3
for i in range(10):
    lp.set_Qs_input_upstream(Qs0*(2-(i+1)/10.))
    lp.evolve_threshold_width_river(1, 1E9)
    plt.plot(lp.x/1000., lp.z, '0.3', linewidth=2)
    plt.pause(.1)
for i in range(10):
    lp.evolve_threshold_width_river(1, 1E9)
    plt.plot(lp.x/1000., lp.z, '0.3', linewidth=2)
    plt.pause(.1)
# Now cut off the water
self.Q -= 1E3
lp.set_Qs_input_upstream(Qs0*1)
for i in range(10):
    lp.evolve_threshold_width_river(1, 1E9)
    plt.plot(lp.x/1000., lp.z, '0.3', linewidth=2)
    plt.pause(.1)
"""
