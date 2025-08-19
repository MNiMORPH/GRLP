# Figure UpliftSubsidence in paper

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
plt.ioff()

import grlp
#reload(grlp)

S0 = 0.01
P_xB = 0.1
z1 = 0

Qamp = 0.5
dt = 3.15E7 * 1E2
nt = int(100)
Bmax = 250.

lp = grlp.LongProfile()
self = lp

self.bcr = z1

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=500, nx=180, x0=10E3)
lp.set_z(S0=-S0, z1=z1)
lp.set_A(k_xA=1.)
lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7/4.*0.7)
lp.set_B(k_xB=Bmax/np.max(lp.x**P_xB), P_xB=P_xB)
lp.set_niter(3)
lp.set_z_bl(z1)
Qs0 = lp.k_Qs * lp.Q[0] * S0**(7/6.)
lp.set_Qs_input_upstream(Qs0)

Q_s_out_1 = []
Q_s_out_2 = []
Q_s_out_3 = []
Q_s_out_4 = []
Q_s_out_5 = []

time_transient_kyr = np.arange(30, 601, 30)

fig = plt.figure(figsize=(13,15))
fig2 = plt.figure(figsize=(5,15))

# Transient 1: Sudden base-level fall
# Plot setup
ax1 = fig.add_subplot(5,2,1)
#ax1.set_xlabel('Downstream distance [km]', fontsize=14, fontweight='bold')
ax1.set_ylabel('Elevation [m]', fontsize=14, fontweight='bold')

ax1sa = fig.add_subplot(5,2,2)
#ax1sa.set_xlabel('Drainage area [km$^2$]', fontsize=14, fontweight='bold')
ax1sa.set_ylabel('Slope [-]', fontsize=14, fontweight='bold')
ax1sa.set_ylim(8E-4,2E-2)

fig2ax1 = fig2.add_subplot(5,1,1)
#fig2ax1.set_xlabel('Time [kyr]', fontsize=14, fontweight='bold')
fig2ax1.set_ylabel('$Q_\mathrm{s,out} / Q_\mathrm{s,in}$', fontsize=14, fontweight='bold')

# Starting case
U = 0.
lp.set_uplift_rate(U/3.15E7)
lp.evolve_threshold_width_river(10, 1E7*3.15E7)
lp.slope_area(verbose=True)
ax1.plot(lp.x/1000., lp.z + 500, color='.5', linewidth=4)
ax1sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=4)

z0 = lp.z.copy()
z_ext_0 = lp.z_ext.copy()

lp.intermittency = 1.

# Calculate
lp.set_z_bl(-100)
for i in range(20):
    lp.evolve_threshold_width_river(1, 3E4*3.15E7)
    lp.compute_Q_s()
    Q_s_out_1.append(lp.Q_s[-1] * lp.intermittency)
    lp.slope_area(verbose=False)
    ax1.plot(lp.x/1000., lp.z + 500, color='.5', linewidth=1, zorder=0)
    ax1sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=1)

_S_end_nouplift = np.abs(lp.S0) * (lp.x[-1]/lp.x[0])**(6/7. \
                                          * (lp.P_xB - lp.P_xA * lp.P_AQ))
Q_s_in_1 = lp.k_Qs * lp.intermittency * lp.Q[-1] * _S_end_nouplift**(7/6.) # m3/s

fig2ax1.plot(time_transient_kyr, Q_s_out_1/Q_s_in_1, 'ko', linewidth=1, zorder=0)

# New equilibrium
lp.evolve_threshold_width_river(10, 1E7*3.15E7)
lp.slope_area(verbose=True)
ax1.plot(lp.x/1000., lp.z + 500, color='0', linewidth=4)
ax1sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '--', color='0', linewidth=4)

# Transient 2: Uplift begins / accelerates

# Plot setup
ax2 = fig.add_subplot(5,2,3)
#ax2.set_xlabel('Downstream distance [km]', fontsize=14, fontweight='bold')
ax2.set_ylabel('Elevation [m]', fontsize=14, fontweight='bold')

ax2sa = fig.add_subplot(5,2,4)
#ax2sa.set_xlabel('Drainage area [km$^2$]', fontsize=14, fontweight='bold')
ax2sa.set_ylabel('Slope [-]', fontsize=14, fontweight='bold')
ax2sa.set_ylim(8E-4,2E-2)

fig2ax2 = fig2.add_subplot(5,1,2)
#fig2ax2.set_xlabel('Time [kyr]', fontsize=14, fontweight='bold')
fig2ax2.set_ylabel('$Q_\mathrm{s,out} / Q_\mathrm{s,in}$', fontsize=14, fontweight='bold')

plt.tight_layout()

# Return to old values
lp.z = z0.copy()
lp.z_ext = z_ext_0.copy()
lp.set_z_bl(z1)
lp.set_Qs_input_upstream(Qs0)
lp.slope_area(verbose=True)
ax2.plot(lp.x/1000., lp.z - lp.z[0] + 500, color='.5', linewidth=4)
ax2sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=4)

# Calculate
U = 1E-3
lp.set_uplift_rate(U/3.15E7)
for i in range(20):
    lp.evolve_threshold_width_river(1, 3E4*3.15E7)
    lp.slope_area(verbose=False)
    ax2.plot(lp.x/1000., lp.z - lp.z[0] + 500, color='.5', linewidth=1)
    ax2sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=1)
    lp.compute_Q_s()
    Q_s_out_2.append(lp.Q_s[-1] * lp.intermittency)

_S_end_nouplift = np.abs(lp.S0) * (lp.x[-1]/lp.x[0])**(6/7. \
                                          * (lp.P_xB - lp.P_xA * lp.P_AQ))
Q_s_in_2 = lp.k_Qs * lp.intermittency * lp.Q[-1] * _S_end_nouplift**(7/6.) # m3/s

fig2ax2.plot(time_transient_kyr, Q_s_out_2/Q_s_in_2, 'ko', linewidth=1, zorder=0)

# New equilibrium
lp.evolve_threshold_width_river(10, 1E14)
lp.slope_area(verbose=True)
ax2.plot(lp.x/1000., lp.z - lp.z[0] + 500, color='0', linewidth=4, zorder=100)
ax2sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='0', linewidth=4)

# Transient 3: Subsidence begins / accelerates

# Plot setup
ax3 = fig.add_subplot(5,2,5)
#ax3.set_xlabel('Downstream distance [km]', fontsize=14, fontweight='bold')
ax3.set_ylabel('Elevation [m]', fontsize=14, fontweight='bold')

ax3sa = fig.add_subplot(5,2,6)
#ax3sa.set_xlabel('Drainage area [km$^2$]', fontsize=14, fontweight='bold')
ax3sa.set_ylabel('Slope [-]', fontsize=14, fontweight='bold')
ax3sa.set_ylim(2E-4,2E-2)

fig2ax3 = fig2.add_subplot(5,1,3)
#fig2ax3.set_xlabel('Time [kyr]', fontsize=14, fontweight='bold')
fig2ax3.set_ylabel('$Q_\mathrm{s,out} / Q_\mathrm{s,in}$', fontsize=14, fontweight='bold')

plt.tight_layout()

# Return to old values
lp.z = z0.copy()
lp.z_ext = z_ext_0.copy()
lp.set_z_bl(z1)
lp.set_Qs_input_upstream(Qs0)
lp.slope_area(verbose=True)
ax3.plot(lp.x/1000., lp.z - lp.z[0] + 500, color='.5', linewidth=4)
ax3sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=4)

# Calculate
U = -5E-4
lp.set_uplift_rate(U/3.15E7)
for i in range(20):
    lp.evolve_threshold_width_river(1, 3E4*3.15E7)
    lp.slope_area(verbose=False)
    ax3.plot(lp.x/1000., lp.z - lp.z[0] + 500, color='.5', linewidth=1)
    ax3sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=1)
    lp.compute_Q_s()
    Q_s_out_3.append(lp.Q_s[-1] * lp.intermittency)

_S_end_nouplift = np.abs(lp.S0) * (lp.x[-1]/lp.x[0])**(6/7. \
                                          * (lp.P_xB - lp.P_xA * lp.P_AQ))
Q_s_in_3 = lp.k_Qs * lp.intermittency * lp.Q[-1] * _S_end_nouplift**(7/6.) # m3/s

fig2ax3.plot(time_transient_kyr, Q_s_out_3/Q_s_in_3, 'ko', linewidth=1, zorder=0)

# New equilibrium
lp.evolve_threshold_width_river(10, 1E14)
lp.slope_area(verbose=True)
ax3.plot(lp.x/1000., lp.z - lp.z[0] + 500, color='0', linewidth=4, zorder=100)
ax3sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='0', linewidth=4)

# Transient 4: Sediment supply doubles

# Plot setup
ax4 = fig.add_subplot(5,2,7)
#ax4.set_xlabel('Downstream distance [km]', fontsize=14, fontweight='bold')
ax4.set_ylabel('Elevation [m]', fontsize=14, fontweight='bold')

ax4sa = fig.add_subplot(5,2,8)
#ax4sa.set_xlabel('Drainage area [km$^2$]', fontsize=14, fontweight='bold')
ax4sa.set_ylabel('Slope [-]', fontsize=14, fontweight='bold')
ax4sa.set_ylim(8E-4,2E-2)

fig2ax4 = fig2.add_subplot(5,1,4)
#fig2ax4.set_xlabel('Time [kyr]', fontsize=14, fontweight='bold')
fig2ax4.set_ylabel('$Q_\mathrm{s,out} / Q_\mathrm{s,in}$', fontsize=14, fontweight='bold')

plt.tight_layout()

# Plot Q_s_out/Q_s_in: equilibrium

# Return to old values
lp.z = z0.copy()
lp.z_ext = z_ext_0.copy()
lp.set_z_bl(z1)
lp.set_Qs_input_upstream(Qs0)
U = 0.
lp.set_uplift_rate(U/3.15E7)
lp.slope_area(verbose=True)
ax4.plot(lp.x/1000., lp.z + 500, color='.5', linewidth=4)
ax4sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=4)

# Calculate
lp.set_Qs_input_upstream(lp.Q_s_0 * 2.)
for i in range(20):
    lp.evolve_threshold_width_river(1, 3E4*3.15E7)
    lp.slope_area(verbose=False)
    ax4.plot(lp.x/1000., lp.z + 500, color='.5', linewidth=1)
    ax4sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=1)
    lp.compute_Q_s()
    Q_s_out_4.append(lp.Q_s[-1] * lp.intermittency)
    
_S_end_nouplift = np.abs(lp.S0) * (lp.x[-1]/lp.x[0])**(6/7. \
                                          * (lp.P_xB - lp.P_xA * lp.P_AQ))
Q_s_in_4 = lp.k_Qs * lp.intermittency * lp.Q[-1] * _S_end_nouplift**(7/6.) # m3/s

fig2ax4.plot(time_transient_kyr, Q_s_out_4/Q_s_in_4, 'ko', linewidth=1, zorder=0)

# New equilibrium
lp.evolve_threshold_width_river(10, 1E14)
lp.slope_area(True)
ax4.plot(lp.x/1000., lp.z + 500, color='0', linewidth=4)
ax4sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='0', linewidth=4)

# Transient 5: Water supply doubles

# Plot setup
ax5 = fig.add_subplot(5,2,9)
ax5.set_xlabel('Downstream distance [km]', fontsize=14, fontweight='bold')
ax5.set_ylabel('Elevation [m]', fontsize=14, fontweight='bold')

ax5sa = fig.add_subplot(5,2,10)
ax5sa.set_xlabel('Drainage area [km$^2$]', fontsize=14, fontweight='bold')
ax5sa.set_ylabel('Slope [-]', fontsize=14, fontweight='bold')
ax5sa.set_ylim(5E-4,2E-2)

fig2ax5 = fig2.add_subplot(5,1,5)
fig2ax5.set_xlabel('Time [kyr]', fontsize=14, fontweight='bold')
fig2ax5.set_ylabel('$Q_\mathrm{s,out} / Q_\mathrm{s,in}$', fontsize=14, fontweight='bold')

plt.tight_layout()

# Plot Q_s_out/Q_s_in: equilibrium

# Return to old values
lp.z = z0.copy()
lp.z_ext = z_ext_0.copy()
lp.set_z_bl(z1)
lp.set_Qs_input_upstream(Qs0)
lp.set_Q(k_xQ=1.433776163432246e-05*2., P_xQ=7/4.*0.7)
U = 0.
lp.set_uplift_rate(U/3.15E7)
lp.slope_area(verbose=True)
ax5.plot(lp.x/1000., lp.z + 500, color='.5', linewidth=4)
ax5sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=4)

# Calculate
for i in range(20):
    lp.evolve_threshold_width_river(1, 3E4*3.15E7)
    lp.slope_area(verbose=False)
    ax5.plot(lp.x/1000., lp.z + 500, color='.5', linewidth=1)
    ax5sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='.5', linewidth=1)
    lp.compute_Q_s()
    Q_s_out_5.append(lp.Q_s[-1] * lp.intermittency)
    
_S_end_nouplift = np.abs(lp.S0) * (lp.x[-1]/lp.x[0])**(6/7. \
                                          * (lp.P_xB - lp.P_xA * lp.P_AQ))
Q_s_in_5 = lp.k_Qs * lp.intermittency * lp.Q[-1] * _S_end_nouplift**(7/6.) # m3/s

fig2ax5.plot(time_transient_kyr, Q_s_out_5/Q_s_in_5, 'ko', linewidth=1, zorder=0)

# New equilibrium
lp.evolve_threshold_width_river(10, 1E14)
lp.slope_area(True)
ax5.plot(lp.x/1000., lp.z + 500, color='0', linewidth=4)
ax5sa.loglog(lp.A[2:-1]/1E6, lp.S[2:-1], '-', color='0', linewidth=4)

# Curve fitting
def plf(x, a, b, c):
    return a*x**b + c

def exf(x, a, b, c):
    return a*np.exp(b*x) + c

popt,pcov = curve_fit( exf, time_transient_kyr[1:], Q_s_out_1[1:]/Q_s_in_1, p0=(1, -0.001, 1) )
fig2ax1.plot(time_transient_kyr, exf(time_transient_kyr, popt[0], popt[1], popt[2]), color='.5', linewidth=2)
print( "1:", np.abs(1/popt[1]) )
popt,pcov = curve_fit( exf, time_transient_kyr[1:], Q_s_out_2[1:]/Q_s_in_2, p0=(-1, -0.001, 0.65) )
fig2ax2.plot(time_transient_kyr, exf(time_transient_kyr, popt[0], popt[1], popt[2]), color='.5', linewidth=2)
print( "2:", np.abs(1/popt[1]) )
popt,pcov = curve_fit( exf, time_transient_kyr[1:], Q_s_out_3[1:]/Q_s_in_3, p0=(1, -0.001, 0.65) )
fig2ax3.plot(time_transient_kyr, exf(time_transient_kyr, popt[0], popt[1], popt[2]), color='.5', linewidth=2)
print( "3:", np.abs(1/popt[1]) )
popt,pcov = curve_fit( exf, time_transient_kyr[1:], Q_s_out_4[1:]/Q_s_in_4, p0=(-1, -0.001, 0.65) )
fig2ax4.plot(time_transient_kyr, exf(time_transient_kyr, popt[0], popt[1], popt[2]), color='.5', linewidth=2)
print( "4:", np.abs(1/popt[1]) )
popt,pcov = curve_fit( exf, time_transient_kyr[1:], Q_s_out_5[1:]/Q_s_in_5, p0=(1, -0.001, 0.65) )
fig2ax5.plot(time_transient_kyr, exf(time_transient_kyr, popt[0], popt[1], popt[2]), color='.5', linewidth=2)
print( "5:", np.abs(1/popt[1]) )

fig.set_tight_layout(True)
fig2.set_tight_layout(True)

#fig.savefig('/home/awickert/Dropbox/Papers/InProgress/LongProfileGravelBedRivers/figures/inProgress/transient_5panel_raw.svg')
#fig.savefig('/home/awickert/Dropbox/Papers/InProgress/LongProfileGravelBedRivers/figures/inProgress/transient_3panel_raw.png')
fig.savefig('transient_5panel_raw.png')
fig.savefig('transient_5panel_raw.pdf')
#fig2.savefig('/home/awickert/Dropbox/Papers/InProgress/LongProfileGravelBedRivers/figures/inProgress/transient_5panel_to_equilibrium_raw.svg')
#fig2.savefig('/home/awickert/Dropbox/Papers/InProgress/LongProfileGravelBedRivers/figures/inProgress/transient_3panel_to_equilibrium_raw.png')
fig2.savefig('transient_5panel_to_equilibrium.png')
fig2.savefig('transient_5panel_to_equilibrium.pdf')

