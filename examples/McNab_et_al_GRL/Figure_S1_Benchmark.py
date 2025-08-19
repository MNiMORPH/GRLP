from grlp import *
from copy import deepcopy
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.signal import find_peaks

def find_lag_times(val, time, scale, threshold=0., can_lead=False, period=False, full=False):


    peak_lags = np.zeros( len(val[0,:]) )
    trough_lags = np.zeros( len(val[0,:]) )

    tps = []

    scl_peaks, __ = find_peaks(scale)
    scl_troughs, __ = find_peaks(-scale)
    scl_tps = np.sort( np.hstack(( scl_peaks, scl_troughs )) )

    for i in range(len(val[0,:])):

        if ( val[:,i].max() - val[:,i].min() ) < threshold:
            peak_lags[i] = np.nan
            trough_lags[i] = np.nan
            continue

        obs_peaks, __ = find_peaks(val[:,i])
        obs_troughs, __ = find_peaks(-val[:,i])
        obs_tps = np.sort( np.hstack(( obs_peaks, obs_troughs )) )
        if not can_lead:
            obs_tps = obs_tps[ np.where( obs_tps >= scl_tps[0] ) ]

        obs_tps_attached = np.zeros( len(obs_tps), dtype=int )
        obs_lag_time = np.zeros( len(obs_tps), dtype=int )

        for j,tp in enumerate(obs_tps):
            if j > len(scl_tps)-1:
                continue
            obs_tps_attached[j] = scl_tps[j]
            obs_lag_time[j] = time[tp] - time[scl_tps[j]]

        peak_lags_i = []
        trough_lags_i = []

        for k,tp in enumerate(obs_tps_attached):

            if obs_lag_time[k] != 0.:
                if any(scl_peaks == tp):
                    peak_lags_i.append( obs_lag_time[k].copy() )
                else:
                    trough_lags_i.append( obs_lag_time[k].copy() )

        if len(peak_lags_i) > 1:
            peak_lags[i] = np.array(peak_lags_i[1:]).mean()
        else:
            peak_lags[i] = np.nan

        if len(trough_lags_i) > 1:
            trough_lags[i] = np.array(trough_lags_i[1:]).mean()
        else:
            trough_lags[i] = np.nan

        tps.append( obs_tps )

    if period:
        for i in range(1,len(val[0,:])):
            if (peak_lags[i] - peak_lags[i-1]) > 0.5*period:
                peak_lags[i] -= period
            if (trough_lags[i] - trough_lags[i-1]) > 0.5*period:
                trough_lags[i] -= period
        for i in range(len(val[0,:])-2,0,-1):
            if (peak_lags[i] - peak_lags[i+1]) > 0.5*period:
                peak_lags[i] -= period
            if (trough_lags[i] - trough_lags[i+1]) > 0.5*period:
                trough_lags[i] -= period

    if full:
        return {'plags': peak_lags, 'tlags': trough_lags, 'obs_tps': tps, 'scl_tps': scl_tps, 'scl_p': scl_peaks, 'scl_t': scl_troughs}
    else:
        return (peak_lags + trough_lags)/2.


# ---- River properties
L = 100.e3
Qw = 10.
B = 98.1202038813591
dx = 1.e3
Qs0 = 0.001

# ---- Set up long profile object
lp = LongProfile()
lp = LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()
lp.set_x(x_ext=np.arange(0., L+dx, dx))
lp.set_Q(Q=Qw)
lp.set_B(B=B)
lp.set_z(S0=(Qs0/(lp.k_Qs*Qw))**(6./7.))
lp.set_Qs_input_upstream(Qs0)
lp.set_z_bl(0.)
lp.set_uplift_rate(0.)
lp.set_niter()
lp.evolve_threshold_width_river(nt=10000, dt=3.15e9)
lp.compute_Q_s()
lp.compute_equilibration_time()

# ---- Record properties
lin_periods = np.logspace(-2., 2., 81) * 3.15e12
lin_gain = np.zeros((len(lin_periods), len(lp.x)))
lin_lag = np.zeros((len(lin_periods), len(lp.x)))
lin_gain_Qs = np.zeros((len(lin_periods), len(lp.x)))
lin_lag_Qs = np.zeros((len(lin_periods), len(lp.x)))
lin_gain_Qs_Qw = np.zeros((len(lin_periods), len(lp.x)))
lin_lag_Qs_Qw = np.zeros((len(lin_periods), len(lp.x)))
for i,p in enumerate(lin_periods):
    lin_gain[i,:] = lp.compute_z_gain(p)
    lin_lag[i,:] = lp.compute_z_lag(p) / p
    lin_gain_Qs[i,:] = lp.compute_Qs_gain(p, A_Qs=0.2)
    lin_lag_Qs[i,:] = lp.compute_Qs_lag(p, A_Qs=0.2) / p
    lin_gain_Qs_Qw[i,:] = lp.compute_Qs_gain(p, A_Q=0.2)
    lin_lag_Qs_Qw[i,:] = lp.compute_Qs_lag(p, A_Q=0.2) / p


# ---- Benchmark
num_periods = np.logspace(-1.5, 1.5, 7) * 3.15e12
num_gain = np.zeros((len(num_periods), len(lp.x)))
num_gain_Qs = np.zeros((len(num_periods), len(lp.x)))
num_lag = np.zeros((len(num_periods), len(lp.x)))
num_lag_Qs = np.zeros((len(num_periods), len(lp.x)))
num_gain_Qs_Qw = np.zeros((len(num_periods), len(lp.x)))
num_lag_Qs_Qw = np.zeros((len(num_periods), len(lp.x)))
A = 0.2
for i,p in enumerate(num_periods):
    # if i>0: break
    print(i)
    time = np.arange(0., 6.*p, p/1000.)
    scale = 1. + A*np.sin(2.*np.pi*time/p)
    z = np.zeros((len(time), len(lp.x)))
    Qs = np.zeros((len(time), len(lp.x)))
    Qs_Qw = np.zeros((len(time), len(lp.x)))
    lp_cp = deepcopy(lp)
    lp_cp_Qw = deepcopy(lp)
    for j,s in enumerate(scale):
        lp_cp.set_Qs_input_upstream(Qs0 * s)
        # lp_cp.set_Q(Qw * s)
        lp_cp.evolve_threshold_width_river(nt=1, dt=p/1000.)
        lp_cp.compute_Q_s()
        z[j,:] = lp_cp.z.copy()
        Qs[j,:] = lp_cp.Q_s.copy()
        
        lp_cp_Qw.set_Q(Qw * s)
        lp_cp_Qw.evolve_threshold_width_river(nt=1, dt=p/1000.)
        lp_cp_Qw.compute_Q_s()
        Qs_Qw[j,:] = lp_cp_Qw.Q_s.copy()

    num_gain[i,:] = (z[4000:,:].max(axis=0) - z[4000:,:].min(axis=0)) / \
        2. / A / ((Qs0/(lp.k_Qs*Qw))**(6./7.)) / (lp.L-lp.x)
    num_lag[i,:] = find_lag_times(z, time, scale, threshold=0.5) / p
    num_gain_Qs[i,:] = (Qs[4000:,:].max(axis=0) - Qs[4000:,:].min(axis=0)) /  (2. * A * Qs0)
    num_lag_Qs[i,:] = find_lag_times(Qs, time, scale, threshold=2.e-5) / p
    num_gain_Qs_Qw[i,:] = (Qs_Qw[4000:,:].max(axis=0) - Qs_Qw[4000:,:].min(axis=0)) /  (2. * A * Qs0)
    num_lag_Qs_Qw[i,:] = find_lag_times(Qs_Qw, time, scale, threshold=5.e-5, can_lead=True) / p



# ---- Plot
fig, axs = plt.subplots(2,3, sharex=True, sharey="row", figsize=(15,10))
xs = np.hstack(( 0, np.arange(9, 90, 10), 98 ))
cm = plt.cm.get_cmap("viridis_r")
cols = cm(lp.x[xs]/lp.L)
for i,x in enumerate(xs):
    if i==10:
        lab1 = "Analytical"
        lab2 = "Numerical"
    else:
        lab1, lab2 = (None,None)
    axs[0,0].plot(lin_periods/lp.equilibration_time, lin_gain[:,x], color=cols[i], label=lab1)
    axs[0,0].scatter(num_periods/lp.equilibration_time, num_gain[:,x], c=[cols[i]]*len(num_periods), label=lab2)
    axs[1,0].plot(lin_periods/lp.equilibration_time, lin_lag[:,x], color=cols[i], label=lab1)
    axs[1,0].scatter(num_periods/lp.equilibration_time, num_lag[:,x], c=[cols[i]]*len(num_periods), label=lab2)

    axs[0,1].plot(lin_periods/lp.equilibration_time, lin_gain_Qs[:,x], color=cols[i], label=lab1)
    axs[0,1].scatter(num_periods/lp.equilibration_time, num_gain_Qs[:,x], c=[cols[i]]*len(num_periods), label=lab2)
    axs[1,1].plot(lin_periods/lp.equilibration_time, lin_lag_Qs[:,x], color=cols[i], label=lab1)
    axs[1,1].scatter(num_periods/lp.equilibration_time, num_lag_Qs[:,x], c=[cols[i]]*len(num_periods), label=lab2)

    axs[0,2].plot(lin_periods/lp.equilibration_time, lin_gain_Qs_Qw[:,x], color=cols[i], label=lab1)
    axs[0,2].scatter(num_periods/lp.equilibration_time, num_gain_Qs_Qw[:,x], c=[cols[i]]*len(num_periods), label=lab2)
    axs[1,2].plot(lin_periods/lp.equilibration_time, lin_lag_Qs_Qw[:,x], color=cols[i], label=lab1)
    axs[1,2].scatter(num_periods/lp.equilibration_time, num_lag_Qs_Qw[:,x], c=[cols[i]]*len(num_periods), label=lab2)

for ax in axs[1,:]:
    ax.set_xlabel(r"$P~/~T_{eq}$")
axs[0,0].set_ylabel(r"$G_z$")
axs[1,0].set_ylabel(r"$\varphi_z~/~P$")
axs[0,0].set_title(r"${\delta}z$")
axs[0,1].set_title(r"${\delta}Q_s$: Varying $Q_{s,0}$")
axs[0,2].set_title(r"${\delta}Q_s$: Varying $Q_w$")
axs[0,0].legend()
# plt.colorbar(cmap[0], label=r"$x~/~L$", ax=axs[0])
for ax in axs:
    for a in ax:
        a.set_xscale("log")
        a.set_aspect(1./a.get_data_ratio())
plt.show()