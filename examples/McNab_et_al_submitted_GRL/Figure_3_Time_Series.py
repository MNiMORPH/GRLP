from grlp import *
from copy import deepcopy
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.signal import find_peaks


# ---- River properties
L = 100.e3
Qw = 10.
B = 98.1202038813591
dx = 1.e3
Qs0 = 0.001


# ---- Set up long profile object
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


# ---- Numerical
num_periods = np.array([0.1, 1., 10.]) * 3.15e12
num_gain = np.zeros((len(num_periods), len(lp.x)))
zs = []
Qss = []
zs_Qw = []
Qss_Qw = []
A = 0.2
for i,p in enumerate(num_periods):
    print(i)
    time = np.arange(0., 6.*p, p/1000.)
    scale = 1. + A*np.sin(2.*np.pi*time/p)
    z = np.zeros((len(time), len(lp.x)))
    Qs = np.zeros((len(time), len(lp.x)))
    z_Qw = np.zeros((len(time), len(lp.x)))
    Qs_Qw = np.zeros((len(time), len(lp.x)))
    lp_cp = deepcopy(lp)
    lp_cp_Qw = deepcopy(lp)
    for j,s in enumerate(scale):
        lp_cp.set_Qs_input_upstream(Qs0 * s)
        lp_cp.evolve_threshold_width_river(nt=1, dt=p/1000.)
        lp_cp.compute_Q_s()
        z[j,:] = lp_cp.z.copy()
        Qs[j,:] = lp_cp.Q_s.copy()
        lp_cp_Qw.set_Q(Qw * s)
        lp_cp_Qw.evolve_threshold_width_river(nt=1, dt=p/1000.)
        lp_cp_Qw.compute_Q_s()
        z_Qw[j,:] = lp_cp_Qw.z.copy()
        Qs_Qw[j,:] = lp_cp_Qw.Q_s.copy()
    zs.append(z)
    Qss.append(Qs)
    zs_Qw.append(z_Qw)
    Qss_Qw.append(Qs_Qw)
    
# ---- Plotting code
def plot_time_series(zs, Qss, scale, periods, Qs0, varying, force_label):
    xs = [0, 20, 40, 60, 80]
    fig, axs = plt.subplots(2,3,sharey="row", sharex="col", figsize=(15,8), gridspec_kw={'height_ratios': [1, 2.5]})
    for i,p in enumerate(num_periods):
        time = np.hstack(( -p, np.arange(0., 6.*p, p/1000.) )) / 3.15e10
        scl = np.hstack(( 1, scale ))
        scl_peaks = np.hstack(( find_peaks(scl)[0], find_peaks(-scl)[0] ))
        Qs_out = np.hstack(( 1., Qss[i][:,-1]/Qs0 ))
        Qs_peaks = np.hstack(( find_peaks(Qs_out)[0], find_peaks(-Qs_out)[0][1:] ))
        for peak in scl_peaks:
            axs[0,i].plot([time[peak]]*2, [0., scl[peak]], ":", c="0.7")
            axs[1,i].plot([time[peak]]*2, [0., 750.], ":", c="0.7")
        axs[0,i].plot()
        axs[0,i].plot(time, scl, c="0.7")
        axs[0,i].scatter(time[scl_peaks], scl[scl_peaks], s=12, c="0.7")
        axs[0,i].plot(time, Qs_out, "--", c="0")
        axs[0,i].scatter(time[Qs_peaks], Qs_out[Qs_peaks], s=12, c="0")
        axs[0,i].set_xlim(-p/3.15e10/4., 3.25*p/3.15e10)
        axs[0,i].set_ylim(0.7, 1.3)
        for x in xs:
            z = np.hstack(( zs[i][0,x], zs[i][:,x] ))
            z_peaks = np.hstack(( find_peaks(z)[0], find_peaks(-z)[0] ))
            axs[1,i].plot(time, z)
            axs[1,i].scatter(time[z_peaks], z[z_peaks], s=12)
        axs[1,i].set_ylim(0.,750.)
    axs[0,0].set_ylabel(force_label + r"$Q_{s,L}/\overline{Q_{s,0}}$")
    axs[1,0].set_ylabel(r"Elevation, $z$ [m]")
    for ax in axs[1,:]:
        ax.set_xlabel("Time [kyr]")
    for i,ax in enumerate(axs[0,:]):
        title = r"$P$ = %.1f x $T_{eq}$" % (num_periods[i]/3.15e12)
        ax.set_title(title)
    fig.suptitle("Varying " + varying + " supply")
    plt.show()

# ---- Plot
plot_time_series(zs, Qss, scale, num_periods, Qs0, "sediment", r"$Q_{s,0}/\overline{Q_{s,0}}$, ")
plot_time_series(zs_Qw, Qss_Qw, scale, num_periods, Qs0, "water", r"$Q_{w}/\overline{Q_{w}}$, ")
