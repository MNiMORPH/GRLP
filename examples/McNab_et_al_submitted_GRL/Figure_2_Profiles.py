from grlp import *
from copy import deepcopy
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap

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
A = 0.2
for i,p in enumerate(num_periods):
    print(i)
    time = np.arange(0., 6.*p, p/1000.)
    scale = 1. + A*np.sin(2.*np.pi*time/p)
    z = np.zeros((len(time), len(lp.x)))
    Qs = np.zeros((len(time), len(lp.x)))
    lp_cp = deepcopy(lp)
    for j,s in enumerate(scale):
        lp_cp.set_Qs_input_upstream(Qs0 * s)
        lp_cp.evolve_threshold_width_river(nt=1, dt=p/1000.)
        lp_cp.compute_Q_s()
        z[j,:] = lp_cp.z.copy()
        Qs[j,:] = lp_cp.Q_s.copy()
    num_gain[i,:] = (z[4000:,:].max(axis=0) - z[4000:,:].min(axis=0)) / \
        2. / (6./7.) / A / ((Qs0/(lp.k_Qs*Qw))**(6./7.)) / (lp.L-lp.x)
    zs.append(z)


# ---- Plot

grays = cm.get_cmap("gray_r", 100)

ts = np.linspace(3000,3250,6,dtype=int)

fig, axs = plt.subplots(3, 3, sharey="row", figsize=(15,9))

for i,p in enumerate(num_periods):
    
    time = np.arange(0., 4.*p, p/1000.)
    axs[0,i].plot(time/3.15e10, scale[:4000], c="gray")
    axs[0,i].scatter(time[ts]/3.15e10, scale[ts], c=ts/len(ts), cmap="gray_r", edgecolors="k")
    
    for j,t in enumerate(ts):
        axs[1,i].plot(lp.x/1000., zs[i][t,:], c=grays(j/len(ts)))
    axs[1,i].plot(lp.x/1000., zs[i][0,:], "--", color="r")

    for j,t in enumerate(ts):
        axs[2,i].plot(lp.x/1000., zs[i][t,:] - zs[i][0,:], c=grays(j/len(ts)))
    axs[2,i].plot([0.,lp.x.max()/1000.], [0,0], "--", color="r")

axs[0,0].set_ylabel(r"$Q_{s,0}~/~\bar{Q_{s,0}}$")
axs[1,0].set_ylabel(r"Elevation, $z$ [m]")
axs[2,0].set_ylabel(r"${\delta}z$ [m]")
for i,ax in enumerate(axs[0,:]):
    ax.set_xlabel("Time [kyr]")
    title = r"$P$ = %.1f x $T_{eq}$" % (num_periods[i]/3.15e12)
    ax.set_title(title)
for ax in axs[1,:]:
    ax.set_xlabel("Downstream distance [km]")
for ax in axs[2,:]:
    ax.set_xlabel("Downstream distance [km]")
fig.tight_layout()
plt.show()
