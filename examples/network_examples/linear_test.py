from grlp import *
from copy import deepcopy

# ---- Set up simple object to get k_Qs later
lp = LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()


# ---- River properties
L = 100.e3
mean_Qw = 10.
mean_Qs = 0.001
B = 98.1202038813591


# ---- Network Object
x = np.arange(0, L, 1.e3)
S0 = (mean_Qs/(lp.k_Qs * mean_Qw))**(6./7.)
net = Network()
net.initialize(
    config_file = None,
    x_bl = L,
    z_bl = 0.,
    S0 = [(mean_Qs/(lp.k_Qs * mean_Qw))**(6./7.)],
    upstream_segment_IDs = [[]],
    downstream_segment_IDs = [[]],
    x = [x],
    z = [(L-x)*S0],
    Q = [np.full(len(x), mean_Qw)],
    B = [np.full(len(x), B)],
    overwrite = False
    )
net.set_niter(3)
net.get_z_lengths()
net.evolve_threshold_width_river_network(nt=1000, dt=3.15e10)
net.list_of_LongProfile_objects[0].compute_equilibration_time()


# ---- Evolve periodic
periods = np.logspace(-2.,2.,7) * net.list_of_LongProfile_objects[0].equilibration_time
gain_Qs = np.zeros((len(periods),len(net.list_of_LongProfile_objects[0].x)))

for i,period in enumerate(periods):
    print(i)
    
    neti = deepcopy(net)
    
    time, dt = np.linspace(0., period*4, 4000, retstep=True)
    
    A_Qs = 0.2
    scale = 1. + A_Qs*np.sin(2.*np.pi*time/period)
    
    z = np.zeros((len(time), len(neti.list_of_LongProfile_objects[0].x)))
    Qs = np.zeros((len(time), len(neti.list_of_LongProfile_objects[0].x)))
    
    for j,s in enumerate(scale):
        S0_j =  S0 * (s**(6./7.))
        neti.update_z_ext_external_upstream( S0 = [S0_j] )
        neti.evolve_threshold_width_river_network(nt=1, dt=dt)
        neti.list_of_LongProfile_objects[0].compute_Q_s()
        z[j,:] = neti.list_of_LongProfile_objects[0].z
        Qs[j,:] = neti.list_of_LongProfile_objects[0].Q_s

    gain_Qs[i,:] = (
        (Qs[2000:,:].max(axis=0) - Qs[2000:,:].min(axis=0)) /
        (2. * A_Qs * mean_Qs)
        )


# ---- Linear
lin_periods = np.logspace(-2.5,2.5,51) * net.list_of_LongProfile_objects[0].equilibration_time
lin_gain_Qs = np.zeros((len(lin_periods),len(net.list_of_LongProfile_objects[0].x)))
for i,period in enumerate(lin_periods):
    lin_gain_Qs[i,:] = net.list_of_LongProfile_objects[0].compute_Qs_gain(period, A_Qs=0.2)


# ---- Plot
Teq = net.list_of_LongProfile_objects[0].equilibration_time
plt.plot(lin_periods/Teq, lin_gain_Qs[:,-1])
plt.plot(periods/Teq, gain_Qs[:,-1], 'o')
plt.xscale("log")
plt.xlabel(r"$P$ / $T_{eq}$")
plt.ylabel(r"$G_{Q_s}$")
plt.show()