from grlp import *
import scipy.stats as sts

# ------------------------- CONTINUOUS CASE ---------------------------------- #

def compute_power_law_coefficient(mean, p, L, x0):
    """
    Calculate power law coefficient, k, needed so that mean(y = k*(x+x0)**p)
    over range x=0 to x=L is equal to given value.
    """
    return mean * L * (p+1.) / ((L+x0)**(p+1.) - x0**(p+1.))

# ---- River properties
x0 = 10.e3 # Length used to set upstream area at inlet (x=0)
L = 100.e3
mean_Qw = 10.
mean_Qs = 0.001
B = 98.1202038813591
p = 2. # inverse hack exponent
dx = 5.e2

# initial set up
lp = LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

# x
x = np.arange(0,L,dx)

# Qw
k_x_Qw = compute_power_law_coefficient(mean_Qw, p, L, x0)
Qw = k_x_Qw*((x+x0)**p)

# Qs
k_x_Qs = compute_power_law_coefficient(mean_Qs, p, L, x0)
ssd = p * k_x_Qs * (x+x0)**(p-1.) / (B * (1. - lp.lambda_p))

# z
S0=(mean_Qs/(lp.k_Qs * mean_Qw))**(6./7.)

net = Network()
net.initialize(
    config_file = None,
    x_bl = L,
    z_bl = 0.,
    S0 = [S0],
    upstream_segment_IDs = [[]],
    downstream_segment_IDs = [[]],
    x = [x],
    z = [(L-x)*S0],
    Q = [Qw],
    B = [np.full(len(x),B)],
    overwrite = False
    )
net.set_niter(3)
net.get_z_lengths()
net.list_of_LongProfile_objects[0].set_source_sink_distributed(ssd)
net.evolve_threshold_width_river_network(nt=100, dt=3.15e11)
for seg in net.list_of_LongProfile_objects:
    seg.compute_Q_s()


plt.plot(net.list_of_LongProfile_objects[0].x/1.e3, net.list_of_LongProfile_objects[0].z)
plt.plot(net.list_of_LongProfile_objects[0].x/1.e3, S0*(L - net.list_of_LongProfile_objects[0].x), ":")
plt.xlabel("Downstream distance [km]")
plt.ylabel("Elevation [m]")
plt.show()

plt.plot(net.list_of_LongProfile_objects[0].x/1.e3, net.list_of_LongProfile_objects[0].S/S0)
plt.xlabel("Downstream distance [km]")
plt.ylabel("Actual/Intended Slope")
plt.show()

# ---------------------------------------------------------------------------- #

# ---------------------------- NETWORK CASE ---------------------------------- #

# ---- Some network properties
effective_rainfall = 1.e3*0.4/3.154e10 # precipitation * runoff coefficient
B = 98.1202038813591 # width
sediment_discharge_ratio = 1.e4 # Qw/Qs ratio to hold fixed throughout network

# ---- Segment lengths and areas
mean_segment_length = 10.e3
mean_segment_length_area_ratio = 300. # mean ratio of segment length to segment area
mean_supply_area = mean_segment_length * mean_segment_length_area_ratio / 2. # mean supply area for channel heads
segment_length = sts.gamma(2., scale=mean_segment_length/2.) # object to draw random segment lengths from - gamma distribution after Shreve

# constant supply areas and area/length ratios
segment_length_area_ratio = mean_segment_length_area_ratio
supply_area = mean_supply_area

# # random supply areas and area/length ratios
# segment_length_area_ratio = sts.norm(loc=mean_segment_length_area_ratio, scale=mean_segment_length_area_ratio/10.) # object to draw area:length ratios from
# supply_area = sts.norm(loc=mean_supply_area, scale=mean_supply_area/10.) # object to draw channel head supply areas from

# ---- Generate random network topology, segment lengths and areas
net_topo = Shreve_Random_Network(
    magnitude=10, 
    segment_length=segment_length,
    segment_length_area_ratio=segment_length_area_ratio,
    supply_area=supply_area
    )

# Define nxs, dxs
min_nxs = 5 # minimum number of points per segment
approx_dx = 5.e2 # approximate dx to aim for
dxs = []
nxs = []
for i,L in enumerate(net_topo.segment_lengths):
    nxs.append(max(min_nxs, int(L/approx_dx))) # nxs set to max of min_nxs and length/approx_dx
    dxs.append(L / nxs[-1])

# ---- Go from areas to discharges
# Discharges to add to channel heads
supply_discharges = [
    area*effective_rainfall for area in net_topo.source_areas
    ]
# Discharges to add along segments
internal_discharges = [
    area*effective_rainfall for area in net_topo.segment_areas
    ]

# ---- Some parameters for use during set up
sources = [i for i in range(len(nxs)) if not net_topo.upstream_segment_IDs[i]]

# ---- Basic lp object to get k_Qs for later
lp = LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

# ---- Loop over segments filling lists for network
x_ls = []
z_ls = []
Q_ls = []
Qs_ssd_ls = []
B_ls = []
for i,nx in enumerate(nxs):
    
    # Set up x domain
    down_IDs = downstream_IDs(net_topo.downstream_segment_IDs, i)[1:]
    down_x = sum([net_topo.segment_lengths[j] for j in down_IDs])
    x0 = - down_x - net_topo.segment_lengths[i]
    x1 = x0 + net_topo.segment_lengths[i]
    x = x0 + np.arange( 0, nx, 1 ) * dxs[i]
    x_ls.append(x)
    
    # set width
    B_ls.append(B)
    
    ############################# DEFINING Q & Qs HERE #########################
    # Set discharges
    if i in sources:
        # For channel heads, min discharge is supply from source area, max
        # discharge is min + amount supplied along stream
        min_discharge = supply_discharges[i]
        max_discharge = supply_discharges[i] + internal_discharges[i]
    else:
        # For internal segments, min is supply from upstream, max is min +
        # amount supplied along stream
        up_IDs = upstream_IDs(net_topo.upstream_segment_IDs, i)
        max_discharge = sum([supply_discharges[j] + internal_discharges[j] for j in up_IDs])
        min_discharge = max_discharge - internal_discharges[i]
    # Q(x) given by linear interpolation between min and max
    Q, dQ = np.linspace(min_discharge, max_discharge, len(x)+1, retstep=True)
    # print(dQ)
    # print(dQ, internal_discharges[i]/net_topo.segment_lengths[i]/effective_rainfall)
    Q_ls.append(Q[:-1])
    # SSD function of dQ/dx
    # Qs_ssd = dQ/sediment_discharge_ratio/dxs[i]/B_ls[i]/(1.-lp.lambda_p)
    Qs_ssd = (max_discharge-min_discharge)/sediment_discharge_ratio/net_topo.segment_lengths[i]/B_ls[i]/(1.-lp.lambda_p)
    Qs_ssd = np.full(len(x), Qs_ssd)    
    Qs_ssd_ls.append( Qs_ssd )
    ############################################################################

    # Set initial z
    # S0 = (1./(lp.k_Qs*sediment_discharge_ratio))**(6./7.)
    S0 = ((Q[0]/sediment_discharge_ratio)/(lp.k_Qs*Q[0]))**(6./7.)
    z = (x.max()-x)*S0
    
    # if not mouth, reset downstream elevation to that of downstream segment
    # needs lists to work backwards from downstream end
    if net_topo.downstream_segment_IDs[i]:
        z += z_ls[net_topo.downstream_segment_IDs[i][0]][0] + dxs[i]*S0
    z_ls.append(z)
    
# ---- Update x coordinates to run from 0 at furthest upstream point, record max
x_min = min([min(x) for x in x_ls])
for i,nx in enumerate(nxs):
    x_ls[i] -= x_min
# AW guess: take max x value and add an increment dx for base-level boundary
x_max = max([max(x) + (x[-1]-x[-2]) for x in x_ls])
# Then add on some vertical distance to make a straight line to base level
dz_for_bl = dxs[0]*S0
for _z in z_ls:
    _z += dz_for_bl
    
# ---- Initialize network object
net = Network()
net.initialize(
    config_file = None,
    x_bl = x_max,
    z_bl = 0.,
    S0 = S0 * np.ones(len(sources)),
    upstream_segment_IDs = net_topo.upstream_segment_IDs,
    downstream_segment_IDs = net_topo.downstream_segment_IDs,
    x = x_ls,
    z = z_ls,
    Q = Q_ls,
    B = B_ls,
    overwrite = False
    )
net.set_niter(3)
net.get_z_lengths()
net.compute_network_properties()

# ---- Set segment source-sink-distributed term
for i,seg in enumerate(net.list_of_LongProfile_objects):
    seg.set_source_sink_distributed(Qs_ssd_ls[i])


# ---- Evolve, check for steady state
dt = 2*3.15e11
for i in range(100):
    net.evolve_threshold_width_river_network(nt=1, dt=dt)
    net.list_of_LongProfile_objects[0].compute_Q_s()
    plt.scatter(i*dt/3.15e10, net.list_of_LongProfile_objects[0].Q_s[-1])
plt.xlabel("Time [kyr]")
plt.ylabel(r"$Q_s$ [m$^3$ s$^{-1}$]")
plt.show()



# ---- Plot planform, Q(x), profile
fig, axs = plt.subplots(1,3,sharex=True)

# Planform
planform = plot_network(net, show=False)
for seg in planform.keys():
    axs[0].plot(planform[seg]['x'], planform[seg]['y'])
axs[0].set_xlabel("Downstream distance [km]")

# Discharge
for seg in net.list_of_LongProfile_objects:
    axs[1].plot(seg.x/1.e3, seg.Q)
axs[1].set_ylabel(r"$Q$ [m$^3$ s$^{-1}$]")
axs[1].set_xlabel("Downstream distance [km]")

# Profile
for seg in net.list_of_LongProfile_objects:
    axs[2].plot(seg.x/1.e3, seg.z)
# Linear version
L = net.list_of_LongProfile_objects[0].x_ext[0][-1]
x = np.linspace(0., L, 100)
z = S0*(L-x)
axs[2].plot(x/1.e3, z, ":")
axs[2].set_xlabel("Downstream distance [km]")
axs[2].set_ylabel("Elevation [m]")

plt.show()


# ---- Plot Q/Qs
plt.plot([0,net.list_of_LongProfile_objects[0].x[-1]/1.e3], [1,1], ":")
for i,seg in enumerate(net.list_of_LongProfile_objects):
    seg.compute_Q_s()
    if i in sources:
        plt.plot(seg.x/1.e3, seg.Q/seg.Q_s/sediment_discharge_ratio, c="black")
    else:
        plt.plot(seg.x/1.e3, seg.Q/seg.Q_s/sediment_discharge_ratio, c="grey")    
plt.xlabel("Downstream distance [km]")
plt.ylabel(r"[ Actual $Q$/$Q_s$ ] / [ Intended $Q$/$Q_s$ ]")
plt.show()

# ---------------------------------------------------------------------------- #
