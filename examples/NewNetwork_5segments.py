import numpy as np
from matplotlib import pyplot as plt
import importlib
#plt.ion()
plt.ioff()

import grlp
importlib.reload(grlp)

dt = 3.15E7*10
_B = 100 # uniform

nseg = 5
numel = [5, 7, 4, 8, 6]

upstream_segment_IDs = [[], [], [0,1], [], [2,3]]
downstream_segment_IDs = [[2], [2], [4], [4], []]

z = []
Q_in_list = [5., 10., 15., 10., 25.]
Q = []
B = []
for i in range(nseg):
    # Start z as all zeros
    z.append( np.zeros( numel[i] ) )
    # Start Q as constant within each segment
    Q.append( Q_in_list[i] * np.ones( numel[i] ) )
    # Uniform valley width: Simpler for example
    B.append( _B * np.ones( numel[i] ) )

# Custom for just this test network
x = [
      1000 * np.array([2, 4, 6.5, 9, 10]),
      1000 * np.array([0, 1, 2, 3, 6, 8, 10.5]),
      1000 * np.array([12, 15, 18, 20]),
      1000 * np.array([2, 6, 8, 12, 14, 16, 18, 20]),
      1000 * np.array([23, 24, 27, 29, 29.5, 30])
    ]

# Base level
x_bl = 1000*32
z_bl = 0

# Upstream boundary condition: 1.5% grade
S0 = 0.015

net = grlp.Network()

net.initialize(
                config_file = None,
                x_bl = x_bl,
                z_bl = z_bl,
                S0 = S0,
                Q_s_0 = None,
                upstream_segment_IDs = upstream_segment_IDs,
                downstream_segment_IDs = downstream_segment_IDs,
                x = x,
                z = z,
                Q = Q,
                B = B,
                overwrite=False
                )

i = 0
for lp in segments:
    lp.set_ID(i)
    lp.set_upstream_segment_IDs(upstream_segment_list[i])
    lp.set_downstream_segment_IDs(downstream_segment_list[i])
    lp.set_intermittency(1)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    # Local or global x -- really just care about diffs for solver
    dx=5000.
    if i > 0:
        nx=6
    else:
        nx=3
    x0=0
    _x = np.arange(x0, x0+dx*nx, dx)
    #_x[-3] += 1E-6
    x_ext = np.hstack((_x[0]-dx, _x, _x[-1]+dx))
    lp.set_x(x_ext=x_ext)
    lp.set_z(S0=-S0, z1=z1)
    lp.set_niter()
    #lp.set_z_bl(z1)
    lp.set_Q(Qlist[i])
    lp.set_B(100.)
    # START HERE
    if len(upstream_segment_list[i]) == 0:
        Qs0 = lp.k_Qs * lp.Q[0] * S0**(7/6.)
        lp.set_Qs_input_upstream(Qs0)
    i += 1
    lp.set_uplift_rate(0)
    

for i in range(len(segments))[::-1]:
    lp = segments[i]
    if len(lp.downstream_segment_IDs) > 0:
        downseg = segments[lp.downstream_segment_IDs[0]]
        lp.z += downseg.z_ext[0] - 7.5
        lp.z_ext += downseg.z_ext[0] - 7.5

def i_empty(_list):
    out = []
    for i in range(len(_list)):
        if len(_list)[i] == 0:
            out.append(i)
    return out

# Works only for 1 downstream stream in each place
for segment in segments[::-1]:
    if len(segment.downstream_segment_IDs) > 0:
        x0 = segment.x[-1]
        downID = segment.downstream_segment_IDs[0]
        segment.x += segments[downID].x_ext[0] - x0
        segment.x_ext += segments[downID].x_ext[0] - x0
for segment in segments:
    segment.x += 7000
    segment.x_ext += 7000

#segments[1].Q[0] = 10.

i = 0
for lp in segments:
    if len(downstream_segment_list[i]) == 0:
        lp.set_z_bl(lp.z_ext[-1])
    i += 1

i = 0
for lp in segments:
    if len(downstream_segment_list[i]) == 0:
        lp.set_z_bl(lp.z_ext[-1])
    i += 1

# External sed input lists
for lp in segments:
    lp.cr_t_incision = []
    lp.cr_z = []
    lp.cr_dz_incision = []

for lp in segments:
    lp.z0 = lp.z.copy()

# Instantiate network object
net = grlp.Network(segments)
net.get_z_lengths()
net.set_niter()
net.build_ID_list()

# More sediment input upstream -- this doesn't account for downstream fining
# or for increases in water input -- but for faking it, isn't too bad of a 
# start.
net.list_of_LongProfile_objects[0].set_Qs_input_upstream(
                                 1.5 * net.list_of_LongProfile_objects[0].Q_s_0)
net.list_of_LongProfile_objects[1].set_Qs_input_upstream(
                                 1 * net.list_of_LongProfile_objects[0].Q_s_0)
net.list_of_LongProfile_objects[3].set_Qs_input_upstream(
                                 1.5 * net.list_of_LongProfile_objects[0].Q_s_0)
net.evolve_threshold_width_river_network(nt=100, dt=dt*100)

_imovie = 0
plt.figure()
for lp in net.list_of_LongProfile_objects:
    for _id in lp.downstream_segment_IDs:
        dsseg = net.list_of_LongProfile_objects[_id]
        _xjoin = [lp.x[-1], dsseg.x[0]]
        _zjoin = [lp.z[-1], dsseg.z[0]]
        plt.plot(_xjoin, _zjoin, 'k-', linewidth=4, alpha=.5)
    if len(lp.downstream_segment_IDs) == 0:
        plt.plot(lp.x_ext[1:], lp.z_ext[1:], '-', linewidth=4, alpha=.5)#, label=lp.)
    else:
        plt.plot(lp.x, lp.z, '-', linewidth=4, alpha=.5)#, label=lp.)
plt.ylim(-500,1500)
#plt.xlim(0,11000)
plt.tight_layout()
plt.text(4000, 1200, "Start.", fontsize=26, fontweight='bold')
plt.draw()
#plt.savefig('/home/awickert/Downloads/RiverNetMovie/RiverNet4seg' + '%03d' %_imovie + '.png')
plt.savefig('/home/andy/Desktop/RiverNetMovie/RiverNet5seg' + '%03d' %_imovie + '.png')
_imovie += 1
#plt.pause(.1)
#plt.legend()

def set_external_sediment_input_parameters( _lp=lp ):
    # Ignore aggradational effects; just simply toy-model test
    # Or maybe aggradation as a simple switch.
    dz_incision = (_lp.z_old - _lp.z)
    dz_incision[dz_incision < 0] = 0.
    _lp.cr_dz_incision.append( dz_incision )
    _lp.cr_t_incision.append( _lp.t )
    _lp.cr_z.append( _lp.z.copy() )

def compute_dzdt_external_sediment_input(_lp = lp, A0=0.001,
                                            decay_time=5*3.15E7):
    # Simple exp decay
    # Consider making A0 a f(Q_s), locally
    input_aggradation_m_per_s = 0.
    if A0 is None:
        #_lp.compute_Q_s()
        A0 = lp.Q_s_0/2.
    for i in range(len(_lp.cr_dz_incision)):
        # Reduce sediment inputs if the river has since aggraded
        dzmod = (_lp.cr_z[i] + _lp.cr_dz_incision[i] - _lp.z)
        dzmod[dzmod < 0] = 0. # If z is higher, nothing happens
        # If incised past this point, still just use prior value
        # (later incision handled later)
        dzmod[dzmod > _lp.cr_dz_incision[i]] = \
            _lp.cr_dz_incision[i][dzmod > _lp.cr_dz_incision[i]]
        # All other values stay the same... no need for another line.
        # Next, let's compute the added sediment from each step of incision
        # That has not yet been re-buried
        # ASSUMING ONE DT PER TIME THAT THIS IS CALLED!!!!!!!!!
        # IN THIS LINE, RIGHT BELOW
        dzdt_incision_i = dzmod / _lp.dt
        _dx_tmp = (_lp.dx_ext[:-1] + _lp.dx_ext[1:])/2.
        Q_s_in_extra = A0 * dzmod \
                        * np.exp( -(_lp.t - _lp.cr_t_incision[i]) / decay_time )
        input_aggradation_m_per_s += Q_s_in_extra / _lp.B / _dx_tmp
    return input_aggradation_m_per_s
    
#lp.set_uplift_rate(input_aggradation_m_per_s)

dt /= 2

net.list_of_LongProfile_objects[-1].set_z_bl(-300)
for ti in range(30):
    for ti_inner in range(20):
        for lp in net.list_of_LongProfile_objects:
            lp.z_old = lp.z.copy()
        net.evolve_threshold_width_river_network(nt=1, dt=dt)
        for lp in net.list_of_LongProfile_objects:
            if lp.ID == 0:
                set_external_sediment_input_parameters(lp)
                lp.dzdt_sed = compute_dzdt_external_sediment_input(lp, A0=0.0001, decay_time=150*3.15E7)#, A0 = 0.003)
                lp.set_source_sink_distributed( -lp.dzdt_sed )
    #print(lp.dzdt_sed)
    print(net.list_of_LongProfile_objects[3].z[0])
    #plt.figure()
    plt.cla()
    for lp in net.list_of_LongProfile_objects:
        z_old = lp.z.copy()
        for _id in lp.downstream_segment_IDs:
            dsseg = net.list_of_LongProfile_objects[_id]
            _xjoin = [lp.x[-1], dsseg.x[0]]
            _zjoin = [lp.z[-1], dsseg.z[0]]
            plt.plot(_xjoin, _zjoin, 'k-', linewidth=4, alpha=.5)
        if len(lp.downstream_segment_IDs) == 0:
            plt.plot(lp.x_ext[1:], lp.z_ext[1:], '-', linewidth=4, alpha=.5)#, label=lp.)
        else:
            plt.plot(lp.x, lp.z, '-', linewidth=4, alpha=.5)#, label=lp.)
    plt.ylim(-500,1500)
    #plt.xlim(0,11000)
    plt.text(-10000, 1100, "Base-level fall: 300 m\n$\Delta t =$ 5 years", fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.draw()
    plt.savefig('/home/andy/Desktop/RiverNetMovie/RiverNet5seg' + '%03d' %_imovie + '.png')
    _imovie += 1
    #plt.pause(0.00001)
    #plt.legend()

# Egregious copy/paste

net.list_of_LongProfile_objects[0].set_Qs_input_upstream(
                                 2 * net.list_of_LongProfile_objects[0].Q_s_0)
for ti in range(60):
    for ti_inner in range(20):
        for lp in net.list_of_LongProfile_objects:
            lp.z_old = lp.z.copy()
        net.evolve_threshold_width_river_network(nt=1, dt=dt)
        for lp in net.list_of_LongProfile_objects:
            if lp.ID == 0:
                set_external_sediment_input_parameters(lp)
                lp.dzdt_sed = compute_dzdt_external_sediment_input(lp, A0=0.0001, decay_time=150*3.15E7)#, A0 = 0.003)
                lp.set_source_sink_distributed( -lp.dzdt_sed )
    #print(lp.dzdt_sed)
    print(net.list_of_LongProfile_objects[3].z[0])
    #plt.figure()
    plt.cla()
    for lp in net.list_of_LongProfile_objects:
        z_old = lp.z.copy()
        for _id in lp.downstream_segment_IDs:
            dsseg = net.list_of_LongProfile_objects[_id]
            _xjoin = [lp.x[-1], dsseg.x[0]]
            _zjoin = [lp.z[-1], dsseg.z[0]]
            plt.plot(_xjoin, _zjoin, 'k-', linewidth=4, alpha=.5)
        if len(lp.downstream_segment_IDs) == 0:
            plt.plot(lp.x_ext[1:], lp.z_ext[1:], '-', linewidth=4, alpha=.5)#, label=lp.)
        else:
            plt.plot(lp.x, lp.z, '-', linewidth=4, alpha=.5)#, label=lp.)
    plt.ylim(-500,1500)
    #plt.xlim(0,11000)
    #plt.text(4000, 150, "Base-level fall: 20 m\nBlue branch still 2x sed\n$\Delta t =$ 2 years", fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.text(-24000, 1100, "Doubling sediment supply:\nblue branch\n$\Delta t =$ 5 years", fontsize=16, fontweight='bold')
    plt.draw()
    plt.savefig('/home/andy/Desktop/RiverNetMovie/RiverNet5seg' + '%03d' %_imovie + '.png')
    _imovie += 1
    #plt.pause(0.00001)
    #plt.legend()


for lp in net.list_of_LongProfile_objects:
    lp.compute_Q_s()
    print(lp.Q_s)

