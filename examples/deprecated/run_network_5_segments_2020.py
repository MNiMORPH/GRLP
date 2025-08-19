import numpy as np
from matplotlib import pyplot as plt
import importlib
#plt.ion()
plt.ioff()


import grlp
importlib.reload(grlp)

S0 = 0.015
P_xB = 0.2
z1 = 0+7.5 # Interesting -- my RHS always has to be 0 (not quite right)

dt = 3.15E7

nseg = 5

segments = []
for i in range(nseg):
    segments.append(grlp.LongProfile())

Qlist = [5., 10., 15., 10., 25.]
upstream_segment_list = [[], [], [0,1], [], [2,3]]
downstream_segment_list = [[2], [2], [4], [4], []]

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
    dx=500.
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
        Qs0 = lp.k_Qs * lp.Q[0] * (1*S0)**(7/6.)
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


#for lp in segments:
#    lp.Q = None
#    lp.set_Q(k_xQ=5E-4, P_xQ=lp.P_xQ)
#    #for id in segment.upstream_segment_IDs:
        
    

#segments[0].dQ[-1] = 10
#segments[1].dQ[-1] = 10
#lp.dQ[0] = Qlist[1] - Qlist[0] # BAND-AID!
# DOWNSTREAM DQ?

net = grlp.Network(segments)
net.get_z_lengths()
net.set_niter()
net.build_ID_list()
#net.set_dQ()
#segments[2].dQ[0] = 0
#segments[0].dQ[-1] = 10
#segments[1].dQ[-1] = 10
net.evolve_threshold_width_river_network(nt=100, dt=dt)

_imovie = 0
plt.figure()
for lp in segments:
    for _id in lp.downstream_segment_IDs:
        dsseg = net.list_of_LongProfile_objects[_id]
        _xjoin = [lp.x[-1], dsseg.x[0]]
        _zjoin = [lp.z[-1], dsseg.z[0]]
        plt.plot(_xjoin, _zjoin, 'k-', linewidth=4, alpha=.5)
    if len(lp.downstream_segment_IDs) == 0:
        plt.plot(lp.x_ext[1:], lp.z_ext[1:], '-', linewidth=4, alpha=.5)#, label=lp.)
    else:
        plt.plot(lp.x, lp.z, '-', linewidth=4, alpha=.5)#, label=lp.)
plt.ylim(-25,200)
plt.xlim(0,11000)
plt.tight_layout()
plt.text(4000, 150, "Start.", fontsize=26, fontweight='bold')
plt.draw()
#plt.savefig('/home/awickert/Downloads/RiverNetMovie/RiverNet4seg' + '%03d' %_imovie + '.png')
_imovie += 1
plt.pause(0.5)
#plt.legend()

net.list_of_LongProfile_objects[0].set_Qs_input_upstream(
                                 2 * net.list_of_LongProfile_objects[0].Q_s_0)

for ti in range(10):
    net.evolve_threshold_width_river_network(nt=10, dt=dt)
    #plt.figure()
    plt.cla()
    for lp in segments:
        for _id in lp.downstream_segment_IDs:
            dsseg = net.list_of_LongProfile_objects[_id]
            _xjoin = [lp.x[-1], dsseg.x[0]]
            _zjoin = [lp.z[-1], dsseg.z[0]]
            plt.plot(_xjoin, _zjoin, 'k-', linewidth=4, alpha=.5)
        if len(lp.downstream_segment_IDs) == 0:
            plt.plot(lp.x_ext[1:], lp.z_ext[1:], '-', linewidth=4, alpha=.5)#, label=lp.)
        else:
            plt.plot(lp.x, lp.z, '-', linewidth=4, alpha=.5)#, label=lp.)
    plt.ylim(-25,200)
    plt.xlim(0,11000)
    plt.tight_layout()
    plt.text(4000, 150, "Doubling sediment supply:\nblue branch\n$\Delta t =$ 10 years", fontsize=16, fontweight='bold')
    plt.draw()
    #plt.savefig('/home/awickert/Downloads/RiverNetMovie/RiverNet4seg' + '%03d' %_imovie + '.png')
    _imovie += 1
    plt.pause(0.5)
    #plt.legend()

net.list_of_LongProfile_objects[-1].set_z_bl(-20)
for ti in range(10):
    net.evolve_threshold_width_river_network(nt=2, dt=dt)
    #plt.figure()
    plt.cla()
    for lp in segments:
        for _id in lp.downstream_segment_IDs:
            dsseg = net.list_of_LongProfile_objects[_id]
            _xjoin = [lp.x[-1], dsseg.x[0]]
            _zjoin = [lp.z[-1], dsseg.z[0]]
            plt.plot(_xjoin, _zjoin, 'k-', linewidth=4, alpha=.5)
        if len(lp.downstream_segment_IDs) == 0:
            plt.plot(lp.x_ext[1:], lp.z_ext[1:], '-', linewidth=4, alpha=.5)#, label=lp.)
        else:
            plt.plot(lp.x, lp.z, '-', linewidth=4, alpha=.5)#, label=lp.)
    plt.ylim(-25,200)
    plt.xlim(0,11000)
    plt.text(4000, 150, "Base-level fall: 20 m\nBlue branch still 2x sed\n$\Delta t =$ 2 years", fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.draw()
    #plt.savefig('/home/awickert/Downloads/RiverNetMovie/RiverNet4seg' + '%03d' %_imovie + '.png')
    _imovie += 1
    plt.pause(0.5)
    #plt.legend()


for lp in segments:
    lp.compute_Q_s()
    print(lp.Q_s)

