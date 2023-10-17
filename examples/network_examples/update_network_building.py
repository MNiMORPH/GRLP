from grlp import *

# ---- River properties
x0 = 10.e3
L = 100.e3
mean_Q = 10.
mean_Qs = 0.001
B = 98.1202038813591

net, net_topo = generate_random_network(10, L, B, mean_Q, mean_Qs, evolve=True)

# Plot
for lp in net.list_of_LongProfile_objects:
    # If not downstream-most segment
    if len( lp.downstream_segment_IDs ) > 0:
        for _id in lp.downstream_segment_IDs:
            dsseg = net.list_of_LongProfile_objects[_id]
            _xjoin = [lp.x[-1], dsseg.x[0]]
            _zjoin = [lp.z[-1], dsseg.z[0]]
            plt.plot(_xjoin, _zjoin, 'k-', linewidth=4, alpha=.5)
    else:
        plt.plot(lp.x_ext[0][-2:], lp.z_ext[0][-2:], 'k-', linewidth=4, alpha=.5)
    plt.plot(lp.x, lp.z, '-', linewidth=4, alpha=.5)#, label=lp.)

plt.show()
