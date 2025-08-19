#! /usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import importlib
#plt.ion()
plt.ioff()

import grlp

#%load_ext autoreload
#%autoreload

importlib.reload(grlp)
#del net

dt = 3.15E7*10
_B = 100 # uniform

nseg = 5
numel = [5, 7, 4, 8, 10]

upstream_segment_IDs = [[], [], [0,1], [], [2,3]]
downstream_segment_IDs = [[2], [2], [4], [4], []]

z = []
#Q_in_list = [5., 5., 10., 5, 15.]
# Test constant 
Q_in_list = [5., 5., 10., 10., 20.]
#Q_in_list = [5., 5., 5., 5, 5.]
Q_in_list = [[2.5, 2.5, 2.5, 5, 5], [7.5, 7.5, 7.5, 7.5, 5, 5, 5], 10, 10., 20.]
Q_in_list = [[2.5, 2.5, 2.5, 5, 5],
             [7.5, 7.5, 7.5, 7.5, 5, 5, 5],
             [10, 11, 13, 14],
             [10., 11, 13, 14, 17, 21, 22, 25],
             [39., 40., 41., 42., 43., 44., 45., 45., 46., 47.]]
Q = []
B = []
print( "" )
print( "**************" )
print( "" )
for i in range(nseg):
    # Start z as all zeros
    z.append( np.zeros( numel[i] ) )
    # Start Q as constant within each segment
    Q.append( Q_in_list[i] * np.ones( numel[i] ) )
    # Uniform valley width: Simpler for example
    B.append( _B * np.ones( numel[i] ) )
    print( numel[i] )
print( "" )
print( "**************" )
print( "" )

# Custom for just this test network
x = [
      1000 * np.array([2, 4, 6.5, 9, 10]),
      1000 * np.array([0, 1, 2, 3, 6, 8, 10.5]),
      1000 * np.array([12, 15, 18, 20]),
      1000 * np.array([2, 6, 8, 12, 14, 16, 18, 20]),
      1000 * np.array([23, 24, 25, 26, 27, 28, 29, 30, 31, 32])
    ]

"""
# Uniform dx test
# Custom for just this test network
x = [
      1000 * np.array([2, 4, 6, 8, 10]),
      1000 * np.array([-2, 0, 2, 4, 6, 8, 10]),
      1000 * np.array([12, 14, 16, 18]),
      1000 * np.array([4, 6, 8, 10, 12, 14, 16, 18]),
      1000 * np.array([20, 22, 24, 26, 28, 30])
    ]
"""

# Base level
x_bl = 1000*37
z_bl = 0

# Upstream boundary condition: 1.5% grade
S0 = [0.015, 0.015, 0.015]

# Instantiate network object
net = grlp.Network()
self = net # For debugging, etc.

# Initialize network by passing x, z, etc. information to model
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

# Should do this above
net.set_niter(1)
net.get_z_lengths()

# For testing
#net.evolve_threshold_width_river_network(nt=100, dt=dt*10)

#"""
# For plotting
# WHEN RUN FOR NT=10, GET BACKWARDS SLOPE ON TRIBUTARY
# THIS IS WHERE WE NEED TO ADD IN CLOSED BASINS AS ANOTHER SEGMENT TYPE
net.evolve_threshold_width_river_network(nt=100, dt=100*dt)
lp = net.list_of_LongProfile_objects[2]

lp.z[1] += 5
lp.z_ext[0][2] += 20
lp.z_ext[1][2] += 20

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

nsteps = 50
for i in range(nsteps):
    net.evolve_threshold_width_river_network(nt=1, dt=1*dt)
    if i % np.ceil(nsteps/10) == 0:
        for lp in net.list_of_LongProfile_objects:
            # If not downstream-most segment
            if len( lp.downstream_segment_IDs ) > 0:
                for _id in lp.downstream_segment_IDs:
                    dsseg = net.list_of_LongProfile_objects[_id]
                    _xjoin = [lp.x[-1], dsseg.x[0]]
                    _zjoin = [lp.z[-1], dsseg.z[0]]
                    plt.plot(_xjoin, _zjoin, 'k-', linewidth=1, alpha=1)
            else:
                plt.plot(lp.x_ext[0][-2:], lp.z_ext[0][-2:], 'k-', linewidth=1, alpha=1)
            plt.plot(lp.x, lp.z, '-', linewidth=1, alpha=1)#, label=lp.)
    if i < (nsteps-1):
        lp.z[2] += 20
        lp.z_ext[0][3] += 30
        lp.z_ext[1][3] += 30

for lp in net.list_of_LongProfile_objects:
    # If not downstream-most segment
    if len( lp.downstream_segment_IDs ) > 0:
        for _id in lp.downstream_segment_IDs:
            dsseg = net.list_of_LongProfile_objects[_id]
            _xjoin = [lp.x[-1], dsseg.x[0]]
            _zjoin = [lp.z[-1], dsseg.z[0]]
            plt.plot(_xjoin, _zjoin, 'k-', linewidth=2, alpha=.5)
    else:
        plt.plot(lp.x_ext[0][-2:], lp.z_ext[0][-2:], 'k-', linewidth=2, alpha=.5)
    plt.plot(lp.x, lp.z, 'k-', linewidth=2, alpha=.5)#, label=lp.)

plt.tight_layout()
plt.show()
#"""

