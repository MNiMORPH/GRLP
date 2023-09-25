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

# Shorter for test
# Custom for just this test network
x = [
      1000 * np.array([2, 4, 6, 8]),
      1000 * np.array([10, 12, 14, 16]),
    ]

# Base level
x_bl = 1000*18
z_bl = 0

# Upstream boundary condition: 1.5% grade
S0 = [0.015]

nseg = len(x)
numel = []
for _x in x:
    numel.append(len(_x))

upstream_segment_IDs = [[], [0]]
downstream_segment_IDs = [[1], []]

z = []
#Q_in_list = [5., 5., 10., 5, 15.]
# Test constant 
#Q_in_list = [5., 5., 10., 5, 15.]
# HM! DOUBLE Q AND GET AN UNEXPECTEDLY LARGE DROP IN SLOPE.
# Expect 0.5**(6/7.) = 0.44545
# Instead, get 0.3572
# This could be the hint behind the rest of the model misfit
Q_in_list = [5., 10.]
# Let's try a linear ramp. Same issue?
# Probing into the discrepancy at the junction
# Looks good!
Q_in_list = [np.arange(6,10), np.arange(10,14)]
# Is the problem more just due to discontinuities in water dsicharge?
Q_in_list = [np.array([2,3,8,9]), np.arange(10,14)]
# Looks good!
# So: discontinuity AT tributary junction is the problem <-- NEXT STEP
# Let's try moving the discontinuity
# Upstream: Looks good.
Q_in_list = [np.array([5,5,5,10]), np.array([10,10,10,10])]
# Downstream: Looks good.
Q_in_list = [np.array([5,5,5,5]), np.array([5,10,10,10])]
# At junction: Looks visually good.
Q_in_list = [np.array([5,5,5,5]), np.array([10,10,10,10])]
# Hm -- let's check slopes, then.
# Back to the upstream example.
Q_in_list = [np.array([5,5,10,10]), np.array([10,10,10,10])]

Q_in_list = [5., 5.]


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

# ITERATIONS CAUSE SEGMENTS DOWNSTREAM OF CONFLUENCE TO *DIVERGE FROM*
# AND BECOME GENTLER IN SLOPE THAN THEY SHOULD BE VIA THEORY
net.set_niter(1)
net.get_z_lengths()

# For testing
#net.evolve_threshold_width_river_network(nt=100, dt=dt*10)

#"""
# For plotting
# WHEN RUN FOR NT=10, GET BACKWARDS SLOPE ON TRIBUTARY
# THIS IS WHERE WE NEED TO ADD IN CLOSED BASINS AS ANOTHER SEGMENT TYPE
net.evolve_threshold_width_river_network(nt=36, dt=100*dt)

# Predict slopes -- without tributary-network inputs
S_predicted = []
S0val = S0[0]
Q0 = Q_in_list[0]
for lp in net.list_of_LongProfile_objects:
    Mean_Q = (lp.Q[:-1] + lp.Q[1:])/2.
    S_predicted.append( S0val * (Q0/Mean_Q)**(6/7.) )

# Print slope calc
print ( "Slopes:" )
_iter = 0
for lp in net.list_of_LongProfile_objects:
    print ("Measured", np.diff(lp.z)/np.diff(lp.x) )
    print ("Predicted", S_predicted[_iter])
    _iter += 1


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
plt.tight_layout()
plt.show()
#"""

