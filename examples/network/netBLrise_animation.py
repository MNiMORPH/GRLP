#! /usr/bin/python3
"""
Animate a river network's response to base-level rise.

A five-segment synthetic network is brought to equilibrium, its base level is
raised by 30 m, and the transient aggradation is animated with
matplotlib.animation.FuncAnimation. The animation is written to a movie file
(netBLrise.gif) in the working directory -- no numbered PNG frames and no
hard-coded output paths. To save an .mp4 instead, use writer='ffmpeg'.
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

import grlp

# ---------------------------------------------------------------------------
# Five-segment synthetic network
# ---------------------------------------------------------------------------
dt = 3.15E7 * 10          # 10 yr, in seconds
_B = 100                  # uniform valley width [m]

nseg = 5
numel = [5, 7, 4, 8, 10]

upstream_segment_IDs   = [[], [], [0, 1], [], [2, 3]]
downstream_segment_IDs = [[2], [2], [4], [4], []]

# Discharge within each segment [m3/s]
Q_in_list = [[2.5, 2.5, 2.5, 5, 5],
             [7.5, 7.5, 7.5, 7.5, 5, 5, 5],
             [10, 11, 13, 14],
             [10., 11, 13, 14, 17, 21, 22, 25],
             [39., 40., 41., 42., 43., 44., 45., 45., 46., 47.]]

# Downvalley node positions of each segment [m]
x = [
      1000 * np.array([2, 4, 6.5, 9, 10]),
      1000 * np.array([0, 1, 2, 3, 6, 8, 10.5]),
      1000 * np.array([12, 15, 18, 20]),
      1000 * np.array([2, 6, 8, 12, 14, 16, 18, 20]),
      1000 * np.array([23, 24, 25, 26, 27, 28, 29, 30, 31, 32])
    ]

z = [np.zeros(numel[i]) for i in range(nseg)]
Q = [np.array(Q_in_list[i], dtype=float) for i in range(nseg)]
B = [_B * np.ones(numel[i]) for i in range(nseg)]

x_bl = 1000 * 33
z_bl = 0
S0 = [0.015, 0.015, 0.015]      # 1.5% grade at the three channel heads

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
                overwrite = False,
                )
net.set_niter(1)
net.get_z_lengths()

# Bring the network to equilibrium, then record the pre-rise profile
net.evolve_threshold_width_river_network(nt=100, dt=100*dt)
for lp in net.list_of_LongProfile_objects:
    lp.x_init = lp.x.copy()
    lp.z_init = lp.z.copy()
z_bl_init = z_bl                 # base level before the rise

# Raise base level by 30 m; the network aggrades in response
net.set_z_bl(z_bl + 30)

# ---------------------------------------------------------------------------
# Animation
# ---------------------------------------------------------------------------
segs = net.list_of_LongProfile_objects

def draw_profiles(ax, use_init, **style):
    """Plot every segment and its downstream connector. use_init selects the
    saved pre-rise reference (x_init, z_init) or the current state; the outlet
    connector runs from the last node to the downstream ghost (base level)."""
    for lp in segs:
        xx = lp.x_init if use_init else lp.x
        zz = lp.z_init if use_init else lp.z
        if len(lp.downstream_segment_IDs) > 0:
            for _id in lp.downstream_segment_IDs:
                ds = segs[_id]
                dsx = ds.x_init if use_init else ds.x
                dsz = ds.z_init if use_init else ds.z
                ax.plot([xx[-1], dsx[0]], [zz[-1], dsz[0]], **style)
        else:
            zbl = z_bl_init if use_init else lp.z_bl
            ax.plot([xx[-1], lp.x_ghost_downstream], [zz[-1], zbl], **style)
        ax.plot(xx, zz, **style)

steps_per_frame = 5
nframes = 120

fig, ax = plt.subplots(figsize=(9, 5))

def update(frame):
    net.evolve_threshold_width_river_network(nt=steps_per_frame, dt=dt)
    ax.cla()
    draw_profiles(ax, use_init=True, color='k', linewidth=2)             # pre-rise
    draw_profiles(ax, use_init=False, color='C0', linewidth=4, alpha=.5)  # current
    ax.set_ylim(-50, 450)
    ax.set_xlabel('Distance downvalley in network [m]', fontsize=14)
    ax.set_ylabel('Elevation above pre-rise outlet [m]', fontsize=14)
    years = round(frame * steps_per_frame * dt / 3.15E7)
    ax.text(0.78, 0.85, '%05d years' % years, transform=ax.transAxes,
            family='monospace')
    fig.tight_layout()

anim = animation.FuncAnimation(fig, update, frames=nframes, interval=50)
anim.save('netBLrise.gif', writer='pillow', fps=20)
print('Wrote netBLrise.gif')
