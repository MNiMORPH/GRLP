#!/usr/bin/env python3

import numpy as np
import random
from grlp import *

def add_segment(nx_ls, down_ls, up_ls, down_ID, nx, trunk_nx_ls=None):

    ID = len(nx_ls)

    nx_ls.append( nx )

    down_ls.append([down_ID])
    up_ls[down_ID].append(ID)
    up_ls.append([])

    if trunk_nx_ls is not None:
        trunk_nx_ls.append( nx )
        return nx_ls, down_ls, up_ls, ID, trunk_nx_ls
    else:
        return nx_ls, down_ls, up_ls, ID

def downstream_IDs(down_ls, i):
    """
    Search list of downstream IDs, return all segments downstream of
    specified point.

    Uses recursive call of _downstream_IDs function.

    """
    IDs = []
    def _downstream_IDs(down_ls, i):
        IDs.append(i)
        if down_ls[i]:
            _downstream_IDs(down_ls, down_ls[i][0])
    _downstream_IDs(down_ls, i)
    return IDs

def upstream_IDs(up_ls, i):
    """
    Search list of upstrean IDs, return all segments upstream of
    specified point.

    Uses recursive call of _upstream_IDs function.

    """
    IDs = []
    def _upstream_IDs(up_ls, i):
        IDs.append(i)
        for j in up_ls[i]:
            _upstream_IDs(up_ls, j)
    _upstream_IDs(up_ls, i)
    return IDs

def add_branch(nx_ls, down_ls, up_ls, down_ID, nx_max):
    seg_nx_max = nx_max - sum([nx_ls[i] for i in downstream_IDs(down_ls, down_ID)])
    if seg_nx_max > 0:
        nx_ls, down_ls, up_ls, ID1 = add_segment(nx_ls, down_ls, up_ls, down_ID, seg_nx_max)
        nx_ls, down_ls, up_ls, ID2 = add_segment(nx_ls, down_ls, up_ls, down_ID, seg_nx_max)
        if random.randint(0,1):
            add_branch(nx_ls, down_ls, up_ls, ID1, nx_max)
            add_branch(nx_ls, down_ls, up_ls, ID2, nx_max)

def network_upstream_IDs(net, i):
    IDs = []
    def _upstream_IDs(net, i):
        for j in net.list_of_LongProfile_objects[i].upstream_segment_IDs:
            IDs.append(j)
            _upstream_IDs(net, j)
    _upstream_IDs(net, i)
    return IDs

def plot_network(net, mouth):

    """
    Create schematic representation of network planfrom.

    """

    DICT = {}

    def plot_segment(net, i, y, y_down=False):
        if net.list_of_LongProfile_objects[i].downstream_segment_IDs:
            x = np.hstack(( net.list_of_LongProfile_objects[i].x_ext[1:]/1000., net.list_of_LongProfile_objects[i].x_ext[-1]/1000. ))
            y = np.hstack(( np.full(len(net.list_of_LongProfile_objects[i].x_ext[1:]), y), y_down))
        else:
            x = net.list_of_LongProfile_objects[i].x_ext[1:]/1000.
            y = np.full(len(net.list_of_LongProfile_objects[i].x_ext[1:]/1000.), y)
        plt.plot(x, y)
        DICT[i] = {}
        DICT[i]['x'] = x
        DICT[i]['y'] = y

    def plot_all_upstream(net, mouth, y_init):
        scl = [-1, 1]
        for i,seg in enumerate(net.list_of_LongProfile_objects[mouth].upstream_segment_IDs):
            up_segs = network_upstream_IDs(net, seg)
            up_heads = [j for j in up_segs if not net.list_of_LongProfile_objects[j].upstream_segment_IDs]
            if len(up_heads) == 0:
                y = y_init + scl[i]
                plot_segment(net, seg, y, y_init)
            else:
                y = y_init + len(up_heads)*scl[i]
                plot_segment(net, seg, y, y_init)
                plot_all_upstream(net, seg, y)

    # ---- Plot mouth
    plot_segment(net, mouth, 0)

    # ---- Plot upstream
    plot_all_upstream(net, 0, 0)

    # ---- Show
    plt.show()

    return DICT


def strahler_order(net):

    check = np.full( len(net.IDs), 0 )
    strahler = np.full( len(net.IDs), 1 )
    O = 1
    while not np.array_equal(strahler, check):
        check = strahler.copy()
        for i in net.IDs:
            up_IDs = network_upstream_IDs(net, i)
            of_order = [j for j in up_IDs if strahler[j] == O]
            if len(of_order) > 1:
                strahler[i] += 1
        O += 1

    return strahler

def build_randomised_network(nx_max):

    """
    Builds lists of segment length, upstream and downstream segment IDs.

    For specified total length and randomised segment length.
    Also randomised tributary branching.

    """

    # ---- Initialise lists
    nx_list = [random.choice(np.arange(10,nx_max-2))]
    trunk_nx_list = [nx_list[0]]
    downstream_segment_list = [[]]
    upstream_segment_list = [[]]
    down_trunk_ID = 0

    # ---- Generate network
    while (nx_max - sum(trunk_nx_list)) > 0:

        seg_nx_max = nx_max-sum(trunk_nx_list)

        poss_nx = np.arange(2, seg_nx_max-2)
        if len(poss_nx) > 1:
            # nx = poss_nx[random.randint(0,len(poss_nx)-1)]
            nx = random.choice(poss_nx)
        else:
            nx = seg_nx_max

        # add trunk and tributary segment
        nx_list, downstream_segment_list, upstream_segment_list, trunk_ID, trunk_nx_list = add_segment(nx_list, downstream_segment_list, upstream_segment_list, down_trunk_ID, nx, trunk_nx_ls=trunk_nx_list)
        nx_list, downstream_segment_list, upstream_segment_list, trib_ID = add_segment(nx_list, downstream_segment_list, upstream_segment_list, down_trunk_ID, nx)
        if sum(trunk_nx_list) == nx_max:
            break

        # add branches to tributary
        if random.randint(0,1):
            add_branch(nx_list, downstream_segment_list, upstream_segment_list, trib_ID, nx_max)

        # update trunk ID
        down_trunk_ID = trunk_ID

    return nx_list, upstream_segment_list, downstream_segment_list

def set_up_network_object(nx_list, nx_max, dx, upstream_segment_list, downstream_segment_list, Q_max, Qs_max, evolve=False):

    """
    Uses lists of segment length, upstream and downstream segment IDs to build
    instance of grlp.Network.

    As wells as configuration lists requires network total discharge and
    sediment supply.

    Optionally evolve for a while (aiming for steady-state).

    """

    # ---- Set up segment LongProfile objects
    segments = []
    heads = [i for i in range(len(nx_list)) if not upstream_segment_list[i]]
    Q_in = Q_max / len(heads)
    Qs_in = Qs_max / len(heads)
    for i,nx in enumerate(nx_list):

        down_IDs = downstream_IDs(downstream_segment_list, i)[1:]
        down_nx = sum([nx_list[i] for i in down_IDs])
        x0 = nx_max - down_nx - nx_list[i]
        x1 = x0 + nx_list[i]
        x = np.arange( (x0-1)*dx, (x1+1)*dx, dx )

        lp = LongProfile()
        lp.set_ID(i)
        lp.set_upstream_segment_IDs(upstream_segment_list[i])
        lp.set_downstream_segment_IDs(downstream_segment_list[i])
        lp.set_intermittency(1)
        lp.basic_constants()
        lp.bedload_lumped_constants()
        lp.set_hydrologic_constants()
        lp.set_x(x_ext=x)

        S0 = (Qs_max / (lp.k_Qs * Q_max))**(6./7.)
        lp.set_z(S0=-S0, z1=0.)
        lp.set_niter()
        lp.set_uplift_rate(0)
        lp.set_B(150.)

        if i in heads:
            lp.set_Q(Q=Q_in)
        else:
            up_IDs = upstream_IDs(upstream_segment_list, i)
            up_heads = [j for j in up_IDs if not upstream_segment_list[j]]
            lp.set_Q(Q=Q_in*len(up_heads))

        if i in heads:
            lp.set_Qs_input_upstream(Qs_in)

        # if lp.downstream_segment_IDs:
            lp.z += segments[lp.downstream_segment_IDs[0]].z_ext[0]
            lp.z_ext += segments[lp.downstream_segment_IDs[0]].z_ext[0]

        if not lp.downstream_segment_IDs:
            lp.set_z_bl(0.)

        segments.append(lp)

    net = Network(segments)
    net.get_z_lengths()
    net.set_niter()
    net.build_ID_list()
    if evolve:
        net.evolve_threshold_width_river_network(nt=1000, dt=3.15e10)

    return net

def build_standardised_network(nx_max, nx_seg, nx_trib=None):
    """
    Builds lists of segment length, upstream and downstream segment IDs.

    For specified total length and segment length. 
    """

    # ---- If tributary length not specified, set to trunk segment length
    if not nx_trib:
        nx_trib = nx_seg

    # ---- Initialise lists
    if nx_seg % 2 == 0:
        nx_list = [nx_seg+1]
    else:
        nx_list = [nx_seg+2]
    trunk_nx_list = [nx_list[0]]
    downstream_segment_list = [[]]
    upstream_segment_list = [[]]
    down_trunk_ID = 0

    # ---- Generate network
    while (nx_max - sum(trunk_nx_list)) > 0:

        nx = min(nx_seg, nx_max-sum(trunk_nx_list))
        nx_tr = min(nx_trib, nx_max-sum(trunk_nx_list))

        # add trunk and tributary segment
        nx_list, downstream_segment_list, upstream_segment_list, trunk_ID, trunk_nx_list = add_segment(nx_list, downstream_segment_list, upstream_segment_list, down_trunk_ID, nx, trunk_nx_ls=trunk_nx_list)
        nx_list, downstream_segment_list, upstream_segment_list, trib_ID = add_segment(nx_list, downstream_segment_list, upstream_segment_list, down_trunk_ID, nx_tr)
        if sum(trunk_nx_list) == nx_max:
            break

        # # add branches to tributary
        # if random.randint(0,1):
        #     add_branch(nx_list, downstream_segment_list, upstream_segment_list, trib_ID, nx_max)

        # update trunk ID
        down_trunk_ID = trunk_ID

    return nx_list, upstream_segment_list, downstream_segment_list

if __name__ == "__main__":

    nx_max = 101
    dx = 1.e3
    Q_max = 16.
    Qs_max = 0.001163
    # nx_list, up_list, down_list = build_randomised_network(nx_max)
    nx_list, up_list, down_list = build_standardised_network(nx_max, 20, 2)
    net = set_up_network_object(nx_list, nx_max, dx, up_list, down_list, Q_max, Qs_max)

    visualisation = plot_network(net, 0)
    print("Strahler order is: %d" % (strahler_order(net).max()))