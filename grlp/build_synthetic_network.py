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
    poss_nx = np.arange(2,seg_nx_max-2)
    if len(poss_nx) > 1:
        nx = (random.choice(poss_nx), random.choice(poss_nx))
    else:
        nx = (seg_nx_max, seg_nx_max)

    if seg_nx_max > 0:
        nx_ls, down_ls, up_ls, ID1 = add_segment(nx_ls, down_ls, up_ls, down_ID, nx[0])
        nx_ls, down_ls, up_ls, ID2 = add_segment(nx_ls, down_ls, up_ls, down_ID, nx[1])
        if random.choice([0,1]):
            add_branch(nx_ls, down_ls, up_ls, ID1, nx_max)
            add_branch(nx_ls, down_ls, up_ls, ID2, nx_max)


def plot_network(net, mouth, show=True):

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
            up_segs = net.find_upstream_IDs(seg)
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
    if show:
        plt.show()
    else:
        plt.close()

    return DICT


def build_randomised_network(nx_max):

    """
    Builds lists of segment length, upstream and downstream segment IDs.

    For specified total length and randomised segment length.
    Also randomised tributary branching.

    """

    # ---- Initialise lists
    nx_list = [random.choice(np.arange(2,nx_max-2))]
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


def set_up_network_object(nx_list, dx, upstream_segment_list, downstream_segment_list, Q_max, Qs_max, evolve=False):
    """
    Uses lists of segment length, upstream and downstream segment IDs to build
    instance of grlp.Network.
    As well as configuration lists requires network total discharge and
    sediment supply.
    Optionally evolve for a while (aiming for steady-state).
    """

    # ---- Some parameters for use during set up
    segments = []
    sources = [i for i in range(len(nx_list)) if not upstream_segment_list[i]]
    Q_in = Q_max / len(sources)
    Qs_in = Qs_max / len(sources)
    x_min = 0

    # ---- Loop over segments setting up LongProfile objects
    for i,nx in enumerate(nx_list):

        # Some basic set up
        lp = LongProfile()
        lp.set_ID(i)
        lp.set_upstream_segment_IDs(upstream_segment_list[i])
        lp.set_downstream_segment_IDs(downstream_segment_list[i])
        lp.set_intermittency(1)
        lp.basic_constants()
        lp.bedload_lumped_constants()
        lp.set_hydrologic_constants()
        lp.set_niter()
        lp.set_uplift_rate(0)

        # Set up x domain
        down_IDs = downstream_IDs(downstream_segment_list, i)[1:]
        down_nx = sum([nx_list[i] for i in down_IDs])
        x0 = - down_nx - nx_list[i]
        x1 = x0 + nx_list[i]
        x = np.arange( (x0-1)*dx, (x1+1)*dx, dx )
        x_min = min(x_min, min(x)+dx)
        lp.set_x(x_ext=x)

        # set width
        lp.set_B(150.)

        # Set initial z
        S0 = (Qs_max / (lp.k_Qs * Q_max))**(6./7.)
        lp.set_z(S0=-S0, z1=0.)

        if i in sources:
            # if segment is a source, set Q and Qs to input values
            lp.set_Q(Q=Q_in)
            lp.set_Qs_input_upstream(Qs_in)
        else:
            # otherwise set Q based on number of upstream sources
            # Qs will be set by input from upstream segments
            up_IDs = upstream_IDs(upstream_segment_list, i)
            num_sources = len([j for j in up_IDs if not upstream_segment_list[j]])
            lp.set_Q(Q=Q_in*num_sources)

        if lp.downstream_segment_IDs:
            # if not mouth, reset downstream elevation to that of downstream segment
            lp.z += segments[lp.downstream_segment_IDs[0]].z_ext[0]
            lp.z_ext += segments[lp.downstream_segment_IDs[0]].z_ext[0]
        else:
            # otherwise set network base level to zero
            lp.set_z_bl(0.)

        # add LongProfile object to segment list
        segments.append(lp)
 
    # ---- Update x coordinates to run from 0 at furthest upstream point
    for i,seg in enumerate(segments):
        segments[i].set_x(x_ext=segments[i].x_ext-x_min)

    # ---- Initialise and set up network object with list of LongProfile objects
    net = Network(segments)
    net.get_z_lengths()
    net.set_niter()
    net.build_ID_list()

    # ---- If requested evolve network, aiming for steady state
    if evolve:
        net.evolve_threshold_width_river_network(nt=1000, dt=3.15e10)

    return net

class Simple_Network:
    """
    Set up simple network with specified total length and segment length.
    Single-segment tributaries added between each trunk segment.
    Optionally specify different length for tributaries.
    """


    def __init__(self, nx_total, nx_trunk_seg, nx_trib_seg=None):
        self.nx_total = nx_total
        self.nx_trunk_seg = nx_trunk_seg
        if nx_trib_seg:
            self.nx_trib_seg = nx_trib_seg
        else:
            self.nx_trib_seg = nx_trunk_seg

        self.upstream_segment_IDs = None
        self.downstream_segment_IDs = None
        self.nxs = None

        self.build_network()

    def add_segment(self, down_ID, trunk=False):
        """
        Add segment to network.
        """

        # ID number for new segment
        ID = len(self.nxs)

        # Check whether trunk or tributary, assign length
        if trunk:
            nx = self.nx_trunk_seg
        else:
            nx = self.nx_trib_seg

        # Make sure length won't leave floating point at end
        # Can't make a segment with only one point...
        # Instead add one to final segment
        if self.nx_total - sum(self.trunk_nxs) <= nx + 1:
            nx = self.nx_total - sum(self.trunk_nxs)

        # Update nx lists
        self.nxs.append(nx)
        if trunk:
            self.trunk_nxs.append(nx)

        # Update topology lists
        self.downstream_segment_IDs.append([down_ID])
        self.upstream_segment_IDs[down_ID].append(ID)
        self.upstream_segment_IDs.append([])

        # Return ID for future use
        return ID

    def build_network(self):
        """
        Generate network topology.
        """

        # initialise lists
        self.upstream_segment_IDs = [[]]
        self.downstream_segment_IDs = [[]]
        self.nxs = [self.nx_trunk_seg]
        self.trunk_nxs = []
        trunk_ID = 0

        # add segments until total length is reached
        while self.nx_total - sum(self.trunk_nxs) > 0:
            __ = self.add_segment(trunk_ID)
            trunk_ID = self.add_segment(trunk_ID, trunk=True)

