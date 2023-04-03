import random
from grlp import *


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


def plot_network(net, show=True):

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

    def plot_all_upstream(net, i, y_init):
        scl = [-1, 1]
        up_IDs = net.list_of_LongProfile_objects[i].upstream_segment_IDs
        up_IDs_up_sources = [len(net.find_upstream_IDs(ID)) for ID in up_IDs]
        up_IDs = np.array(up_IDs)[np.array(up_IDs_up_sources).argsort()][::-1]
        for i,ID in enumerate(up_IDs):
            up_segs = net.find_upstream_IDs(ID)
            up_heads = [j for j in up_segs if not net.list_of_LongProfile_objects[j].upstream_segment_IDs]
            if len(up_heads) == 0:
                y = y_init + scl[i]
                plot_segment(net, ID, y, y_init)
            else:
                y = y_init + len(up_heads)*scl[i]
                plot_segment(net, ID, y, y_init)
                plot_all_upstream(net, ID, y)

    # ---- Find mouth
    mouth = [seg.ID for seg in net.list_of_LongProfile_objects if not seg.downstream_segment_IDs][0]

    # ---- Plot mouth
    plot_segment(net, mouth, 0)

    # ---- Plot upstream
    plot_all_upstream(net, 0, 0)

    # ---- Show
    if show:
        plt.yticks(labels=[], ticks=[])
        plt.xlabel("Downstream distance [km]")
        plt.show()
    else:
        plt.close()

    return DICT


def set_up_network_object(
    nx_list, dx, upstream_segment_list, downstream_segment_list, Q_in, Qs_in, B, evolve=False):
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
        x = np.arange( (x0-1), (x1+1), 1. ) * dx
        lp.set_x(x_ext=x)

        # set width
        lp.set_B(B)

        # Set initial z
        S0 = (Qs_in / (lp.k_Qs * Q_in))**(6./7.)
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
    x_min = min([min(lp.x_ext)+dx for lp in segments])
    for i,seg in enumerate(segments):
        segments[i].set_x(x_ext=segments[i].x_ext-x_min)

    # ---- Initialise and set up network object with list of LongProfile objects
    net = Network(segments)
    net.get_z_lengths()
    net.set_niter()
    net.build_ID_list()
    net.compute_network_properties()

    # ---- If requested evolve network, aiming for steady state
    if evolve:
        net.evolve_threshold_width_river_network(nt=1000, dt=3.15e10)
        for seg in net.list_of_LongProfile_objects: seg.compute_Q_s()

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


class Shreve_Random_Network:
    """
    Create random network of specified magnitude.
    Follows algorithm described in Shreve (1974, Water Resource Res.).
    Creates lists of upstream/downstream segment IDs for us in GRLP.
    """

    def __init__(self, magnitude):
        self.magnitude = magnitude
        self.links = None
        self.upstream_segment_IDs = None
        self.downstream_segment_IDs = None
        self.build_network_topology()
        self.build_lists()

    @property
    def probability_external(self):
        """
        Calculate liklihood that next link is external.
        """
        try:
            p = (self.N_internal - self.N_external)
            p *= (self.magnitude - self.N_external)
            p /= (self.N_internal - self.N_external + 1.)
            p /= (2.*self.magnitude - 2. - self.N_internal - self.N_external)
        except ZeroDivisionError:
            p = np.inf
        return p

    @property
    def N_external(self):
        """
        Number of external (i.e. source) links.
        """
        return sum(self.links)

    @property
    def N_internal(self):
        """
        Number of internal links.
        """
        return len(self.links) - sum(self.links)

    def build_network_topology(self):
        """
        Build network topology as binary string.
        Ones represent external links, zeros internal links.
        Algorithm from Shreve (1974, Water Resource Res.).
        """

        # Initialise some properties
        self.links = []
        k = 1
        last = 1

        # Start looping
        while True:

            # Generate random number between 0 and 1.
            # If greater than liklihood of external link, add internal link.
            if random.random() > self.probability_external:
                self.links.append(0)
                if last:
                    k += 1
                else:
                    last = 0

            # If less than liklihood of external link, add external link.
            else:
                self.links.append(1)
                if not last:
                    last = 1
                if last:
                    if k > 1:
                        k -= 1
                    else:
                        break
        
        # Remove extra initial link
        self.links.remove(0)

    def build_lists(self):
        """
        Build lists for use in GRLP.
        """

        # Initialise topology lists
        self.upstream_segment_IDs = [[]]
        self.downstream_segment_IDs = [[]]

        # Initialise some other useful properties
        down_segs = []
        down_seg = 0
        seg = 1

        # Loop through links
        for i,l in enumerate(self.links):

            # Update lists
            self.upstream_segment_IDs[down_seg].append(seg)
            self.upstream_segment_IDs.append([])
            self.downstream_segment_IDs.append([down_seg])
            
            # If more to come, continue through network
            if i < len(self.links)-1:
            
                # If internal link, continue upstream
                if not l:
                    down_segs.append(down_seg)
                    down_seg = seg

                # If external link, work back downstream to last free internal link
                else:
                    while len(self.upstream_segment_IDs[down_seg]) > 1:
                        down_seg = down_segs.pop(-1)
                seg += 1
