import random
import copy
from scipy.optimize import minimize
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


def plot_network_deprecated(net, show=True):

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

def get_simple_network_setup_params(
    upstream_segment_list,
    downstream_segment_list,
    L,
    mean_Q,
    mean_Qs,
    min_nxs=5,
    approx_dx=5.e2):
    """
    Find spatial discretisation and input discharges for given network
    geometry and desired length, approximate discretisation, mean discharges.
    """
    
    # ---- Spatial discretisation
    
    # Number of links
    num_links = len(upstream_segment_list)
    
    # Find sources
    sources = [
        i for i,up_ids in enumerate(upstream_segment_list) if len(up_ids)==0]
    
    # Find maximum topological length,
    # i.e. number of downstream segments to outlet
    max_topo_length = max([
        len(downstream_IDs(downstream_segment_list, i))
        for i in range(num_links)])
        
    # Find length of each link so that total length equals L
    link_length = L / max_topo_length
    
    # Find how many nodes for each link, set dx
    link_n = max(min_nxs, int(link_length/approx_dx))
    nxs = [link_n for i in range(num_links)]
    dx = link_length / link_n
    
    # ---- Sediment & water discharge
    
    # Find number of sources upstream of each point
    up_sources = []
    for i in range(num_links):
        count = 0
        up_IDs = upstream_IDs(upstream_segment_list, i)
        for ID in up_IDs:
            if len(upstream_IDs(upstream_segment_list, ID)) == 1:
                count += 1
        up_sources.append(count)
        
    # Find input sediment and water discharge to give specified means
    Q_in = mean_Q / np.mean(up_sources)
    Qs_in = mean_Qs / np.mean(up_sources)
    
    # Return
    return nxs, dx, Q_in, Qs_in

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
    
    # ---- Basic lp object to get k_Qs for later
    lp = LongProfile()
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()

    # ---- Loop over segments filling lists for network
    x_ls = []
    z_ls = []
    Q_ls = []
    B_ls = []
    for i,nx in enumerate(nx_list):
        
        # Set up x domain
        down_IDs = downstream_IDs(downstream_segment_list, i)[1:]
        down_nx = sum([nx_list[j] for j in down_IDs])
        x0 = - down_nx - nx_list[i]
        x1 = x0 + nx_list[i]
        x = np.arange( x0, x1, 1. ) * dx
        x_ls.append(x)
        
        # set width
        B_ls.append(B)
        
        # Set initial z
        S0 = (Qs_in/(lp.k_Qs*Q_in))**(6./7.)
        z = (x.max()-x)*S0
        
        # if not mouth, reset downstream elevation to that of downstream segment
        # needs lists to work backwards from downstream end
        if downstream_segment_list[i]:
            z += z_ls[downstream_segment_list[i][0]][0] + dx*S0
        z_ls.append(z)
        
        # Set discharge
        if i in sources:
            # if segment is a source, set Q input values
            Q_ls.append(np.full(len(x), Q_in))
        else:
            # otherwise set based on number of upstream sources
            up_IDs = upstream_IDs(upstream_segment_list, i)
            num_sources = len([j for j in up_IDs if not upstream_segment_list[j]])
            Q_ls.append(np.full(len(x), Q_in*num_sources))
            
    # ---- Update x coordinates to run from 0 at furthest upstream point, record max
    x_min = min([min(x) for x in x_ls])
    for i,nx in enumerate(nx_list):
        x_ls[i] -= x_min
    # AW guess: take max x value and add an increment dx for base-level boundary
    x_max = max([max(x) for x in x_ls]) + dx
    # Then add on some vertical distance to make a straight line to base level
    dz_for_bl = dx*S0
    for _z in z_ls:
        _z += dz_for_bl

    # ---- Initialize network object
    net = Network()
    net.initialize(
        config_file = None,
        x_bl = x_max,
        z_bl = 0.,
        S0 = S0 * np.ones(len(nx_list)),
        upstream_segment_IDs = upstream_segment_list,
        downstream_segment_IDs = downstream_segment_list,
        x = x_ls,
        z = z_ls,
        Q = Q_ls,
        B = B_ls,
        overwrite = False
        )
    net.set_niter(3)
    net.get_z_lengths()

    # ---- If requested evolve network, aiming for steady state
    if evolve:
        net.evolve_threshold_width_river_network(nt=1000, dt=3.15e10)
    
    # ---- Compute Qs
    for seg in net.list_of_LongProfile_objects: seg.compute_Q_s()

    return net

def generate_random_network(magnitude, length, width, mean_Q, mean_Qs, evolve=False, topology=None):
    """
    Generate a random network with given magnitude, length, width, and mean
    discharges.
    """
    
    # Get random network topology
    net_topo = Shreve_Random_Network(magnitude=magnitude, topology=topology)
    
    # Get setup parameters
    nxs, dx, Q_in, Qs_in = get_simple_network_setup_params(
        net_topo.upstream_segment_IDs,
        net_topo.downstream_segment_IDs,
        length,
        mean_Q,
        mean_Qs)
        
    # Set up the object
    net = set_up_network_object(
        nxs,
        dx,
        net_topo.upstream_segment_IDs,
        net_topo.downstream_segment_IDs,
        Q_in,
        Qs_in,
        width,
        evolve)
        
    # Return
    return net, net_topo

def plot_network(net, show=True):
    """
    Generate a plotable network planform from a network object.
    """

    def check_for_segment_conflicts(ID, segs_by_topo_length, net, ys):
        """
        Check for segments that overlap with the given segment.
        """
        
        topo_length = len(net.find_downstream_IDs(ID))
        
        for nearby_ID in segs_by_topo_length[topo_length]:
            if nearby_ID != ID:
                if ys[ID] == ys[nearby_ID]:
                    return nearby_ID
                                
        for nearby_ID in segs_by_topo_length[topo_length+1]:
            if (
                nearby_ID != ID and
                net.list_of_LongProfile_objects[nearby_ID]. \
                downstream_segment_IDs):
                down_ID = (
                    net.list_of_LongProfile_objects[nearby_ID]. \
                    downstream_segment_IDs[0])
                if (
                    ys[ID] >= min(ys[nearby_ID], ys[down_ID]) and 
                    ys[ID] <= max(ys[nearby_ID], ys[down_ID])
                    ):
                    return nearby_ID
                    
        for nearby_ID in segs_by_topo_length[topo_length-1]:
            if nearby_ID != ID:
                if ys[ID] == ys[nearby_ID]:
                    return nearby_ID
                    
        return False

    def create_planform(net, ys):
        """
        Generate the final x and y coordinates to plot.
        """
        
        planform = {}
        for i,seg in enumerate(net.list_of_LongProfile_objects):
            if not seg.downstream_segment_IDs:
                x = np.hstack(( seg.x, seg.x_ext[0][-1], seg.x_ext[0][-1] ))/1000.
                y = np.hstack(( np.full(len(seg.x),ys[i]), ys[i], ys[i] ))
            else:
                x = np.hstack(( seg.x, seg.x_ext[0][-1], seg.x_ext[0][-1] ))/1000.
                y = np.hstack(( 
                    np.full(len(seg.x),ys[i]), 
                    ys[i], 
                    ys[seg.downstream_segment_IDs[0]]
                    ))
            planform[i] = {'x': x, 'y': y}
        return planform
        
    def plot_planform(planform):
        """
        Plot the planform.
        """
        
        for i in planform:
            plt.plot(planform[i]['x'], planform[i]['y'])
        plt.show()
            
    # ---- Organise segments by distance upstream (topological length)
    # Used later to check for conflicts between segments.
    segs_by_topo_length = {0: [], net.max_topological_length+2: []}
    for i in range(1,net.max_topological_length+2):
        segs_by_topo_length[i] = []
    for i,seg in enumerate(net.list_of_LongProfile_objects):
        topo_length = len(net.find_downstream_IDs(seg.ID))
        segs_by_topo_length[topo_length].append(seg.ID)
    
    # ---- Set up arrays to fill
    ys = np.full( len(net.list_of_LongProfile_objects), np.nan )
    sides = np.full( len(net.list_of_LongProfile_objects), 0)
    up_sides = np.full( len(net.list_of_LongProfile_objects), -1)
    connections = [
        [np.nan,np.nan] for i in range(len(net.list_of_LongProfile_objects))
        ]
    
    # ---- Loop over segments building planform
    for i,seg in enumerate(net.list_of_LongProfile_objects):
        
        # ---- Check if outlet
        if not seg.downstream_segment_IDs:
            ys[i] = 0.
            connections[i][1] = copy.copy(ys[i])
        
        # ---- Otherwise, add segment on to downstream one
        else:
            
            # Some info about the segment
            down_ID = seg.downstream_segment_IDs[0]
            topo_length = len(net.find_downstream_IDs(seg.ID))

            # Add segment based on relationship to downstream segment
            # Record what side of downstream segment the segment is on
            # Update "up_sides" so that next segment goes on the other side
            ys[i] = ys[down_ID] + up_sides[down_ID]
            sides[i] = copy.copy(up_sides[down_ID])
            up_sides[down_ID] *= -1
            
            # Check for conflict
            conflicting_id = check_for_segment_conflicts(
                seg.ID, segs_by_topo_length, net, ys)
            
            # If there is a conflict, move downstream until reaching a segment
            # with the right direction to fix the conflict
            if conflicting_id:
                seg_to_adjust = seg.ID
                while sides[conflicting_id] != sides[seg_to_adjust]:
                    seg_to_adjust = (
                        net.list_of_LongProfile_objects[seg_to_adjust]. \
                        downstream_segment_IDs[0])
            
            # Move everything upstream of that segment out the way until the
            # conflict is addressed
            while check_for_segment_conflicts(
                seg.ID, segs_by_topo_length, net, ys):
                up_IDs_down_ID = net.find_upstream_IDs(seg_to_adjust)
                ys[up_IDs_down_ID] += sides[seg_to_adjust]

    # ---- Get everything starting from zero and positive
    ys -= ys.min() - 1.

    # ---- Create final planform
    planform = create_planform(net, ys)
    if show:
        plot_planform(planform)

    return planform

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

    def __init__(self, magnitude, topology=None):
        self.magnitude = magnitude
        self.links = topology
        self.upstream_segment_IDs = None
        self.downstream_segment_IDs = None
        if not self.links:
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

def power_law(x, k, p):
    """
    Simple power law: y = k*(x^p).
    """
    return k * (x**p)

def power_law_misfit(pars, x, y):
    """
    Compute misfit between some data x,y and power law with given parameters.
    
    k = pars[0]
    p = pars[1]
    x = data x
    y = data y
    """
    k = pars[0]
    p = pars[1]
    modely = power_law(x, k, p)
    misfit = np.mean( np.sqrt( (modely - y)**2. ) )
    return misfit
    
def optimize_power_law(x, y):
    """
    Find optimal power law parameters for data x, y.
    """
    x_scale = x / max(x)
    y_scale = y / max(y)
    fit = minimize(
        power_law_misfit, 
        [1.,1.], 
        args=(x_scale,y_scale), 
        bounds=[[0,None],[0,None]])
    k_scale = fit.x[0]
    p = fit.x[1]
    k = k_scale / (max(x)**p) * max(y)
    return k, p
    
def find_network_hack_parameters(net):
    """
    Find optimal Hack parameters for a given network object.
    
    First get distances downstream from source for each network segment. Then
    optimise power law describing increasing discharge downstream.
    """
    d = []
    Q = []
    for lp in net.list_of_LongProfile_objects:
        upstream_IDs = net.find_upstream_IDs(lp.ID)
        ds = []
        for up_id in upstream_IDs:
            ds.append(min(net.list_of_LongProfile_objects[up_id].x))
        d.append(np.max(lp.x) - min(ds))
        Q.append(np.mean(lp.Q))
    k, p = optimize_power_law(d, Q)
    return {'k': k, 'p': p, 'd': d, 'Q': Q}

def find_network_hack_parameters_non_dim(net):
    """
    Find optimal Hack parameters for a given network object.
    
    First get distances downstream from source for each network segment. Then
    optimise power law describing increasing discharge downstream.
    """
    
    topo_lengths = []
    upstream_sources = []
    upstream_segments = []
    for seg in net.list_of_LongProfile_objects:
            
        up_IDs = net.find_upstream_IDs(seg.ID)
        upstream_segments.append(len(up_IDs))
        up_sources = [ID for ID in up_IDs if ID in net.sources]
        upstream_sources.append(len(up_sources))
        
        topo_length = 0
        for s in up_sources:
            topo_length_i = 0
            down_IDs = net.find_downstream_IDs(s)
            for down_ID in down_IDs:
                topo_length_i += 1
                if down_ID == seg.ID:
                    break
            topo_length = max(topo_length, topo_length_i)
        topo_lengths.append(topo_length)

    k_src, p_src = optimize_power_law(np.array(topo_lengths), np.array(upstream_sources))
    k_seg, p_seg = optimize_power_law(np.array(topo_lengths), np.array(upstream_segments))
    return {'k_src': k_src, 'p_src': p_src, 'k_seg': k_seg, 'p_seg': p_seg, 'lengths': topo_lengths, 'sources': upstream_sources, 'segments': upstream_segments}
