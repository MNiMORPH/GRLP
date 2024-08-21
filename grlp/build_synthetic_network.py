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

# def get_simple_network_setup_params(
#     upstream_segment_list,
#     downstream_segment_list,
#     L,
#     mean_Q,
#     mean_Qs,
#     min_nxs=5,
#     approx_dx=5.e2):
#     """
#     Find spatial discretisation and input discharges for given network
#     geometry and desired length, approximate discretisation, mean discharges.
#     """
# 
#     # ---- Spatial discretisation
# 
#     # Number of links
#     num_links = len(upstream_segment_list)
# 
#     # Find sources
#     sources = [
#         i for i,up_ids in enumerate(upstream_segment_list) if len(up_ids)==0]
# 
#     # Find maximum topological length,
#     # i.e. number of downstream segments to outlet
#     max_topo_length = max([
#         len(downstream_IDs(downstream_segment_list, i))
#         for i in range(num_links)])
# 
#     # Find length of each link so that total length equals L
#     link_length = L / max_topo_length
# 
#     # Find how many nodes for each link, set dx
#     link_n = max(min_nxs, int(link_length/approx_dx))
#     nxs = [link_n for i in range(num_links)]
#     dx = link_length / link_n
# 
#     # ---- Sediment & water discharge
# 
#     # Find number of sources upstream of each point
#     up_sources = []
#     for i in range(num_links):
#         count = 0
#         up_IDs = upstream_IDs(upstream_segment_list, i)
#         for ID in up_IDs:
#             if len(upstream_IDs(upstream_segment_list, ID)) == 1:
#                 count += 1
#         up_sources.append(count)
# 
#     # Find input sediment and water discharge to give specified means
#     Q_in = mean_Q / np.mean(up_sources)
#     Qs_in = mean_Qs / np.mean(up_sources)
# 
#     # Return
#     return nxs, dx, Q_in, Qs_in

# def set_up_network_object(nx_list, dxs, segment_lengths, upstream_segment_list,
#     downstream_segment_list, supply_discharge_list, internal_discharge_list,
#     sediment_discharge_ratio, width_list, evolve=False):
#     """
#     Uses lists of segment length, upstream and downstream segment IDs to build
#     instance of grlp.Network.
#     As well as configuration lists requires network total discharge and
#     sediment supply.
#     Optionally evolve for a while (aiming for steady-state).
#     """
# 
#     # ---- Some parameters for use during set up
#     sources = [i for i in range(len(nx_list)) if not upstream_segment_list[i]]
# 
#     # ---- Basic lp object to get k_Qs for later
#     lp = LongProfile()
#     lp.basic_constants()
#     lp.bedload_lumped_constants()
#     lp.set_hydrologic_constants()
# 
#     # ---- Loop over segments filling lists for network
#     x_ls = []
#     z_ls = []
#     Q_ls = []
#     Qs_ssd_ls = []
#     B_ls = []
#     for i,nx in enumerate(nx_list):
# 
#         # Set up x domain
#         down_IDs = downstream_IDs(downstream_segment_list, i)[1:]
#         down_x = sum([segment_lengths[j] for j in down_IDs])
#         x0 = - down_x - segment_lengths[i]
#         x1 = x0 + segment_lengths[i]
#         x = x0 + np.arange( 0, nx_list[i], 1 ) * dxs[i]
#         x_ls.append(x)
# 
#         # Set discharges
#         if i in sources:
#             min_discharge = supply_discharge_list[i]
#             max_discharge = supply_discharge_list[i] + internal_discharge_list[i]
#         else:
#             up_IDs = upstream_IDs(upstream_segment_list, i)
#             max_discharge = sum([supply_discharge_list[j] + internal_discharge_list[j] for j in up_IDs])
#             min_discharge = max_discharge - internal_discharge_list[i]
#         Q, dQ = np.linspace(min_discharge, max_discharge, len(x), retstep=True)
#         Q_ls.append(Q)
#         mean_discharge = (min_discharge + max_discharge)/2.
#         B = np.linspace(
#             width_list[i]*min_discharge/mean_discharge,
#             width_list[i]*max_discharge/mean_discharge,
#             len(x)
#             )
#         B_ls.append(B)
#         Qs_ssd = dQ/sediment_discharge_ratio/dxs[i]/B/(1.-lp.lambda_p)
#         Qs_ssd = np.full(len(x), Qs_ssd)
#         Qs_ssd_ls.append( Qs_ssd )
# 
#         # Set initial z
#         S0 = (1./(lp.k_Qs*sediment_discharge_ratio))**(6./7.)
#         z = (x.max()-x)*S0
# 
#         # if not mouth, reset downstream elevation to that of downstream segment
#         # needs lists to work backwards from downstream end
#         if downstream_segment_list[i]:
#             z += z_ls[downstream_segment_list[i][0]][0] + dxs[i]*S0
#         z_ls.append(z)
# 
#     # ---- Update x coordinates to run from 0 at furthest upstream point, record max
#     x_min = min([min(x) for x in x_ls])
#     for i,nx in enumerate(nx_list):
#         x_ls[i] -= x_min
#     # AW guess: take max x value and add an increment dx for base-level boundary
#     x_max = max([max(x) + (x[-1]-x[-2]) for x in x_ls])
#     # Then add on some vertical distance to make a straight line to base level
#     dz_for_bl = dxs[0]*S0
#     for _z in z_ls:
#         _z += dz_for_bl
# 
#     # ---- Initialize network object
#     net = Network()
#     net.initialize(
#         config_file = None,
#         x_bl = x_max,
#         z_bl = 0.,
#         S0 = S0 * np.ones(len(nx_list)),
#         upstream_segment_IDs = upstream_segment_list,
#         downstream_segment_IDs = downstream_segment_list,
#         x = x_ls,
#         z = z_ls,
#         Q = Q_ls,
#         B = B_ls,
#         overwrite = False
#         )
#     net.set_niter(3)
#     net.get_z_lengths()
# 
#     # ---- Set segment source-sink-distributed term
#     for i,seg in enumerate(net.list_of_LongProfile_objects):
#         seg.set_source_sink_distributed(Qs_ssd_ls[i])
# 
#     # ---- If requested evolve network, aiming for steady state
#     if evolve:
#         net.evolve_threshold_width_river_network(nt=100, dt=3.15e12)
# 
#     # ---- Compute Qs
#     for seg in net.list_of_LongProfile_objects: seg.compute_Q_s()
# 
#     return net

def generate_x_domain(segment_lengths, downstream_segments, min_nxs, approx_dx):

    # Generate x list
    x_ls = []
    for i,L in enumerate(segment_lengths):
        nx = max(min_nxs, int(L/approx_dx))
        dx = L / nx
        down_IDs = downstream_IDs(downstream_segments, i)[1:]
        down_x = sum([segment_lengths[j] for j in down_IDs])
        x0 = - down_x - segment_lengths[i]
        x1 = x0 + segment_lengths[i]
        x = x0 + np.arange( 0, nx, 1 ) * dx
        x_ls.append(x)
    
    # Update x coordinates to run from 0 at furthest upstream point, record max
    x_min = min([min(x) for x in x_ls])
    for i,L in enumerate(segment_lengths):
        x_ls[i] -= x_min
        
    # Take max x value and add an increment dx for base-level boundary
    x_max = max([max(x) + (x[-1]-x[-2]) for x in x_ls])

    return x_ls, x_max


def generate_discharges(supply_discharges, internal_discharges,
    upstream_segments, xs):
    """
    Produce list of discharge arrays for set up of network object, based on
    given lists of upstream supply discharges, and internal supply discharges.
    """
    
    # ---- Set up lists
    Q_ls = []
    dQ_ls = []
    
    # ---- Loop over segments
    for i,L in enumerate(supply_discharges):
        
        # If an inlet segment, minimum discharge is that supplied, maximum
        # discharge is that plus the internally supplied discharge.
        if not upstream_segments[i]:
            min_discharge = supply_discharges[i]
            max_discharge = supply_discharges[i] + internal_discharges[i]
            
        # Otherwise, we need to sum the discharges upstream.
        # Then, the minimum discharge is the maximum minus that supplied
        # internally to that segment.
        else:
            up_IDs = upstream_IDs(upstream_segments, i)
            max_discharge = sum(
                [supply_discharges[j] + internal_discharges[j] for j in up_IDs]
                )
            min_discharge = max_discharge - internal_discharges[i]
            
        # Interpolate linearly between minimum and maximum discharge.
        Q, dQ = np.linspace(
            min_discharge, max_discharge, len(xs[i])+1, retstep=True
            )
            
        # Save the discharge array and corresponding rate of increase.
        Q_ls.append(Q[:-1])
        dQ_ls.append(dQ)
        
    return Q_ls, dQ_ls


def generate_ssds(discharges, sediment_discharge_ratio, xs, widths,
    upstream_segments, downstream_segments, lambda_p):
    """
    Produce list of source-sink-distributed arrays for set up of network object
    based on given lists of discharge, valley width. Sets
    source-sink-distributed term to cancel out effects on slope of variation in
    discharge.
    """
    
    # ---- Prepare list for output
    Qs_ssd_ls = []
    
    # ---- Loop over segments
    for i,Q in enumerate(discharges):

        # Compute dQ and dx
        dQ = np.diff(Q)[0]
        dx = np.diff(xs[i])[0]
        
        # Compute basic ssd term - fill array
        Qs_ssd = dQ/sediment_discharge_ratio/dx/widths[i]/(1.-lambda_p)
        Qs_ssd = np.full(len(xs[i]), Qs_ssd)
        
        # We need to do some corrections!
        
        # If an inlet segment, half ssd at the upstream end, since we only have
        # half a cell there.
        if not upstream_segments[i]:
            Qs_ssd[0] /= 2.
            
        # Otherwise, we need to adjust the value at the upstream end to account
        # for different widths of the upstream segments. We average all the
        # widths, weighted by dx.
        else:
            Bs = [widths[i][0]]
            dxs = [dx]
            for up_id in upstream_segments[i]:
                Bs.append(widths[up_id][-1])
                dxs.append(np.diff(xs[up_id])[-1])
            B_mean = np.sum(np.array(Bs)*np.array(dxs))/np.sum(dxs)
            Qs_ssd[0] *= widths[i][0] / B_mean
            
        # If an outlet segment, half ssd at the downstream end, since we only
        # have half a cell there.
        if not downstream_segments[i]:
            Qs_ssd[-1] /= 2.
        
        # Save the list.
        Qs_ssd_ls.append( Qs_ssd )
    
    return Qs_ssd_ls
    
def generate_variable_widths(mean_width, mean_discharge, discharges):
    
    widths = []
    for i,Q in enumerate(discharges):
        widths.append( Q * mean_width / mean_discharge )
    
    min_width = np.concatenate(widths).min()
    max_width = np.concatenate(widths).max()
    width_range = max_width - min_width
    
    return widths

def generate_zs(xs, x_max, S0, downstream_segments):
    
    zs = []
    for i,x in enumerate(xs):
        
        z = (x.max()-x)*S0

        # if not mouth, reset downstream elevation to that of downstream segment
        # needs lists to work backwards from downstream end
        if downstream_segments[i]:
            dx = x[1] - x[0]
            z += zs[downstream_segments[i][0]][0] + dx*S0
            
        zs.append(z)

    # Then add on some vertical distance to make a straight line to base level        
    for z in zs:
        dz_for_bl = (x_max-xs[0][-1])*S0
        z += dz_for_bl

    return zs


def generate_random_network(magnitude=None, max_length=None, segment_lengths=None,
    segment_length=None, internal_discharges=None, supply_discharges=None, 
    segment_length_area_ratio=None, supply_area=None, approx_dx=1.e2,
    min_nxs=5, mean_discharge=None, effective_rainfall=1.,
    sediment_discharge_ratio=1.e4, mean_width=100., variable_width=False,
    topology=None, evolve=False):
    
    if not max_length and not segment_length and not segment_lengths:
        print(
            "Error: " +
            "you must specify max_length, segment_length or segment_lengths. " +
            "Exiting.")
        return None
        
    if not mean_discharge and not effective_rainfall and not supply_discharges:
        print(
            "Error: " +
            "you must specify mean_discharge or effective_rainfall or " +
            "supply_discharges. Exiting.")
        return None

    # ---- Basic lp object to get properties later
    lp = LongProfile()
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()


    # ---- Generate network topology
    net_topo = Shreve_Random_Network(
        magnitude=magnitude, 
        segment_length=segment_length,
        segment_length_area_ratio=segment_length_area_ratio,
        supply_area=supply_area,
        max_length=max_length,
        topology=topology
        )


    # ---- Set up x-domain
    
    # If segment lengths are not provided, set to those from topology
    if not segment_lengths:
        segment_lengths = net_topo.segment_lengths
    
    # Generate x list
    x_ls, x_max = generate_x_domain(
        segment_lengths,
        net_topo.downstream_segment_IDs,
        min_nxs,
        approx_dx)


    # ---- Set channel head supply areas
    # If not present in topology object, set to one for channel heads and zero
    # for internal segments
    if net_topo.source_areas:
        supply_areas = net_topo.source_areas
    else:
        supply_areas = []
        for i in range(len(net_topo.upstream_segment_IDs)):
            if len(upstream_IDs(net_topo.upstream_segment_IDs, i)) == 1:
                supply_areas.append(1.)
            else:
                supply_areas.append(0.)
    
    
    # ---- Set internal supply areas
    # If not present in topology object, set to zero for all segments
    if net_topo.segment_areas:
        internal_areas = net_topo.segment_areas
    else:
        internal_areas = [0 for i in net_topo.upstream_segment_IDs]
    
    
    # ---- If mean discharge is provided, need to find effective rainfall
    if mean_discharge:
    
        # Find total length, for normalising
        total_length = sum(segment_lengths)
    
        # Find area upstream of each segment
        # Weighted by segment length relative to total length
        scl_areas = []
        for i in range(len(net_topo.upstream_segment_IDs)):
            area = 0
            up_IDs = upstream_IDs(net_topo.upstream_segment_IDs, i)
            for ID in up_IDs:
                area += supply_areas[ID] + internal_areas[ID]
            area -= internal_areas[i]/2.
            scl_areas.append(area * segment_lengths[i] / total_length)
        mean_area = np.sum(scl_areas)
            
        # Find effective rainfall to give specified mean discharge
        effective_rainfall = mean_discharge / mean_area


    # ---- Set supply and internal discharges
    # If not provided, set using areas and effective rainfall
    if not supply_discharges:
        supply_discharges = np.array(supply_areas) * effective_rainfall
    if not internal_discharges:
        internal_discharges = np.array(internal_areas) * effective_rainfall


    # ---- Generate discharge list
    sources = [
        i for i in range(len(net_topo.upstream_segment_IDs))
            if not net_topo.upstream_segment_IDs[i]
        ]
    Q_ls = generate_discharges(
        supply_discharges,
        internal_discharges,
        sources,
        net_topo.upstream_segment_IDs,
        x_ls
        )


    # ---- Set widths
    if variable_width:
        B_ls = generate_variable_widths(mean_width, mean_discharge, Q_ls)
    else:
        B_ls = [np.full(len(x), mean_width) for x in x_ls]


    # ---- Generate ssd list
    Qs_ssd_ls = generate_ssds(Q_ls, sediment_discharge_ratio, x_ls, B_ls, lp.lambda_p)


    # ---- Set initial z
    S0 = (1./(lp.k_Qs*sediment_discharge_ratio))**(6./7.)
    z_ls = generate_zs(x_ls, x_max, S0, net_topo.downstream_segment_IDs)


    # ---- Initialize network object
    net = Network()
    net.initialize(
        config_file = None,
        x_bl = x_max,
        z_bl = 0.,
        S0 = S0 * np.ones(len(sources)),
        upstream_segment_IDs = net_topo.upstream_segment_IDs,
        downstream_segment_IDs = net_topo.downstream_segment_IDs,
        x = x_ls,
        z = z_ls,
        Q = Q_ls,
        B = B_ls,
        overwrite = False
        )
    net.set_niter(3)
    net.get_z_lengths()
    
    # ---- Set segment source-sink-distributed term
    for i,seg in enumerate(net.list_of_LongProfile_objects):
        seg.set_source_sink_distributed(Qs_ssd_ls[i])

    # ---- If requested evolve network, aiming for steady state
    if evolve:
        net.evolve_threshold_width_river_network(nt=100, dt=3.15e12)
    
    # ---- Compute Qs
    for seg in net.list_of_LongProfile_objects: seg.compute_Q_s()
    
    return net, net_topo


def plot_network(net, show=True):
    """
    Generate a plotable network planform from a network object.
    """

    def check_for_segment_conflicts(ID, segs_by_topo_length, net, ys):
        """
        Check for segments that overlap with the given segment.
        """
        
        seg = net.list_of_LongProfile_objects[ID]
        x_max = seg.x_ext[0].max()
        x_min = seg.x_ext[0].min()
        y = ys[ID]

        for other_topo_length in segs_by_topo_length.keys():
            for other_ID in segs_by_topo_length[other_topo_length]:
                if other_ID != ID:
                    
                    other_seg = net.list_of_LongProfile_objects[other_ID]
                    other_x_max = other_seg.x_ext[0].max()
                    other_x_min = other_seg.x_ext[0].min()
                    other_y = ys[other_ID]
                    
                    if (
                        (other_x_min <= x_min <= other_x_max) or 
                        (other_x_min <= x_max <= other_x_max)
                        ):
                        if y == other_y:
                            return other_ID
                    
                    if other_seg.downstream_segment_IDs:
                        down_y = ys[other_seg.downstream_segment_IDs[0]]
                        if (
                            (x_min <= other_x_max <= x_max) and
                            (min(other_y,down_y) <= y <= max(other_y,down_y))
                            ):
                            return other_ID
                    
                    # # Don't think I need this bit...
                    # if seg.downstream_segment_IDs:
                    #     down_ID = seg.downstream_segment_IDs[0]
                    #     if down_ID != other_ID:
                    #         down_y = ys[down_ID]
                    #         if (
                    #             (other_x_min <= x_max <= other_x_max) and 
                    #             (min(y,down_y) <= other_y <= max(y,down_y))
                    #             ):
                    #             return other_ID
                            
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
                if seg.x.max() >= \
                    net.list_of_LongProfile_objects[conflicting_id].x.max():
                    while sides[conflicting_id] != sides[seg_to_adjust]:
                        seg_to_adjust = (
                            net.list_of_LongProfile_objects[seg_to_adjust]. \
                            downstream_segment_IDs[0])
                else:
                    while sides[seg.ID] == sides[seg_to_adjust]:
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

# class Simple_Network:
#     """
#     Set up simple network with specified total length and segment length.
#     Single-segment tributaries added between each trunk segment.
#     Optionally specify different length for tributaries.
#     """
# 
# 
#     def __init__(self, nx_total, nx_trunk_seg, nx_trib_seg=None):
#         self.nx_total = nx_total
#         self.nx_trunk_seg = nx_trunk_seg
#         if nx_trib_seg:
#             self.nx_trib_seg = nx_trib_seg
#         else:
#             self.nx_trib_seg = nx_trunk_seg
# 
#         self.upstream_segment_IDs = None
#         self.downstream_segment_IDs = None
#         self.nxs = None
# 
#         self.build_network()
# 
#     def add_segment(self, down_ID, trunk=False):
#         """
#         Add segment to network.
#         """
# 
#         # ID number for new segment
#         ID = len(self.nxs)
# 
#         # Check whether trunk or tributary, assign length
#         if trunk:
#             nx = self.nx_trunk_seg
#         else:
#             nx = self.nx_trib_seg
# 
#         # Make sure length won't leave floating point at end
#         # Can't make a segment with only one point...
#         # Instead add one to final segment
#         if self.nx_total - sum(self.trunk_nxs) <= nx + 1:
#             nx = self.nx_total - sum(self.trunk_nxs)
# 
#         # Update nx lists
#         self.nxs.append(nx)
#         if trunk:
#             self.trunk_nxs.append(nx)
# 
#         # Update topology lists
#         self.downstream_segment_IDs.append([down_ID])
#         self.upstream_segment_IDs[down_ID].append(ID)
#         self.upstream_segment_IDs.append([])
# 
#         # Return ID for future use
#         return ID
# 
#     def build_network(self):
#         """
#         Generate network topology.
#         """
# 
#         # initialise lists
#         self.upstream_segment_IDs = [[]]
#         self.downstream_segment_IDs = [[]]
#         self.nxs = [self.nx_trunk_seg]
#         self.trunk_nxs = []
#         trunk_ID = 0
# 
#         # add segments until total length is reached
#         while self.nx_total - sum(self.trunk_nxs) > 0:
#             __ = self.add_segment(trunk_ID)
#             trunk_ID = self.add_segment(trunk_ID, trunk=True)


class Shreve_Random_Network:
    """
    Create random network of specified magnitude.
    Follows algorithm described in Shreve (1974, Water Resource Res.).
    Creates lists of upstream/downstream segment IDs for us in GRLP.
    Optionaly assigns segment lengths and drainage areas with constant values
    or from scipy.stats random number generators.
    """

    def __init__(self, magnitude, segment_length=None,
        segment_length_area_ratio=None, supply_area=None, max_length=None,
        segment_lengths=None, segment_length_area_ratios=None,
        source_areas=None, segment_areas=None, topology=None):
        
        self.magnitude = magnitude
        self.links = topology
        self.upstream_segment_IDs = None
        self.downstream_segment_IDs = None
        self.segment_length = segment_length
        self.segment_length_area_ratio = segment_length_area_ratio
        self.supply_area = supply_area
        self.max_length = max_length
        self.segment_lengths = segment_lengths
        self.segment_length_area_ratios = segment_length_area_ratios
        self.source_areas = None
        self.segment_areas = None
        if not self.links:
            self.build_network_topology()
        self.build_lists()
        if not segment_lengths:
            if segment_length or max_length:
                self.set_segment_lengths()
        if segment_length_area_ratio or segment_length_area_ratios:
            self.set_segment_areas()

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

                # If external link, work back downstream to last free internal
                # link
                else:
                    while len(self.upstream_segment_IDs[down_seg]) > 1:
                        down_seg = down_segs.pop(-1)
                seg += 1

    def set_segment_lengths(self):
        
        if self.segment_length:
            self.segment_lengths = []
            for i in range(len(self.upstream_segment_IDs)):
                try:
                    segment_length = float(self.segment_length)
                except TypeError:
                    segment_length = self.segment_length.rvs(size=1)[0]
                self.segment_lengths.append(segment_length)
                
        else:
            self.segment_lengths = np.ones(len(self.upstream_segment_IDs))
            
        if self.max_length:
            
            # Find maximum length
            lengths = []
            for i in range(len(self.upstream_segment_IDs)):
                length = 0
                down_IDs = downstream_IDs(self.downstream_segment_IDs, i)
                for j in down_IDs:
                    length += self.segment_lengths[j]
                lengths.append(length)
            max_length_pre_scale = max(lengths)
            
            # Rescale segment lengths
            scale = self.max_length / max_length_pre_scale
            self.segment_lengths = [
                length*scale for length in self.segment_lengths
                ]
                
    def set_segment_areas(self):
        
        self.source_areas = []
        for i in range(len(self.upstream_segment_IDs)):
            if not self.upstream_segment_IDs[i]:
                try:
                    source_area = float(self.supply_area)
                except TypeError:
                    source_area = self.supply_area.rvs(size=1)[0]
                self.source_areas.append(source_area)
            else:
                self.source_areas.append(0.)

        if not self.segment_length_area_ratios:
            self.segment_length_area_ratios = []
            for i in range(len(self.upstream_segment_IDs)):
                try:
                    self.segment_length_area_ratios.append(
                        float(self.segment_length_area_ratio)
                        )
                except TypeError:
                    self.segment_length_area_ratios.append(
                        self.segment_length_area_ratio.rvs(size=1)[0]
                        )
                        
        self.segment_areas = []
        for i in range(len(self.upstream_segment_IDs)):        
            self.segment_areas.append(
                self.segment_lengths[i] * 
                self.segment_length_area_ratios[i]
                )


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
