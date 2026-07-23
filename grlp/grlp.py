import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse import spdiags, identity, block_diag
from scipy import sparse
from scipy.sparse.linalg import spsolve, isolve
from scipy.stats import linregress
import networkx as nx
import warnings
import sys
import copy
from scipy.optimize import minimize

from . import solver


def _power_law(x, k, p):
    """
    Simple power law: y = k*(x^p).
    """
    return k * (x**p)


def _power_law_misfit(pars, x, y):
    """
    Compute misfit between some data x,y and power law with given parameters.

    k = pars[0]
    p = pars[1]
    x = data x
    y = data y
    """
    k = pars[0]
    p = pars[1]
    modely = _power_law(x, k, p)
    misfit = np.mean( np.sqrt( (modely - y)**2. ) )
    return misfit


def _optimize_power_law(x, y):
    """
    Find optimal power law parameters for data x, y.
    """
    x_scale = x / max(x)
    y_scale = y / max(y)
    fit = minimize(
        _power_law_misfit,
        [1.,1.],
        args=(x_scale,y_scale),
        bounds=[[0,None],[0,None]])
    k_scale = fit.x[0]
    p = fit.x[1]
    k = k_scale / (max(x)**p) * max(y)
    return k, p


class LongProfile(object):
    """
    Gravel-bed river long-profile solution builder and solver
    """
    def __init__(self):
        self.z = None
        self.x = None
        self.A = None
        self.Q = None
        self.B = None
        self.D = None # Grain size; needed only to resolve width and depth
        self.b = None # Width and depth need not be resoloved to compute
        self.h = None # long-profile evolution
        self.S0 = None # S0 for Q_s_0 where there is a defined boundary input
        self.dx_2cell = None
        self.x_ghost_upstream = None # boundary ghost-node position, upstream
        self.x_ghost_downstream = None # and downstream: the two scalars that
                                       # replace the padded x_ext ends
        self.Q_s_0 = None
        self.Q_ghost_upstream = None
        self.Q_ghost_downstream = None
        self.A_ghost_upstream = None # boundary ghost-node drainage area, held
        self.A_ghost_downstream = None # so an imported-A boundary survives a
                                       # later conversion to discharge
        self.z_bl = None
        self.ssd = 0. # distributed sources or sinks
        self.sinuosity = 1.
        self.intermittency = 1.
        self.t = 0
        self.upstream_segment_IDs = []
        self.downstream_segment_IDs = []
        self.ID = None
        self.downstream_fining_subsidence_equivalent = 0.
        self.gravel_fractional_loss_per_km = None
        #self.downstream_dx = None # not necessary if x_ext given
        #self.basic_constants()
        self.L = None
        # The extra water added at upstream tributary junctions,
        # spread evenly across all branches [a kludge]
        # Setting default to 0
        self.dQ_up_jcn = 0.

    def set_ID(self, ID):
        """
        Set the ID of this segment
        """
        self.ID = ID

    def set_upstream_segment_IDs(self, upstream_segment_IDs):
        """
        Set a list of ID numbers assigned to upstream river segments
        Requires list or None input
        """
        self.upstream_segment_IDs = upstream_segment_IDs

    def set_downstream_segment_IDs(self, downstream_segment_IDs):
        """
        Set a list of ID numbers assigned to downstream river segments
        Requires list or None input
        """
        self.downstream_segment_IDs = downstream_segment_IDs

    #def set_downstream_dx(self, downstream_dx)
    #    """
    #    Downstream dx, if applicable, for linking together segments in a
    #    network. This could be part of x_ext
    #    """
    #    self.downstream_dx = downstream_dx

    def basic_constants(self):
        self.lambda_p = 0.35
        self.rho_s = 2650.
        self.rho = 1000.
        self.g = 9.805
        self.epsilon = 0.2 # Parker channel criterion
        self.tau_star_c = 0.0495
        self.phi = 3.97 # coefficient for Wong and Parker MPM

    def bedload_lumped_constants(self):
        self.k_qs = self.phi * ((self.rho_s - self.rho)/self.rho)**.5 \
                    * self.g**.5 * self.epsilon**1.5 * self.tau_star_c**1.5
        self.k_b = 0.17 * self.g**(-.5) \
                   * ((self.rho_s - self.rho)/self.rho)**(-5/3.) \
                   * (1+self.epsilon)**(-5/3.) * self.tau_star_c**(-5/3.)
        self.k_Qs = self.k_b * self.k_qs

    def set_hydrologic_constants(self, P_xA=7/4., P_AQ=0.7, P_xQ=None):
        self.P_xA = P_xA # inverse Hack exponent
        self.P_AQ = P_AQ # drainage area -- discharge exponent
        if P_xQ:
            warnings.warn("P_xQ may be inconsistent with P_xA and P_AQ")
            self.P_xQ = P_xQ
        else:
            self.P_xQ = P_xA * P_AQ

    def set_intermittency(self, I):
        self.intermittency = I

    def set_x(self, x=None, x_ext=None, dx=None, nx=None, x0=None,
              verbose=True):
        """
        Set x directly or calculate it.
        Pass one of three options:
        x alone
        x_ext alone (this will also define x)
        dx, nx, and x0
        """
        if x is not None:
            if verbose:
                warnings.warn("\n"+
                      "Passing x alone leaves boundary conditions undefined."+
                      "\n"+
                      "The ghost-node positions default to a linear "+
                      "extrapolation of the end cells;"+
                      "\n"+
                      "pass x_ext to set them explicitly, and be sure to "+
                      "define the base level and upstream supply."+
                      "\n"
                      )
            self.x = np.array(x)
            self.dx = np.diff(self.x)
            self.dx_2cell = self.x[2:] - self.x[:-2]
            # Ghost-node positions by linear extrapolation (x alone leaves the
            # boundary otherwise undefined; this is the sensible default)
            self.x_ghost_upstream = 2*self.x[0] - self.x[1]
            self.x_ghost_downstream = 2*self.x[-1] - self.x[-2]
        if x_ext is not None:
            if x is not None:
                warnings.warn("\n"+
                                "x_ext will overwrite x values just set by "+
                                "passing x")
            x_ext = np.array(x_ext)
            self.x = x_ext[1:-1]
            # Preserve the supplied ghost positions exactly (they need not be
            # a linear extrapolation of the interior grid)
            self.x_ghost_upstream = x_ext[0]
            self.x_ghost_downstream = x_ext[-1]
            self.dx_2cell = self.x[2:] - self.x[:-2]
            self.dx = np.diff(self.x)
        elif (dx is not None) and (nx is not None) and (x0 is not None):
            self.x = np.arange(x0, x0+dx*nx, dx)
            self.x_ghost_upstream = x0 - dx # = 2*x[0] - x[1]
            self.x_ghost_downstream = self.x[-1] + dx # = 2*x[-1] - x[-2]
            self.dx = dx * np.ones(len(self.x) - 1)
            self.dx_2cell = np.ones(len(self.x) - 1)
        elif x is None:
            sys.exit("Need x OR x_ext OR (dx, nx, x0)")
        self.nx = len(self.x)
        self.L = self.x_ghost_downstream - self.x_ghost_upstream
        if (nx is not None) and (nx != self.nx):
            warnings.warn("Choosing x length instead of supplied nx")

    def set_z(self, z=None, z_ext=None, S0=None, z1=0):
        """
        Set z directly or calculate it

        ::

            S0 = initial slope (negative for flow from left to right)
                 unlike in the paper, this is a dz/dx value down the valley,
                 so we account for sinuosity as well at the upstream boundary.
            z1 = elevation value at RHS

        Code works for single segments; the Network class manages z values
        on its own.
        """
        if z is not None:
            self.z = z
        elif z_ext is not None:
            if self.z is not None:
                self.z = z_ext[1:-1]
        elif self.x.any() and (S0 is not None):
            self.z = self.x * S0 + (z1 - self.x[-1] * S0)
        else:
            sys.exit("Error defining variable")
        #self.dz = self.z_ext[2:] - self.z_ext[:-2] # dz over 2*dx!

    def set_A(self, A=None, A_ext=None, k_xA=None, P_xA=None):
        """
        Set A directly or calculate it
        """
        if A is not None:
            self.A = A
            self.A_ghost_upstream = 2*A[0] - A[1]
            self.A_ghost_downstream = 2*A[-1] - A[-2]
        elif A_ext is not None:
            A_ext = np.array(A_ext)
            self.A = A_ext[1:-1]
            self.A_ghost_upstream = A_ext[0]
            self.A_ghost_downstream = A_ext[-1]
        elif self.x.any():
            self.k_xA = k_xA
            if P_xA:
                self.P_xA = P_xA
            self.A = self.k_xA * self.x**self.P_xA
            self.A_ghost_upstream = self.k_xA * self.x_ghost_upstream**self.P_xA
            self.A_ghost_downstream = self.k_xA \
                                       * self.x_ghost_downstream**self.P_xA
        else:
            sys.exit("Error defining variable")

    def set_Q(self, Q=None, Q_ext=None, q_R=None, A_R=None, P_AQ=None,
              k_xQ=None, P_xQ=None, update_Qs_input=True):
        """
        Set Q directly or calculate it
        q_R = storm rainfall rate [m/hr]

        Set only for a 1D array: Q is handled externally for river networks
        """
        if k_xQ is not None:
            self.k_xQ = k_xQ
        if P_xQ is not None:
            self.P_xQ = P_xQ
        if Q is not None:
            # Check if it is a scalar or an array
            if hasattr(Q, "__iter__"):
                self.Q = Q
            else:
                # Assuming "x" is known already
                self.Q = Q * np.ones(self.x.shape)
            Q_ext = np.hstack( (2*self.Q[0]-self.Q[1],
                                self.Q,
                                2*self.Q[-1]-self.Q[-2]) )
        elif Q_ext is not None:
            self.Q = Q_ext[1:-1]
        elif q_R and A_R:
            if P_AQ:
                self.P_AQ = P_AQ
            q_R = q_R/3600. # to m/s
            A_ext = np.hstack(( self.A_ghost_upstream, self.A,
                                self.A_ghost_downstream ))
            Q_ext = q_R * np.minimum(A_R, A_ext) \
                    * (A_ext/np.minimum(A_R,A_ext))**self.P_AQ
            self.Q = Q_ext[1:-1]
        elif self.x.any() and k_xQ and P_xQ:
            self.Q = k_xQ * self.x**P_xQ
            _x_ext = np.hstack(( self.x_ghost_upstream, self.x,
                                 self.x_ghost_downstream ))
            Q_ext = k_xQ * _x_ext**P_xQ
        else:
            sys.exit("Error defining variable")
        # [was VERY helpful; currently deprecated]
        # dQ_ext_upwind over the cell and the cell upstream of it
        # Therefore, the same cell in which Q increases
        # also experiences the nonzero dQ/dx (23.09.23)
        #
        # This is an update over Eq. D3 in Wickert & Schildgen (2019),
        # who applied a central difference, which then spread the dQ out
        # over two cells.
        #
        # Update x2: The impulsive increase in Q and dQ may not interact well
        # with dz/dx and d2z/dx2, which are calculated across the target cell
        # and both its neighbors. Therefore, I believe that we should go back
        # to the 2-cell approach and then "smear" Q by averaging it also over
        # the target cell and its neighbors.
        #
        # [old material:]
        # This then combines with the 1/4 factor in the coefficients
        # for the stencil that results from (2*dx)**2
        # Explicit boundary ghost discharges: analytic when Q is a power law,
        # linear extrapolation when Q is an array. The walking solver uses
        # these at a channel head / river mouth (linear fallback if unset).
        self.Q_ghost_upstream = Q_ext[0]
        self.Q_ghost_downstream = Q_ext[-1]
        self.dQ_ext_2cell = Q_ext[2:] - Q_ext[:-2]
        #dQ_ext_upwind = Q_ext[1:-1] - Q_ext[:-2]
        # Keep sediment supply tied to water supply, except
        # by changing S_0, to only turn one knob for one change (Q/Qs)
        if update_Qs_input:
            if self.Q_s_0:
                self.set_Qs_input_upstream(self.Q_s_0)

    def set_B(self, B=None, k_xB=None, P_xB=None):
        """
        Set B directly or calculate it: B = k_xB * x**P_xB
        """
        if B is not None:
            # Check if it is a scalar or an array
            if hasattr(B, "__iter__"):
                self.B = B
            else:
                # Assuming "x" is known already
                self.B = B * np.ones(self.x.shape)
        elif k_xB and self.x.any():
            self.B = k_xB * self.x**P_xB
            self.k_xB = k_xB
            self.P_xB = P_xB

    def set_uplift_rate(self, U):
        """
        Uplift rate if positive -- or equivalently, rate of base-level fall
        Subsidence (or base-level rise) accomplished by negative uplift
        This can also be thought as a general source (uplift) or sink
        (subsidence) term -- for example, from landslide inputs (source) or
        downstream fining (sink).
        SI units (m/s)
        """
        # Keeping sign positive now and including as adding to river
        # instead of dropping base level
        self.U = U # not sure this is the best -- flipping the sign

    def set_source_sink_distributed(self, ssd):
        self.ssd = ssd

    def set_Sternberg_gravel_loss(self, gravel_fractional_loss_per_km=None ):
        """
        Based on Dingle et al. (2017).
        """
        """
        distance_downstream_from_boundary = self.x - self.x[0]
        gravel_input = self.Q_s_0 * \
                       np.exp( -gravel_fractional_loss_per_km/1000.
                                * distance_downstream_from_boundary)
        # Q_s_0 may cause it to break; perhaps just consider Q_s[0]
        self.downstream_fining_subsidence_equivalent = np.hstack((
                0, np.diff(gravel_input) / ( (1-self.lambda_p) * self.B[1:] )
                ))
        """
        # Use finite difference
        # Though this is now separated from the implicit solution
        # Must include in semi-implicit by iterating through calls to this
        if gravel_fractional_loss_per_km is not None:
            self.gravel_fractional_loss_per_km = gravel_fractional_loss_per_km
        elif self.gravel_fractional_loss_per_km is not None:
            pass
        else:
            raise ValueError('You must define gravel_fractional_loss_per_km.')
        self.compute_Q_s()
        # gravel_fractional_loss_per_km is per kilometer; divide by 1000 to get
        # the per-meter coefficient used against the (metre-based) geometry.
        self.downstream_fining_subsidence_equivalent = \
                - self.gravel_fractional_loss_per_km / 1000. * self.Q_s \
                / ( (1-self.lambda_p) * self.B )

    def set_niter(self, niter):
        self.niter = niter

    def set_Qs_input_upstream(self, Q_s_0):
        """
        S0, the boundary-condition slope, is set as a function of Q_s_0.
        Note that here I use S in a different sense than in the paper:
        sinuosity is external to S here, meaning that it has to be included
        explicitly in the equation to compute S0. This is so S0 impacts
        dz/dx|boundary directly, instead of needing to be modulated each time.
        """
        self.Q_s_0 = Q_s_0
        # Q[0] is centerpoint of S?
        self.S0 = np.sign(self.Q[0]) * self.sinuosity * \
                      ( np.abs(Q_s_0) /
                        ( self.k_Qs
                              * np.abs(self.Q[0])) )**(6/7.)
        # Give upstream cell the same width as the first cell in domain
        self.z_ghost_upstream = self.z[0] + self.S0 * self.dx[0]

    def set_S0(self, S0):
        """
        Set the upstream boundary by its slope S0 rather than by the sediment
        supply Q_s_0 -- convenient when the slope is known independently, e.g.
        measured from a DEM. The equivalent Q_s_0 is derived from S0 and the
        current discharge (the exact inverse of the relation in
        set_Qs_input_upstream) and applied through that setter, so the boundary
        remains a sediment supply: a later change in Q adjusts the boundary
        slope, as a real physical forcing should.
        S0 follows the boundary convention (sinuosity external, positive for a
        river descending downstream); its sign is taken from the discharge.
        """
        Q_s_0 = self.k_Qs * np.abs(self.Q[0]) \
                    * ( np.abs(S0) / self.sinuosity )**(7/6.)
        self.set_Qs_input_upstream(Q_s_0)

    def set_z_bl(self, z_bl):
        """
        Set the right-hand Dirichlet boundary conditions, i.e. the base level,
        given in the variable "z_bl" (elevation, base level)

        For 1D single-segment mode, not network.
        """
        self.z_bl = z_bl

    def set_x_bl(self, x_bl):
        self.x_bl = x_bl
        self.x_ghost_downstream = self.x_bl

    def evolve_threshold_width_river(self, nt=1, dt=3.15E7):
        """
        Solve the triadiagonal matrix through time, with a given
        number of time steps (nt) and time-step length (dt)

        A single segment is solved as a one-edge network by the unified walking
        solver: this LongProfile is wrapped in a Network of one segment and
        stepped there (Network._evolve_by_walking). No padded ghost arrays are
        maintained; the diagnostics read the node-based fields directly.
        """
        if (len(self.upstream_segment_IDs) > 0) or \
           (len(self.downstream_segment_IDs) > 0):
            warnings.warn("Unset boundary conditions for river segment"+
                          "in network.\n"+
                          "Local solution on segment will not be sensible.")
        self.nt = nt
        self.ID = 0
        net = Network( [self] )
        net.t = self.t
        net.set_niter(self.niter)
        net.get_z_lengths()
        net.evolve_threshold_width_river_network(nt, dt)
        self.Qs_internal = 1/(1-self.lambda_p) * np.cumsum(self.dz_dt) \
                            * self.B + self.Q_s_0

    def build_LHS_coeff_C0(self, dt=3.15E7):
        """
        Build the LHS coefficient for the tridiagonal matrix.
        This is the "C0" coefficient, which is likely to be constant and
        uniform unless there are dynamic changes in width (epsilon_0 in
        k_Qs), sinuosity, or intermittency, in space and/or through time

        See eq. D3. "1/4" subsumed into "build matrices".
        For C1 (other function), Q/B included as well.
        """
        self.dt = dt # Needed to build C0, C1
        self.C0 = self.k_Qs * self.intermittency \
                    / ((1-self.lambda_p) * self.sinuosity**(7/6.)) \
                    * self.dt

    def analytical_threshold_width(self, P_xQ=None, x0=None, x1=None,
                                   z0=None, z1=None):
        """
        Analytical: no uplift
        """
        if x0 is None:
            x0 = self.x[0]
        if x1 is None:
            x1 = self.x[-1]
        if z0 is None:
            z0 = self.z[0]
        if z1 is None:
            z1 = self.z[-1]
        if P_xQ is None:
            P_xQ = self.P_xQ
        #print P_xQ
        #self.zanalytical2 = (z1 - z0) * (self.x**e - x0**e)/(x1**e - x0**e) + z0
        self.P_a = 1 - 6*P_xQ/7. # beta
        self.k_a = 1/(x1**self.P_a - x0**self.P_a) * (z1 - z0) # alpha
        self.c_a = z0 - x0**self.P_a/(x1**self.P_a - x0**self.P_a) * (z1 - z0) # gamma
        self.zanalytical = self.k_a * self.x**self.P_a + self.c_a
        return self.zanalytical

    def analytical_threshold_width_perturbation(self, P_xQ=None, x0=None, x1=None,
                                   z0=None, z1=None, U=None):
        """
        DEPRECATED and known to be incorrect.

        This early ("First attempt at perturbation theory") closed-form attempt
        at a steady-state solution with uplift does not reproduce the numerical
        model: its constants of integration blow up and catastrophically
        cancel, returning unphysical elevations. It is retained only for the
        historical record.

        Use analytical_threshold_width_uplift instead, which gives the correct
        steady-state solution with uplift (or base-level fall).
        """
        warnings.warn(
            "analytical_threshold_width_perturbation is deprecated and "
            "incorrect; use analytical_threshold_width_uplift instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if x0 is None:
            x0 = self.x[0]
        if x1 is None:
            x1 = self.x[-1]
        if z0 is None:
            z0 = self.z[0]
        if z1 is None:
            z1 = self.z[-1]
        if P_xQ is None:
            P_xQ = self.P_xQ
        if U is None:
            U = self.U
        # Base state coefficients (no perturbation)
        #self.P_a = 1 - 6*P_xB/7. # beta
        #self.k_a = (z1 - z0)/(x1**self.P_a - x0**self.P_a) # alpha
        #self.c_a = z0 - x0**self.P_a/(x1**self.P_a - x0**self.P_a) * (z1 - z0) # gamma
        # Coefficients
        K = self.k_Qs * self.sinuosity * self.intermittency \
            / (1 - self.lambda_p) \
            * abs(self.k_a * self.P_a)**(1/6.) \
            * self.k_xQ / self.k_xB
        P = self.P_xQ
        #print(P)
        # Constants of integration
        #c1 = self.U * (x0**(P+2) - x1**(P+2)) / (K*(P-2)*(self.P_a + P - 2) \
        #     + (x1**self.P_a - x0**self.P_a) / self.P_a
        #     - z1
        c1 = ( self.U * (x0**(P+2) - x1**(P+2)) / (K*(P-2)*(self.P_a + P - 2)) \
               + z0 - z1 ) \
             / ( (x0**self.P_a - x1**self.P_a) / self.P_a )
        c2 = - (c1 * x1**self.P_a)/self.P_a \
             + (U * x1**(2-P))/(K * (P-2) * (self.P_a + P - 2)) + z1
        self.zanalytical = c1 * self.x**self.P_a / self.P_a \
            - self.U * self.x**(2-P) / (K * (P-2) * (self.P_a + P - 2)) \
            + c2

    #def analytical_threshold_width_perturbation_2(self):
    #    self.analytical_threshold_width()

    def analytical_threshold_width_uplift(self, nx_fine=20001):
        """
        Semi-analytical steady-state long profile *with uplift* (equivalently,
        distributed base-level fall), for the transport-limited threshold-width
        gravel river.

        Uplift is applied to the valley surface at rate U (set_uplift_rate).
        At steady state (dz/dt = 0), mass conservation

            (1 - lambda_p) * B * dz/dt = -dQ_s/dx + (1 - lambda_p) * B * U

        forces the long-term-averaged bedload discharge to grow downstream as
        uplifted material is exported:

            Q_s(x) = I * Q_s_0  +  (1 - lambda_p) * U * \\int_{x0}^{x} B dx'

        where I is the intermittency and Q_s_0 is the (during-flood) sediment
        supply set at the upstream boundary. Note the boundary flux carried by
        the long-term average is I * Q_s_0, consistent with compute_Q_s.

        The threshold-width transport law Q_s = k_Qs * I * Q * (S/sinuosity)**(7/6)
        then gives the valley slope, and the bed is recovered by integrating
        up from base level:

            S(x) = sinuosity * ( Q_s(x) / (k_Qs * I * Q(x)) )**(6/7)
            z(x) = z_bl + \\int_{x}^{x_bl} S dx'

        Both integrals are evaluated by the trapezoidal rule on a fine grid of
        nx_fine points (there is no elementary closed form for general P_xQ and
        P_xB), then sampled at the model nodes. Because the reference grid is
        independent of the model grid, comparing the numerical solution to this
        one exhibits first-order convergence under model grid refinement.

        Requires the power-law discharge and width parameters (k_xQ, P_xQ,
        k_xB, P_xB), the uplift rate U, the intermittency, the upstream
        sediment supply (Q_s_0, via set_Qs_input_upstream), and base level.

        Sets and returns self.zanalytical (elevations at self.x).
        """
        for attr in ("k_xQ", "P_xQ", "k_xB", "P_xB", "U", "Q_s_0"):
            if getattr(self, attr, None) is None:
                raise ValueError(
                    "analytical_threshold_width_uplift requires the power-law "
                    "parameters (k_xQ, P_xQ, k_xB, P_xB), the uplift rate, and "
                    "the upstream sediment supply Q_s_0 to be set; missing: "
                    + attr + "."
                )
        lambda_p = self.lambda_p
        I = self.intermittency
        # Fine, model-independent reference grid spanning the boundary nodes.
        x0 = self.x_ghost_upstream
        x_bl = self.x_ghost_downstream
        z_bl = self.z_bl
        xf = np.linspace(x0, x_bl, nx_fine)
        Bf = self.k_xB * xf**self.P_xB
        Qf = self.k_xQ * xf**self.P_xQ
        # Cumulative trapezoidal integral of B from x0 (numpy-only, so no
        # dependence on the moving scipy.integrate cumtrapz/cumulative_trapezoid
        # name).
        intB = np.concatenate((
            [0.],
            np.cumsum(0.5 * (Bf[1:] + Bf[:-1]) * np.diff(xf))
        ))
        Qs = I * self.Q_s_0 + (1 - lambda_p) * self.U * intB
        Sf = self.sinuosity * (Qs / (self.k_Qs * I * Qf))**(6/7.)
        # z by integrating slope upstream from base level:
        # z(x) = z_bl + \int_x^{x_bl} S dx'. Integrate on the reversed grid.
        dz_up = np.concatenate((
            [0.],
            np.cumsum(0.5 * (Sf[1:] + Sf[:-1]) * np.diff(xf))
        ))
        zf = z_bl + (dz_up[-1] - dz_up)
        self.zanalytical = np.interp(self.x, xf, zf)
        return self.zanalytical

    def compute_Q_s(self):
        """
        Slope and sediment discharge at each node. A single segment is solved
        as a one-edge network (exactly as evolve_threshold_width_river does),
        so its slopes and sediment discharges come from walking to each node's
        real neighbour -- the channel-head ghost from the boundary slope S0, the
        outlet ghost from base level -- rather than from a padded z_ext array.
        Sets self.S and self.Q_s.
        """
        Network( [self] ).compute_Q_s()

    def compute_channel_width(self):
        if self.D is not None:
            self.b = 0.17 / ( self.g**.5
                              * ((self.rho_s - self.rho)/self.rho)**(5/3.)
                              * (1+self.epsilon)**(5/3.)
                              * (self.tau_star_c**(5/3.)) ) \
                            * self.Q * self.S**(7/6.) / self.D**1.5
        else:
            raise ValueError('Set grain size to compute channel width.')

    def compute_flow_depth(self):
        if self.D is not None:
            self.h = (self.rho_s - self.rho)/self.rho * (1+self.epsilon) \
                        * self.tau_star_c * self.D / self.S
        else:
            raise ValueError('Set grain size to compute channel depth.')

    def slope_area(self, verbose=False):
        # Slope from the unified walk (a one-edge network). compute_Q_s returns
        # the down-valley channel slope (divided by sinuosity); the slope-area
        # analysis uses the topographic gradient dz/dx, so multiply sinuosity
        # back in.
        self.compute_Q_s()
        self.S = self.S * self.sinuosity
        logS = np.log10(self.S)
        logA = np.log10(self.A)
        out = linregress(logA[1:-1], logS[1:-1]) # remove edge effects
        self.theta = -out.slope
        self.ks = 10.**out.intercept
        self.thetaR2 = out.rvalue**2.
        if verbose:
            print("Concavity = ", self.theta)
            print("k_s = ", self.ks)
            print("R2 = ", out.rvalue**2.)

    def compute_diffusivity(self):
        """
        Compute diffusivity at each point along valley.
        From linearized version of threshold width equation.
        """
        self.compute_Q_s()
        self.diffusivity = (7./6.) * self.k_Qs * self.intermittency * self.Q * self.S**(1./6.) \
            / self.sinuosity**(7./6.) / self.B / (1. - self.lambda_p)

    def compute_length(self):
        """
        Compute total segment length, spanning the boundary ghost nodes.
        """
        self.L = self.x_ghost_downstream - self.x_ghost_upstream

    def compute_equilibration_time(self):
        """
        Compute valley equilibration time (sensu Paola et al., 1992).
        From scaling of linearized version of threshold width equation.
        """
        self.compute_diffusivity()
        if not self.L:
            self.compute_length()
        self.equilibration_time = self.L**2. / self.diffusivity.mean()

    def compute_e_folding_time(self, n):
        """
        Compute valley e-folding times as function of wavenumber.
        From solving linearized version of threshold width equation.
        """
        self.compute_equilibration_time()
        return 1./self.diffusivity.mean()/self.compute_wavenumber(n)**2.

    def compute_wavenumber(self, n):
        """
        Compute wavenumber for series solutions to linearized version of
        threshold width equation.
        """
        if not self.L:
            self.compute_length()
        return (2*n + 1) * np.pi / 2. / self.L

    def compute_series_coefficient(self, n, period):
        """
        Compute coefficient for series solutions to linearized version of
        threshold width equation.
        """
        if not self.L:
            self.compute_length()
        return 4 * np.pi * period / self.L / \
            ((period**2.) * (self.diffusivity.mean()**2.) *
                (self.compute_wavenumber(n)**4.) + (4.*np.pi**2.))

    def compute_z_series_terms(self, period, nsum):
        """
        Compute amplitudes of cos(2*pi*t/P) and sin(2*pi*t/P) terms in
        solutions to linearized version of threshold width equation.
        """
        cos_term = 0
        sin_term = 0
        for n in range(nsum):
            coeff_cos = (np.cos(self.compute_wavenumber(n) * self.x) *
                            self.compute_series_coefficient(n, period))
            cos_term += coeff_cos*self.diffusivity.mean()
            sin_term += coeff_cos*2.*np.pi/period/self.compute_wavenumber(n)**2.
        return cos_term, sin_term

    def compute_Qs_series_terms(self, period, nsum):
        """
        Compute amplitudes of cos(2*pi*t/P) and sin(2*pi*t/P) terms in
        solutions to linearized version of threshold width equation.
        """
        cos_term = 0
        sin_term = 0
        for n in range(nsum):
            coeff_sin = (np.sin(self.compute_wavenumber(n) * self.x) *
                            self.compute_series_coefficient(n, period))
            cos_term += coeff_sin*self.diffusivity.mean()*self.compute_wavenumber(n)
            sin_term += coeff_sin*2.*np.pi/period/self.compute_wavenumber(n)
        return cos_term, sin_term

    def compute_z_gain(self, period, nsum=100):
        """
        Compute gain (relative amplitude) between periodic forcing in sediment
        or water supply and response in valley elevation.
        From solving linearized version of threshold width equation.
        """
        if not self.L:
            self.compute_length()
        self.compute_diffusivity()
        cos_term, sin_term = self.compute_z_series_terms(period, nsum)
        return (6./7.) * \
            np.sqrt( ( (self.L - self.x) - sin_term)**2. + cos_term**2. ) \
            / (self.L - self.x)

    def compute_Qs_gain(self, period, A_Qs=0., A_Q=0., nsum=100):
        """
        Compute gain (relative amplitude) between periodic forcing in sediment
        or water supply and response in valley elevation.
        From solving linearized version of threshold width equation.
        """
        self.compute_diffusivity()
        cos_term, sin_term = self.compute_Qs_series_terms(period, nsum)
        return np.sqrt( (A_Qs/(A_Qs - A_Q) - sin_term)**2. + cos_term**2. ) \


    def compute_z_lag(self, period, nsum=100):
        """
        Compute lag time between periodic forcing in sediment or water supply
        and response in valley elevation.
        From solving linearized version of threshold width equation.
        """
        # Basic calculation
        if not self.L:
            self.compute_length()
        self.compute_diffusivity()
        cos_term, sin_term = self.compute_z_series_terms(period, nsum)
        lag = -(period/(2.*np.pi)) * \
                np.arctan( -cos_term / ((self.L - self.x) - sin_term) )

        # Search for and correct any cycle skipping.
        # Arctan function can only resolve -0.25 < lag/period < 0.25
        for i in range(1,len(self.x)):
            if lag[i] < lag[i-1] and lag[i-1] - lag[i] > period/4.:
                lag[i:] += period/2.
            if i != len(self.x)-1:
                if lag[i+1] > lag[i] and lag[i+1] - lag[i] > period/4.:
                    lag[:i+1] += period/2.
                    
        # At least for parameter ranges we are interested in, we expect that
        # lag at the inlet should not be too big. So we can also catch
        # cycle-skipping that effects the whole length of the valley by
        # imposing the condition that lag at the inlet should not exceed
        # 0.5*period.
        while lag[0] > 0.5*period:
            lag -= 0.5*period

        return lag

    def compute_Qs_lag(self, period, A_Qs=0., A_Q=0., nsum=1000):
        """
        Compute lag time between periodic forcing in sediment or water supply
        and response in bedload sediment discharge.
        From solving linearized version of threshold width equation.
        """
        # Basic calculation
        self.compute_diffusivity()
        cos_term, sin_term = self.compute_Qs_series_terms(period, nsum)
        lag = -(period/(2.*np.pi)) * \
                np.arctan( -cos_term / ((A_Qs/(A_Qs - A_Q) - sin_term))  )

        # Search for and correct any cycle skipping.
        # Arctan function can only resolve -0.25 < lag/period < 0.25
        for i in range(1,len(self.x)):
            if lag[i] < lag[i-1] and lag[i-1] - lag[i] > period/4.:
                lag[i:] += period/2.
            if i != len(self.x)-1:
                if lag[i+1] > lag[i] and lag[i+1] - lag[i] > period/4.:
                    lag[:i+1] += period/2.

        # At least for parameter ranges we are interested in, we expect that
        # lag at the inlet should not be too big. So we can also catch
        # cycle-skipping that effects the whole length of the valley by
        # imposing the condition that lag at the inlet should not exceed
        # 0.5*period.
        while lag[0] > 0.5*period:
            lag -= 0.5*period

        return lag


class Network(object):
    """
    Gravel-bed river long-profile solution builder and solver
    """

    def __init__(self, list_of_LongProfile_objects=None):
        """
        Instantiate the Network object with a list of Long Profile objects.
        Other funcitons will iterate over this network and its connectivity
        to create a full tridiagonal matrix solver.
        """
        self.list_of_LongProfile_objects = list_of_LongProfile_objects
        self.t = 0
        self.Q_s_0 = None
        self.S0 = None

    def build_ID_list(self):
        self.IDs = []
        for lp in self.list_of_LongProfile_objects:
            # IDs
            self.IDs.append(lp.ID)
        self.IDs = np.array(self.IDs)

    def get_z_lengths(self):
        self.list_of_segment_lengths = []
        for lp in self.list_of_LongProfile_objects:
            self.list_of_segment_lengths.append(len(lp.z))

    def compute_Q_s(self):
        """
        Sediment discharge and slope at each point in the network, computed by
        walking the topology to each node's real neighbours rather than reading
        maintained ghost arrays. Sets S and Q_s on each segment, using the same
        slope / sediment-discharge relationship as LongProfile.compute_Q_s.

        At a confluence, the head node has one upstream neighbour per incoming
        tributary; as in the single-segment case, S and Q_s there are the
        average over the tributaries.
        """
        for lp in self.list_of_LongProfile_objects:
            # Downstream ghost: base level at the outlet, else the first node
            # of the downstream segment
            if len(lp.downstream_segment_IDs) == 0:
                z_down = lp.z_bl
                if lp.x_ghost_downstream is not None:
                    x_down = lp.x_ghost_downstream
                else:
                    x_down = 2*lp.x[-1] - lp.x[-2]
            else:
                downseg = self.list_of_LongProfile_objects[
                              lp.downstream_segment_IDs[0]]
                z_down = downseg.z[0]
                x_down = downseg.x[0]
            # Upstream ghost(s): the boundary slope S0 at a channel head, else
            # the last node of each incoming tributary
            z_up = []
            x_up = []
            if len(lp.upstream_segment_IDs) == 0:
                _xg = 2*lp.x[0] - lp.x[1]
                x_up.append( _xg )
                z_up.append( lp.z[0] + lp.S0 * ( lp.x[0] - _xg ) )
            else:
                for upseg_ID in lp.upstream_segment_IDs:
                    upseg = self.list_of_LongProfile_objects[upseg_ID]
                    z_up.append( upseg.z[-1] )
                    x_up.append( upseg.x[-1] )
            # Assemble one ghost-padded profile per upstream neighbour and apply
            # the single-segment slope / Q_s relationship, then average (the
            # profiles differ only at the head node)
            S = []
            Q_s = []
            for _zu, _xu in zip(z_up, x_up):
                _z = np.hstack(( _zu, lp.z, z_down ))
                _x = np.hstack(( _xu, lp.x, x_down ))
                _dx = _x[2:] - _x[:-2]
                S.append( np.abs( (_z[2:] - _z[:-2]) / _dx) / lp.sinuosity )
                Q_s.append(
                    -np.sign( _z[2:] - _z[:-2] ) \
                    * lp.k_Qs * lp.intermittency * lp.Q * S[-1]**(7/6.)
                    )
            lp.S = np.mean(S, axis=0)
            lp.Q_s = np.mean(Q_s, axis=0)

    def set_niter(self, niter):
        # MAKE UNIFORM IN BASE CLASS
        self.niter = niter


    def set_x_bl(self, x_bl=None):
        """
        Set the downstream base-level node's x-position on the single river
        mouth: move the outlet ghost to x_bl (node-based; defaults to one cell
        beyond the last node). Pairs with set_z_bl (elevation) to place base
        level at (x_bl, z_bl); the walker reads x_ghost_downstream at the outlet.

        Expects a single channel-mouth segment (as set_z_bl does).

        ::

            !!!!!
            MAYBE I SHOULD CALL z0 --> z_bl
        """
        if len(self.list_of_channel_mouth_segment_IDs) == 1:
            ID = self.list_of_channel_mouth_segment_IDs[0]
        else:
            sys.exit( ">1 channel-mouth-segment ID listed.\n"+
                      "Simulation not set up to manage >1 river mouth.\n"+
                      "Exiting" )
        lp = self.list_of_LongProfile_objects[ID]
        # Default: one cell beyond the mouth (the linear-extrapolation ghost)
        if x_bl is None:
            x_bl = lp.x[-1] + lp.dx[-1]
        lp.set_x_bl( x_bl )

        # We should have some code to account for changes in both x and z
        # with base-level change, and remeshes the downstream-most segment,
        # as needed

    def set_Qs_input_upstream(self, S0=None, Q_s_0=None):
        """
        Set the upstream sediment-supply boundary at each channel head, from a
        prescribed boundary slope S0 or an input sediment discharge Q_s_0. Each
        may be a scalar (applied at every head) or an iterable (one value per
        head, in channel-head-ID order). Delegates per head to the single-segment
        LongProfile.set_Qs_input_upstream for the Q_s_0 case.
        """
        heads = self.list_of_channel_head_segment_IDs
        # Broadcast a scalar to every head, or take one value per head
        def _per_head(val):
            try:
                iter(val)
                return list(val)
            except TypeError:
                return [val] * len(heads)
        if Q_s_0 is not None:
            self.Q_s_0 = Q_s_0
            for ID, _Qs0 in zip(heads, _per_head(Q_s_0)):
                self.list_of_LongProfile_objects[ID].set_Qs_input_upstream(_Qs0)
        elif S0 is not None:
            self.S0 = S0
            for ID, _S0 in zip(heads, _per_head(S0)):
                lp = self.list_of_LongProfile_objects[ID]
                lp.S0 = _S0
                lp.z_ghost_upstream = lp.z[0] + lp.S0 * lp.dx[0]

    def set_z_bl (self, z0):
        """
        Set downstream boundary (ultimate base level, singular): External

        This function will set only the downstream-most boundary condition.

        It expects a list of length (1) for the class variable:
        self.list_of_channel_mouth_segment_IDs.
        This assumption will have to be relaxed if the code ever be updated
        to allow multiple river mouths.

        Args:
            z0 (float): Base-level elevation. Sets z_bl for the mouth seg

        Returns:
            None
        """
        # !!!!! MAYBE I SHOULD CALL z0 --> z_bl
        if len(self.list_of_channel_mouth_segment_IDs) == 1:
            ID = self.list_of_channel_mouth_segment_IDs[0]
        else:
            sys.exit( ">1 channel-mouth-segment ID listed.\n"+
                      "Simulation not set up to manage >1 river mouth.\n"+
                      "Exiting" )
        # SET DOWNSTREAM BOUNDARY (ULTIMATE BASE LEVEL, SINGULAR): EXTERNAL
        lp = self.list_of_LongProfile_objects[ID]
        lp.z_bl = z0

        # We should have some code to account for changes in both x and z
        # with base-level change, and remeshes the downstream-most segment,
        # as needed

    def create_list_of_channel_head_segment_IDs(self):
        """
        Finds all segments that do not have any upstream tributary segments.

        This is similar to "find_sources", but does not assume that Q_s_0
        has been set.

        Therefore, it is set by topology rather than boundary conditions.
        """
        # Source nodes (no upstream edges); each has one out-edge, the head
        # segment. Sorted by ID so the order matches sediment-supply inputs.
        head_ids = [
            data["segment_id"]
            for node in self.graph.nodes if self.graph.in_degree(node) == 0
            for _, _, data in self.graph.out_edges(node, data=True)
        ]
        self.list_of_channel_head_segment_IDs = sorted(head_ids)

    def create_list_of_channel_mouth_segment_IDs(self):
        """
        Create a list of segments that are channel mouths.

        The length of this list should be 1 (convergent network), but it is
        remaining a list in case of future work including distributary networks.

        Though in principle possible, GRLP need not be run with multiple
        tributary networks (and therefore mutliple mouths).
        It seems cleaner (to me, Wickert) to run each tributary network
        as a separate instance of GRLP.
        """
        # Outlet nodes (no downstream edges); each has one in-edge, the mouth
        # segment. A convergent network has exactly one.
        mouth_ids = [
            data["segment_id"]
            for node in self.graph.nodes if self.graph.out_degree(node) == 0
            for _, _, data in self.graph.in_edges(node, data=True)
        ]
        self.list_of_channel_mouth_segment_IDs = sorted(mouth_ids)

        if len(self.list_of_channel_mouth_segment_IDs) == 1:
            self.channel_mouth_segment_ID = \
                self.list_of_channel_mouth_segment_IDs[0]
        elif len(self.list_of_channel_mouth_segment_IDs) == 0:
            sys.exit("Ahmm... why are there no river mouths?")
        else:
            sys.exit("Ahmm... why are there multiple river mouths?")

    def build_graph(self):
        """
        Build a NetworkX directed-graph representation of the network topology.

        Edges are river segments -- each edge carries its ``segment_id`` and the
        ``LongProfile`` object itself -- and nodes are the junctions between
        them: one ``("source", i)`` per channel head, one ``("jcn", c)`` per
        confluence (named by the segment ``c`` flowing out of it), and a single
        ``("outlet",)``. A one-segment network is a two-node, one-edge graph.

        Flow direction is downstream (edges point from the upstream node to the
        downstream node), so NetworkX ancestor/descendant queries map directly
        onto upstream/downstream reaches.

        This is the emerging source of truth for topology; the per-segment
        upstream/downstream ID lists are retained during the transition.
        Built directly from those ID lists (an empty upstream list marks a
        source, an empty downstream list marks the outlet).
        """
        def upstream_node(seg_id):
            lp = self.list_of_LongProfile_objects[seg_id]
            if not lp.upstream_segment_IDs:
                return ("source", seg_id)
            else:
                return ("jcn", seg_id)

        def downstream_node(seg_id):
            lp = self.list_of_LongProfile_objects[seg_id]
            if not lp.downstream_segment_IDs:
                return ("outlet",)
            return ("jcn", lp.downstream_segment_IDs[0])

        G = nx.DiGraph()
        self._edge_of_segment = {}
        for lp in self.list_of_LongProfile_objects:
            u, v = upstream_node(lp.ID), downstream_node(lp.ID)
            G.add_edge(u, v, segment_id=lp.ID, segment=lp)
            self._edge_of_segment[lp.ID] = (u, v)
        self.graph = G
        return self.graph

    def update_Q(self, Q=None):
        """
        Set discharge within each segment.
        Use the discharge provided during initialize()
        unless a new Q be provided.
        """
        # Update Q (list of arrays)
        _idx = 0
        # Then pass the information to each long-profile object
        for lp in self.list_of_LongProfile_objects:
            #print( _idx )
            lp.Q = Q[_idx]
            #print( len(lp.Q) )
            _idx += 1



    def set_intermittency(self, intermittency):
        """
        Set the flow intermittency value within the channel network
        """
        _isscalar = False
        try:
            iter( intermittency )
        except:
            _isscalar = True
        if _isscalar:
            for lp in self.list_of_LongProfile_objects:
                lp.set_intermittency(intermittency)
        else:
            i = 0
            for lp in self.list_of_LongProfile_objects:
                lp.set_intermittency(intermittency[i])
                i += 1

    def compute_land_areas_around_confluences(self):
        """
        Approximate area at confluence as the sum of lengths
        from the midpoints between nodes above and below the confluence,
        each multiplied by its respective valley width
        I assume that the width above the confluence is basically constant
        until the confluence, and that the same goes for the
        width below
        This is a first attempt
        """
        for lp in self.list_of_LongProfile_objects:
            # COULD CLEAN THIS UP
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if len(lp.upstream_segment_IDs) > 0:
                land_areas_above_confluence = []
                for ID in lp.upstream_segment_IDs:
                    upseg = self.list_of_LongProfile_objects[ID]
                    half_dx = ( lp.x[0] - upseg.x[-1] ) / 2.
                    representative_width = upseg.B[-1]
                    land_areas_above_confluence.append( half_dx *
                                                        representative_width )
                # Assuming a convergent network
                half_dx = ( lp.x[1] - lp.x[0] ) / 2.
                representative_width = lp.B[0]
                land_area_below_confluence = half_dx * representative_width

                lp.land_area_around_confluence = np.sum((
                                          np.sum(land_areas_above_confluence),
                                          land_area_below_confluence ))

            else:
                #print("Ahhhhmmmm.... why do you have no dx_ext?")
                lp.land_area_around_confluence = None

            ##print(lp.land_area_around_confluence)
            #raise ValueError("WHY OH WHY")


    def initialize(self,
                    config_file = None,
                    x_bl = None,
                    z_bl = None,
                    S0 = None,
                    Q_s_0 = None,
                    upstream_segment_IDs = None,
                    downstream_segment_IDs = None,
                    x = None,
                    z = None,
                    Q = None,
                    dQ = None,
                    B = None,
                    overwrite=False
                    ):
        #print( locals.keys() )
        """
        Run only once, at beginning of program.
        """

        # FOR NOW, RECORD S0, Q_s_0
        # Used in loop over solver
        self.S0 = S0
        self.Q_s_0 = Q_s_0

        #########################################
        # FIRST, CHECK IF A CONFIG FILE EXISTS. #
        # ENSURE NO DUPLICITY WITH PASSED VARS. #
        #########################################

        """
        # SOME SKETCHUP OF LOOPING OVER INTERNAL FCN VARIABLES, BUT THESE
        # UNFORTUNATELY DO NOT WORK WITHIN A CLASS.
        # FIGURE IT OUT LATER.

        def print_args(i,j,k):
            x = None
            for x in locals():
                print(locals()[x])

        item = None # preallocate so it isn't later added to keys()
        for item in locals():
            print( locals()[item] )

        """

        if config_file is not None:
            sys.exit("Code not yet set to work with a config file.")

        # FOR NOW, JUST LET CODE FAIL IF WE LEAVE TOO MANY THINGS AS NONE

        ############################################################
        # SECOND, IF LONG-PROFILE OBJECTS PASSED IN __INIT__, NOTE #
        ############################################################

        _build_segments = True
        if self.list_of_LongProfile_objects is not None:
            if overwrite:
                print("Overwriting prior network segments.")
            else:
                _build_segments = False

        ######################################
        # THIRD, INPUT AND BUILD THE NETWORK #
        ######################################

        # Required information:
        # x
        # z
        # Q
        # B
        # upstream_segment_IDs
        # downstream_segment_IDs
        if _build_segments:
            nseg = len(x)
            segments = []
            for i in range(nseg):
                segments.append( LongProfile() )
        # Class var; clunkier name
        self.list_of_LongProfile_objects = segments

        i = 0
        for lp in segments:
            # IDs and network-ID connections
            lp.set_ID(i)
            lp.set_upstream_segment_IDs( upstream_segment_IDs[i] )
            lp.set_downstream_segment_IDs( downstream_segment_IDs[i] )
            # x, z, Q
            # !!!!!!!!!!!!!!!!! NEEDS NETWROK INFO / MAYBE NOT NECESSARY
            # BECAUSE OF HOW IT REQUIRES X_EXT
            # lp.set_x( x = x[i], verbose=False )
            # TURN INTO FUNCTION ????? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            lp.x = x[i]
            lp.dx = np.diff(x[i])
                # LET'S CHANGE THE SIGN CONVENTION FOR S0??
                # !!!!!!!!!!!!
                # NOTING HERE BUT IT IS SET IN MULTIPLE OTHER PLACES
            # Setting variables directly.
            # "Setter" functions also calculate derived products
            # in a way that works for a single profile, but not for
            # a network.
            # STREAMLINE FUNCTIONS LATER -- UPDATE ONLY NAMED VAR?
            # OTHER FCN TO DO THE REST?
            lp.z = z[i]
            # !!!!!!!!!!!!!!!!!!!
            # Need to manage dQ_ext_2cell in network
            # TO DO HERE
            # lp.set_Q( Q = Q[i] )
            # The other GRLP-ey stuff
            lp.set_intermittency( 1 ) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # These all check out -- no problems with network.
            lp.basic_constants()
            lp.bedload_lumped_constants()
            # REMOVE??? PROBABLY. THIS SHOULD BE SET EXTERNALLY FOR A NETWORK.
            # lp.set_hydrologic_constants()
            #lp.set_z_bl(z1)
            lp.set_B( B = B[i] )
            # COULD WRITE A FUNCTION AROUND THIS
            # BUT I REALLY WANT TO REWRITE MORE IN TERMS OF SOURCE/SINK
            # DO SOMETHING HERE !!!!!
            lp.set_uplift_rate( 0 ) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            # FM:
            # Setting dQ of upstream segments. Used when dealing with junctions.
            # If dQ for each segment is provided, we look to the list for the
            # values of any upstream segments.
            # Otherwise, set to zero.
            lp.dQ_up_jcn = []
            for up_id in upstream_segment_IDs[i]:
                if dQ is not None:
                    lp.dQ_up_jcn.append( dQ[up_id] )
                else:
                    lp.dQ_up_jcn.append( 0. )
              
            i += 1

        # Generate list of all segment IDs to store within this Network object
        self.build_ID_list()

        ####################################
        #  THIRD: SET UP THE NETWORK X,Z   #
        # INTERIOR AND BOUNDARY CONDITIONS #
        ####################################

        # Required information:
        # x_bl
        # z_bl
        # Only one among:
        #   S0
        #   Q_s_0

        # Build the NetworkX topology graph (edges = segments, nodes =
        # junctions). Emerging source of truth for topology.
        self.build_graph()

        # Identify channel head and mouth segments (derived from the graph)
        self.create_list_of_channel_head_segment_IDs()
        self.create_list_of_channel_mouth_segment_IDs()

        # Set discharge within each segment
        self.update_Q( Q )

        # Boundary conditions: upstream sediment supply (S0 or Q_s_0) at the
        # channel heads, and downstream base level (z_bl) at the mouth. The
        # solver walks the topology and needs no padded x_ext / z_ext / Q_ext.
        self.set_Qs_input_upstream( S0 = S0, Q_s_0 = Q_s_0 )            # b.c.
        self.set_z_bl( z_bl )                                           # b.c.
        self.set_x_bl( x_bl )                                           # b.c.

        # Land area around each conflunece: Special case to help with dz/dt
        self.compute_land_areas_around_confluences()

        """
        # DEBUG
        lp = self.list_of_LongProfile_objects[0]
        i = 0
        print( lp.z_ext )
        print( lp.x_ext )
        print( np.abs( (lp.z_ext[i][2:] - lp.z_ext[i][:-2]) \
                 / lp.dx_ext_2cell[i] )**(1/6.) )
        """

    def evolve_threshold_width_river_network(self, nt=1, dt=3.15E7):
        """
        Solve the triadiagonal matrix through time, with a given
        number of time steps (nt) and time-step length (dt)
        """
        self.nt = nt
        self.dt = dt
        # The numerical engine lives in the solver module: it walks the topology
        # to assemble one global sparse system (multi-tributary confluences use
        # the conservative three-node junction cell) and Picard-iterates.
        return solver.evolve(self, nt, dt)

    def find_downstream_IDs(self, ID):
        """
        Return all segment IDs downstream of (and including) the specified
        segment, in order from that segment to the outlet.

        Walks the topology graph's unique out-edge chain.
        Added: FM, 03/2021. Migrated to the NetworkX graph.
        """
        IDs = [ID]
        node = self._edge_of_segment[ID][1]  # downstream node of this segment
        while self.graph.out_degree(node) > 0:
            (_, next_node, data), = list(self.graph.out_edges(node, data=True))
            IDs.append(data["segment_id"])
            node = next_node
        return IDs

    def find_upstream_IDs(self, ID):
        """
        Return all segment IDs upstream of (and including) the specified
        segment, by depth-first walk up the topology graph's in-edges.
        Added: FM, 03/2021. Migrated to the NetworkX graph.
        """
        IDs = []
        def _upstream_IDs(seg_id):
            IDs.append(seg_id)
            upstream_node = self._edge_of_segment[seg_id][0]
            for _, _, data in self.graph.in_edges(upstream_node, data=True):
                _upstream_IDs(data["segment_id"])
        _upstream_IDs(ID)
        return IDs

    def compute_absolute_lengths(self):
        """
        Return mean distance from source to mouth.
        Added: FM, 03/2021.
        """
        
        # First get downstream distance at outlet.
        x_max = (
            self.list_of_LongProfile_objects[self.channel_mouth_segment_ID].x[-1]
            )
            
        # Get average path lengths from heads to outlet
        Ls = []
        for i in self.list_of_channel_head_segment_IDs:
            Ls.append( x_max - self.list_of_LongProfile_objects[i].x[0] )
        self.mean_head_length = np.mean(Ls)
        
        # Get average path lengths from all points to outlet, weighted by dx.
        Ls = np.array([])
        dxs = np.array([])
        for seg in self.list_of_LongProfile_objects:
            Ls = np.append( Ls, x_max - seg.x )
            dxs = np.append( dxs, np.append( np.diff(seg.x),
                                             seg.x[-1] - seg.x[-2] ) )
        self.mean_length = np.sum(Ls*dxs)/np.sum(dxs)

    def compute_tokunaga_metrics(self):
        """
        Compute Tokunaga's metrics for the network.
        
        Assumes various things are already calculated, so only really works
        inside self.compute_network_properties().
        
        Put together following examples in Tarboton (1996, J. Hydrology) and
        Pelletier and Turcotte (2000). Seems to work at least for those. 
        """
        
        # Count numbers of streams of order i that join streams of order j
        N = np.zeros(( self.orders[-1], self.orders[-1]+1 ))
        for i in self.orders[:-1]:
            for stream in self.streams_by_order[i]:
                for ID in stream:
                    seg = self.list_of_LongProfile_objects[ID]
                    downID = seg.downstream_segment_IDs[0]
                    if downID not in stream:
                        down_seg = self.list_of_LongProfile_objects[downID]
                        adjacentID = [
                            id
                            for id in down_seg.upstream_segment_IDs
                            if id != ID
                            ][0]
                        adjacent_order = self.segment_orders[adjacentID]+1
                        N[i-1,adjacent_order-1] += 1

        # Get averages by dividing by number of streams j
        T = np.zeros(( self.orders[-1], self.orders[-1] ))
        for i in self.orders[1:]:
            T[:i-1,i-1] = N[:i-1,i] / self.order_counts[i]

        # Tokunaga's e_k - Average number of streams i flowing into streams
        # of i+k.
        # i.e. e_1 is the average of the numbers of order 1 streams flowing
        # into order 2 streams, order 2 streams flowing into order 3 streams,
        # etc.
        # Average is weighted by the numbers of receiving streams.
        e_k = np.zeros(max(self.orders))
        for k in range(1,max(self.orders)):
            weights = [c for c in self.order_counts.values()][k:]
            e_k[k] = sum(T.diagonal(k)*weights) / sum(weights)
                
        # Ratios of e_k / e_k-1
        # Again, average is weighted by number of receiver streams involved in
        # computing e_k and e_k-1.
        with np.errstate(invalid="ignore", divide="ignore"):
            K = e_k[2:] / e_k[1:-1]
        if len(K) > 0:
            weights = np.array([
                self.order_counts[o]+self.order_counts[o-1]
                    for o in self.orders[2:]
                ])
            K_mean = sum(K*weights)/sum(weights)
        else:
            K_mean = np.nan
        
        self.tokunaga = {'counts': N, 'average_counts': T, 'e_k': e_k,
            'K': K, 'K_mean': K_mean}

    def compute_topological_lengths(self):
        """
        Compute various metrics related to segment topological lengths.
        """
        
        # First compute topological length for each segment.
        # Number of segments to the outlet, including the segment itself.
        self.segment_topogical_lengths = []
        for seg in self.list_of_LongProfile_objects:
            self.segment_topogical_lengths.append(
                len(self.find_downstream_IDs(seg.ID))
                )
                
        # Get maximum and mean
        self.max_topological_length = max(self.segment_topogical_lengths)
        self.mean_topological_length = np.mean(self.segment_topogical_lengths)
        
        # Get mean for channel heads only
        self.head_topological_lengths = []
        for i in self.list_of_channel_head_segment_IDs:
            self.head_topological_lengths.append(
                self.segment_topogical_lengths[i]
                )
        self.mean_head_topological_length = np.mean(
            self.head_topological_lengths
            )
            
    def compute_topological_widths(self):
        """
        Compute topological widths, defined as the number of segments at a
        given topological distance from the outlet.
        """
        
        # First organise segments into distances from the outlet
        self.segs_by_topological_length = {}
        for i,seg in enumerate(self.list_of_LongProfile_objects):
            topo_length = len(self.find_downstream_IDs(seg.ID))
            if topo_length in self.segs_by_topological_length.keys():
                self.segs_by_topological_length[topo_length].append(seg.ID)
            else:
                self.segs_by_topological_length[topo_length] = [seg.ID]
        
        # Count numbers of segments at each distance from outlet
        self.topological_widths = {}
        for topo_length in self.segs_by_topological_length.keys():
            self.topological_widths[topo_length] = len(
                self.segs_by_topological_length[topo_length]
                )
        
        # Get maximum and mean
        self.max_topological_width = max(self.topological_widths.values())
        self.mean_topological_width = np.mean(
            list(self.topological_widths.values())
            )

        
    def compute_strahler_orders(self):
        """
        Assign segments Strahler orders; construct Strahler streams.
        """
        
        def _step_down(i):
            """
            Recursive function to step through network, assigning Strahler
            orders.
            """
            up_orders = [
                self.segment_orders[j]
                    for j in 
                        self.list_of_LongProfile_objects[i].upstream_segment_IDs
                ]
            if up_orders:
                self.segment_orders[i] = max(up_orders)
                if len(np.where(up_orders == max(up_orders))[0]) > 1:
                    self.segment_orders[i] += 1
            for j in self.list_of_LongProfile_objects[i].downstream_segment_IDs:
                _step_down(j)

        # compute strahler orders, working down from each source
        self.segment_orders = np.ones(len(self.IDs), dtype=int)
        for i in self.list_of_channel_head_segment_IDs:
            _step_down(i)

        # prepare dictionary to store Strahler streams
        self.streams_by_order = {}
        for order in np.unique(self.segment_orders):
            self.streams_by_order[order] = []
            
        # organise segments into streams            
        for seg in self.list_of_LongProfile_objects:
            seg_order = self.segment_orders[seg.ID]
            
            # identify neighbouring segments
            up_IDs = self.find_upstream_IDs(seg.ID)
            down_IDs = self.find_downstream_IDs(seg.ID)
            adjacent_IDs = np.hstack((up_IDs, down_IDs))
        
            # if no streams of that order yet, initiate list
            if not self.streams_by_order[seg_order]:
                self.streams_by_order[seg_order].append([seg.ID])
            
            # otherwise, check if the appropriate stream is already initiated
            # if so, append; if not initiate the list.
            else:    
                for i,stream in enumerate(self.streams_by_order[seg_order]):
                    if any(ID in stream for ID in adjacent_IDs):
                        self.streams_by_order[seg_order][i].append(seg.ID)
                if not seg.ID in np.hstack(self.streams_by_order[seg_order]):
                    self.streams_by_order[seg_order].append([seg.ID])

    def compute_horton_ratios(self):
        """
        Compute Horton's bifurcation, length and area(discharge) ratios.
        
        Horton proposes that the number of streams per stream order decreases
        by a roughly constant ratio as stream order increases. This ratio
        is termed the "bifurcation ratio". We can measue it using linear
        regression in linear-logarithmic space.
        
        Similarly, average stream length and average stream area increase
        steadily with increasing stream order. We compute them in the same
        way. We don't really have area, so we use discharge instead, to create
        a "discharge ratio".
        """
        
        # ---- BIFURCATION RATIO
        
        # count number of streams in each order
        self.order_counts = {}
        for order in self.streams_by_order.keys():
            self.order_counts[order] = len(self.streams_by_order[order])

        # bifurcation ratio
        self.orders = list(self.order_counts.keys())
        fit = np.polyfit(
            self.orders,
            np.log10([self.order_counts[o] for o in self.orders]),
            1)
        self.bifurcation_ratio = 10.**(-fit[0])
        self.bifurcation_scale = 10.**(fit[1] -fit[0])

        # ---- LENGTH RATIO

        # compute average stream lengths in each order
        self.stream_lengths = {}
        self.order_lengths = {}
        for o in self.orders:
            self.stream_lengths[o] = []
            for stream in self.streams_by_order[o]:
                l = 0.
                for segID in stream:
                    l += (
                        self.list_of_LongProfile_objects[segID].x.max() -
                        self.list_of_LongProfile_objects[segID].x.min()
                        )
                self.stream_lengths[o].append(l)
            self.order_lengths[o] = np.mean([l for l in self.stream_lengths[o]])

        # compute length ratio
        fit = np.polyfit(
            self.orders,
            np.log10([self.order_lengths[o] for o in self.orders]),
            1)
        self.length_ratio = 10.**(fit[0])
        self.length_scale = 10.**(fit[0] + fit[1])

        # ---- DISCHARGE RATIO

        # compute stream discharges
        self.stream_discharges = {}
        self.order_discharges = {}
        for o in self.orders:
            self.stream_discharges[o] = []
            for stream in self.streams_by_order[o]:
                self.stream_discharges[o] = [
                    self.list_of_LongProfile_objects[segID].Q.mean()
                    for segID in stream
                    ]
            self.order_discharges[o] = np.mean(
                [q for q in self.stream_discharges[o]]
                )

        # compute discharge ratio
        fit = np.polyfit(
            self.orders,
            np.log10([self.order_discharges[o] for o in self.orders]),
            1)
        self.discharge_ratio = 10.**fit[0]
        self.discharge_scale = 10.**(fit[0] + fit[1])
        
    def compute_jarvis_E(self):
        """
        Compute Jarvis's (1972, WRR) E metric.
        """
        
        heads_sum = 0
        interior_sum = 0
        
        for i,seg in enumerate(self.list_of_LongProfile_objects):
            if seg.ID in self.list_of_channel_head_segment_IDs:
                heads_sum += len(self.find_downstream_IDs(seg.ID))
            else:
                n_upstream_heads = [
                    i
                    for i in self.find_upstream_IDs(seg.ID)
                    if i in self.list_of_channel_head_segment_IDs
                    ]
                interior_sum += (
                    len(n_upstream_heads) * 
                    len(self.find_downstream_IDs(seg.ID))
                    )
                
        self.jarvis_E = interior_sum / heads_sum    

        
    def compute_mean_diffusivity(self):
        """
        Compute mean diffusivity (and related properties).
        
        We weight by dx, to avoid bias towards more densely sampled parts of
        the network. 
        """
        
        # slope and sediment discharge, by walking the topology (sets seg.S)
        self.compute_Q_s()
        
        # get dxs
        dxs = np.hstack(
            [np.append( np.diff(seg.x), seg.x[-1] - seg.x[-2] )
             for seg in self.list_of_LongProfile_objects]
            )
            
        # discharge
        Q_stack = np.hstack([seg.Q for seg in self.list_of_LongProfile_objects])
        self.mean_Q = (Q_stack * dxs).sum() / dxs.sum()
        
        # width
        B_stack = np.hstack([seg.B for seg in self.list_of_LongProfile_objects])
        self.mean_B = (B_stack * dxs).sum() / dxs.sum()
        
        # slope
        S_stack = np.hstack([seg.S for seg in self.list_of_LongProfile_objects])
        self.mean_S = (S_stack * dxs).sum() / dxs.sum()
        
        # diffusivity
        for seg in self.list_of_LongProfile_objects:
            seg.diffusivity = (7./6.) * seg.k_Qs * seg.intermittency * seg.Q * seg.S**(1./6.) \
                / seg.sinuosity**(7./6.) / seg.B / (1. - seg.lambda_p)
        diff_stack = np.hstack(
            [seg.diffusivity for seg in self.list_of_LongProfile_objects])
        self.mean_diffusivity = (diff_stack * dxs).sum() / dxs.sum()

    def compute_network_properties(self):
        """
        Compute various network properties.
        """
        self.compute_topological_lengths()
        self.compute_topological_widths()
        self.compute_absolute_lengths()
        self.compute_strahler_orders()
        self.compute_horton_ratios()
        self.compute_tokunaga_metrics()
        self.compute_jarvis_E()
        self.compute_mean_diffusivity()

    def find_hack_parameters(self):
        """
        Find optimal Hack parameters for this network.

        First get distances downstream from source for each network segment.
        Then optimise power law describing increasing discharge downstream.
        """
        d = []
        Q = []
        for lp in self.list_of_LongProfile_objects:
            up_IDs = self.find_upstream_IDs(lp.ID)
            ds = []
            for up_id in up_IDs:
                ds.append(min(self.list_of_LongProfile_objects[up_id].x))
            d.append(np.max(lp.x) - min(ds))
            Q.append(np.mean(lp.Q))
        k, p = _optimize_power_law(d, Q)
        return {'k': k, 'p': p, 'd': d, 'Q': Q}

    def find_hack_parameters_non_dim(self):
        """
        Non-dimensional Hack parameters for this network: power laws relating
        topological length to the numbers of upstream sources and segments.
        """
        self.create_list_of_channel_head_segment_IDs()
        sources = self.list_of_channel_head_segment_IDs
        topo_lengths = []
        upstream_sources = []
        upstream_segments = []
        for seg in self.list_of_LongProfile_objects:

            up_IDs = self.find_upstream_IDs(seg.ID)
            upstream_segments.append(len(up_IDs))
            up_sources = [ID for ID in up_IDs if ID in sources]
            upstream_sources.append(len(up_sources))

            topo_length = 0
            for s in up_sources:
                topo_length_i = 0
                down_IDs = self.find_downstream_IDs(s)
                for down_ID in down_IDs:
                    topo_length_i += 1
                    if down_ID == seg.ID:
                        break
                topo_length = max(topo_length, topo_length_i)
            topo_lengths.append(topo_length)

        k_src, p_src = _optimize_power_law(
            np.array(topo_lengths), np.array(upstream_sources))
        k_seg, p_seg = _optimize_power_law(
            np.array(topo_lengths), np.array(upstream_segments))
        return {'k_src': k_src, 'p_src': p_src, 'k_seg': k_seg, 'p_seg': p_seg,
                'lengths': topo_lengths, 'sources': upstream_sources,
                'segments': upstream_segments}

    def plot(self, show=True):
        """
        Build a planform plot of the network and (optionally) show it.

        Segments are laid out on alternating sides of the trunk, ordered by
        topological length, with overlaps resolved; returns the planform
        coordinates. Works on any Network, synthetic or DEM-derived. (Formerly
        the free function build_synthetic_network.plot_network.)
        """

        def check_for_segment_conflicts(ID, segs_by_topo_length, net, ys):
            """
            Check for segments that overlap with the given segment.
            """

            seg = net.list_of_LongProfile_objects[ID]
            x_max = seg.x.max()
            x_min = seg.x.min()
            y = ys[ID]

            for other_topo_length in segs_by_topo_length.keys():
                for other_ID in segs_by_topo_length[other_topo_length]:
                    if other_ID != ID:

                        other_seg = net.list_of_LongProfile_objects[other_ID]
                        other_x_max = other_seg.x.max()
                        other_x_min = other_seg.x.min()
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
                # Downstream connection point: the first node of the downstream
                # segment, or -- at the outlet -- the base-level ghost position
                if seg.downstream_segment_IDs:
                    x_down = net.list_of_LongProfile_objects[
                                 seg.downstream_segment_IDs[0]].x[0]
                elif seg.x_ghost_downstream is not None:
                    x_down = seg.x_ghost_downstream
                else:
                    x_down = 2*seg.x[-1] - seg.x[-2]
                if not seg.downstream_segment_IDs:
                    x = np.hstack(( seg.x, x_down, x_down ))/1000.
                    y = np.hstack(( np.full(len(seg.x),ys[i]), ys[i], ys[i] ))
                else:
                    x = np.hstack(( seg.x, x_down, x_down ))/1000.
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

        # ---- Topological lengths (sets self.max_topological_length)
        self.compute_topological_lengths()

        # ---- Organise segments by distance upstream (topological length)
        # Used later to check for conflicts between segments.
        segs_by_topo_length = {0: [], self.max_topological_length+2: []}
        for i in range(1,self.max_topological_length+2):
            segs_by_topo_length[i] = []
        for i,seg in enumerate(self.list_of_LongProfile_objects):
            topo_length = len(self.find_downstream_IDs(seg.ID))
            segs_by_topo_length[topo_length].append(seg.ID)

        # ---- Set up arrays to fill
        ys = np.full( len(self.list_of_LongProfile_objects), np.nan )
        sides = np.full( len(self.list_of_LongProfile_objects), 0)
        up_sides = np.full( len(self.list_of_LongProfile_objects), -1)
        connections = [
            [np.nan,np.nan] for i in range(len(self.list_of_LongProfile_objects))
            ]

        # ---- Loop over segments building planform
        for i,seg in enumerate(self.list_of_LongProfile_objects):

            # ---- Check if outlet
            if not seg.downstream_segment_IDs:
                ys[i] = 0.
                connections[i][1] = copy.copy(ys[i])

            # ---- Otherwise, add segment on to downstream one
            else:

                # Some info about the segment
                down_ID = seg.downstream_segment_IDs[0]
                topo_length = len(self.find_downstream_IDs(seg.ID))

                # Add segment based on relationship to downstream segment
                # Record what side of downstream segment the segment is on
                # Update "up_sides" so that next segment goes on the other side
                ys[i] = ys[down_ID] + up_sides[down_ID]
                sides[i] = copy.copy(up_sides[down_ID])
                up_sides[down_ID] *= -1

                # Check for conflict
                conflicting_id = check_for_segment_conflicts(
                    seg.ID, segs_by_topo_length, self, ys)

                # If there is a conflict, move downstream until reaching a segment
                # with the right direction to fix the conflict
                if conflicting_id:
                    seg_to_adjust = seg.ID
                    if seg.x.max() >= \
                        self.list_of_LongProfile_objects[conflicting_id].x.max():
                        while sides[conflicting_id] != sides[seg_to_adjust]:
                            seg_to_adjust = (
                                self.list_of_LongProfile_objects[seg_to_adjust]. \
                                downstream_segment_IDs[0])
                    else:
                        while sides[seg.ID] == sides[seg_to_adjust]:
                            seg_to_adjust = (
                                self.list_of_LongProfile_objects[seg_to_adjust]. \
                                downstream_segment_IDs[0])

                # Move everything upstream of that segment out the way until the
                # conflict is addressed
                while check_for_segment_conflicts(
                    seg.ID, segs_by_topo_length, self, ys):
                    up_IDs_down_ID = self.find_upstream_IDs(seg_to_adjust)
                    ys[up_IDs_down_ID] += sides[seg_to_adjust]

        # ---- Get everything starting from zero and positive
        ys -= ys.min() - 1.

        # ---- Create final planform
        planform = create_planform(self, ys)
        if show:
            plot_planform(planform)

        return planform





"""
# Test whether the lists have populated
for lp in self.list_of_LongProfile_objects: print(lp.x_ext)
print("")
for lp in self.list_of_LongProfile_objects: print(lp.z_ext)
"""
