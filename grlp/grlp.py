import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse import spdiags, identity, block_diag
from scipy import sparse
from scipy.sparse.linalg import spsolve, isolve
from scipy.stats import linregress
import warnings
import sys

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
        self.dx_ext = None
        self.dx_2cell = None
        self.Q_s_0 = None
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

    def set_x(self, x=None, x_ext=None, dx=None, nx=None, x0=None):
        """
        Set x directly or calculate it.
        Pass one of three options:
        x alone
        x_ext alone (this will also define x)
        dx, nx, and x0
        """
        if x is not None:
            # This doesn't have enough information to work consistently
            # Needs ext
            self.x = np.array(x)
            self.dx = np.diff(self.x)
            self.dx_2cell = self.x[2:] - self.x[:-2]
        elif x_ext is not None:
            self.x_ext = np.array(x_ext)
            self.x = x_ext[1:-1]
            self.dx_ext = np.diff(self.x_ext)
            self.dx_ext_2cell = self.x_ext[2:] - self.x_ext[:-2]
            self.dx_2cell = self.x[2:] - self.x[:-2]
            self.dx = np.diff(self.x)
        elif (dx is not None) and (nx is not None) and (x0 is not None):
            self.x = np.arange(x0, x0+dx*nx, dx)
            self.x_ext = np.arange(x0-dx, x0+dx*(nx+1), dx)
            self.dx = dx * np.ones(len(self.x) - 1)
            self.dx_ext = dx * np.ones(len(self.x) + 1)
            self.dx_2cell = np.ones(len(self.x) - 1)
            self.dx_ext_2cell = self.x_ext[2:] - self.x_ext[:-2]
        else:
            sys.exit("Need x OR x_ext OR (dx, nx, x0)")
        self.nx = len(self.x)
        self.L = self.x_ext[-1] - self.x_ext[0]
        if (nx is not None) and (nx != self.nx):
            warnings.warn("Choosing x length instead of supplied nx")

    def set_z(self, z=None, z_ext=None, S0=None, z1=0):
        """
        Set z directly or calculate it
        S0 = initial slope (negative for flow from left to right)
             unlike in the paper, this is a dz/dx value down the valley,
             so we account for sinuosity as well at the upstream boundary.
        z1 = elevation value at RHS
        """
        if z is not None:
            self.z = z
            self.z_ext = np.hstack((2*z[0]-z[1], z, 2*z[-1]-z[-2]))
        elif z_ext is not None:
            self.z_ext = z_ext
            self.z = z_ext[1:-1]
        elif self.x.any() and self.x_ext.any() and (S0 is not None):
            self.z = self.x * S0 + (z1 - self.x[-1] * S0)
            self.z_ext = self.x_ext * S0 + (z1 - self.x[-1] * S0)
            #print self.z_ext
        else:
            sys.exit("Error defining variable")
        #self.dz = self.z_ext[2:] - self.z_ext[:-2] # dz over 2*dx!

    def set_A(self, A=None, A_ext=None, k_xA=None, P_xA=None):
        """
        Set A directly or calculate it
        """
        if A is not None:
            self.A = A
            self.A_ext = np.hstack((2*A[0]-A[1], A, 2*A[-1]-A[-2]))
        elif A_ext is not None:
            self.A_ext = A_ext
            self.A = self.A_ext[1:-1]
        elif self.x.any() and self.x_ext.any():
            self.k_xA = k_xA
            if P_xA:
                self.P_xA = P_xA
            self.A_ext = self.k_xA * self.x_ext**self.P_xA
            self.A = self.k_xA * self.x**self.P_xA
        else:
            sys.exit("Error defining variable")

    def set_Q(self, Q=None, Q_ext=None, q_R=None, A_R=None, P_AQ=None,
              k_xQ=None, P_xQ=None, update_Qs_input=True):
        """
        Set Q directly or calculate it
        q_R = storm rainfall rate [m/hr]
        """
        if k_xQ is not None:
            self.k_xQ = k_xQ
        if Q is not None:
            # Check if it is a scalar or an array
            if hasattr(Q, "__iter__"):
                self.Q = Q
            else:
                # Assuming "x" is known already
                self.Q = Q * np.ones(self.x.shape)
            # Have to be able to pass Q_ext, created with adjacencies
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Q_ext = np.hstack( (2*self.Q[0]-self.Q[1],
                                self.Q,
                                2*self.Q[-1]-self.Q[-2]) )
        elif Q_ext is not None:
            self.Q = Q_ext[1:-1]
        elif q_R and A_R:
            if P_AQ:
                self.P_AQ = P_AQ
            q_R = q_R/3600. # to m/s
            Q_ext = q_R * np.minimum(A_R, self.A_ext) \
                    * (self.A_ext/np.minimum(A_R,self.A_ext))**self.P_AQ
            self.Q = Q_ext[1:-1]
        elif self.x.any() and self.x_ext.any() and k_xQ and P_xQ:
            self.Q = k_xQ * self.x**P_xQ
            Q_ext = k_xQ * self.x_ext**P_xQ
        else:
            sys.exit("Error defining variable")
        # dQ over 2*dx!
        # See Eq. D3 in Wickert & Schildgen (2019)
        # This then combines with the 1/4 factor in the coefficients
        # for the stencil that results from (2*dx)**2
        self.dQ = Q_ext[2:] - Q_ext[:-2]
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
        elif k_xB and self.x.any() and self.x_ext.any():
            self.B = k_xB * self.x**P_xB
            self.k_xB = k_xB
            self.P_xB = P_xB

    def set_uplift_rate(self, U):
        """
        Uplift rate if positive -- or equivalently, rate of base-level fall
        Subsidence (or base-level rise) accomplished by negative uplift
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
        self.downstream_fining_subsidence_equivalent = \
                - self.gravel_fractional_loss_per_km * self.Q_s \
                / ( (1-self.lambda_p) * self.B )

    def set_niter(self, niter=3):
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
        self.S0 = - np.sign(self.Q[0]) * self.sinuosity * \
                      ( np.abs(Q_s_0) / 
                        ( self.k_Qs * self.intermittency 
                              * np.abs(self.Q[0])) )**(6/7.)
        # Give upstream cell the same width as the first cell in domain
        self.z_ext[0] = self.z[0] - self.S0 * self.dx_ext[0]

    def update_z_ext_0(self):
        # Give upstream cell the same width as the first cell in domain
        self.z_ext[0] = self.z[0] - self.S0 * self.dx_ext[0]

    def compute_coefficient_time_varying(self):
        if self.S0 is not None:
            self.update_z_ext_0()
        self.dzdx_0_16 = np.abs( (self.z_ext[2:] - self.z_ext[:-2]) \
                         / self.dx_ext_2cell )**(1/6.)
        self.C1 = self.C0 * self.dzdx_0_16 * self.Q / self.B
        # Handling C1 for networked rivers
        # Need to link the two segments without skipping the channel head
        # DOESN'T SEEM TO CHANGE ANYTHING!
        # Looks right when both are 0! Any accidental inclusion of its own
        # ghost-node Qs,in?
        if len(self.downstream_segment_IDs) > 0:
            self.C1[-1] = self.C0[-1] \
                          * (np.abs(self.z_ext[-2] - self.z_ext[-1]) \
                                   /self.dx[-1])**(1/6.) \
                          * self.Q[-1] / self.B[-1]
        # This one matters! The above doesn't!!!! (Maybe.)
        # WORK HERE. If turns to 0, fixed. But why? Stays at initial profile?
        if len(self.upstream_segment_IDs) > 0:
            self.C1[0] = self.C0[0] \
                          * (np.abs(self.z_ext[1] - self.z_ext[0]) \
                                   /self.dx[0])**(1/6.) \
                          * self.Q[0] / self.B[0]

    def set_z_bl(self, z_bl):
        """
        Set the right-hand Dirichlet boundary conditions, i.e. the base level,
        given in the variable "z_bl" (elevation, base level)
        """
        self.z_bl = z_bl
        self.z_ext[-1] = self.z_bl

    def set_bcr_Dirichlet(self):
        self.bcr = self.z_bl * ( self.C1[-1] * 7/3. \
                       * (1/self.dx_ext[-2] + 1/self.dx_ext[-1])/2. \
                       + self.dQ[-1]/self.Q[-1] )

    def set_bcl_Neumann_RHS(self):
        """
        Boundary condition on the left (conventionally upstream) side of the
        domain.

        This is for the RHS of the equation as a result of the ghost-node
        approach for the Neumann upstream boundary condition with a prescribed
        transport slope.

        This equals 2*dx * S_0 * left_coefficients
        (2*dx is replaced with the x_ext{i+1} - x_ext{i-1} for the irregular
        grid case)
        """
        # Give upstream cell the same width as the first cell in domain
        # 2*dx * S_0 * left_coefficients
        self.bcl = self.dx_ext_2cell[0] * self.S0 * \
                            - self.C1[0] * ( 7/3./self.dx_ext[0]
                            - self.dQ[0]/self.Q[0]/self.dx_ext_2cell[0] )

    def set_bcl_Neumann_LHS(self):
        """
        Boundary condition on the left (conventionally upstream) side of the
        domain.

        This changes the right diagonal on the LHS of the equation using a
        ghost-node approach by defining a boundary slope that is calculated
        as a function of input water-to-sediment supply ratio.

        LHS = coeff_right at 0 + coeff_left at 0, with appropriate dx
              for boundary (already supplied)
        """
        self.right[0] = -self.C1[0] * 7/3. \
                         * (1/self.dx_ext[0] + 1/self.dx_ext[1])

    def evolve_threshold_width_river(self, nt=1, dt=3.15E7):
        """
        Solve the triadiagonal matrix through time, with a given
        number of time steps (nt) and time-step length (dt)
        """
        if (len(self.upstream_segment_IDs) > 0) or \
           (len(self.downstream_segment_IDs) > 0):
            warnings.warn("Unset boundary conditions for river segment"+
                          "in network.\n"+
                          "Local solution on segment will not be sensible.")
        self.nt = nt
        self.build_LHS_coeff_C0(dt)
        self.set_z_bl(self.z_bl)
        for ti in range(int(self.nt)):
            self.zold = self.z.copy()
            for i in range(self.niter):
                # If I want to keep this, will have to add to the networked
                # river too
                if self.gravel_fractional_loss_per_km is not None:
                    self.set_Sternberg_gravel_loss()
                self.build_matrices()
                self.z_ext[1:-1] = spsolve(self.LHSmatrix, self.RHS)
                #print self.bcl
            self.t += self.dt
            self.z = self.z_ext[1:-1].copy()
            self.dz_dt = (self.z - self.zold)/self.dt
            self.Qs_internal = 1/(1-self.lambda_p) * np.cumsum(self.dz_dt) \
                                * self.B + self.Q_s_0
            if self.S0 is not None:
                self.update_z_ext_0()

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
                    * self.dt / self.dx_ext_2cell

    def build_matrices(self):
        """
        Build the tridiagonal matrix (LHS) and the RHS matrix for the solution
        """
        self.compute_coefficient_time_varying()
        self.left = -self.C1 * ( (7/3.)/self.dx_ext[:-1]
                        - self.dQ/self.Q/self.dx_ext_2cell )
        self.center = -self.C1 * ( (7/3.) \
                              * (-1/self.dx_ext[:-1] \
                                 -1/self.dx_ext[1:]) ) \
                                 + 1.
        self.right = -self.C1 * ( (7/3.)/self.dx_ext[1:] # REALLY?
                        + self.dQ/self.Q/self.dx_ext_2cell )
        # Apply boundary conditions if the segment is at the edges of the
        # network (both if there is only one segment!)
        if len(self.upstream_segment_IDs) == 0:
            #print self.dx_ext_2cell
            self.set_bcl_Neumann_LHS()
            self.set_bcl_Neumann_RHS()
        else:
            self.bcl = 0. # no b.c.-related changes
        if len(self.downstream_segment_IDs) == 0:
            self.set_bcr_Dirichlet()
        else:
            self.bcr = 0. # no b.c.-related changes
        self.left = np.roll(self.left, -1)
        self.right = np.roll(self.right, 1)
        self.diagonals = np.vstack((self.left, self.center, self.right))
        self.offsets = np.array([-1, 0, 1])
        self.LHSmatrix = spdiags(self.diagonals, self.offsets, len(self.z),
                            len(self.z), format='csr')
        self.RHS = np.hstack(( self.bcl+self.z[0],
                               self.z[1:-1],
                               self.bcr+self.z[-1])) \
                               + self.ssd * self.dt \
                               + self.downstream_fining_subsidence_equivalent \
                                      *self.dt \
                               + self.U * self.dt

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
        print(P)
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

    def compute_Q_s(self):
        self.S = np.abs( (self.z_ext[2:] - self.z_ext[:-2]) /
                         (self.dx_ext_2cell) ) / self.sinuosity
        self.Q_s = -np.sign( self.z_ext[2:] - self.z_ext[:-2] ) \
                   * self.k_Qs * self.intermittency * self.Q * self.S**(7/6.)

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
        self.S = np.abs( (self.z_ext[2:] - self.z_ext[:-2]) \
                         / self.dx_ext_2cell )
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

    def compute_equilibration_time(self):
        """
        Compute valley equilibration time (sensu Paola et al., 1992).
        From scaling of linearized version of threshold width equation.
        """
        self.compute_diffusivity()
        self.equilibration_time = self.L**2. / self.diffusivity.mean()

    def compute_e_folding_time(self, n):
        """
        Compute valley e-folding times as function of wavenumber.
        From solving linearized version of threshold width equation.
        """
        self.compute_equilibration_time()
        return 1./self.diffusivity.mean()/self.wavenumber(n)**2.

    def compute_wavenumber(self, n):
        """
        Compute wavenumber for series solutions to linearized version of
        threshold width equation.
        """
        return (2*n + 1) * np.pi / 2. / self.L

    def compute_series_coefficient(self, n, period):
        """
        Compute coefficient for series solutions to linearized version of
        threshold width equation.
        """
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

        return lag


class Network(object):
    """
    Gravel-bed river long-profile solution builder and solver
    """

    def __init__(self, list_of_LongProfile_objects):
        """
        Instantiate the Network object with a list of Long Profile objects.
        Other funcitons will iterate over this network and its connectivity
        to create a full tridiagonal matrix solver.
        """
        self.list_of_LongProfile_objects = list_of_LongProfile_objects
        self.t = 0

    def build_ID_list(self):
        self.IDs = []
        for lp in self.list_of_LongProfile_objects:
            # IDs
            self.IDs.append(lp.ID)
        self.IDs = np.array(self.IDs)

    def build_block_diagonal_matrix_core(self):
        self.block_start_absolute = []
        self.block_end_absolute = []
        self.sparse_matrices = []
        #self.dx_downstream = [] # Should be in input, at least for now
        for lp in self.list_of_LongProfile_objects:
            # Absolute start and end list
            if len(self.block_start_absolute) > 0:
                self.block_start_absolute.append \
                     (self.block_end_absolute[-1])
            else:
                self.block_start_absolute.append(0)
            if len(self.block_end_absolute) > 0:
                self.block_end_absolute.append \
                     (self.block_end_absolute[-1] + lp.LHSmatrix.shape[0])
            else:
                self.block_end_absolute.append(lp.LHSmatrix.shape[0])
            # n-diagonal matrices
            self.sparse_matrices.append(lp.LHSmatrix)
        self.LHSblock_matrix = sparse.lil_matrix(
                                          block_diag(self.sparse_matrices) )
        self.block_start_absolute = np.array(self.block_start_absolute)
        self.block_end_absolute = np.array(self.block_end_absolute) - 1

    def add_block_diagonal_matrix_upstream_boundary_conditions(self):
        for lp in self.list_of_LongProfile_objects:
            for ID in lp.upstream_segment_IDs:
                # Space to edit
                col = self.block_end_absolute[self.IDs == ID][0]
                row = self.block_start_absolute[self.IDs == lp.ID][0]
                # Matrix entry, assuming net aligns with ids
                upseg = self.list_of_LongProfile_objects[ID]
                C0 = upseg.k_Qs * upseg.intermittency \
                        / ((1-upseg.lambda_p) * upseg.sinuosity**(7/6.)) \
                        * self.dt / (2 * lp.dx_ext[0])
                #C0 = upseg.C0[-1] # Should be consistent
                dzdx_0_16 = ( np.abs(lp.z_ext[1] - lp.z_ext[0])
                              / (lp.dx_ext[0]))**(1/6.)
                C1 = C0 * dzdx_0_16 * upseg.Q[-1] / lp.B[0]
                left_new = -C1 * 7/6. * 2 / lp.dx_ext[0]
                self.LHSblock_matrix[row, col] = left_new

    def add_block_diagonal_matrix_downstream_boundary_conditions(self):
        for lp in self.list_of_LongProfile_objects:
            for ID in lp.downstream_segment_IDs:
                # Space to edit
                col = self.block_start_absolute[self.IDs == ID][0]
                row = self.block_end_absolute[self.IDs == lp.ID][0]
                # Matrix entry, assuming net aligns with ids
                downseg = self.list_of_LongProfile_objects[ID]
                C0 = downseg.k_Qs * downseg.intermittency \
                        / ((1-downseg.lambda_p) * downseg.sinuosity**(7/6.)) \
                        * self.dt / (2 * lp.dx_ext[-1])
                dzdx_0_16 = ( np.abs(lp.z_ext[-2] - lp.z_ext[-1])
                              / (lp.dx_ext[0]))**(1/6.)
                C1 = C0 * dzdx_0_16 * lp.Q[-1] / downseg.B[0]
                right_new = -C1 * 7/6. * 2 / lp.dx_ext[-1]
                self.LHSblock_matrix[row, col] = right_new


    """
    def get_z_all(self):
        self.zold = []
        self.list_of_segment_lengths = []
        for lp in self.list_of_LongProfile_objects:
            self.zold.append(lp.z.copy())
            self.list_of_segment_lengths.append(len(lp.z))
        #self.zold = np.array(self.zold)
    """

    def get_z_lengths(self):
        self.list_of_segment_lengths = []
        for lp in self.list_of_LongProfile_objects:
            self.list_of_segment_lengths.append(len(lp.z))

    def stack_RHS_vector(self):
        self.RHS = []
        for lp in self.list_of_LongProfile_objects:
            self.RHS += list(lp.RHS)
        self.RHS = np.array(self.RHS)

    def set_niter(self, niter=3):
        # MAKE UNIFORM IN BASE CLASS
        self.niter = niter

    def update_zext(self):
        # Should just do this less ad-hoc
        for lp in self.list_of_LongProfile_objects:
            for ID in lp.downstream_segment_IDs:
                lp_downstream = np.array(self.list_of_LongProfile_objects) \
                                [self.IDs == ID][0]
                lp.z_ext[-1] = lp_downstream.z_ext[1]
            for ID in lp.upstream_segment_IDs:
                lp_upstream = np.array(self.list_of_LongProfile_objects) \
                                [self.IDs == ID][0]
                lp.z_ext[0] = lp_upstream.z_ext[-2]

    # Newly added
    def update_xext(self):
        # Should just do this less ad-hoc
        for lp in self.list_of_LongProfile_objects:
            for ID in lp.downstream_segment_IDs:
                lp_downstream = np.array(self.list_of_LongProfile_objects) \
                                [self.IDs == ID][0]
                lp.x_ext[-1] = lp_downstream.x_ext[1]
            for ID in lp.upstream_segment_IDs:
                lp_upstream = np.array(self.list_of_LongProfile_objects) \
                                [self.IDs == ID][0]
                lp.x_ext[0] = lp_upstream.x_ext[-2]
            # Update derived dx values 
            lp.dx_ext = np.diff(lp.x_ext)
            lp.dx_ext_2cell = lp.x_ext[2:] - lp.x_ext[:-2]


    def evolve_threshold_width_river_network(self, nt=1, dt=3.15E7):
        """
        Solve the triadiagonal matrix through time, with a given
        number of time steps (nt) and time-step length (dt)
        """
        # self.dt is decided earlier
        self.nt = nt
        self.dt = dt
        self.update_xext()
        self.update_zext()
        for ti in range(int(self.nt)):
            for lp in self.list_of_LongProfile_objects:
                lp.zold = lp.z.copy()
                lp.build_LHS_coeff_C0(dt=self.dt)
                lp.compute_coefficient_time_varying()
            for lp in self.list_of_LongProfile_objects:
                #print lp.C1
                lp.build_matrices()
            self.build_block_diagonal_matrix_core()
            self.add_block_diagonal_matrix_upstream_boundary_conditions()
            self.add_block_diagonal_matrix_downstream_boundary_conditions()
            # b.c. for no links
            """
            for lp in self.list_of_LongProfile_objects:
                if len(lp.upstream_segment_IDs) == 0:
                    lp.set_bcl_Neumann_LHS()
                    lp.set_bcl_Neumann_RHS()
                if len(lp.downstream_segment_IDs) == 0:
                    lp.set_bcr_Dirichlet()
            """
            self.stack_RHS_vector()

            for i in range(self.niter):
                self.update_zext()
                for lp in self.list_of_LongProfile_objects:
                    # Update coefficient for all: elements may call to others
                    # within the net
                    lp.compute_coefficient_time_varying()
                for lp in self.list_of_LongProfile_objects:
                    lp.build_matrices()
                # Update semi-implicit on boundaries
                # Commenting these two out helps solution!
                # Don't understand why. Perhaps error in code for them?
                #self.add_block_diagonal_matrix_upstream_boundary_conditions()
                #self.add_block_diagonal_matrix_downstream_boundary_conditions()
                out = spsolve(sparse.csr_matrix(self.LHSblock_matrix), self.RHS)
                i = 0
                idx = 0
                for lp in self.list_of_LongProfile_objects:
                    lp.z_ext[1:-1] = \
                                    out[idx:idx+self.list_of_segment_lengths[i]]
                    idx += +self.list_of_segment_lengths[i]
                    i += 1
            self.update_zext()
            self.t += self.dt # Update each lp z? Should make a global class
                              # that these both inherit from
            for lp in self.list_of_LongProfile_objects:
                lp.t = self.t
            i = 0
            idx = 0
            for lp in self.list_of_LongProfile_objects:
                lp.z = lp.z_ext[1:-1].copy()
                lp.dz_dt = (lp.z - lp.zold)/self.dt
                #lp.Qs_internal = 1/(1-lp.lambda_p) * np.cumsum(lp.dz_dt)*lp.B \
                #                 + lp.Q_s_0
                if lp.S0 is not None:
                    lp.update_z_ext_0()

    def find_downstream_IDs(self, ID):
        """
        Search list of downstream IDs, return all segments downstream of
        specified point.
        Uses recursive call of _downstream_IDs function.
        Added: FM, 03/2021.
        """
        IDs = []
        def _downstream_IDs(i):
            IDs.append(i)
            down_IDs = self.list_of_LongProfile_objects[i].downstream_segment_IDs
            if down_IDs:
                _downstream_IDs(down_IDs[0])
        _downstream_IDs(ID)
        return IDs

    def find_upstream_IDs(self, ID):
        """
        Search list of upstream IDs, return all segments upstream of
        specified point.
        Uses recursive call of _upstream_IDs function.
        Added: FM, 03/2021.
        """
        IDs = []
        def _upstream_IDs(i):
            IDs.append(i)
            up_IDs = self.list_of_LongProfile_objects[i].upstream_segment_IDs
            for j in up_IDs:
                _upstream_IDs(j)
        _upstream_IDs(ID)
        return IDs

    def find_sources(self):
        """
        Find network sources (heads).
        """
        self.sources = [i for i in self.IDs if self.list_of_LongProfile_objects[i].Q_s_0]

    def compute_mean_discharge(self):
        """
        Return mean discharge throughout network.
        Calculated "pathwise", so trunk segments counted multiple times.
        Interested in mean path from source to mouth.
        Added: FM, 03/2021.
        """
        Q_arr = np.array([])
        heads = [i for i in self.IDs if self.list_of_LongProfile_objects[i].Q_s_0]
        for i in self.sources:
            downstream_path = self.find_downstream_IDs(i)
            for j in downstream_path:
                Q_arr = np.hstack(( Q_arr, self.list_of_LongProfile_objects[j].Q))
        self.mean_discharge =  Q_arr.mean()

    def compute_mean_downstream_distance(self):
        """
        Return mean distance from source to mouth.
        Added: FM, 03/2021.
        """
        x_max = [self.list_of_LongProfile_objects[i].x[-1] for i in self.IDs if not self.list_of_LongProfile_objects[i].downstream_segment_IDs][0]
        x_arr = np.array([])
        for i in self.sources:
            x_arr = np.hstack(( x_arr, self.list_of_LongProfile_objects[i].x[0] ))
        # self.mean_downstream_distance = np.sqrt(((x_max - x_arr)**2).mean())
        self.mean_downstream_distance = (x_max - x_arr).mean()

    def compute_network_properties(self):
        """
        Compute various network properties.
        Added: FM, 03/2021.
        """

        self.find_sources()
        self.compute_mean_downstream_distance()
        self.compute_mean_discharge()

        # recursive function to step through network
        def _step_down(i):
            self.topological_length += 1
            up_orders = [self.segment_orders[j] for j in self.list_of_LongProfile_objects[i].upstream_segment_IDs]
            if up_orders:
                self.segment_orders[i] = max(up_orders)
                if len(np.where(up_orders == max(up_orders))[0]) > 1:
                    self.segment_orders[i] += 1
            for j in self.list_of_LongProfile_objects[i].downstream_segment_IDs:
                _step_down(j)

        # compute strahler orders, working down from each source
        self.segment_orders = np.zeros(len(self.IDs), dtype=int)
        self.max_topological_length = False
        for i in self.sources:
            self.topological_length = -1
            _step_down(i)
            self.max_topological_length = max(self.max_topological_length, self.topological_length)
        del self.topological_length

        # organise segments into streams
        self.streams_by_order = {}
        for order in np.unique(self.segment_orders):
            self.streams_by_order[order+1] = []
        for seg in self.list_of_LongProfile_objects:
            up_IDs = self.find_upstream_IDs(seg.ID)
            down_IDs = self.find_downstream_IDs(seg.ID)
            if not self.streams_by_order[self.segment_orders[seg.ID]+1]:
                self.streams_by_order[self.segment_orders[seg.ID]+1].append([seg.ID])
            else:
                for i,stream in enumerate(self.streams_by_order[self.segment_orders[seg.ID]+1]):
                    if any(ID in stream for ID in up_IDs) or any(ID in stream for ID in down_IDs):
                        self.streams_by_order[self.segment_orders[seg.ID]+1][i].append(seg.ID)
                if not seg.ID in np.hstack(self.streams_by_order[self.segment_orders[seg.ID]+1]):
                    self.streams_by_order[self.segment_orders[seg.ID]+1].append([seg.ID])

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
        
        # compute stream lengths
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
            self.order_discharges[o] = np.mean([q for q in self.stream_discharges[o]])
            
        # compute discharge ratio
        fit = np.polyfit(
            self.orders,
            np.log10([self.order_discharges[o] for o in self.orders]),
            1)
        self.discharge_ratio = 10.**fit[0]
        self.discharge_scale = 10.**(fit[0] + fit[1])    

        # # count number of streams in each order
        # self.order_counts = np.zeros(max(self.segment_orders)+1)
        # lengths = np.zeros(len(self.order_counts))
        # discharges = np.zeros(len(self.order_counts))
        # for i,seg in enumerate(self.list_of_LongProfile_objects):
        #     up_os = self.segment_orders[seg.upstream_segment_IDs]
        #     if self.segment_orders[i] not in up_os:
        #         self.order_counts[self.segment_orders[i]] += 1
        #     lengths[self.segment_orders[i]] += seg.x.max() - seg.x.min()
        #     discharges[self.segment_orders[i]] += seg.Q.mean()
        # self.order_lengths = lengths / self.order_counts / 1.e3
        # self.order_discharges = discharges / self.order_counts
        # 
        # # compute bifurcation ratios directly and estimate with log-fit
        # self.bifurcation_ratios = np.full(len(self.order_counts), np.nan)
        # for i in range(len(self.order_counts)-1):
        #     self.bifurcation_ratios[i] = (self.order_counts[i] / 
        #                                   self.order_counts[i+1])
        # fit = np.polyfit(
        #     np.arange(1,len(self.order_counts)+1,1),
        #     np.log10(self.order_counts),
        #     1)
        # self.bifurcation_ratio = 10.**(-fit[0])
        # self.bifurcation_intercept = fit[1]
        # 
        # # compute length ratios
        # self.length_ratios = np.full(len(self.order_counts), np.nan)
        # for i in range(1,len(self.order_counts)):
        #     self.length_ratios[i] = (self.order_lengths[i] / 
        #                              self.order_lengths[i-1])
        # fit = np.polyfit(
        #     np.arange(1,len(self.order_lengths)+1,1),
        #     np.log10(self.order_lengths),
        #     1)
        # self.length_ratio = 10.**(fit[0])
        # self.length_intercept = fit[1]
        # 
        # # compute discharge ratios
        # fit = np.polyfit(
        #     np.arange(1,len(self.order_discharges)+1,1),
        #     np.log10(self.order_discharges),
        #     1)
        # self.discharge_ratio = 10.**(fit[0])
