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
        self.x_ext = None
        self.dx_ext = None
        self.dx_2cell = None
        self.dx_ext_2cell = None
        self.z_ext = None
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
                print("Passing x alone leaves boundary conditions undefined."+
                      "\n"+
                      "Be sure to define these in order to set:"+
                      "\n"+
                      "x_ext"+
                      "\n"+
                      "dx_ext"+
                      "\n"+
                      "dx_ext_2cell"+
                      "\n"
                      )
            self.x = np.array(x)
            self.dx = np.diff(self.x)
            self.dx_2cell = self.x[2:] - self.x[:-2]
            # This doesn't have enough information to work consistently
            # Needs ext
            #self.x_ext = np.hstack( [ [self.x[0] - self.dx[0]],
            #                          self.x,
            #                          [self.x[-1] + self.dx[-1]] ] )
            #self.dx_ext = np.diff(self.x_ext)
            #self.dx_ext_2cell = self.x_ext[2:] - self.x_ext[:-2]
            # Apply nan and expect that the user will
            self.x_ext = np.hstack( [ np.nan,
                                      self.x,
                                      np.nan ] )
            self.dx_ext = np.diff(self.x_ext)
            self.dx_ext_2cell = self.x_ext[2:] - self.x_ext[:-2]
            # LCR
            self.dx_ext_2cell__left = self.dx_ext[:-1] - self.dx_ext_2cell
            self.dx_ext_2cell__cent = self.dx_ext[:-1] - self.dx_ext[1:]
            self.dx_ext_2cell__right = self.dx_ext[1:] - self.dx_ext_2cell
        if x_ext is not None:
            self.x_ext = np.array(x_ext)
            self.x = x_ext[1:-1]
            self.dx_ext = np.diff(self.x_ext)
            self.dx_ext_2cell = self.x_ext[2:] - self.x_ext[:-2]
            self.dx_2cell = self.x[2:] - self.x[:-2]
            self.dx = np.diff(self.x)
            # LCR
            self.dx_ext_2cell__left = self.dx_ext[:-1] - self.dx_ext_2cell
            self.dx_ext_2cell__cent = self.dx_ext[:-1] - self.dx_ext[1:]
            self.dx_ext_2cell__right = self.dx_ext[1:] - self.dx_ext_2cell
        elif (dx is not None) and (nx is not None) and (x0 is not None):
            self.x = np.arange(x0, x0+dx*nx, dx)
            self.x_ext = np.arange(x0-dx, x0+dx*(nx+1), dx)
            self.dx = dx * np.ones(len(self.x) - 1)
            self.dx_ext = dx * np.ones(len(self.x) + 1)
            self.dx_2cell = np.ones(len(self.x) - 1)
            self.dx_ext_2cell = self.x_ext[2:] - self.x_ext[:-2]
            # LCR. Though could go to a simpler grid with dx here, really
            self.dx_ext_2cell__left = self.dx_ext[:-1] - self.dx_ext_2cell
            self.dx_ext_2cell__cent = self.dx_ext[:-1] - self.dx_ext[1:]
            self.dx_ext_2cell__right = self.dx_ext[1:] - self.dx_ext_2cell
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

        Code works for single segments; the Network class manages z values
        on its own.
        """
        if z is not None:
            self.z = z
            if self.z_ext is None:
                self.z_ext = np.hstack((2*z[0]-z[1], z, 2*z[-1]-z[-2]))
        elif z_ext is not None:
            self.z_ext = z_ext
            if self.z is not None:
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

        Set only for a 1D array: Q is handled externally for river networks
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
        elif k_xB and self.x.any() and self.x_ext.any():
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
        self.downstream_fining_subsidence_equivalent = \
                - self.gravel_fractional_loss_per_km * self.Q_s \
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
        self.z_ext[0] = self.z[0] + self.S0 * self.dx_ext[0]

    def update_z_ext_0(self):
        """
        Give upstream cell the same width as the first cell in domain.
        Used only for GRLP alone (non-networked mode)

        Q1: z_ext and boundary conditions

        Q2: z_ext and C1
            self.C1 = self.C0 * dzdx_0_16 * self.Q / self.B
        And C0 has all local variables.
        """
        # Only one segment: towards applying boundary condition upstream
        self.z_ext[0] = self.z[0] + self.S0 * self.dx_ext[0]

    def compute_coefficient_time_varying(self):
        if self.S0 is not None:
            self.update_z_ext_0()
        # !!!C0!!!
        # NOT YET UPDATED
        # KEEPING self.dx_ext_2cell INSTEAD OF USING MORE PRECISE OPTION
        # Just used for 1-seg modeling
        dzdx_0_16 = np.abs( (self.z_ext[2:] - self.z_ext[:-2]) \
                         / self.dx_ext_2cell )**(1/6.)
        self.C1 = self.C0 * dzdx_0_16 * self.Q / self.B

    def network__compute_coefficient_time_varying(self):
        i=0
        # !!!!!!!!!!!!!!!!!!!!!
        # COUNT ON SCIENTIST ALREADY UPDATING Z_EXT
        #if self.S0 is not None:
        #    self.update_z_ext_0() # <-- Ursprüngliche Quelle
        # !!!C0!!!
        # NOT YET UPDATED
        # KEEPING self.dx_ext_2cell INSTEAD OF USING MORE PRECISE OPTION
        # But this is fine so long as all gradients may be approximated
        # to be linear.
        dzdx_0_16_list = []
        for _idx in range(len( self.z_ext) ):
            # ADD FUNCTIONALITY TO LOOP OVER Q AND WEIGHT BY IT
            dzdx_0_16_list.append(
                      np.abs( (self.z_ext[_idx][2:] - self.z_ext[_idx][:-2])
                                / self.dx_ext_2cell[_idx] )**(1/6.)
            )
        # FM: necessary? we compute again below

        C1_list = []
        for _idx in range(len( self.z_ext) ):
            # dzdx_2cell ? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # FM: could use the list from above?
            dzdx_0_16 = np.abs( (self.z_ext[_idx][2:] - self.z_ext[_idx][:-2]) \
                        / self.dx_ext_2cell[_idx] )**(1/6.)
            C1_list.append( self.C0 * dzdx_0_16 * self.Q / self.B )
        # ADW, 2024.08.10, rev. 2024.08.16
        # Previously, I summed the C1_list.
        # This created the problem that Fergus encountered, in which sediment
        # transport rates beyond the junction were doubled.
        self.C1 = np.mean(C1_list, axis=0)
        # By using "mean", we ensure that everything below the junction is
        # proper
        # The upstream-most C1 will be correct when there is just one stream
        # entering this lower one.
        # And when there are multiple streams, it is not used: the upstream
        # inputs are set separately.
        # Demonstrate this and make sure that it fails if I'm wrong.
        if len(self.upstream_segment_IDs) > 1:
            self.C1[0] = np.nan

    def set_z_bl(self, z_bl):
        """
        Set the right-hand Dirichlet boundary conditions, i.e. the base level,
        given in the variable "z_bl" (elevation, base level)

        For 1D single-segment mode, not network.
        """
        self.z_bl = z_bl
        self.z_ext[-1] = self.z_bl

    def set_x_bl(self, x_bl):
        self.x_bl = x_bl
        self.x_ext[-1] = self.x_bl

    def set_bcr_Dirichlet(self):
        # !!!C0!!!
        # UPDATED BUT JUST USING dx_ext_2cell
        #self.bcr = self.z_bl * ( self.C1[-1] * 7/3. \
        # z_bl not set. Hm.
        # Possibly because it is midway through the network
        # Maybe I really do need to update how I pass things here...
        # !!!!!!!!!!!!!!!! JUST MAKE SOMETHING RUN
        ##self.z_bl = 0

        #sys.exit("ERROR: UNSUPPORTED DQ_EXT_2CELL")
        if type(self.x_ext) is list:
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HACk to make it work for now
            # !!!! NO LONGER BEING USED FOR NETWORK; CHANGE BACK FOR SINGLE CHANNEL
            # OR JUST REMOVE
            # IS USED, FOR RIVER MOUTH. BUT MAYBE IT DOESN'T MATTER?
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.bcr = self.z_bl * self.C1[-1] / self.dx_ext_2cell[0][-1] * (
                         7/3. * ( 1/(self.x_ext[0][-1] - self.x_ext[0][-2]) )
                         +
                         self.dQ_ext_2cell[0][-1]/self.Q[-1]
                         * 1/self.dx_ext_2cell[0][-1]
                         )
                      # !!!!!!!!!!!!!!!!!
                      # dQ_ext_2cell --> dQ_ext_2cell[0]. Expecting list! HACK.
        else:
            # Why do I average over the last two dx here?
            self.bcr = self.z_bl * self.C1[-1] / self.dx_ext_2cell[-1] * (
                           7/3. * ( 1/self.dx_ext[-2] + 1/self.dx_ext[-1])/2.
                           + self.dQ_ext_2cell[-1]/self.Q[-1]
                           * 1/ self.dx_ext_2cell[-1]
                           )
                           #+ self.dQ_ext_2cell[-1]/self.Q[-1] )

         # I HAVE NOT CHECKED WHY PREV CODE DIDN'T HAVE A DX IN THE LAST LINE
         # WITH DQ/Q
         # AND NOW I'VE JUST INCLUDED THE SAME DQ, BUT IN JUST ONE CELL
         # INSTEAD OF TWO.
         # (WHEN I CHANGED TO "UPWIND". JUST CHANGED BACK, AND STILL DON'T KNOW.
         # IN ANY CASE, THE TESTS INCLUDE Z_BL = 0, SO BCR = 0 AND THIS CAN'T
         # BE THE SOURCE OF AN ERROR.

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
        # !!!C0!!!
        # UPDATED BUT JUST USING dx_ext_2cell
        #self.bcl = self.dx_ext_2cell[0] * self.S0 * \
        #                    -self.C1[0] * ( 7/3./self.dx_ext[0]
        #                    - self.dQ_ext_2cell[0]/self.Q[0]/self.dx_ext_2cell[0] )
        # !!!C0!!!
        # Probably not so easy to update as just updating C1
        # BECAUSE IT IS CHANGING THE RHS
        if type(self.x_ext) is list:
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HACk to make it work for now
            # !!!! NO LONGER BEING USED FOR NETWORK; CHANGE BACK FOR SINGLE CHANNEL
            # OR JUST REMOVE
            # YES! IS USED, WHEN THERE ARE NO UPSTREAM SEGMENTS
            # SHOULD IT STILL BE, OR SHOULD NEW METHOD BE USED?
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #self.bcl = self.dx_ext_2cell[0][0] * -self.S0 * \
            #                    self.C1[0] / self.dx_ext_2cell[0][0] \
            #                    * ( 7/3./self.dx_ext[0][0]
            #                    - self.dQ_ext_2cell[0][0]/self.Q[0]/self.dx_ext_2cell[0][0] )
            # RESTORED OLD SIGN CONVENTION HERE: POSITIVE EVEN THOUGH
            # IT REASONABLY SHOULD BE FLIPPED
            # CONSIDER BREAKING THIS AT SOME LATER POINT.
            self.bcl = self.S0 * \
                                self.C1[0] \
                                * ( 7/3./self.dx_ext[0][0]
                                - self.dQ_ext_2cell[0][0]/self.Q[0]/self.dx_ext_2cell[0][0] )
                                # !!!!!!!!!!!!!!!!!
                                # dQ_ext_2cell --> dQ_ext_2cell[0]. Expecting list! HACK.
        else:
            #self.bcl = self.dx_ext_2cell[0] * -self.S0 * \
            #                    -self.C1[0] / self.dx_ext_2cell[0] \
            #                    * ( 7/3./self.dx_ext[0]
            #                    - self.dQ_ext_2cell[0]/self.Q[0]/self.dx_ext_2cell[0] )
            # RESTORED OLD SIGN CONVENTION HERE: POSITIVE EVEN THOUGH
            # IT REASONABLY SHOULD BE FLIPPED
            # CONSIDER BREAKING THIS AT SOME LATER POINT.
            self.bcl = self.S0 * \
                                self.C1[0] \
                                * ( 7/3./self.dx_ext[0]
                                - self.dQ_ext_2cell[0]/self.Q[0]/self.dx_ext_2cell[0] )

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
        # !!!C0!!!
        # UPDATED BUT JUST USING dx_ext_2cell
        #self.right[0] = -self.C1[0] * 7/3. \
        if type(self.x_ext) is list:
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HACK TEST
            # NO LONGER IN USE!
            # raise ValueError('Ooh!')
            # YES! IS USED, WHEN THERE ARE NO UPSTREAM SEGMENTS
            # SHOULD IT STILL BE, OR SHOULD NEW METHOD BE USED?
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.right[0] = -self.C1[0] / self.dx_ext_2cell[0][0] * 7/3. \
                             * (1/self.dx_ext[0][0] + 1/self.dx_ext[0][1])
        else:
            self.right[0] = -self.C1[0] / self.dx_ext_2cell[0] * 7/3. \
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
                    * self.dt

    def build_matrices(self):
        """
        Build the tridiagonal matrix (LHS) and the RHS matrix for the solution
        """
        self.compute_coefficient_time_varying()
        # !!!C0!!!
        # UPDATED WITH STRAIGHT self.dx_ext_2cell
        self.left = -self.C1 / self.dx_ext_2cell \
                        * ( (7/3.)/self.dx_ext[:-1]
                        - self.dQ_ext_2cell/self.Q/self.dx_ext_2cell )
        self.center = -self.C1 / self.dx_ext_2cell \
                              * ( (7/3.)
                              * (-1/self.dx_ext[:-1]
                                 -1/self.dx_ext[1:]) ) \
                                 + 1.
        self.right = -self.C1 / self.dx_ext_2cell \
                              * ( (7/3.)/self.dx_ext[1:] # REALLY?
                                  + self.dQ_ext_2cell/self.Q/self.dx_ext_2cell )
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

    def network__build_matrix_inner(self):
        """
        Build the tridiagonal matrix (LHS) and the RHS matrix for the solution
        """
        self.network__compute_coefficient_time_varying()
        # !!!C0!!!
        # UPDATED WITH STRAIGHT self.dx_ext_2cell
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # NEXT TO DO! SET THESE UP FOR NETWORK OR SINGLE SEGMENT
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # THIS IS A TOTAL HACK, BUT IT SHOULD ACTUALLY BUILD SUCH
        # INNER MATRICES BECAUSE IT NEGLECTS THE BOUNDARY
        # CONDITIONS' IMPACTS, OR SO I GUESS.
        # !!!!!!!!!!!!!!!!!! TO CHECK !!!!!!!!!!!!!!!!!!!!!!!

        # For all nodes that are fully internal to the segment
        # 0s should be fine: internal members okay, and all external
        # should be shared downstream (yep) or split out later for upstream

        # THESE LOOK OKAY
        #self.C0 = self.k_Qs * self.intermittency \
        #            / ((1-self.lambda_p) * self.sinuosity**(7/6.)) \
        #            * self.dt
        #self.C1 = self.C0 * dzdx_0_16 * self.Q / self.B
        # NEXT TEST

        self.left = -self.C1 / self.dx_ext_2cell[0] \
                        * ( (7/3.)/self.dx_ext[0][:-1]
                        - self.dQ_ext_2cell[0]/self.Q/self.dx_ext_2cell[0] )
        self.center = -self.C1 / self.dx_ext_2cell[0] \
                              * ( (7/3.)
                              * (-1/self.dx_ext[0][:-1]
                                 -1/self.dx_ext[0][1:]) ) \
                                 + 1.
        self.right = -self.C1 / self.dx_ext_2cell[0] \
                              * ( (7/3.)/self.dx_ext[0][1:] # REALLY?
                                  + self.dQ_ext_2cell[0]/self.Q/self.dx_ext_2cell[0] )
                                  # !!!!!!!!!!!!!!!!!
                                  # dQ --> dQ[0]. Expecting list! HACK.
        # Far-left "self.center" depends on upstream boundary conditions
        # This is solved by discretizing the gradients and discharges abnove
        # and below the upstream-most node (i.e., tributary junction)
        # for each tributary and for the river reach immediately downstream
        # of the junction.
        # This applies only if there are multiple tributaries. Otherwise,
        # the existing "self.center" should be correct.
        # Although it might be set up differently (straight across junction)
        # Actually, could use this in general, though don't have to.
        # as in, for any number of tribs.
        # Maybe that would be cleaner, even if it provides a different
        # way of making boundary calculations

        """
        # dQs/dx = self.C0 * d/dx (Q S^(7/6))
        self.C0 = self.k_Qs * self.intermittency \
                    / ((1-self.lambda_p) * self.sinuosity**(7/6.)) \
                    * self.dt
        self.C1 = self.C0 * dzdx_0_16 * self.Q / self.B
        """

        # Perhaps add a separate function in the future to do this
        # only once

        # START HERE AFTER MAKING CONFLUENCE AREA FUNCTION

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!
        # ADD NONLINEARITY!
        #if len(self.dx_ext) > 1: # Or even 1
        # Tributary contributions to center: Write to work for any
        # integer > 0 number of tributaries

        if len(self.downstream_segment_IDs) > 0:
            # SHOULD DO ANOTHER CHECK: ACTUALLY TRIB JCN?
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # Update new trib at upseg[-1]
            # Could do with len z_ext or something like this
            # Index -1: Before roll
            self.left[-1] = -self.C1[-1] / self.dx_ext_2cell[0][-1] \
                            * ( (7/3.)/self.dx_ext[0][-2]
                            - (self.Q[-1]-self.Q[-2])/self.Q[-1]/self.dx[-1] )
            self.right[-1] = -self.C1[-1] / self.dx_ext_2cell[0][-1] \
                                  * ( (7/3.)/self.dx_ext[0][-1] # REALLY?
                                      + (self.Q[-1]-self.Q[-2])/self.Q[-1]/self.dx[-1] )


        if len(self.upstream_segment_IDs) > 0:
            # Update new trib at downseg[0]
            _trib_cent = 0.
            self.upseg_trib_coeffs = []
            ##print("ID", self.ID)
            for _tribi in range( len(self.upstream_segment_IDs) ):
                # Slope for nonlinear portion
                # This indexing works assuming only 1 downstream segment
                # (i.e., convergent network)
                # or in general if we keep the inner stride as the upstream
                # segments
                dzdx_0_16 = ( np.abs( self.z_ext[_tribi][0] -
                                      self.z_ext[_tribi][1] )
                              / self.dx_ext[_tribi][0] )**(1/6.)
                #dzdx_0_16 = 1 # DEBUG TEST
                #print ( self.dx_ext )

                # All of C1 except for C0 is included here
                # It varies based on trib junction
                # _trib_coeff = 1 * \
                # 1E0 to play with coefficients and check them
                # AW: 2024.08.10; re-implemented cleanly into trunk 2024.08.16
                # I added dQ_dx to the discharges to account for along-stream
                # additions of water at the tributary junctions.
                # Previously, this was always 0 -- meaning that it was just
                # the two tributaries coming together.
                # Currently, a kludge: user sets self.dQ_up_jcn for each
                # segment. Later could integrate upstream and downstream
                # discharges to calculate internally.
                # Dividing this discharge by len(self.Q_ext) because in the
                # current implementation, the added water is spread evenly
                # across all tributaries. This can also change in an updated
                # version of the code.
                _trib_coeff = dzdx_0_16 * 1E0 * \
                              ( (self.Q_ext[_tribi][0] + \
                                    self.dQ_up_jcn[_tribi] / len(self.Q_ext) ) /
                                self.dx_ext[_tribi][0] ) \
                              / self.land_area_around_confluence
                """
                print("_tribi", _tribi)
                print("TC", _trib_coeff)
                print(dzdx_0_16)
                print(self.Q_ext[_tribi][0])
                print(self.dx_ext[_tribi][0])
                print("")
                """
                #_trib_coeff *= -1
                _trib_cent += _trib_coeff
                #print("T", self.upstream_segment_IDs[_tribi],
                #            self.C0 * _trib_coeff)
                # Svae this value to transmit to upstream parts of matrix
                self.upseg_trib_coeffs.append( self.C0 * _trib_coeff )
            #print(_tribi) # Yep: 2 streams

            # All of C1 except for C0 is included here
            # This is the value for downstream of the confluence
            # Using z_ext so it updates while iterating
            # _TRIBI, 0, ... DOESN'T REALLY MATTER
            dzdx_0_16 = ( np.abs( self.z_ext[0][1] -
                                  self.z_ext[0][2] )
                          / self.dx[0] )**(1/6.)
            #dzdx_0_16 = 1 # DEBUG TEST
            #self.DEBUG_dzdx_0_16__downstream = dzdx_0_16.copy()
            # Positive for right, negative for center
            #_mainstem_cent_right = 1 * \
            # 1E0 to play with coefficients and check them
            _mainstem_cent_right = dzdx_0_16 * 1E0 * \
                                   (self.Q[0] + self.Q[1])/2. / self.dx[0] \
                                   / self.land_area_around_confluence

                                   # Replacing this with just exactly the Q
                                   # at the confluence
                                   #(self.Q[0] + self.Q[1])/2. / self.dx[0] \

            # Changes based on tributaries having more or less Q each!
            #print("")
            #print("TRIB CENT", _trib_cent)
            #print("")
            self.center[0] = self.C0 * ( _trib_cent + _mainstem_cent_right ) \
                              + 1

            """
                    C0 = upseg.k_Qs * upseg.intermittency \
                            / ((1-upseg.lambda_p) * upseg.sinuosity**(7/6.)) \
                            * self.dt
                    # From upseg, could be either dx_ext, so why not 0 : )
                    dzdx_0_16 = ( np.abs( lp.z_ext[_relative_id][0] -
                                          lp.z_ext[_relative_id][1] )
                                  / lp.dx_ext[_relative_id][0] )**(1/6.)
                    C1 = C0 * dzdx_0_16 * upseg.Q[-1] / lp.B[0]
                    # Slight hack but will work in convergent network
                    #left_new = C0 / upseg.dx_ext[0][-1]
                    left_new = C1 / lp.dx_ext[_relative_id][0] \
                                / (lp.dx_ext_2cell[_relative_id][0]/2.)
                    self.LHSblock_matrix[row, col] = left_new
                    _relative_id += 1
            """


            # Right needs to be changed too
            # This should be positive... but I get something that looks right
            # when I make it negative
            # ???????????????????????????????????????????
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.right[0] = - self.C0 * _mainstem_cent_right
            #print("R", self.right[0])
            # But gradient keeps decreasing and should be constant
            # Is sediment being transmitted across properly via
            # slope * discharge?

            # Left will be handled among boundary conditions

            """
            print("! 0 IF WE CONSERVE MASS !")
            print( self.right[0] + self.center[0] - 1 + np.sum(self.upseg_trib_coeffs) )

            print("L-R balance!")
            print( np.sum(self.upseg_trib_coeffs) - self.right[0] )
            """

        # As long as the network is convergent and based on Dirichlet boundary
        # conditions on the downstream end, the single downstream segment
        # should suffice with no modifications

        # Apply boundary conditions if the segment is at the edges of the
        # network (both if there is only one segment!)
        # Change this only if we decide to use the new junction approach
        # for everything
        if len(self.upstream_segment_IDs) == 0:
            #print self.dx_ext_2cell
            self.set_bcl_Neumann_LHS()
            self.set_bcl_Neumann_RHS()
        else:
            self.bcl = 0. # no b.c.-related changes
        if len(self.downstream_segment_IDs) == 0:
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # For river mouth.
            # Keep this????
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.set_bcr_Dirichlet()
        else:
            self.bcr = 0. # no b.c.-related changes
        self.left = np.roll(self.left, -1)
        #print(self.right)
        self.right = np.roll(self.right, 1)
        #print(self.right)
        self.diagonals = np.vstack((self.left, self.center, self.right))
        self.offsets = np.array([-1, 0, 1])
        self.LHSmatrix = spdiags(self.diagonals, self.offsets, len(self.z),
                            len(self.z), format='csr')
        #print("Diag")
        #print(self.diagonals)
        #print("Array)
        #plt.figure(); plt.imshow(np.log10(self.LHSmatrix.todense()), interpolation='nearest')
        #plt.show()

        self.RHS = np.hstack(( self.bcl+self.zold[0],
                               self.zold[1:-1],
                               self.bcr+self.zold[-1])) \
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

    def compute_Q_s(self):
        S = []
        Q_s = []
        # Ensure that this function works even if there is no list involved
        if type(self.z_ext) is np.ndarray:
            z_ext = [ self.z_ext ]
        else:
            dx_ext_2cell = [ dx_ext_2cell ]
        if type(self.dx_ext_2cell) is np.ndarray:
            dx_ext_2cell = list(self.dx_ext_2cell)
        else:
            dx_ext_2cell = self.dx_ext_2cell
        # Next, loop through slopes and sediment discharges
        for _z, _dx in zip(z_ext, dx_ext_2cell):
            S.append( np.abs( (_z[2:] - _z[:-2]) / _dx) / self.sinuosity )
            Q_s.append(
                -np.sign( _z[2:] - _z[:-2] ) \
                * self.k_Qs * self.intermittency * self.Q * S[-1]**(7/6.)
                )
        self.S = np.mean(S, axis=0)
        self.Q_s = np.mean(Q_s, axis=0)

        # # old, non-network way
        # self.S = np.abs( (self.z_ext[0][2:] - self.z_ext[0][:-2]) /
        #                  (self.dx_ext_2cell[0]) ) / self.sinuosity
        # self.Q_s = -np.sign( self.z_ext[0][2:] - self.z_ext[0][:-2] ) \
        #            * self.k_Qs * self.intermittency * self.Q * self.S**(7/6.)

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

    def compute_length(self):
        """
        Compute total segment length.
        Average over external arrays for each tributary.
        """
        Ls = []
        for x in self.x_ext:
            Ls.append(x.max() - x.min())
        self.L = np.mean(Ls, axis=0)

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
        return 1./self.diffusivity.mean()/self.wavenumber(n)**2.

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

    def map_block_diagonal_matrix_blocks(self):
        """
        Obtain the (row, column) indices for the cells at the beginnings
        and ends of blocks in the block-diagonal matrix.

        These are "absolute" -- that is, within the full block-diagonal
        matrix, rather than simply being indexed within a single block.

        Because the block-diagonal coefficient matrix and each of its blocks
        are square, we need only record one index for each. The index repeats.
        """
        self.block_start_absolute = []
        self.block_end_absolute = []

        for lp in self.list_of_LongProfile_objects:
            # If a start_absolute has already been defined, then the next
            # end_absolute is the cell immediately preceding.
            if len(self.block_start_absolute) > 0:
                self.block_start_absolute.append \
                     (self.block_end_absolute[-1])
            # Otherwise, we are at the very beginning and it is time to define
            # the first start_absolute at 0
            else:
                self.block_start_absolute.append(0)
            # If a start_absolute has already been defined, stride to the
            # next block and mark the location of its first cell.
            # Each LHSmatrix is square, as is the overall block-diagonal
            # matrix.
            if len(self.block_end_absolute) > 0:
                self.block_end_absolute.append \
                     (self.block_end_absolute[-1] + lp.LHSmatrix.shape[0])
            # Finally, if we are at the very start (before any block-end
            # position has been recorded), then note the position of the
            # end of the first block.
            else:
                self.block_end_absolute.append(lp.LHSmatrix.shape[0])
            # In the above definitions for block_end_absolute, we ignore
            # 0-indexing, which we correct below.
        # Create numpy arrays of the start, end indices, correcting the
        # offset at the end.
        self.block_start_absolute = np.array(self.block_start_absolute)
        self.block_end_absolute = np.array(self.block_end_absolute) - 1

    def create_block_diagonal_matrix_with_internal_tridiagonals(self):
        """
        Add the internal set of tridiagonal-matrix values to each of the
        blocks within the block-diagonal matrix structure.

        Currently, this assumes that a LHSmatrix has already been developed
        for each block. Considering that this has been built using a
        function orginally built for a non-network structure,
        this may be worth revisiting!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        """
        sparse_blocks = []
        for lp in self.list_of_LongProfile_objects:
            # n-diagonal matrices
            sparse_blocks.append(lp.LHSmatrix)
        # Assemble
        self.LHSblock_matrix = sparse.lil_matrix(
                                          block_diag(sparse_blocks) )

    def add_block_diagonal_matrix_upstream_boundary_conditions(self):
        """
        Add internal upstream boundary conditions
        """
        for lp in self.list_of_LongProfile_objects:
            lp.later_trib_coeffs = [] # DEBUG
            if len(lp.upstream_segment_IDs) == 0:
                pass
            else:
                _relative_id = 0
                for ID in lp.upstream_segment_IDs:
                    # Space to edit
                    col = self.block_end_absolute[self.IDs == ID][0]
                    row = self.block_start_absolute[self.IDs == lp.ID][0]
                    # Matrix entry, assuming net aligns with ids
                    upseg = self.list_of_LongProfile_objects[ID]
                    # OLD
                    #C0 = upseg.k_Qs * upseg.intermittency \
                    #        / ((1-upseg.lambda_p) * upseg.sinuosity**(7/6.)) \
                    #        * self.dt / (2 * lp.dx_ext[0])
                    # !!!C0!!!
                    # I've updated dx_ext_2cell to acknowledge boundaries
                    # This should therefore be used here.
                    """
                    # Commented out 19 Sept 2023
                    C0 = upseg.k_Qs * upseg.intermittency \
                            / ((1-upseg.lambda_p) * upseg.sinuosity**(7/6.)) \
                            * self.dt / lp.dx_ext_2cell[_relative_id][0]
                    # From earlier
                    #C0 = upseg.C0[-1] # Should be consistent
                    # !!!!!!!!!!!!!!!!!!!! [0] z_ext
                    dzdx_0_16 = ( np.abs(lp.z_ext[0][1] - lp.z_ext[0][0])
                                  / (lp.dx_ext[_relative_id][0]))**(1/6.)
                    C1 = C0 * dzdx_0_16 * upseg.Q[-1] / lp.B[0]
                    left_new = -C1 * 7/6. * 2 / lp.dx_ext[_relative_id][0]
                    self.LHSblock_matrix[row, col] = left_new
                    _relative_id += 1
                    """

                    """
                    C0 = upseg.k_Qs * upseg.intermittency \
                            / ((1-upseg.lambda_p) * upseg.sinuosity**(7/6.)) \
                            * self.dt
                    # From upseg, could be either dx_ext, so why not 0 : )
                    dzdx_0_16 = ( np.abs( lp.z_ext[_relative_id][0] -
                                          lp.z_ext[_relative_id][1] )
                                  / lp.dx_ext[_relative_id][0] )**(1/6.)
                    C1 = C0 * dzdx_0_16 * upseg.Q[-1] / \
                                    lp.land_area_around_confluence
                    # Slight hack but will work in convergent network
                    #left_new = C0 / upseg.dx_ext[0][-1]
                    left_new = C1 / lp.dx_ext[_relative_id][0]
                    self.LHSblock_matrix[row, col] = left_new
                    #if _relative_id == 0:
                    print("L", upseg.ID, left_new)
                    _relative_id += 1
                    """

                    """
                    # THIS BASICALLY SEEMED TO WORK, BUT I WONDER...
                    # I keep getting different values from the calcs above
                    # in some places
                    # but not others
                    # and I don't know why
                    # and it makes a total mess
                    # So here.
                    self.LHSblock_matrix[row, col] = \
                                    lp.upseg_trib_coeffs[_relative_id]
                    print("L", upseg.ID, lp.upseg_trib_coeffs[_relative_id])
                    _relative_id += 1
                    """

                    # Revisiting this question
                    # in case I should calculate locally
                    # for some reason

                    # Same as above.
                    C0 = upseg.k_Qs * upseg.intermittency \
                            / ((1-upseg.lambda_p) * upseg.sinuosity**(7/6.)) \
                            * self.dt

                    # From earlier
                    #C0 = upseg.C0[-1] # Should be consistent
                    # !!!!!!!!!!!!!!!!!!!! [0] z_ext

                    # Same as above
                    # This is specific for the upstream reach
                    dzdx_0_16 = ( np.abs( lp.z_ext[_relative_id][1]
                                          - lp.z_ext[_relative_id][0] )
                                  / (lp.dx_ext[_relative_id][0]) )**(1/6.)

                    # Half value -- because we are using upseg.Q?
                    # But perhaps that is as it should be, and might
                    # solve some of our problem.
                    # But actually, the above code also uses
                    # that upstream Q....
                    # ... though it really doesn't use a "C1" as stated
                    #C1 = C0 * dzdx_0_16 * upseg.Q[-1] / lp.B[0]

                    # Try something like above -- still the same?
                    _trib_coeff = -dzdx_0_16 * \
                                  ( lp.Q_ext[_relative_id][0] /
                                    lp.dx_ext[_relative_id][0] ) \
                                  / lp.land_area_around_confluence

                    # Very slightly different from upstream values !!!???WHY?
                    # Hm, in fact, they seem the same, upon running the code
                    # So perhaps no reason to do this
                    # Just use values from above function, in array?
                    self.LHSblock_matrix[row, col] = \
                                    C0 * _trib_coeff
                    ##print("L", upseg.ID, lp.upseg_trib_coeffs[_relative_id],
                    ##            C0 * _trib_coeff)

                    lp.later_trib_coeffs.append(C0 * _trib_coeff) # DEBUG

                    #print(upseg.left)
                    # THESE ARE EXACTLY IDENTICAL TO THOSE CALCULATED ABOVE.
                    # SO: NOT THE ANSWER FOR FIXING THE ITERATION PROBELM
                    # (THEY REMAIN IDENTICAL THROUGH ITERATIONS)
                    # BUT: POSSIBLY WORTH USING TO NOT SPEND TIME
                    # CALCULATING THE ABOVE.
                    self.LHSblock_matrix[row, col] = -lp.upseg_trib_coeffs[_relative_id]

                    _relative_id += 1

                    """
                    C1 = C0 * dzdx_0_16 * upseg.Q[-1] \
                          / lp.land_area_around_confluence


                    left_new = -C1 * 7/6. * 2 / lp.dx_ext[_relative_id][0]
                    self.LHSblock_matrix[row, col] = left_new
                    _relative_id += 1
                    """

                    """
                    C0 = upseg.k_Qs * upseg.intermittency \
                            / ((1-upseg.lambda_p) * upseg.sinuosity**(7/6.)) \
                            * self.dt
                    # From upseg, could be either dx_ext, so why not 0 : )
                    dzdx_0_16 = ( np.abs( lp.z_ext[_relative_id][0] -
                                          lp.z_ext[_relative_id][1] )
                                  / lp.dx_ext[_relative_id][0] )**(1/6.)
                    C1 = C0 * dzdx_0_16 * upseg.Q[-1] / \
                                    lp.land_area_around_confluence
                    # Slight hack but will work in convergent network
                    #left_new = C0 / upseg.dx_ext[0][-1]
                    left_new = C1 / lp.dx_ext[_relative_id][0]
                    self.LHSblock_matrix[row, col] = left_new
                    #if _relative_id == 0:
                    print("L", upseg.ID, left_new)
                    _relative_id += 1
                    """

                    """
                    _trib_coeff = dzdx_0_16 * 1E0 * \
                                  ( self.Q_ext[_tribi][0] /
                                    self.dx_ext[_tribi][0] ) \
                                  / self.land_area_around_confluence
                    _trib_cent += _trib_coeff
                    print("T", self.upstream_segment_IDs[_tribi],
                                self.C0 * _trib_coeff)
                    # Svae this value to transmit to upstream parts of matrix
                    self.upseg_trib_coeffs.append( self.C0 * _trib_coeff )
                    """

    def add_block_diagonal_matrix_downstream_boundary_conditions(self):
        """
        Add internal downstream boundary conditions
        """
        for lp in self.list_of_LongProfile_objects:
            for ID in lp.downstream_segment_IDs:
                # Space to edit
                col = self.block_start_absolute[self.IDs == ID][0]
                row = self.block_end_absolute[self.IDs == lp.ID][0]
                """
                # Matrix entry, assuming net aligns with ids
                downseg = self.list_of_LongProfile_objects[ID]
                # OLD
                #C0 = downseg.k_Qs * downseg.intermittency \
                #        / ((1-downseg.lambda_p) * downseg.sinuosity**(7/6.)) \
                #        * self.dt / (2 * lp.dx_ext[-1])
                # !!!C0!!!
                # I've updated dx_ext_2cell to acknowledge boundaries
                # This should therefore be used here.
                # !!!!!!!!!!!!!!!!! [0]
                C0 = downseg.k_Qs * downseg.intermittency \
                        / ((1-downseg.lambda_p) * downseg.sinuosity**(7/6.)) \
                        * self.dt
                        # We've removed this bit
                        # / lp.dx_ext_2cell[0][-1]
                # !!!!!!!!!!!!!!!!!!!!!! [0] z_ext, dx_ext
                # Hmph. The other part of this solution will be using 2*dx.
                # I think I should have a consistent equation across this
                # boundary.
                #dzdx_0_16 = ( np.abs(lp.z_ext[0][-2] - lp.z_ext[0][-1])
                #              / (lp.dx_ext[0][-1]))**(1/6.)
                # This should just be the same as the end of the standard
                # calculation
                # For a convergent network, can just choose any part
                # of dx_ext_2cell; 0 is the most sensible.
                # Then, naturally, the end.
                dzdx_0_16 = np.abs( (lp.z_ext[0][-3] - lp.z_ext[0][-1]) \
                               / lp.dx_ext_2cell[0][-1] )**(1/6.)
                # In fact, should everytihng be the same?
                # C1 values are the same here
                # But of course, are smaller than those from downstream
                # by a factor of 3.84ish
                C1 = C0 * dzdx_0_16 * lp.Q[-1] / lp.B[-1]

                dQ_term = 0
                for _iter_i in range(len(lp.x_ext)):
                    dQ_term += ( 1 / lp.Q[-1] ) \
                               * lp.dQ_ext_upwind[_iter_i][-1] \
                                / lp.dx_ext[_iter_i][-1]
                right_new = -C1 / lp.dx_ext_2cell[0][-1] \
                              * ( (7/3.)/lp.dx_ext[0][-1] # REALLY?
                                  + dQ_term)
                """
                """
                right_new_noC1 = 1 / lp.dx_ext_2cell[0][-1] \
                              * ( (7/3.)/lp.dx_ext_2cell[0][-1] # REALLY?
                                  + lp.dQ_ext_2cell[0][-1]/lp.Q[-1]
                                    / lp.dx_ext_2cell[0][-1] )
                """

                # dQ/dx term is miniscule! Doesn't really matter.
                # Turns out, it does when I have more significant
                # differences at the tributary junction.
                #right_new = -C1 / lp.dx_ext[0][-1] \
                #              * (7/3.)/lp.dx_ext_2cell[0][-1]

                # GOOD FOR 2 SEGS
                # BUT NOT ACCOUNTING PROPERLY WHEN 1 OF THOSE SEGS
                # HAS TRIBUTARIES
                # HAVE NOT LOOKED INTO EXACTLY WHY (ABOVE)
                # JUST INTUITED AND CHANGED TO LP.RIGHT[0]

                # AH, PROBABLY DQ_TERM
                # WELL, NOT JUST THAT. BUT COULD SORT THIS,
                # OR JUST ACCEPT THAT I DID IT CORRECTLY ABOVE.
                #####self.LHSblock_matrix[row, col] = right_new

                self.LHSblock_matrix[row, col] = lp.right[0]

                """
                # Notes from above
                self.C0 = self.k_Qs * self.intermittency \
                            / ((1-self.lambda_p) * self.sinuosity**(7/6.)) \
                            * self.dt

                self.left = -self.C1 / self.dx_ext_2cell \
                                * ( (7/3.)/self.dx_ext[:-1]
                                - self.dQ_ext_2cell/self.Q/self.dx_ext_2cell )
                self.center = -self.C1 / self.dx_ext_2cell \
                                      * ( (7/3.)
                                      * (-1/self.dx_ext[:-1]
                                         -1/self.dx_ext[1:]) ) \
                                         + 1.
                self.right = -self.C1 / self.dx_ext_2cell \
                                      * ( (7/3.)/self.dx_ext[1:] # REALLY?
                                          + self.dQ_ext_2cell/self.Q/self.dx_ext_2cell )
                """

                # Perhaps I need the other half of the "handshake" across
                # matrices to also be the d/dx (dQs/dx) form, with 1-cell
                # rather than 2-cell calculations.

                # The upstream term can balance it better, but isn't actually
                # adjusted for local slope. Hence right_new[0] sort of works
                """
                dzdx_0_16 = np.abs( (lp.z[-1] - downseg.z[0]) \
                               / lp.dx_ext[0][-1] )**(1/6.)
                C1 = C0 * dzdx_0_16 * (lp.Q[-1]+downseg.Q[0])/2. / lp.B[-1]
                right_new = -C1 / lp.dx_ext[0][-1] \
                              * ( (7/3.)/lp.dx_ext[0][-1] # REALLY?
                                  + (downseg.Q[0] - lp.Q[-1])/lp.Q[-1]
                                    / lp.dx_ext[0][-1] )
                self.LHSblock_matrix[row, col] = right_new
                """

                #print("OLD,NEW:", lp.right[0], right_new*2)
                #self.LHSblock_matrix[row, col] = lp.right[0]

                #self.LHSblock_matrix[row, col] = right_new*2



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

    def set_niter(self, niter):
        # MAKE UNIFORM IN BASE CLASS
        self.niter = niter

    def create_x_ext_lists(self):
        """
        ##########################################################
        # GENERATE LISTS OF x_ext: 1 FOR EACH INCOMING TRIBUTARY #
        ##########################################################

        Set up "network" lists of "_ext" variables: one per upstream-linked
        segment and a minimum of 1 if no links are present
        Currently building this for convergent networks only
        """
        # Pad x_ext with nans
        _nan1 = np.array([np.nan])
        # Loop through long profiles (segments) in network
        for lp in self.list_of_LongProfile_objects:
            lp.x_ext = []
            x_inner = np.concatenate( [_nan1, lp.x, _nan1] )
            for _iter_i in range(np.max( (1, len(lp.upstream_segment_IDs)) )):
                lp.x_ext.append(x_inner.copy())


    def update_x_ext_internal(self):
        """
        ###################################################
        # POPULATE x_ext LISTS WITH VALUES FROM NEIGHBORS #
        ###################################################

        x_ext[0] of downstream segment set to x[-1] of upstream segment.

        The upstream-most segments have x_ext set by another function, and
        to default to a spacing that is the same as that between x[0] and x[1].

        x_ext[-1] of upstream segment set to x[0] of downstream segment.

        The final downstream segment has x_ext[-1] for the lateral position
        of the base-level cell.
        """

        # Order:
        # Inner (close to each other): upstream.
        # Outer (strides from each other): downstream.
        # Right now, this is moot: Downstream = 1
        for lp in self.list_of_LongProfile_objects:
            _idx = 0
            # SET UPSTREAM BOUNDARIES: INTERNAL
            # Here, max so the downstream-most segment gets looped through
            # too, even if it has no downseg ID
            for i_downseg in range(np.max((1, len(lp.downstream_segment_IDs)))):
                for upseg_ID in lp.upstream_segment_IDs:
                    upseg = self.list_of_LongProfile_objects[upseg_ID]
                    lp.x_ext[_idx][0] = upseg.x[-1]
                    _idx += 1
            # SET DOWNSTREAM BOUNDARIES: INTERNAL
            _idx = 0
            for downseg_ID in lp.downstream_segment_IDs:
                # For each downstream ID, must update for each upstream
                # ID
                downseg = self.list_of_LongProfile_objects[downseg_ID]
                # Min = 1 so downseg still updated for headwaters segments
                for i_upseg in range(np.max((1, len(lp.upstream_segment_IDs)))):
                    lp.x_ext[_idx][-1] = downseg.x[0]
                    _idx += 1

    def update_x_ext_external_upstream(self):
        """
        Update x_ext at external upstream boundaries.

        This, together with z, provides sediment inputs as these locations.

        By default, this is set to be the same as the dx between the
        first and the second cells: The upstream boundary condition is
        set by the slope, and so z_ext[0] can be used with a fixed x_ext[0]
        to set this -- and then, there are fewer moving parts / places
        with possible error.

        In fact, I am writing this function to make sure that z_ext[0] is
        *always* just at the same spacing as that between x[0] and x[1].
        Let's keep this simple!

        But nonetheless keeping this in its own function to highlight the
        conceptual difference and in case we want to change this functionality
        later.
        """
        # SET UPSTREAM BOUNDARIES: EXTERNAL
        for ID in self.list_of_channel_head_segment_IDs:
            lp = self.list_of_LongProfile_objects[ID]
            for x_ext_array in lp.x_ext:
                dx0 = lp.x[1] - lp.x[0]
                x_ext_array[0] = lp.x[0] - dx0

    def set_x_bl(self, x_bl):
        """
        Alias for `update_x_ext_external_downstream`.
        !!!!!
        MAYBE I SHOULD CALL z0 --> z_bl
        """
        update_x_ext_external_downstream( x_bl )

    def update_x_ext_external_downstream(self, x_base_level=None):
        """
        Set downstream boundary (ultimate base level, singular): External

        This function will set only the downstream-most boundary condition.

        It expects a list of length (1) for the class variable:
        self.list_of_channel_mouth_segment_IDs.
        This assumption will have to be relaxed if the code ever be updated
        to allow multiple river mouths.

        Args:
            x0 (float): Base-level downvalley position (mouth seg x_ext[-1])

        Returns:
            None
        """
        if len(self.list_of_channel_mouth_segment_IDs) == 1:
            ID = self.list_of_channel_mouth_segment_IDs[0]
        else:
            sys.exit( ">1 channel-mouth-segment ID listed.\n"+
                      "Simulation not set up to manage >1 river mouth.\n"+
                      "Exiting" )

        # SET DOWNSTREAM BOUNDARY (ULTIMATE BASE LEVEL, SINGULAR): EXTERNAL

        # Flag for whether x should be calculated internally
        _calcx = False

        # If provided, use x position
        if x_base_level is not None:
            _x_bl = x_base_level
        # Otherwise, default to one dx beyond river-mouth position
        else:
            _calcx = True

        lp = self.list_of_LongProfile_objects[ID]

        # Loop over each 1D array within the list: each trib connection
        for _x_ext in lp.x_ext:
            # This should be the same each time, but nonetheless,
            # calculating it locally for each array makes me more comfortable
            # because I am not assuming that one equals they other.
            # Though if they are unequal... potential big problems!
            if _calcx:
                _x_bl = lp.x[-1] + lp.dx[-1]
            _x_ext[-1] = _x_bl

        # We should have some code to account for changes in both x and z
        # with base-level change, and remeshes the downstream-most segment,
        # as needed

    def create_z_ext_lists(self):
        """
        ##########################################################
        # GENERATE LISTS OF z_ext: 1 FOR EACH INCOMING TRIBUTARY #
        ##########################################################

        Set up "network" lists of "_ext" variables: one per upstream-linked
        segment and a minimum of 1 if no links are present
        Currently building this for convergent networks only
        """
        """
        # Pad x_ext with nans
        _nan1 = np.array([np.nan])
        # Loop through long profiles (segments) in network
        for lp in self.list_of_LongProfile_objects:
            lp.z_ext = np.max( (1, len(lp.upstream_segment_IDs)) ) * \
                [ np.concatenate( [_nan1, lp.z, _nan1] ) ]
            print( "" )
            print( lp.ID )
            print( "" )
            print( "Z_EXT" )
            print( lp.z_ext[:] )
            print( "" )
        """
        # z_ext messup happens after this.

        # HUH! I WONDER WHAT I WAS THINKING WHEN I WROTE THAT.
        # BUT I SEEM TO HAVE AN ANSWER.

        _nan1 = np.array([np.nan])
        for lp in self.list_of_LongProfile_objects:
            lp.z_ext = []
            z_inner = np.concatenate( [_nan1, lp.z, _nan1] )
            for _iter_i in range(np.max( (1, len(lp.upstream_segment_IDs)) )):
                lp.z_ext.append(z_inner.copy())

    def update_z_ext_internal(self):
        """
        ###################################################
        # POPULATE x_ext LISTS WITH VALUES FROM NEIGHBORS #
        ###################################################

        z_ext[0] of downstream segment set to z[-1] of upstream segment.

        The upstream-most segments have z_ext set based on sediment-supply
        boundary conditions, rather than being set here.

        z_ext[-1] of upstream segment set to z[0] of downstream segment.
        This simulates the "internal base level" communicated
        among tributaries in the network.

        The final downstream segment has z_ext[-1] set as base level.
        """

        # Order:
        # Inner (close to each other): upstream.
        # Outer (strides from each other): downstream.
        # Right now, this is moot: Downstream = 1
        for lp in self.list_of_LongProfile_objects:
            # SET UPSTREAM BOUNDARIES: INTERNAL
            _idx = 0
            # Here, max so the downstream-most segment gets looped through
            # too, even if it has no downseg ID
            for i_downseg in range(np.max((1, len(lp.downstream_segment_IDs)))):
                for upseg_ID in lp.upstream_segment_IDs:
                    upseg = self.list_of_LongProfile_objects[upseg_ID]
                    lp.z_ext[_idx][0] = upseg.z_ext[0][-2]
                    _idx += 1
            # SET DOWNSTREAM BOUNDARIES: INTERNAL
            _idx = 0
            for downseg_ID in lp.downstream_segment_IDs:
                # For each downstream ID, must update for each upstream
                # ID
                downseg = self.list_of_LongProfile_objects[downseg_ID]
                # Min = 1 so downseg still updated for headwaters segments
                for i_upseg in range(np.max((1, len(lp.upstream_segment_IDs)))):
                    lp.z_ext[_idx][-1] = downseg.z_ext[0][1]
                    _idx += 1

    def update_z_ext_external_upstream(self, S0=None, Q_s_0=None):
        """
        Update z_ext at external upstream boundaries.

        This provides sediment inputs as these locations.

        If the value is iterable, it will provide upstream boundary conditions
        in the same order as that of the provided headwater segments.

        If it is a scalar, it will provide the same value for each segment.

        If self.S0 and/or self.Q_s_0 have already been set, you can run this
        function without passing any variables.

        Anything that you do pass here will overwrite previously-set values
        for these variables.

        Note: This overlaps somewhat with lp.set_Qs_input_upstream(Q_s_0).
        However, it is more flexible (S0 or Q_s_0) and expects z_ext to
        be a single array within a list (as opposed to an array outside of a
        list).
        """


        # THIS IS RATHER MESSY AT HANDLING INTERNAL VS EXTERNAL S0, Q_S_0
        # WORKS BETTER FOR NOW, BUT SHOULD REWRITE
        if self.S0 is not None and self.Q_s_0 is not None:
            # Use Q_s_0
            pass

        elif S0 is not None and Q_s_0 is not None:
            sys.exit( "Choose only one of S0, Q_s_0.\n"+
                      "(Q_s_0 is used to generate S0.)" )

        if S0 is None and Q_s_0 is None:
            if self.S0 is not None and self.Q_s_0 is not None:
                warnings.warn( "\nUnclear whether to update Q_s_0 and S0 "+
                               "based on input S0 or input Q_s_0.\n"+
                               "Leaving function without updating values." )
            return

        # Before starting, set a bool as a scalar check.
        _is_scalar = False

        # And a flag for using Q_s_0
        _use_Q_s_0 = False

        #########################################
        # IF Q_s_0 IS USED, FIRST CONVERT TO S0 #
        #########################################

        # First, check on whether it exists already. Set or just use this
        if Q_s_0 is not None:
            # Set the flag
            _use_Q_s_0 = True
            # Set the internal variable
            self.Q_s_0 = Q_s_0
            # Scalar or array?
            try:
                iter(Q_s_0)
                Q_s_0 = np.array(Q_s_0).squeeze() # Use a numpy array
            except:
                _is_scalar=True

        elif self.Q_s_0 is not None:
            # Set the flag
            _use_Q_s_0 = True
            Q_s_0 = self.Q_s_0 # set the local variable (ease of use here)
            try:
                iter(Q_s_0)
                Q_s_0 = np.array(Q_s_0).squeeze() # Enforce numpy array
                                                  # (should be one already)
            except:
                _is_scalar=True

        # Second, if array, check length
        if not _is_scalar and _use_Q_s_0:
            if len(Q_s_0) != len(self.list_of_channel_head_segment_IDs):
                sys.exit( "Q_s_0 array length is "+str(len(Q_s_0))+'.\n'
                          "Number of headwater segments is "+
                            str(len(self.list_of_channel_head_segment_IDs))+
                            '.\n'+
                          "These should be equal or Q_s_0 should be "+
                          "passed as a scalar\n"+
                          "(same value everywhere)."
                 )

        # Third, set self.Q_s_0 for all the segments
        # This is done whether Q_s_0 be scalar or array type
        _idx = 0
        if _use_Q_s_0:
            for ID in self.list_of_channel_head_segment_IDs:
                lp = self.list_of_LongProfile_objects[ID]
                # If Q_s_0 is a scalar value, apply it everywhere
                if _is_scalar:
                    lp.Q_s_0 = Q_s_0
                # Otherwise, iterate over the supplied Q_s_0
                else:
                    lp.Q_s_0 = Q_s_0[_idx]
                _idx += 1

        # Fourth, compute the S0 values
        if _use_Q_s_0:
            _Q0 = []
            _sinuosity = []
            _intermittency = []
            _k_Qs = []
            for ID in self.list_of_channel_head_segment_IDs:
                lp = self.list_of_LongProfile_objects[ID]
                # We are using Q[0] instead of (Q_ext[0] + Q_ext[1])/2
                # Q_ext[1] = Q[0], by definitions
                # Therefore, I will set all Q_ext[0] = Q[0].
                # This will be superfluous for the code, but at least
                # consistent internally
                _Q0.append(lp.Q[0])
                _sinuosity.append(lp.sinuosity)
                _intermittency.append(lp.intermittency)
                _k_Qs.append(lp.k_Qs)
            _Q0 = np.array(_Q0)
            _sinuosity = np.array(_sinuosity)
            _intermittency = np.array(_intermittency)
            _k_Qs = np.array(_k_Qs)
            # Note: Negative S0 if sloping downstream.
            # This is reverse to the usuaal (backwards) sign convention.
            S0 = - np.sign(_Q0) * _sinuosity * \
              ( np.abs(Q_s_0) /
                ( _k_Qs * _intermittency
                      * np.abs(_Q0)) )**(6/7.)

        ################################################
        #        IF S0 BE PROVIDED, JUST USE IT        #
        # OTHERWISE, THIS USES THE ABOVE-CALCULATED S0 #
        ################################################

        if S0 is not None:
            self.S0 = S0
        else:
            S0 = self.S0

        ##print("SCALAR?", _is_scalar)
        ##print( S0 )
        ##print( np.atleast_1d(np.array(S0).squeeze) )

        # S0 might be iterable even if Q_s_0 be not
        try:
            iter(S0)
            # Enforce numpy array
            # And ensure that it isn't squeezed down to 0D
            # if there is just 1 S0 value given
            S0 = np.atleast_1d(np.array(S0).squeeze())
        except:
            _is_scalar=True

        ##print("SCALAR?", _is_scalar)

        # FIFTH: Set S0 and z_ext[0]
        _idx = 0
        for ID in self.list_of_channel_head_segment_IDs:
            lp = self.list_of_LongProfile_objects[ID]
            lp.S0 = S0[_idx]
            # if _is_scalar:
            #     lp.S0 = S0
            # else:
            #     lp.S0 = S0[_idx]
            _idx += 1
            # Hard-coding: Expecting only one segment in list
            # Because this is just for the channel-head segments
            ##print("ID", lp.ID)
            lp.z_ext[0][0] = lp.z[0] + lp.S0 * lp.dx[0]

    def set_z_bl (self, z0):
        """
        Alias for `update_z_ext_external_downstream`.
        !!!!!
        MAYBE I SHOULD CALL z0 --> z_bl
        """
        update_z_ext_external_downstream( z0 )

    def update_z_ext_external_downstream(self, z0):
        """
        Set downstream boundary (ultimate base level, singular): External

        This function will set only the downstream-most boundary condition.

        It expects a list of length (1) for the class variable:
        self.list_of_channel_mouth_segment_IDs.
        This assumption will have to be relaxed if the code ever be updated
        to allow multiple river mouths.

        Args:
            z0 (float): Base-level elevation. Sets z_ext[-1] for the mouth seg

        Returns:
            None
        """
        if len(self.list_of_channel_mouth_segment_IDs) == 1:
            ID = self.list_of_channel_mouth_segment_IDs[0]
        else:
            sys.exit( ">1 channel-mouth-segment ID listed.\n"+
                      "Simulation not set up to manage >1 river mouth.\n"+
                      "Exiting" )
        # SET DOWNSTREAM BOUNDARY (ULTIMATE BASE LEVEL, SINGULAR): EXTERNAL
        lp = self.list_of_LongProfile_objects[ID]
        # Loop over each 1D array within the list: each trib connection
        for _z_ext in lp.z_ext:
            _z_ext[-1] = z0
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
        self.list_of_channel_head_segment_IDs = []
        for lp in self.list_of_LongProfile_objects:
            if not lp.upstream_segment_IDs:
                self.list_of_channel_head_segment_IDs.append(lp.ID)

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
        self.list_of_channel_mouth_segment_IDs = []
        for lp in self.list_of_LongProfile_objects:
            if not lp.downstream_segment_IDs:
                self.list_of_channel_mouth_segment_IDs.append(lp.ID)

        if len(self.list_of_channel_mouth_segment_IDs) == 1:
            self.channel_mouth_segment_ID = \
                self.list_of_channel_mouth_segment_IDs[0]
        elif len(self.list_of_channel_mouth_segment_IDs) == 0:
            sys.exit("Ahmm... why are there no river mouths?")
        else:
            sys.exit("Ahmm... why are there multiple river mouths?")

    def update_dx_ext(self):
        """
        Create dx_ext arrays -- one for each upstream link -- from x_ext.
        """
        for lp in self.list_of_LongProfile_objects:
            lp.dx_ext = []
            for x_ext in lp.x_ext:
                lp.dx_ext.append( np.diff(x_ext) )

    def update_dx_2cell(self):
        """
        Create dx_2cell arrays: One, internal to each segment.

        This function assumes that distance increases from left to right.

        Look here in case sign errors are encountered, but it also feels
        safer to keep this as such to make sure that we don't miss such errors.
        It is also possible that GRLP will run anyway because of the abs()
        involved with the slope calculation.
        """
        for lp in self.list_of_LongProfile_objects:
            lp.dx_2cell = lp.x[2:] - lp.x[:-2]

    def update_dx_ext_2cell(self):
        """
        Create dx_ext arrays -- one for each upstream link -- from x_ext.

        This function assumes that distance increases from left to right.

        Look here in case sign errors are encountered, but it also feels
        safer to keep this as such to make sure that we don't miss such errors.
        It is also possible that GRLP will run anyway because of the abs()
        involved with the slope calculation.
        """
        for lp in self.list_of_LongProfile_objects:
            lp.dx_ext_2cell = []
            for x_ext in lp.x_ext:
                lp.dx_ext_2cell.append( x_ext[2:] - x_ext[:-2] )

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

    def create_Q_ext_lists(self):
        """
        ##########################################################
        # GENERATE LISTS OF Q_ext: 1 FOR EACH INCOMING TRIBUTARY #
        ##########################################################

        Set up "network" lists of "_ext" variables: one per upstream-linked
        segment and a minimum of 1 if no links are present
        Currently building this for convergent networks only

        Run after "update_Q"
        """

        """
        # Pad Q_ext with nans
        _nan1 = np.array([np.nan])
        # Loop through long profiles (segments) in network
        for lp in self.list_of_LongProfile_objects:
            lp.Q_ext = np.max( (1, len(lp.upstream_segment_IDs)) ) * \
                [ np.concatenate( [_nan1, lp.Q, _nan1] ).copy() ]
        """
        # ^ Even wtih .copy(), this caused the two arrays to be linked in
        # memory. Therefore, the discharge of the second would always
        # be chosen. This casued extreme weirdness with slopes
        # at confluences -- it would seem as if there were more or less
        # discharge coming from the tributaries when compared
        # to the mainstem

        _nan1 = np.array([np.nan])
        for lp in self.list_of_LongProfile_objects:
            lp.Q_ext = []
            Q_inner = np.concatenate( [_nan1, lp.Q, _nan1] )
            for _iter_i in range(np.max( (1, len(lp.upstream_segment_IDs)) )):
                lp.Q_ext.append(Q_inner.copy())



    def update_Q_ext_from_Q(self):
        """
        Run in order after "update_Q()" and "create_Q_ext_lists()".

        This sets the [1:-1] (i.e., non-boundary) values for each Q_ext array
        within each Q_ext list
        """
        for lp in self.list_of_LongProfile_objects:
            # List of arrays
            for Q_ext_array in lp.Q_ext:
                Q_ext_array[1:-1] = lp.Q

    def update_Q_ext_internal(self):
        """
        ###################################################
        # POPULATE Q_ext LISTS WITH VALUES FROM NEIGHBORS #
        ###################################################

        Q_ext[0] of downstream segment set to Q[-1] of upstream segment.
        This is done for each Q_ext array within the list of arrays,
        corresponding to each tributary junction.

        The upstream-most segments have Q_ext set by another function;
        this becomes part of the broader upstream boundary condition
        (including how Q_s_0 is managed).
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REVISIT
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        Q_ext[-1] of upstream segment set to Q[0] of downstream segment.

        The final downstream segment has Q_ext[-1] = Q[-1], set by another
        function:
        Seems reasonable that at the very mouth, we have no new water inputs,
        but the flow does need to stay continuous and exit the domain.
        !!!!!!!!!!!!!!!!!!!!
        DO THIS
        !!!!!!!!!!!!!!!!!!!!
        lp.Q_ext[
        Run this after update_Q, to make sure that it is using
        the most recent discharge values
        """
        # Order:
        # Inner (close to each other): upstream.
        # Outer (strides from each other): downstream.
        # Right now, this is moot: Downstream = 1
        for lp in self.list_of_LongProfile_objects:
            _idx = 0
            # SET UPSTREAM BOUNDARIES: INTERNAL
            # Here, max so the downstream-most segment gets looped through
            # too, even if it has no downseg ID
            for i_downseg in range(np.max((1, len(lp.downstream_segment_IDs)))):
                for upseg_ID in lp.upstream_segment_IDs:
                    upseg = self.list_of_LongProfile_objects[upseg_ID]
                    lp.x_ext[_idx][0] = upseg.x[-1]
                    _idx += 1
            # SET DOWNSTREAM BOUNDARIES: INTERNAL
            _idx = 0
            for downseg_ID in lp.downstream_segment_IDs:
                # For each downstream ID, must update for each upstream
                # ID
                downseg = self.list_of_LongProfile_objects[downseg_ID]
                # Min = 1 so downseg still updated for headwaters segments
                for i_upseg in range(np.max((1, len(lp.upstream_segment_IDs)))):
                    lp.x_ext[_idx][-1] = downseg.x[0]
                    _idx += 1

        # Order:
        # Inner (close to each other): upstream.
        # Outer (strides from each other): downstream.
        # Right now, this is moot: Downstream = 1

        for lp in self.list_of_LongProfile_objects:
            _idx = 0
            # SET UPSTREAM BOUNDARIES: INTERNAL
            # Here, max so the downstream-most segment gets looped through
            # too, even if it has no downseg ID
            for i_downseg in range(np.max((1, len(lp.downstream_segment_IDs)))):
                for upseg_ID in lp.upstream_segment_IDs:
                    upseg = self.list_of_LongProfile_objects[upseg_ID]
                    ##print( upseg )
                    lp.Q_ext[_idx][0] = upseg.Q[-1]
                    ##print("!!!!!!!!!!!!!!!!!!!!!!!!")
                    ##print(_idx)
                    ##print("Q_ext", lp.Q_ext)
                    ##print("!!!!!!!!!!!!!!!!!!!!!!!!")
                    _idx += 1
            # SET DOWNSTREAM BOUNDARIES: INTERNAL
            _idx = 0
            for downseg_ID in lp.downstream_segment_IDs:
                # For each downstream ID, must update for each upstream
                # ID
                downseg = self.list_of_LongProfile_objects[downseg_ID]
                # Min = 1 so downseg still updated for headwaters segments
                for i_upseg in range(np.max((1, len(lp.upstream_segment_IDs)))):
                    lp.Q_ext[_idx][-1] = downseg.Q[0]
                    _idx += 1

    def update_Q_ext_external_upstream(self):
        """
        Update Q_ext at external upstream boundaries.

        Based on how S0 is defined and the need for a Qs:Qw ratio to set S0,
        the upstream slope (and sediment-supply) boundary condition,
        Q at all upstream boundaries is simply set to be identical to the
        same value as it is internally.

        Functionally, this doesn't matter: slopes are based on Q[0]
        rather than (Q_ext[0] + Q_ext[1])/2, where Q_ext[1] = Q[0].
        """
        # SET UPSTREAM BOUNDARIES: EXTERNAL
        for ID in self.list_of_channel_head_segment_IDs:
            lp = self.list_of_LongProfile_objects[ID]
            for Q_ext_array in lp.Q_ext:
                Q_ext_array[0] = lp.Q[0]

    def update_Q_ext_external_downstream(self):
        """
        Set discharge at downstream boundary (ultimate base level, singular)

        It expects a list of length (1) for the class variable:
        self.list_of_channel_mouth_segment_IDs.
        This assumption will have to be relaxed if the code ever be updated
        to allow multiple river mouths.

        Here, we just assume that the downstream-boundary water discharge (Q)
        is identical to the one just above -- no new tributaries join as it
        enters the ocean, lake, river, basin, etc.
        """
        if len(self.list_of_channel_mouth_segment_IDs) == 1:
            ID = self.list_of_channel_mouth_segment_IDs[0]
        else:
            sys.exit( ">1 channel-mouth-segment ID listed.\n"+
                      "Simulation not set up to manage >1 river mouth.\n"+
                      "Exiting" )
        lp = self.list_of_LongProfile_objects[ID]
        # Assume that there is just one river mouth
        # Could easily make this become a loop
        # Nope: Multiple arrays if there are also upstream segments
        # Just loop over them all
        for Q_ext_array in lp.Q_ext:
                Q_ext_array[-1] = lp.Q[-1]

    def update_dQ_ext_2cell(self):
        """
        Use segment adjacencies to set changes in discharge down segments.
        For a convergent network:

        Q_ext[-1] of the upstream segment should see a large-ish increase
        because the discharge after the tributary junction will can be
        significantly higher than that above.

        Q_ext[0] of the downstream segment should see only a modest difference
        overall (i.e., when summing the upstream tributaries), but we in fact
        require these to be separated into a dQ_ext_upwind for each river-segment
        combination (here, typically two tributaries joining into one
        downstream river segment). This is needed to properly weight
        slopes for the C1 coefficient.

        # NOW BACK TO 2CELL, FOR CONGRUENCE WITH SLOPE CALCULATIONS.
        # NEED TO SMEAR Q FOR IT TO WORK, THOUGH
        """
        for lp in self.list_of_LongProfile_objects:
            lp.dQ_ext_2cell = []
            for Q_ext_array in lp.Q_ext:
                lp.dQ_ext_2cell.append( (Q_ext_array[2:] - Q_ext_array[:-2]) )

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
                lp.set_intermittency = intermittency
        else:
            i = 0
            for lp in self.list_of_LongProfile_objects:
                lp.set_intermittency = intermittency[i]
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

        # Identify channel head and mouth segments
        self.create_list_of_channel_head_segment_IDs()
        self.create_list_of_channel_mouth_segment_IDs()

        # Generate arrays of x, including networked links
        self.create_x_ext_lists()
        self.update_x_ext_internal()
        self.update_x_ext_external_upstream()                           # b.c.
        self.update_x_ext_external_downstream( x_bl )                   # b.c.

        # From these, generate arrays of dx
        self.update_dx_ext()
        self.update_dx_2cell()
        self.update_dx_ext_2cell()

        # Generate arrays of Q based on externally provided (user-set) values
        self.update_Q( Q )
        self.create_Q_ext_lists()
        self.update_Q_ext_from_Q()
        self.update_Q_ext_internal()
        self.update_Q_ext_external_upstream()  # b.c., Q_ext[0] = Q[0]
        self.update_Q_ext_external_downstream()   # b.c., Q_ext[-1] = Q[-1]
        self.update_dQ_ext_2cell()

        # Generate arrays of z based on externally provided (user-set) values
        self.create_z_ext_lists()
        self.update_z_ext_internal()
        self.update_z_ext_external_upstream( S0 = S0, Q_s_0 = Q_s_0 )  # b.c.
        self.update_z_ext_external_downstream( z_bl )                   # b.c.

        # # Generate arrays of Q based on externally provided (user-set) values
        # self.update_Q( Q )
        # self.create_Q_ext_lists()
        # self.update_Q_ext_from_Q()
        # self.update_Q_ext_internal()
        # self.update_Q_ext_external_upstream()  # b.c., Q_ext[0] = Q[0]
        # self.update_Q_ext_external_downstream()   # b.c., Q_ext[-1] = Q[-1]
        # self.update_dQ_ext_2cell()

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
        self.update_z_ext_internal() # FM: set again later inside iteration?
        # self.dt is decided earlier

        for ti in range(int(self.nt)):
            for lp in self.list_of_LongProfile_objects:
                lp.build_LHS_coeff_C0(dt=self.dt)
                lp.zold = lp.z.copy() # DEBUG OR NOT? KEEP OR NOT? < -- KEEP! FOR ITERATOR.
                lp.network__compute_coefficient_time_varying() # <-- Quelle der Problem
                #lp.zold = lp.z.copy() # DEBUG OR NOT? KEEP OR NOT?
            """
            # Zap from here: Only inside iteration
            for lp in self.list_of_LongProfile_objects:
                lp.network__build_matrix_inner()
            self.map_block_diagonal_matrix_blocks()
            self.create_block_diagonal_matrix_with_internal_tridiagonals()
            self.add_block_diagonal_matrix_upstream_boundary_conditions()
            self.add_block_diagonal_matrix_downstream_boundary_conditions()
            """
            # b.c. for no links
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            """
            for lp in self.list_of_LongProfile_objects:
                if len(lp.upstream_segment_IDs) == 0:
                    lp.set_bcl_Neumann_LHS()
                    lp.set_bcl_Neumann_RHS()
                if len(lp.downstream_segment_IDs) == 0:
                    lp.set_bcr_Dirichlet()
            """
            # Leave to iterations?
            #self.stack_RHS_vector()

            for _iter_i in range(self.niter):
                self.update_z_ext_internal()
                for lp in self.list_of_LongProfile_objects:
                    # Update coefficient for all: elements may call to others
                    # within the net
                    lp.network__compute_coefficient_time_varying()
                for lp in self.list_of_LongProfile_objects:
                    lp.network__build_matrix_inner()
                # Try adding all of these here??
                self.map_block_diagonal_matrix_blocks()
                self.create_block_diagonal_matrix_with_internal_tridiagonals()
                self.add_block_diagonal_matrix_upstream_boundary_conditions()
                self.add_block_diagonal_matrix_downstream_boundary_conditions()
                self.stack_RHS_vector()
                # Returning z to z_old if in the iteration loops
                # Shouldnt be necessary on first iteration, but will
                # add this check after making sure that other things
                # are working
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # Removed reference to this: dzdx_0_16 now only from z_ext
                # and friends. Can update as much as we want before the end :)
                #for lp in self.list_of_LongProfile_objects:
                #    lp.z = lp.zold.copy()
                # Update semi-implicit on boundaries
                # Commenting these two out helps solution!
                # Don't understand why. Perhaps error in code for them?
                #self.add_block_diagonal_matrix_upstream_boundary_conditions()
                #self.add_block_diagonal_matrix_downstream_boundary_conditions()
                out = spsolve(sparse.csr_matrix(self.LHSblock_matrix), self.RHS)

                lp_i = 0
                idx = 0
                for lp in self.list_of_LongProfile_objects:
                    # IMPORTANT: NEED TO UPDATE Z FIRST
                    # Z_EXT VALUES BASED ON Z VALUES, WHEN
                    # INTERNAL VALUES ARE UPDATED FROM UPSTREAM
                    # Hm: Might not want to update z within iterator.

                    # !!!!!!!!!!!!!!!!!!!!!!!
                    # ITERATOR NOT WORKING UNTIL THIS IS FIXED
                    #

                    # NEED TO UPDATE lp.z: OTHER FUNCTIONS ARE BASED ON THIS
                    lp.z = out[idx:idx+self.list_of_segment_lengths[lp_i]]

                    # Not sure if this be the way to go:
                    #for _tribi in range(len(lp.z_ext)):
                    #    lp.z_ext[_tribi][1:-1] = lp.z

                    # Or maybe this is it:
                    self.update_z_ext_internal()
                    self.update_z_ext_external_upstream( S0 = self.S0,
                                                          Q_s_0 = self.Q_s_0 )
                    # This probably needn't be run, but whatever
                    #self.update_z_ext_external_downstream( z_bl )

                    idx += +self.list_of_segment_lengths[lp_i]
                    lp_i += 1

                    # Update boundaries in later functions

            # Simplify; Hard-code single iteration
            i = 0
            idx = 0
            for lp in self.list_of_LongProfile_objects:
                lp.z = out[idx:idx+self.list_of_segment_lengths[i]]
                for _tribi in range(len(lp.z_ext)):
                    lp.z_ext[_tribi][1:-1] = lp.z
                idx += +self.list_of_segment_lengths[i]
                i += 1
            self.update_z_ext_internal()
            # Update upstream boundary condition: Elevation may change to keep
            # slope constant
            self.update_z_ext_external_upstream( S0 = self.S0, Q_s_0 = self.Q_s_0 )

            self.t += self.dt # Update each lp z? Should make a global class
                              # that these both inherit from
            for lp in self.list_of_LongProfile_objects:
                lp.t = self.t

            for lp in self.list_of_LongProfile_objects:
                lp.dz_dt = (lp.z - lp.zold)/self.dt
                #lp.Qs_internal = 1/(1-lp.lambda_p) * np.cumsum(lp.dz_dt)*lp.B \
                #                 + lp.Q_s_0

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
            dxs = np.append( dxs, seg.dx_ext[0][1:] )
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
        
        # get dxs
        dxs = np.hstack(
            [seg.dx_ext[0][1:] for seg in self.list_of_LongProfile_objects]
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
        for seg in self.list_of_LongProfile_objects: seg.compute_diffusivity()
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





"""
# Test whether the lists have populated
for lp in self.list_of_LongProfile_objects: print(lp.x_ext)
print("")
for lp in self.list_of_LongProfile_objects: print(lp.z_ext)
"""
