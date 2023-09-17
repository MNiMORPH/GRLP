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
        elif x_ext is not None:
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
        self.S0 = np.sign(self.Q[0]) * self.sinuosity * \
                      ( np.abs(Q_s_0) / 
                        ( self.k_Qs * self.intermittency 
                              * np.abs(self.Q[0])) )**(6/7.)
        # Give upstream cell the same width as the first cell in domain
        self.z_ext[0] = self.z[0] + self.S0 * self.dx_ext[0]
        sys.exit(4)

    def update_z_ext_0(self):
        # Give upstream cell the same width as the first cell in domain
        # UPDATE TO WORK WITH LIST OF ARRAYS
        # Before fixing, think about what z_ext and C1 are really doing
        # and how best to set them
        # Does z_ext set up the boundary conditions in networked mode?
        
        # Q1: z_ext and boundary conditions
        
        # Q2: z_ext and C1
        #     self.C1 = self.C0 * self.dzdx_0_16 * self.Q / self.B
        # And C0 has all local variables.
        # So what we really need is the channel slope at the local site
        # How do we weight it? Equally? By the discharge from each river?
        # Yes! Q_s = k_Qs * Q * S^(7/6), and we are differentiating
        # d Q_s / d x
        # Therefore, S for the coefficient will be set via:
        # Qs = k_Qs * Q * S^{7/6}
        # START HERE!!!!!!!!!!!!
        if type( self.z_ext ) is np.ndarray:
            # Only one segment: towards applying boundary condition upstream
            self.z_ext[0] = self.z[0] + self.S0 * self.dx_ext[0]
        elif type( self.z_ext ) is list:
            for connected_array in self.z_ext:
            # Wait a minute! We don't necessarily want the upstream boundary
            # conditions to be enforced everywhere; we have internal segments.
                if len(self.upstream_segment_IDs) == 0:
                    # WHY ARE WE DOING THIS?
                    # ALREADY ADDRESSED WITHIN NETWORK?
                    # I GUESS IT DOESN'T HURT, BUT REVISIT AND STREAMLINE
                    # WHEN THERE ARE FEWER MOVING PARTS
                    # If upstream-most, then apply boundary condition
                    # This should run just once in the loop and then exit
                    # DX_EXT SHOULD BE LIST OF ARRAYS, SO NEED [0][0]
                    # TO GET FIRST (ONLY) ARRAY, AND THEN ITS FIRST ELEMENT
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    connected_array[0] = self.z[0] + self.S0 * self.dx_ext[0][0]
                else:
                    # Otherwise, obtain upstream neighbors' elevations
                    # PERHAPS NEED TO GET INFO FROM NETWORK!
                    # We can't, buuuut: This is already wrapped in:
                    # "if lp.S0 is not None"
                    # So: Let's throw an error here.
                    print("S0 was specified for a non-headwaters segment. Exiting.")
                    sys.exit(1)
    
    def compute_coefficient_time_varying(self):
        i=0
        print( ((self.z_ext[i][2:] - self.z_ext[i][:-2]) \
                         / self.dx_ext_2cell[i] )**(1/6.)
             )
        print("FCN YO!")
        if self.S0 is not None:
            self.update_z_ext_0() # <-- Ursprüngliche Quelle
        # DEBUG
        i=0
        #print( ((self.z_ext[i][2:] - self.z_ext[i][:-2]) \
        #                 / self.dx_ext_2cell[i] )**(1/6.)
        #     )
        #print("FCN Whoooooa!")
        # !!!C0!!!
        # NOT YET UPDATED
        # KEEPING self.dx_ext_2cell INSTEAD OF USING MORE PRECISE OPTION
        # But this is fine so long as all gradients may be approximated
        # to be linear.
        # TRY: not network
        # EXCEPT: has a network
        # PROBABLY NOT THE BEST WAY TO DO THIS
        try:
            self.dzdx_0_16 = np.abs( (self.z_ext[2:] - self.z_ext[:-2]) \
                             / self.dx_ext_2cell )**(1/6.)
        except:
            # Inefficient: repeats whole segment calc when just the b.c.
            # change is needed
            dzdx_0_16_list = []
            for i in range(len( self.z_ext) ):
                # ADD FUNCTIONALITY TO LOOP OVER Q AND WEIGHT BY IT
                dzdx_0_16_list.append( 
                          np.abs( (self.z_ext[i][2:] - self.z_ext[i][:-2]) \
                         / self.dx_ext_2cell[i] )**(1/6.)
                )
                # UMMMMM WE NEED WEIGHTING BY Q
                # I GUESS I WILL WRITE A NEW FUNC TO PASS IT
            # PLACEHOLDER TO SEE IF THINGS JUST WORK, THOUGH
            # INSTEAD OF STRAIGHT MEAN, GIVE DISCHARGE WEIGHTING
            self.dzdx_0_16 = np.mean(dzdx_0_16_list)
                
                
        # UPDATE THIS TO INCLUDE THE DIV BY self.dx_ext_2cell THAT USED
        # TO BE IN C0?
        # OR JUST RETURN CODE TO HOW IT USED TO BE, SO self.dx_ext_2cell
        # ONCE MORE BE PART OF C0?
        # CURRENTLY, THE self.dx_ext_2cell IS APPLIED EXTERNALLY TO C1
        self.C1 = self.C0 * self.dzdx_0_16 * self.Q / self.B
        # Handling C1 for networked rivers
        # Need to link the two segments without skipping the channel head
        # DOESN'T SEEM TO CHANGE ANYTHING!
        # Looks right when both are 0! Any accidental inclusion of its own
        # ghost-node Qs,in?
        if len(self.downstream_segment_IDs) > 0:
            # !!!C0!!!
            # LEAVING AS-IS.
            # ADDRESSING CHANGES TO C0 IN LATER INSTANCES, NOT HERE
            # but now need things to work with lists of arrays]
            # LIKELY AGAIN NEED DISCHARGE WEIGHTING HERE; HACKING TO 
            # SEE IF I CAN GET IT TO WORK
            # JUST SETTING EVERYTHING TO PICK THE FIRST ARRAY IN THE LIST
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.C1[-1] = self.C0 \
                          * (np.abs(self.z_ext[0][-2] - self.z_ext[0][-1]) \
                                   /self.dx[-1])**(1/6.) \
                          * self.Q[-1] / self.B[-1]
        # This compute_coefficient_timeone matters! The above doesn't!!!! (Maybe.)
        # WORK HERE. If turns to 0, fixed. But why? Stays at initial profile?
        if len(self.upstream_segment_IDs) > 0:
            # !!!C0!!!
            # LEAVING AS-IS.
            # ADDRESSING CHANGES TO C0 IN LATER INSTANCES, NOT HERE
            # NEED DISCHARGE WEIGHTING HERE, BUT HACKING TO SEE IF I CAN
            # MAKE IT WORK
            # JUST SETTING EVERYTHING TO PICK THE FIRST ARRAY IN THE LIST 
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.C1[0] = self.C0 \
                          * (np.abs(self.z_ext[0][1] - self.z_ext[0][0]) \
                                   /self.dx[0])**(1/6.) \
                          * self.Q[0] / self.B[0]

    def set_z_bl(self, z_bl):
        """
        Set the right-hand Dirichlet boundary conditions, i.e. the base level,
        given in the variable "z_bl" (elevation, base level)
        """
        self.z_bl = z_bl
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # UPDATE
        self.z_ext[0][-1] = self.z_bl
        
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
        self.z_bl = 0
        self.bcr = self.z_bl * ( self.C1[-1] / self.dx_ext_2cell[0][-1] * 7/3. \
                       * (1/self.dx_ext[0][-2] + 1/self.dx_ext[0][-1])/2. \
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
        # !!!C0!!!
        # UPDATED BUT JUST USING dx_ext_2cell
        #self.bcl = self.dx_ext_2cell[0] * self.S0 * \
        #                    -self.C1[0] * ( 7/3./self.dx_ext[0]
        #                    - self.dQ[0]/self.Q[0]/self.dx_ext_2cell[0] )
        # !!!C0!!!
        # Probably not so easy to update as just updating C1
        # BECAUSE IT IS CHANGING THE RHS
        self.bcl = self.dx_ext_2cell[0][0] * self.S0 * \
                            self.C1[0] / self.dx_ext_2cell[0][0] \
                            * ( 7/3./self.dx_ext[0][0]
                            - self.dQ[0]/self.Q[0]/self.dx_ext_2cell[0][0] )

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
        self.right[0] = -self.C1[0] / self.dx_ext_2cell[0][0] * 7/3. \
                         * (1/self.dx_ext[0][0] + 1/self.dx_ext[0][1])

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
        self.left = -self.C1 / self.dx_ext_2cell[0] \
                        * ( (7/3.)/self.dx_ext[0][:-1]
                        - self.dQ/self.Q/self.dx_ext_2cell[0] )
        self.center = -self.C1 / self.dx_ext_2cell[0] \
                              * ( (7/3.)
                              * (-1/self.dx_ext[0][:-1]
                                 -1/self.dx_ext[0][1:]) ) \
                                 + 1.
        self.right = -self.C1 / self.dx_ext_2cell[0] \
                              * ( (7/3.)/self.dx_ext[0][1:] # REALLY?
                                  + self.dQ/self.Q/self.dx_ext_2cell[0] )
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
            if len(lp.upstream_segment_IDs) == 0:
                pass
            elif len(lp.upstream_segment_IDs) == 1:
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
                    C0 = upseg.k_Qs * upseg.intermittency \
                            / ((1-upseg.lambda_p) * upseg.sinuosity**(7/6.)) \
                            * self.dt / lp.dx_ext_2cell[0]
                    # From earlier
                    #C0 = upseg.C0[-1] # Should be consistent
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    dzdx_0_16 = ( np.abs(lp.z_ext[0][1] - lp.z_ext[0][0])
                                  / (lp.dx_ext[0][0]))**(1/6.)
                    C1 = C0 * dzdx_0_16 * upseg.Q[-1] / lp.B[0]
                    left_new = -C1 * 7/6. * 2 / lp.dx_ext[0][0]
                    self.LHSblock_matrix[row, col] = left_new
            elif len(lp.upstream_segment_IDs) > 1:
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


    def add_block_diagonal_matrix_downstream_boundary_conditions(self):
        for lp in self.list_of_LongProfile_objects:
            for ID in lp.downstream_segment_IDs:
                # Space to edit
                col = self.block_start_absolute[self.IDs == ID][0]
                row = self.block_end_absolute[self.IDs == lp.ID][0]
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
                        * self.dt / lp.dx_ext_2cell[0][-1]
                # !!!!!!!!!!!!!!!!!!!!!! [0] z_ext, dx_ext
                dzdx_0_16 = ( np.abs(lp.z_ext[0][-2] - lp.z_ext[0][-1])
                              / (lp.dx_ext[0][0]))**(1/6.)
                C1 = C0 * dzdx_0_16 * lp.Q[-1] / downseg.B[0]
                right_new = -C1 * 7/6. * 2 / lp.dx_ext[0][-1]
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
            lp.x_ext = np.max( (1, len(lp.upstream_segment_IDs)) ) * \
                [ np.concatenate( [_nan1, lp.x, _nan1] ) ]

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
        for lp in self.list_of_LongProfile_objects:
            _idx = 0
            # SET UPSTREAM BOUNDARIES: INTERNAL
            for upseg_ID in lp.upstream_segment_IDs:
                upseg = self.list_of_LongProfile_objects[upseg_ID]
                lp.x_ext[_idx][0] = upseg.x[-1]
                _idx += 1
            # SET DOWNSTREAM BOUNDARIES: INTERNAL
            _idx = 0
            for downseg_ID in lp.downstream_segment_IDs:
                downseg = self.list_of_LongProfile_objects[downseg_ID]
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
        # z_ext messup happens after this.

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
        for lp in self.list_of_LongProfile_objects:
            # SET UPSTREAM BOUNDARIES: INTERNAL
            _idx = 0
            for upseg_ID in lp.upstream_segment_IDs:
                upseg = self.list_of_LongProfile_objects[upseg_ID]
                lp.z_ext[_idx][0] = upseg.z[-1]
                _idx += 1
            # SET DOWNSTREAM BOUNDARIES: INTERNAL
            _idx = 0
            for downseg_ID in lp.downstream_segment_IDs:
                downseg = self.list_of_LongProfile_objects[downseg_ID]
                lp.z_ext[_idx][-1] = downseg.z[0]
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
        
        if S0 is not None and Q_s_0 is not None:
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

        # S0 might be iterable even if Q_s_0 be not
        try:
            iter(S0)
            S0 = np.array(S0).squeeze() # Enforce numpy array
        except:
            _is_scalar=True
        
        # FIFTH: Set S0 and z_ext[0]
        _idx = 0
        for ID in self.list_of_channel_head_segment_IDs:
            lp = self.list_of_LongProfile_objects[ID]
            if _is_scalar:
                lp.S0 = S0
            else:
                lp.S0 = S0[_idx]
            _idx += 1
            # Hard-coding: Expecting only one segment in list
            # Because this is just for the channel-head segments
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
                    B = None,
                    overwrite=False
                    ):
        #print( locals.keys() )
        """
        Run only once, at beginning of program.
        """
        
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
            lp.set_x( x = x[i], verbose=False )
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
            # Need to manage dQ in network
            # TO DO HERE
            lp.set_Q( Q = Q[i] )
            # The other GRLP-ey stuff
            lp.set_intermittency( 1 ) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # These all check out -- no problems with network.
            lp.basic_constants()
            lp.bedload_lumped_constants()
            lp.set_hydrologic_constants()
            lp.set_niter()
            #lp.set_z_bl(z1)
            lp.set_B( B = B[i] )
            # COULD WRITE A FUNCTION AROUND THIS
            # BUT I REALLY WANT TO REWRITE MORE IN TERMS OF SOURCE/SINK
            # DO SOMETHING HERE !!!!!
            lp.set_uplift_rate( 0 ) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        
        # Generate arrays of z based on external values
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " )
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " )
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " )
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " )
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " )
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " )
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " )
        self.create_z_ext_lists()
        self.update_z_ext_internal()
        self.update_z_ext_external_upstream( S0 = S0, Q_s_0 = Q_s_0 )  # b.c.
        self.update_z_ext_external_downstream( z_bl )                   # b.c.

        # DEBUG
        lp = self.list_of_LongProfile_objects[0]
        i = 0
        print( lp.z_ext )
        print( lp.x_ext )
        print( np.abs( (lp.z_ext[i][2:] - lp.z_ext[i][:-2]) \
                 / lp.dx_ext_2cell[i] )**(1/6.) )

    def evolve_threshold_width_river_network(self, nt=1, dt=3.15E7):
        """
        Solve the triadiagonal matrix through time, with a given
        number of time steps (nt) and time-step length (dt)
        """
        self.nt = nt
        self.dt = dt
        self.update_z_ext_internal()
        # self.dt is decided earlier
        for ti in range(int(self.nt)):
            for lp in self.list_of_LongProfile_objects:
                # Debug
                i=0
                print( np.abs( (lp.z_ext[i][2:] - lp.z_ext[i][:-2]) \
                         / lp.dx_ext_2cell[i] )**(1/6.) )
                lp.zold = lp.z.copy()
                print( "HEY!" )
                lp.build_LHS_coeff_C0(dt=self.dt)
                print( np.abs( (lp.z_ext[i][2:] - lp.z_ext[i][:-2]) \
                         / lp.dx_ext_2cell[i] )**(1/6.) )
                lp.zold = lp.z.copy()
                print( "Ho!" )
                lp.compute_coefficient_time_varying() # <-- Quelle der Problem
                print( np.abs( (lp.z_ext[i][2:] - lp.z_ext[i][:-2]) \
                         / lp.dx_ext_2cell[i] )**(1/6.) )
                lp.zold = lp.z.copy()
                print( "Whooooa!" )
                print( lp )
            for lp in self.list_of_LongProfile_objects:
                #print lp.C1
                lp.build_matrices()
                # Debug
                i=0
                print( np.abs( (lp.z_ext[i][2:] - lp.z_ext[i][:-2]) \
                         / lp.dx_ext_2cell[i] )**(1/6.) )
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
                self.update_z_ext_internal()
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
                    # ??????????????????????????????????
                    # THIS SEEMS POTENTIALLY PROBLEMATIC
                    # ??????????????????????????????????
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    lp.z_ext[0][1:-1] = \
                                    out[idx:idx+self.list_of_segment_lengths[i]]
                    idx += +self.list_of_segment_lengths[i]
                    i += 1
            self.update_z_ext_internal()
            self.t += self.dt # Update each lp z? Should make a global class
                              # that these both inherit from
            for lp in self.list_of_LongProfile_objects:
                lp.t = self.t
            i = 0
            idx = 0
            for lp in self.list_of_LongProfile_objects:
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                lp.z = lp.z_ext[0][1:-1].copy()
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

"""
# Test whether the lists have populated
for lp in self.list_of_LongProfile_objects: print(lp.x_ext)
print("")
for lp in self.list_of_LongProfile_objects: print(lp.z_ext)
"""

