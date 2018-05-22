import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse import spdiags, identity
from scipy.sparse.linalg import spsolve, isolve
from scipy.stats import linregress
import warnings
import sys

class LongProfile(object):

    def __init__(self):
        self.z = None
        self.x = None
        self.A = None
        self.Q = None
        self.B = None
        self.dx_ext = None
        self.dx_2cell = None
        self.Q_s_0 = None
        self.sinuosity = 1.
        self.intermittency = 0.01
        self.t = 0
        #self.basic_constants()

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
        self.I = I
                
    def set_x(self, x=None, x_ext=None, dx=None, nx=None, x0=None):
        """
        Set x directly or calculate it.
        Pass one of three options:
        x alone
        x_ext alone (this will also define x)
        dx, nx, and x0
        """
        if x is not None:
            self.x = np.array(x)
            diff = np.diff(self.x)
            dx_mean = np.mean(diff)
            if (diff == dx_mean).all():
                self.dx = dx_mean
                self.dx_isscalar = False
            else:
                #sys.exit("Uniform x spacing required")
                self.dx = diff
                self.dx_2cell = self.x[2:] - self.x[:-2]
                self.dx_isscalar = True
        elif x_ext is not None:
            self.x_ext = np.array(x_ext)
            self.x = x_ext[1:-1]
            diff = np.diff(self.x_ext)
            dx_mean = np.mean(diff)
            if (diff == dx_mean).all():
                self.dx_ext = dx_mean
                self.dx_isscalar = True
            else:
                #sys.exit("Uniform x spacing required")
                self.dx_ext = diff
                self.dx_ext_2cell = self.x_ext[2:] - self.x_ext[:-2]
                self.dx_2cell = self.x[2:] - self.x[:-2]
                self.dx = np.diff(self.x)
                self.dx_isscalar = False
        elif (dx is not None) and (nx is not None) and (x0 is not None):
            self.x = np.arange(x0, x0+dx*nx, dx)
            self.dx = dx
            self.dx_isscalar = True
            self.x_ext = np.hstack((self.x[0]-dx, self.x, self.x[-1]+dx))
        else:
            sys.exit("Need x OR x_ext OR (dx, nx, x0)")
        self.nx = len(self.x)
        if (nx is not None) and (nx != self.nx):
            warnings.warn("Choosing x length instead of supplied nx")
            
    def set_z(self, z=None, z_ext=None, S0=None, z1=0):
        """
        Set z directly or calculate it
        S0 = initial slope (negative for flow from left to right)
        z1 = elevation value at RHS
        """
        if z:
            self.z = z
            self.z_ext = np.hstack((2*z[0]-z[1], z, 2*z[-1]-z[-2]))
        elif z_ext:
            self.z_ext = z_ext
            self.z = z_ext[1:-1]
        elif self.x.any() and self.x_ext.any() and S0:
            self.z = self.x * S0 + (z1 - self.x[-1] * S0)
            self.z_ext = self.x_ext * S0 + (z1 - self.x_ext[-1] * S0)
        else:
            sys.exit("Error defining variable")
        #self.dz = self.z_ext[2:] - self.z_ext[:-2] # dz over 2*dx!
            
    def set_A(self, A=None, A_ext=None, k_xA=None, P_xA=None):
        """
        Set A directly or calculate it
        """
        if A:
            self.A = A
            self.A_ext = np.hstack((2*A[0]-A[1], A, 2*A[-1]-A[-2]))
            self.dA = self.A_ext[2:] - self.A_ext[:-2]
        elif A_ext:
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
        self.dA = self.A_ext[2:] - self.A_ext[:-2] # dA over 2*dx!
                                                   # Not sure if used.

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
            Q_ext = np.hstack((2*self.Q[0]-self.Q[1], self.Q, 2*self.Q[-1]-self.Q[-2]))
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
        self.dQ = Q_ext[2:] - Q_ext[:-2] # dQ over 2*dx!
        # Keep sediment supply tied to water supply, except
        # by changing S_0, to only turn one knob for one change (Q/Qs)
        if update_Qs_input:
            if self.Q_s_0:
                self.set_Qs_input_upstream(self.Q_s_0)

    def set_B(self, B=None, B_ext=None, k_xB=None, P_xB=None):
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
            B_ext = np.hstack((2*self.B[0]-self.B[1], self.B, 2*self.B[-1]-self.B[-2]))
        elif B_ext is not None:
            self.B = B_ext[1:-1]
        elif k_xB and self.x.any() and self.x_ext.any():
            self.B = k_xB * self.x**P_xB
            B_ext = k_xB * self.x_ext**P_xB
            self.k_xB = k_xB
            self.P_xB = P_xB
        self.dB = B_ext[2:] - B_ext[:-2] # dB over 2*dx!
        
    def set_uplift_rate(self, U):
        """
        Uplift rate if positive -- or equivalently, rate of base-level fall
        Subsidence (or base-level rise) accomplished by negative uplift
        """
        self.U = -U # not sure this is the best -- flipping the sign

    def set_niter(self, niter=3):
        self.niter = niter
        
    def set_Qs_input_upstream(self, Q_s_0):
        self.Q_s_0 = Q_s_0
        # Q[0] is centerpoint of S?
        self.S0 = -((1/self.k_Qs) * (Q_s_0/self.Q[0]))**(6/7.)
        if self.dx_isscalar:
            self.z_ext[0] = self.z[0] - self.S0 * self.dx
            #self.z_ext[0]
        else:
            # Give upstream cell the same width as the first cell in domain
            self.z_ext[0] = self.z[0] - self.S0 * self.dx_ext[0]
        
    def update_z_ext_0(self):
        if self.dx_isscalar:
            self.z_ext[0] = self.z[0] - self.S0 * self.dx
        else:
            # Give upstream cell the same width as the first cell in domain
            self.z_ext[0] = self.z[0] - self.S0 * self.dx_ext[0]

    def compute_coefficient_time_varying(self):
        self.update_z_ext_0()
        if self.dx_isscalar:
            self.dzdt_0_16 = np.abs( (self.z_ext[2:] - self.z_ext[:-2]) \
                             / (2*self.dx) )**(1/6.)
        else:
            self.dzdt_0_16 = np.abs( (self.z_ext[2:] - self.z_ext[:-2]) \
                             / self.dx_ext_2cell )**(1/6.)
        self.C1 = self.C0 * self.dzdt_0_16 * self.Q / self.B

    def set_z_bl(self, z_bl):
        self.z_bl = z_bl
        self.z_ext[-1] = self.z_bl

    def set_bcr_Dirichlet(self):
        #self.bcr_value = bcr
        if self.dx_isscalar:
            self.bcr = self.z_bl * ( self.C1[-1] * 7/6. \
                                     + self.dQ[-1]/self.Q[-1]/4. \
                                     - self.dB[-1]/self.B[-1]/4. )
                                     #+ self.z[-1]
        else:
            self.bcr = self.z_bl * ( self.C1[-1] * 7/3. \
                           * (-1/self.dx_ext[-2] - 1/self.dx_ext[-1]) \
                           + self.dQ[-1]/self.Q[-1]/self.dx_ext_2cell[0] \
                           - self.dB[-1]/self.B[-1]/self.dx_ext_2cell[0] )
        
    def set_bcl_Neumann_RHS(self):
        """
        Set the LHS boundary condition
        """
        #self.bcl = self.z[0] + 2*self.dx*self.S0*self.left[0]
        #self.bcl = 2*self.dx*self.S0*self.left[0]
        if self.dx_isscalar:
            self.bcl = -2 * self.dx * self.S0 * \
                       self.C1[0] * ( 7/6. - self.dQ[0]/self.Q[0]/4. \
                                           + self.dB[0]/self.B[0]/4.)
        else:
            # Give upstream cell the same width as the first cell in domain
            self.bcl = -self.dx_ext_2cell[0] * self.S0 * \
                                  self.C1[0] * ( 7/3./self.dx_ext[0]
                                - self.dQ[0]/self.Q[0]/self.dx_ext_2cell[0]
                                + self.dB[0]/self.B[0]/self.dx_ext_2cell[0] )# \
                                #* 107.60426753093884
                                            
    def set_bcl_Neumann_LHS(self):
        """
        from ghost node approach
        """
        if self.dx_isscalar:
            #self.right[0] = -2 * self.C1[0] * 7/6. * self.Q[0]
            self.right[0] = -2 * self.C1[0] * 7/6.
            #self.right[0] = self.left[0] + self.right[0] # should be the same as the above
        else:
            # PROBLEM POSSIBLY HERE -- CHECK THIS GHOST NODE APPROACH RIGOROUSLY
            #
            self.right[0] = -self.C1[0] * ( (7/3.) \
                                        * (-1/self.dx_ext[0] \
                                           -1/self.dx_ext[1]) )
            
    
    def evolve_threshold_width_river(self, nt=1, dt=3.15E7):
        self.dt = dt
        self.nt = nt
        if self.dx_isscalar:
            self.C0 = self.k_Qs/(1-self.lambda_p) * self.sinuosity \
                      * self.intermittency * self.dt / self.dx**2
        else:
            # 250 = patch for dx = 500
            self.C0 = self.k_Qs/(1-self.lambda_p) * self.sinuosity \
                      * self.intermittency * self.dt / self.dx_ext_2cell
        for ti in range(int(self.nt)):
            self.zold = self.z.copy()
            self.set_z_bl(self.z_bl + self.U * self.dt)
            for i in range(self.niter):
                self.compute_coefficient_time_varying()
                if self.dx_isscalar:
                    self.left = -self.C1 * ( (7/6.) - self.dQ/self.Q/4. \
                                + self.dB/self.B/4.)
                    self.center = self.C1 * 2 * ( (7/6.) ) + 1.
                    self.right = -self.C1 * ( (7/6.) + self.dQ/self.Q/4. \
                                 - self.dB/self.B/4. )
                else:
                    self.left = -self.C1 * ( (7/3.)/self.dx_ext[:-1]
                                    - self.dQ/self.Q/self.dx_ext_2cell \
                                    + self.dB/self.B/self.dx_ext_2cell)
                    self.center = -self.C1 * ( (7/3.) \
                                          * (-1/self.dx_ext[:-1] \
                                             -1/self.dx_ext[1:]) ) \
                                             + 1.
                    self.right = -self.C1 * ( (7/3.)/self.dx_ext[1:]
                                    + self.dQ/self.Q/self.dx_ext_2cell \
                                    - self.dB/self.B/self.dx_ext_2cell )
                self.set_bcl_Neumann_LHS()
                self.set_bcl_Neumann_RHS()
                self.set_bcr_Dirichlet()
                self.left = np.roll(self.left, -1)
                self.right = np.roll(self.right, 1)
                # More boundary conditions: RHS Dirichlet
                #self.bcr = self.center[-1] * self.bcr_value \
                #           + self.left[-1] * self.bcr_value
                # TO DO: Fix Dirichlet b.c. here
                diagonals = np.vstack((self.left, self.center, self.right))
                offsets = np.array([-1, 0, 1])
                LHSmatrix = spdiags(diagonals, offsets, len(self.z), 
                                    len(self.z), format='csr')
                RHS = np.hstack((self.bcl+self.z[0], self.z[1:-1], self.bcr+self.z[-1]))
                #print np.mean(self.z)
                self.z_ext[1:-1] = spsolve(LHSmatrix, RHS)
                print self.bcl
            self.t += self.dt
            #print np.mean(self.z)
            self.z = self.z_ext[1:-1].copy()
            self.dz_dt = (self.z - self.zold)/self.dt
            #S = np.diff(self.z_ext)[1:] + np.diff(self.z_ext)[:-1]
            #self.b = 2.61 * self.Q * S**(7/6.)/
            self.Qs_internal = 1/(1-self.lambda_p) * np.cumsum(self.dz_dt)*self.B + self.Q_s_0
            self.update_z_ext_0()
        print self.bcl, self.bcr, self.right[0], self.C1[0], self.z[0]
        print np.mean(self.left/1E10), np.mean(self.center/1E10), np.mean(self.right/1E10)
        print np.mean(self.left/1E10)/np.mean(self.center/1E10)
        print self.bcl/self.C1[0], self.dzdt_0_16[0]
        print self.bcl
        print self.z[0]
        #print self.dzdt_0_16
    
    def analytical_threshold_width(self, P_xB=None, P_xQ=None, x0=None, x1=None, 
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
        if P_xB is None:
            P_xB = self.P_xB
        #print P_xB
        #print P_xQ
        #e = 1 + 6*(P_xB - P_xQ)/7.
        #self.zanalytical2 = (z1 - z0) * (self.x**e - x0**e)/(x1**e - x0**e) + z0
        self.P_a = 1 + 6*(P_xB - P_xQ)/7. # beta
        self.k_a = 1/(x1**self.P_a - x0**self.P_a) * (z1 - z0) # alpha
        self.c_a = z0 - x0**self.P_a/(x1**self.P_a - x0**self.P_a) * (z1 - z0) # gamma
        self.zanalytical = self.k_a * self.x**self.P_a + self.c_a
        return self.zanalytical
        
    def analytical_threshold_width_perturbation(self, P_xB=None, P_xQ=None, x0=None, x1=None, 
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
        if P_xB is None:
            P_xB = self.P_xB
        if U is None:
            U = self.U
        # Base state coefficients (no perturbation)
        #self.P_a = 1 + 6*(P_xB - P_xQ)/7. # beta
        #self.k_a = (z1 - z0)/(x1**self.P_a - x0**self.P_a) # alpha
        #self.c_a = z0 - x0**self.P_a/(x1**self.P_a - x0**self.P_a) * (z1 - z0) # gamma
        # Coefficients
        K = self.k_Qs * self.sinuosity * self.intermittency \
            / (1 - self.lambda_p) \
            * abs(self.k_a * self.P_a)**(1/6.) \
            * self.k_xQ / self.k_xB
        P = self.P_xQ - self.P_xB + (self.P_xB - 1.)/6.
        print P
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
        self.S = np.abs( (self.z_ext[2:] - self.z_ext[:-2]) / (2*self.dx) ) \
                 / self.sinuosity
        self.Q_s = self.k_Qs * self.intermittency * self.Q * self.S**(7/6.)

    def slope_area(self, verbose=False):
        self.S = np.abs( (self.z_ext[2:] - self.z_ext[:-2]) / (2*self.dx) )
        logS = np.log10(self.S)
        logA = np.log10(self.A)
        out = linregress(logA[1:-1], logS[1:-1]) # remove edge effects
        self.theta = -out.slope
        self.ks = 10.**out.intercept
        self.thetaR2 = out.rvalue**2.
        if verbose:
            print "Concavity = ", self.theta
            print "k_s = ", self.ks
            print "R2 = ", out.rvalue**2.
        
