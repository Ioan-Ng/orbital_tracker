import numpy as np
from scipy.integrate import solve_ivp
import Tools as tl
#Our class where we will parse orbit data, either in keplearian element form or r and v vectors.
#UNITS
#everything will be in meters, m
class Orbit_Propagator:
    def __init__(self, initial_state, t_span, keplearian_data, central_body,degrees,ta_check):
        self.ta_check = ta_check #if using true if using true anom but false if we are using mean anom
        self.degrees = degrees #if the angles come in degrees we want to change them to radians when converting to vectors
        self.central_body = central_body
        self.mu = self.central_body["mu"]
        self.t_span = t_span #desribes t0 and t_end e.g (0,100)
        #we first want to check if we need to convert our state variables from keplerian to vecotr form
        if keplearian_data:
            self.r0, self.v0 = tl.kpToVector(initial_state, self.mu, self.degrees,self.ta_check)
        #otherwise we can just keep our parsed data as is
        else:
            self.r0 = initial_state[:3]
            self.v0 = initial_state[3:]

        #now we must set up our initial y0 for our differntial equation
        self.r0 = np.array(self.r0)
        self.v0 = np.array(self.v0)
        y0 = self.r0.tolist() +self.v0.tolist()

    def propagate_orbit(self):
        sol = solve_ivp(self.universal_gravitation_diff_eq,self.t_span,self.y0,method = "RK45",t_eval = np.linspace(self.t_span[0], self.t_span[1],10e4))
        #get the numerical results from the differnetial equation
        
        #make them into vector form for pyvista plotting
        self.xs = self.sol.y[0]    
        self.ys = self.sol.y[1]
        self.zs = self.sol.y[2]

        self.rs = np.vstack((self.xs,self.ys,self.zs)).T
    def universal_gravitation_diff_eq(self, current_state,t):
        rx,ry,rz,vx,vy,vz = current_state

        #get the position vector r, we will dentoe the underline as r_
        r_ = np.array([rx,ry,rz])
        r_norm = np.linalg.norm(r_)

        #governing equation
        ax,ay,az = -self.mu*r_*(r_norm**3)

        return[vx,vy,vz,ax,ay,az]

