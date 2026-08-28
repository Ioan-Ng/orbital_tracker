import numpy as np
from scipy.integrate import solve_ivp
import Tools as tl
#Our class where we will parse orbit data, either in keplearian element form or r and v vectors.
class Orbit_Propagator():
    def __init__(self, initial_state, t_span, keplearian_data, central_body):
        self.central_body = central_body
        self.mu = self.central_body["mu"]
        self.t_span = t_span #desribes t0 and t_end e.g (0,100)
        #we first want to check if we need to convert our state variables from keplerian to vecotr form
        if keplearian_data:
            self.r0, self.v0 = tl.kpToVector(self.mu)
        #otherwise we can just keep our parsed data as is
        else:
            self.r0 = initial_state[:3]
            self.v0 = initial_state[3:]

        #now we must set up our initial y0 for our differntial equation

        y0 = self.r0.tolist() +self.v0.tolist()

    def propagate_orbit(self):
        sol = solve_ivp(self.universal_gravitation_diff_eq,self.t_span,method = "RK45",t_eval = np.linspace(self.t_span[0], self.t_span[1],3,10e4))
        #get the numerical results from the differnetial equation
        self.rs = sol.y
    def universal_gravitation_diff_eq(self, current_state,t):
        rx,ry,rz,vx,vy,vz = current_state

        #get the position vector r, we will dentoe the underline as r_
        r_ = np.array([rx,ry,rz])
        r_norm = np.linalg.norm(r_)

        #governing equation
        ax,ay,az = -self.mu*r_*(r_norm**3)

        return[vx,vy,vz,ax,ay,az]