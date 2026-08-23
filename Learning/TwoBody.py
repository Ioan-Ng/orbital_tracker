import numpy as np
from scipy.integrate import solve_ivp
import tools as tools 

class TwoBodyPropagator:
    def __init__(self,state, t_span,ma,mb):
        self.r0 = state[:3]
        self.v0 = state[3:]
        self.ma = ma
        self.mb = mb
        self.t_span = t_span
        self.y0 = np.concatenate([self.r0, self.v0])

    def diff_eq(self,t,y):
        G = 6.674e-20
        rx,ry,rz,vx,vy,vz = y

        r = np.array([rx,ry,rz])
        r_mag = np.linalg.norm(r)

        ax,ay,az = -r*G*(self.ma+self.mb)/(r_mag**3)
        return[vx,vy,vz,ax,ay,az]

    def propagate_orbit(self):

        
        self.sol = solve_ivp(self.diff_eq, self.t_span, self.y0, method = "DOP853", t_eval=np.linspace(0, self.t_span[1], 86400)) 
        self.xs = self.sol.y[0]    
        self.ys = self.sol.y[1]
        self.zs = self.sol.y[2]
        self.rs = np.vstack((self.xs,self.ys,self.zs)).T