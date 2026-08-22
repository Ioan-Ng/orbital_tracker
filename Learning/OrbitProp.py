import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import planet_data as pd
from mpl_toolkits.mplot3d import Axes3D
import tools as tools 
class OrbitPropagator:
    def __init__(self, state0, t_span,coes = False, cb = pd.earth):
        if coes:
            self.r0,self.v0 = tools.coesToRV(state0, deg = False,mu = cb["mu"])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]
        
        
        self.t_span = t_span
        self.cb = cb

        self.y0 = self.r0.tolist() + self.v0.tolist()
        
    def propagate_orbit(self):

        
        self.sol = solve_ivp(self.diff_eq, self.t_span, self.y0, method = "DOP853", t_eval=np.linspace(0, self.t_span[1], 86400)) 
        self.xs = self.sol.y[0]    
        self.ys = self.sol.y[1]
        self.zs = self.sol.y[2]

        self.rs = np.vstack((self.xs,self.ys,self.zs)).T
    def diff_eq(self,t,y):
        rx,ry,rz,vx,vy,vz = y

        r = np.array([rx,ry,rz])
        r_mag = np.linalg.norm(r)
        ax,ay,az = -r*self.cb['mu']/(r_mag**3)
        return[vx,vy,vz,ax,ay,az]#
    
    
