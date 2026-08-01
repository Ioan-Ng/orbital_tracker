import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import planet_data as pd
from mpl_toolkits.mplot3d import Axes3D
class OrbitPropagator:
    def __init__(self, r0, v0, t_span, dt, cb =pd.sun):
        self.r0 = r0
        self.v0 = v0
        self.t_span = t_span
        self.dt = dt
        self.cb = cb
        
        def propagate_orbit(self):

            self.y0 = self.r0 + self.v0
            self.sol = solve_ivp(self.diff_eq, self.t_span, self.y0, method = "DOP853",args=(mu,), t_eval=np.linspace(0, 1*24*60*60, 1000000)) 
        
        def diff_eq(self,t,y):
            rx,ry,rz,vx,vy,vz = y

            r = np.array([rx,ry,rz])
            r_mag = np.linalg.norm(r)

            ax,ay,ax = -r*self.cb['mu']/(r_mag**3)
            return[vx,vy,vz,ax,ay,az]#
        
        def plot(self, show_plot = False, save_plot = False, Title ="Orbit"):
            fig = plt.figure(figsize=(18,6))
            ax = fig.add_subplot(111, projection = "3d")
            
            #plot trajectory
            ax.plot(self.rs[:,0], self.rs[:,1], self.rs[:,2], "b", label = "Trajectory") #the w is the colour, white
            ax.plot([self.rs[0, 0]], [self.rs[0, 1]], [self.rs[0, 2]], "wo", label="Initial Position")
            ax.plot([self.rs[len(self.rs)-1, 0]], [self.rs[len(self.rs)-1, 1]], [self.rs[len(self.rs)-1, 2]], "go", label="Final Position")

            #plot the central body, the sphere plotting fucntion from https://stackoverflow.com/questions/11140163/plotting-a-3d-cube-a-sphere-and-a-vector
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x = self.cd['radius']*np.cos(u)*np.sin(v)
            y = self.cd['radius']*np.sin(u)*np.sin(v)
            z = self.cd['radius']*np.cos(v)
            ax.plot_surface(x, y, z, color="r")

            #plot the x y z axis
            l = self.cd['radius']*2
            
            x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
            u,v,w = [[1,0,0],[0,1,0],[0,0,1]] #just making the arrows one unit
            ax.quiver(x,y,z,u,v,w, color = "b") #u,v,w is where the arrows of x y z end
            max_val = np.max(np.abs(self.rs))
            ax.set_xlim([-max_val,max_val])
            ax.set_ylim([-max_val,max_val])
            ax.set_zlim([-max_val,max_val])

            ax.set_xlabel(["X (km)"])
            ax.set_ylabel(["Y (km)"])
            ax.set_zlabel(["Z (km)"])

            ax.set_aspect("equal")
            ax.set_title(Title)


            plt.legend()
            
            if show_plot:
                plt.show()
            if save_plot:
                plt.savefig(title+ ".png",dpi = 300)