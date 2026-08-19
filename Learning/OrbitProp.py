import numpy as np
import pyvista as pv

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import planet_data as pd
from mpl_toolkits.mplot3d import Axes3D
import tools as tools 
class OrbitPropagator:
    def __init__(self, state0, t_span,coes = False, cb = pd.earth):
        if coes:
            self.r0,self.v0 = tools.coesToRV(state0, deg = True,mu = cb["mu"])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]
        
        
        self.t_span = t_span
        self.cb = cb

        self.y0 = self.r0.tolist() + self.v0.tolist()
        
    def propagate_orbit(self):

        
        self.sol = solve_ivp(self.diff_eq, self.t_span, self.y0, method = "DOP853", t_eval=np.linspace(0, self.t_span[1], 1000000)) 
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
    


    def plot(self, show_plot=True, save_plot=False, Title="Orbit"):
        # Initialize high-performance GPU plotter
        plotter = pv.Plotter(off_screen=not show_plot)
        plotter.title = Title

        # 1. Trajectory line mesh
        trajectory = pv.PolyData(self.rs)
        # Add scalar array for point order to draw a continuous line path
        trajectory.lines = np.hstack([[len(self.rs)] + list(range(len(self.rs)))])
        plotter.add_mesh(trajectory, color="blue", line_width=3, label="Trajectory")

        # 2. Initial and Final Position Markers
        plotter.add_points(np.array([self.rs[0]]), color="white", point_size=12, render_points_as_spheres=True, label="Initial! Position")
        plotter.add_points(np.array([self.rs[-1]]), color="green", point_size=12, render_points_as_spheres=True, label="Final Position")

        # 3. Central Body (GPU-rendered parametric sphere)
        central_body = pv.Sphere(radius=self.cb['radius'], center=(0, 0, 0), theta_resolution=60, phi_resolution=60)
        plotter.add_mesh(central_body, color="cornflowerblue", show_edges=False, smooth_shading=True)

        # 4. Coordinate Axes (Replaces ax.quiver with clean 3D arrows)
        axis_length = self.cb['radius'] * 2.0
        plotter.add_axes_at_origin(x_color="red", y_color="green", z_color="blue", line_width=2, labels_off=False)

        # 5. Labels & Visual Environment
        plotter.add_legend(bcolor=(0.1, 0.1, 0.1, 0.5), face=None)
        plotter.show_grid(xlabel="X (km)", ylabel="Y (km)", zlabel="Z (km)")
        plotter.set_background("black")  # Ideal for space/orbital views

        # Equal aspect ratio scaling (built into PyVista by default)
        max_val = np.max(np.abs(self.rs))
        plotter.camera.clipping_range = (0.1, max_val * 10)

        # Handle Saving/Displaying
        if save_plot:
            plotter.screenshot(f"{Title}.png")
            
        if show_plot:
            plotter.show()
        else:
            plotter.close()