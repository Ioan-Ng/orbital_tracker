import numpy as np
import matplotlib.pyplot as plt
import planet_data as pd
from mpl_toolkits.mplot3d import Axes3D

def plot(rs,labels,show_plot = True, save_plot = False, Title ="Multiple Orbits",cb = pd.sun):
        fig = plt.figure(figsize=(18,6))
        ax = fig.add_subplot(111, projection = "3d")
        max = 0 
        #plot trajectory
        n  = 0 
        for r in rs:
            ax.plot(r[:,0], r[:,1], r[:,2], "b", label = labels[n]) #the w is the colour, white
            n += 1
            ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], "wo", label="Initial Position")
            ax.plot([r[len(r)-1, 0]], [r[len(r)-1, 1]], [r[len(r)-1, 2]], "go", label="Final    Position")
            max_val = np.max(np.abs(r))
            if max_val >= max:
                max = max_val
        #plot the central body, the sphere plotting fucntion from https://stackoverflow.com/questions/11140163/plotting-a-3d-cube-a-sphere-and-a-vector
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = cb['radius']*np.cos(u)*np.sin(v)
        y = cb['radius']*np.sin(u)*np.sin(v)
        z = cb['radius']*np.cos(v)
        ax.plot_surface(x, y, z, color="r")

        #plot the x y z axis
        l = cb['radius']*2
        
        x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
        u,v,w = [[1,0,0],[0,1,0],[0,0,1]] #just making the arrows one unit
        ax.quiver(x,y,z,u,v,w, color = "b") #u,v,w is where the arrows of x y z end
        
        ax.set_xlim([-max,max])
        ax.set_ylim([-max,max])
        ax.set_zlim([-max,max])

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