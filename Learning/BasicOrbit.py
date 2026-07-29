import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
#this fuctnion is just a tempate for getting the values we need to parse into the ode solver that is part of scipy libary

def plot(rs):
    fig = plt.figure(figsize=(18,6))
    ax = fig.add_subplot(111, projection = "3d")
    
    #plot trajectory
    ax.plot(rs[:,0], rs[:,1], rs[:,2], "b", label = "Trajectory") #the w is the colour, white
    ax.plot([rs[0, 0]], [rs[0, 1]], [rs[0, 2]], "wo", label="Initial Position")

    #plot the central body, the sphere plotting fucntion from https://stackoverflow.com/questions/11140163/plotting-a-3d-cube-a-sphere-and-a-vector
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = r_mag*np.cos(u)*np.sin(v)
    y = r_mag*np.sin(u)*np.sin(v)
    z = r_mag*np.cos(v)
    ax.plot_surface(x, y, z, color="r")

    #plot the x y z axis
    l = r_mag*2
    
    x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
    u,v,w = [[1,0,0],[0,1,0],[0,0,1]] #just making the arrows one unit
    ax.quiver(x,y,z,u,v,w, color = "b") #u,v,w is where the arrows of x y z end
    max_val = np.max(np.abs(rs))
    ax.set_xlim([-max_val,max_val])
    ax.set_ylim([-max_val,max_val])
    ax.set_zlim([-max_val,max_val])

    ax.set_xlabel(["X (km)"])
    ax.set_ylabel(["Y (km)"])
    ax.set_zlabel(["Z (km)"])

    ax.set_aspect("equal")
    ax.set_title("Model of diff eq")


    plt.legend()
    plt.show()

def diff_eq(t,y,mu):
    #t is time, y is the state, it gives us the postiton and velocity at t and mu is the gravity thing
    rx,ry,rz,vx,vy,vz = y
    r = np.array([rx,ry,rz]) #r is our postion vecotr

    r_norm = np.linalg.norm(r) #r_norm is the magnitude of our position vector
    
    #n2l
    ax,ay,az = -r*mu/(r_norm**3)

    return [vx,vy,vz, ax,ay,az]

if __name__ == "__main__":
    #initial conditions
    mu = 398600.4418 #km^3/s^2
    r_mag = 7000 #km
    v_mag = np.sqrt(mu/r_mag) #km/s

    #initial position and velocity vectors
    r0 = [r_mag+1000, 0, 0] #initial position
    v0 = [0,v_mag,0] #initial velocity


    time_span = 8000 #seconds
    dt = 1 #time step in seconds


    #number of steps
    n_steps = int(np.ceil(time_span/dt))
    
    #initalise our arrays
    ys = np.zeros((n_steps, 6)) #pre allocate memory for efficiency 
    ts = np.zeros((n_steps,1))

    #intialise
    y0 = r0 + v0
    ys[0] = np.array(y0)
    step = 1 # since step = 0 is our inital values

    solver = ode(diff_eq)
    solver.set_integrator('lsoda')
    solver.set_initial_value(y0,0)
    solver.set_f_params(mu) #extra fucntion parameters that don't rely on the previous results, ie constants


    while solver.successful() and step<n_steps:
        solver.integrate(solver.t + dt)
        ts[step] = solver.t
        ys[step] = solver.y
        step += 1
        
    rs = ys[:,:3]
    plot(rs)

