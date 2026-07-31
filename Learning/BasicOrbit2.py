
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D


def plot(rs):
    fig = plt.figure(figsize=(18,6))
    ax = fig.add_subplot(111, projection = "3d")
    
    #plot trajectory
    ax.plot(rs[:,0], rs[:,1], rs[:,2], "b", label = "Trajectory") #the w is the colour, white
    ax.plot([rs[0, 0]], [rs[0, 1]], [rs[0, 2]], "wo", label="Initial Position")
    ax.plot([rs[len(rs)-1, 0]], [rs[len(rs)-1, 1]], [rs[len(rs)-1, 2]], "go", label="Final Position")

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



#initial conditions
mu = 398600.4418 #km^3/s^2
r_mag = 6371 #km
v_mag = np.sqrt(mu/r_mag) #km/s

#initial position and velocity vectors
r0 = [r_mag+2000, 0, 0] #initial position
v0 = [0,v_mag,0] #initial velocity


time_span = (0,1*24*60*60) #seconds

y0 = r0 + v0

sol = solve_ivp(diff_eq, time_span, y0, method = "DOP853",args=(mu,), t_eval=np.linspace(0, 1*24*60*60, 1000000)) 
print(sol)
ts = sol.t       
xs = sol.y[0]    
ys = sol.y[1]
zs = sol.y[2]

rs = np.vstack((xs,ys,zs)).T

plt.plot(ts,xs, "r")
plt.ylabel("displacement")
plt.xlabel("time")
plt.show()
plot(rs)