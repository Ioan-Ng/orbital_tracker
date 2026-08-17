import numpy as np
import math as m
import matplotlib.pyplot as plt
import planet_data as pd
from mpl_toolkits.mplot3d import Axes3D

def plot(rs,labels,show_plot = True, save_plot = False, Title ="Multiple Orbits",cb = pd.earth):
        fig = plt.figure(figsize=(18,6))
        ax = fig.add_subplot(111, projection = "3d")
        max = 0 
        #plot trajectory
        n  = 0 
        for r in rs:
            ax.plot(r[:,0], r[:,1], r[:,2], "r", label = labels[n],color='red', linewidth=2.5, zorder=10, clip_on=False) #the w is the colour, white
            n += 1
            ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], "wo", label="Initial Position")
            ax.plot([r[len(r)-1, 0]], [r[len(r)-1, 1]], [r[len(r)-1, 2]], "go", label="Final    Position")
            max_val = np.max(np.abs(r))
            if max_val >= max:
                max = max_val
        #plot the central body, the sphere plotting fucntion from https://stackoverflow.com/questions/11140163/plotting-a-3d-cube-a-sphere-and-a-vector
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = cb['radius'] * np.outer(np.cos(u), np.sin(v))
        y = cb['radius'] * np.outer(np.sin(u), np.sin(v))
        z = cb['radius'] * np.outer(np.ones(np.size(u)), np.cos(v))
        # Force the orbit line to always render above the planet regardless of camera view
        # Transparent sphere option
        ax.plot_surface(x, y, z, color='blue', alpha=0.5, zorder=1)

  
        #plot the x y z axis
        l = cb['radius']*2
        
        x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
        u,v,w = [[1,0,0],[0,1,0],[0,0,1]] #just making the arrows one unit
        ax.quiver(x,y,z,u,v,w, color = "r", zorder = 4) #u,v,w is where the arrows of x y z end
        
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

def coesToRV(coes, deg = False, mu = pd.earth["mu"]):
    
    if deg:
        a,e,i,ta,aop,raan = coes #major axis, eccentricity, inclincation, true anaomoly, argument of periapsis , right ascension of the acsending node,
        i = np.degrees(i)
        ta = np.degrees(ta)
        aop = np.degrees(aop)
        raan = np.degrees(raan)
    else:
        a,e,i,ta,aop,raan = coes
    E = ecc_anomaly([ta,e], "tae")

    r_norm = a*(1-e**2)/(1+e*np.cos(ta))

    r_perif = r_norm*np.array([m.cos(ta),m.sin(ta),0])
    v_perif = m.sqrt(mu*a)/r_norm*np.array([-m.sin(E),m.cos(E)*m.sqrt(1-e**2),0])


    perif2eci = np.transpose(eci2perif(raan,aop,i))
    r = np.dot(perif2eci,r_perif)
    v = np.dot(perif2eci,v_perif)

    return r,v

def ecc_anomaly(arr,method, tol = 1e-8):
    if method == "newton":
        Me,e = arr
        if Me< np.pi/2: E0 = Me+e/2
        else: E0 = Me-e
        for n in range(400):
            ratio(E0-e*np.sin(E0)-Me)/(1-e*np.cos(E0))
            if abs(ratio) < tol:
                if n == 0: return E0
                else: return E1
            else:
                E1 = E0-ratio
                E0=E1
        return False
    elif method == "tae":
        ta,e = arr
        return 2*m.atan(np.sqrt((1-e)/1+e))*m.tan(ta/2)
    else:
        print("invalid")

def eci2perif(raan, aop, i):
    row0 = [-m.sin(raan)*m.cos(i)*m.sin(aop) + m.cos(raan)*m.cos(aop),m.cos(raan)*m.cos(i)*m.sin(aop)+m.sin(raan)*m.cos(aop),m.sin(i)*m.sin(aop)]
    row1 = [-m.sin(raan)*m.cos(i)*m.cos(aop) - m.cos(raan)*m.sin(aop),m.cos(raan)*m.cos(i)*m.cos(aop)-m.sin(raan)*m.sin(aop),m.sin(i)*m.cos(aop)]
    row2 = [m.sin(raan)*m.sin(i), -m.cos(raan)*m.sin(i),m.cos(i)]
    return np.array([row0,row1,row2])