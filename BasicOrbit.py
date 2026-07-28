import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


def diff_eq(t,y,mu):
    #t is time, y is the state, it gives us the postiton and velocity at t and mu is the gravity thing
    rx,ry,rz,vx,vy,vz = y
    r = np.array([rx,ry,rz]) #r is our postion vecotr

    r_norm = np.linalg.norm(r) #r_norm is the magnitude of our position vector
    
    #n2l
    ax,ay,az = -r*mu/(r_norm**3)

    return [vx,vy,vz, ax,ay,az]
