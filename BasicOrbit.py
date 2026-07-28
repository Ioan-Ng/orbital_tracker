import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

#this fuctnion is just a tempate for getting the values we need to parse into the ode solver that is part of scipy libary


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
    r0 = [r_mag, 0, 0] #initial position
    v0 = [0,v_mag,0] #initial velocity


    time_span = 1000 #seconds
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
    