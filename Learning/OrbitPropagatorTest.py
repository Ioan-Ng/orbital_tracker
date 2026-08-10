import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
import planet_data as pd
from sys import path
import tools as tools
from OrbitProp import OrbitPropagator as OP
cb = pd.sun

if __name__ == "__main__":
    r_magA = cb["radius"]
    v_magA = np.sqrt(cb["mu"]/r_magA)

    r0A = [r_magA*1.1, 0, 0] #initial position
    v0A = [0,v_magA,0] #initial velocity
    t_span = (0,1*24*60*60)#seconds

    r_magB =cb["radius"]
    v_magB = np.sqrt(cb["mu"]/r_magB)
    r0B = [r_magB*1.5, 0, r_magB*1.5] #initial position
    v0B = [0,v_magB*0.8,0]
    
    opA = OP(r0A, v0A, t_span, cb)
    opA.propagate_orbit()

    opB = OP(r0B,v0B,t_span,cb)
    opB.propagate_orbit()
    labels = ["A","B"]
    tools.plot([opA.rs,opB.rs],labels)