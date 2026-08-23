import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
import planet_data as pd
from sys import path
import tools as tools
import pyvista as pv
from TwoBody import TwoBodyPropagator as tb
from OrbitProp import OrbitPropagator as OP
cb = pd.earth

if __name__ == "__main__":
    # r_magA = cb["radius"]
    # v_magA = np.sqrt(cb["mu"]/r_magA)

    # r0A = [r_magA*1.1, 0, 0] #initial position
    # v0A = [0,v_magA,0] #initial velocity
    t_span = (0,2.36e6)#seconds

    # r_magB =cb["radius"]+1000
    # v_magB = np.sqrt(cb["mu"]/r_magB)*1.35
    # r0B = [r_magB, 0, 0] #initial position
    # v0B = [0,v_magB,0.4]
    
    # r_magC = cb["radius"]
    # v_magC = np.sqrt(cb["mu"]/r_magC)

    # r0C = [r_magC*1.2, 120, 0] #initial position
    # v0C = [0,v_magC*0.9,v_magC*0.9] #initial velocity

    # opA = OP(r0A, v0A, t_span, cb)
    # opA.propagate_orbit()

    # opB = OP(r0B,v0B,t_span,cb)
    # opB.propagate_orbit()

    # opC = OP(r0C, v0C, t_span,cb)
    # opC.propagate_orbit()
    # labels = ["A","B","C"]
    # tools.plot([opA.rs,opB.rs,opC.rs],labels)
    # c0  =[cb["radius"] + 414,0.0006189,51.6393,0.0,234.1955,105.6372]
    # op0 = OP(c0, t_span,coes = True)
    # op0.propagate_orbit()
    # tools.plot([op0.rs], "a")
    neutron = pd.neutron_star
    moon = pd.moon
    r0 = [10000.0, 0.0, 0.0]  # Total distance 10,000 km
    v0 = [0.0, 6096.7, 0.0]  # km/s
    state0 = r0 + v0
    t_span = (0, 31)

    op = tb(state0,t_span,2.7846e30,2.7846e30,neutron,neutron)
    op.propagate_orbit()
    print(op.rs)
    tools.plotTwoBody(op.rs,neutron,neutron)