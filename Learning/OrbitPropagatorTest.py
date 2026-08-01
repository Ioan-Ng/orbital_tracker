import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
import planet_data as pd
from sys import path

from OrbitProp import OrbitPropagator as OP
cb = pd.sun

if __name__ == "__main__":
    r_mag = cb["radius"]
    v_mag = np.sqrt(cb["mu"]/r_mag)

    r0 = [r_mag*1.4, 0, 0] #initial position
    v0 = [0,v_mag,0] #initial velocity
    t_span = (0,1*24*60*600)#seconds

    op = OP(r0, v0, t_span, cb)
    op.propagate_orbit()
    op.plot()