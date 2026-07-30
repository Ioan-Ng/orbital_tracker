import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import figure

def horizontal_spring_ode(t,y,omega):
    x, v = y
    dxdt = v

    #our ode is
    dvdt = -(omega**2)*x
    return[dxdt,dvdt]

omega = 2
t_span = (0, 10)          
initial_ys = [10.0, 0.0]    
sol = solve_ivp(horizontal_spring_ode, t_span, initial_ys, method = "RK45",args=(omega,), t_eval=np.linspace(0, 10, 100))
print(sol)
ts = sol.t       
xs = sol.y[0]    
vs = sol.y[1]    

figure(figsize=(10, 4))
plt.plot(ts,xs, "r")
plt.ylabel("displacement")
plt.xlabel("time")
plt.show()
