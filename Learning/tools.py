import numpy as np
import math as m
import matplotlib.pyplot as plt
import planet_data as pd
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
import pandas as pds
from OrbitProp import OrbitPropagator as OP
def plot(
    rs,
    labels,
    show_plot=True,
    save_plot=False,
    Title="Multiple Orbits",
    cb=pd.earth,
):
    plotter = pv.Plotter()

    # 1. Central Body
    sphere = pv.Sphere(radius=cb["radius"], center=(0, 0, 0))
    plotter.add_mesh(
        sphere, color="blue", opacity=.4, label="Central Body", show_edges=False
    )

    max_val = cb["radius"]

    # 2. Plot Trajectories
    colors = [
    "red", "blue", "green", "orange", "yellow", "purple", "brown", "pink",
    "cyan", "magenta", "lime", "teal", "indigo", "violet", "gold", 
    "coral", "maroon", "navy", "olive", "turquoise"
]
    k = 0
    for n, r in enumerate(rs):
        
        points = np.asarray(r, dtype=np.float64)

        # FIXED: Use pv.MultipleLines to connect continuous points
        trajectory = pv.MultipleLines(points)

        # Plot line
        k += 1
        plotter.add_mesh(
            trajectory,
            color=colors[k%20],
            line_width=3,
            label=labels[n] if n < len(labels) else f"Orbit {n+1}",
        )

        # Initial Position
        start_pt = pv.PolyData(points[0:1])
        plotter.add_mesh(
            start_pt,
            color="white",
            point_size=12,
            render_points_as_spheres=True,
            label="Initial Position" if n == 0 else None,
        )

        # Final Position
        end_pt = pv.PolyData(points[-1:])
        plotter.add_mesh(
            end_pt,
            color="green",
            point_size=12,
            render_points_as_spheres=True,
            label="Final Position" if n == 0 else None,
        )

        curr_max = np.max(np.abs(points))
        if curr_max > max_val:
            max_val = curr_max

    # 3. Setup Scene
    plotter.add_axes_at_origin()
    plotter.add_title(Title)
    plotter.show_grid(
        xtitle="X (km)",
        ytitle="Y (km)",
        ztitle="Z (km)",
        bounds=[-max_val, max_val, -max_val, max_val, -max_val, max_val],
    )
    plotter.add_legend(size =(0.2,0.4),loc = "upper right")

    if save_plot:
        plotter.screenshot(f"{Title}.png")

    if show_plot:
        plotter.show()
    else:
        plotter.close()
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
    p = a * (1 - e**2)
    v_perif = np.sqrt(mu / p) * np.array([-np.sin(ta), e + np.cos(ta), 0])
   

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
            ratio = (E0-e*np.sin(E0)-Me)/(1-e*np.cos(E0))
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

def true_anomaly(arr):
    E,e = arr
    return 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))


def orbitsPropagate(file, cb):
    rs =[]
    labels = []
    d2r = 2*np.pi/360
    mu = cb["mu"]
    data = pds.read_csv(file)
    for i in range(len(data)-1):
        print(i)
        line = data.loc[i]
        
        mean_motion = line.MEAN_MOTION
        eccentricity = line.ECCENTRICITY
        inc = line.INCLINATION*d2r
        raan = line.RA_OF_ASC_NODE*d2r
        mean_anomaly = line.MEAN_ANOMALY*d2r
        periapsis = line.ARG_OF_PERICENTER*d2r
        epoch = line.EPOCH 
        a = (mu/((2*np.pi*mean_motion)/86400)**2)**(1/3)
        t_span = (0,1*24*60*60)
        
        E = ecc_anomaly([mean_anomaly,eccentricity],"newton")
        trueAnomaly = true_anomaly([E, eccentricity])
        r = a * (1 - eccentricity * np.cos(E))
        
        c  = [a,eccentricity,inc,trueAnomaly,periapsis,raan]
        op = OP(c,t_span,coes = True)
        op.propagate_orbit()
        rs.append(op.rs)  
        labels.append(getattr(line, 'OBJECT_NAME', f'Sat {id}'))
    plot(rs,labels)