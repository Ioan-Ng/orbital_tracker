from OrbitPropagator import Orbit_Propagator as OP
import Tools as tools
import PlanetData as pd
import pandas as pds
import numpy as np
file = "SpaceStationData.csv"

cb = pd.earth


mu = cb["mu"]
data = pds.read_csv(file)
d2r = 2*np.pi/360

for i in range(len(data)):
    line = data.loc[i]
    
    mean_motion = line.MEAN_MOTION
    eccentricity = line.ECCENTRICITY
    inc = line.INCLINATION*d2r
    raan = line.RA_OF_ASC_NODE*d2r
    mean_anomaly = line.MEAN_ANOMALY*d2r
    periapsis = line.ARG_OF_PERICENTER*d2r
    epoch = line.EPOCH 
    a = (mu/((2*np.pi*mean_motion)/86400)**2)**(1/3)
  
    
    E = tools.ecc_anomaly(mean_anomaly,eccentricity,ta_check=False)
    trueAnomaly = tools.true_anomaly([E, eccentricity])
    r = a * (1 - eccentricity * np.cos(E))
    
    c  = [a,eccentricity,inc,trueAnomaly,periapsis,raan]
    orbit = OP(c,t_span = (0,90*60),keplearian_data = False,central_body = cb, degrees=True,ta_check = True)
    orbit.propagate_orbit()
    rs =orbit.rs


tools.animate_plot(rs, t_span=(0,90*60), central_body = pd.earth)