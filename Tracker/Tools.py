import numpy as np
import math as m

#main fucntion to convert our state elemetns into vecotrs
def kpToVector(intial_state,mu, degrees, ta_check):
    #ta_check says true if we are using true anomaly and will be false if we are using mean anomaly
    if degrees:
        a,e,i,ta,aop,raan = intial_state #major axis, eccentricity, inclincation, true anaomoly, argument of periapsis , right ascension of the acsending node,
        i = np.radians(i)
        ta = np.radians(ta)
        aop = np.radians(aop)
        raan = np.radians(raan)
    else:
        a,e,i,ta,aop,raan = intial_state

    E = ecc_anomaly(ta,e,ta_check)
    r_norm = a*(1-e**2)/(1+e*np.cos(ta))

    r_perif = r_norm*np.array([m.cos(ta),m.sin(ta),0])
    p = a * (1 - e**2)
    v_perif = np.sqrt(mu / p) * np.array([-np.sin(ta), e + np.cos(ta), 0])
   
    #still not to sure why this works, i need to look into it
    eci_t0_perif = np.transpose(eci_t0_perif(raan,aop,i))
    r = np.dot(eci_t0_perif,r_perif)
    v = np.dot(eci_t0_perif,v_perif)

    return r,v
def ecc_anomaly(ta,e, ta_check):
    tol = 1e-8
    #if we are using mean anomlay then we have to use newtwon rapson mtehod to solve our ecc anom
    if not ta_check:
        mean_anomaly = ta
        #e0 and e1 are the things we use to iterate for newton rapson
        if mean_anomaly < np.pi/2:
            E0 = mean_anomaly + e/2
        else: 
            E0 = mean_anomaly - e
        for n in range(400):
            ratio = (E0 - e*np.sin(E0) - mean_anomaly)/(1 - e*np.cos(E0))
            #if we are close enough within range then we are happy, our tolerance is 1e-8 but we can chnang this if we want to ve more precise
            if abs(ratio) < tol:
                if n == 0: 
                    return E0
                else:
                    return E1
            else:
                E1 = E0-ratio
                E0 = E1
        return False
    elif ta_check :
        #we we do have the tru anomaly then we can just calculate directly using formula 
        return 2*m.atan(np.sqrt((1-e)/1+e))*m.tan(ta/2)
    else:
        #if something goofy is parsed
        print("invalid method")

def eci_t0_perif(raan, aop, i):
    #don't know too well what is going on here but 
    row0 = [-m.sin(raan)*m.cos(i)*m.sin(aop) + m.cos(raan)*m.cos(aop),m.cos(raan)*m.cos(i)*m.sin(aop)+m.sin(raan)*m.cos(aop),m.sin(i)*m.sin(aop)]
    row1 = [-m.sin(raan)*m.cos(i)*m.cos(aop) - m.cos(raan)*m.sin(aop),m.cos(raan)*m.cos(i)*m.cos(aop)-m.sin(raan)*m.sin(aop),m.sin(i)*m.cos(aop)]
    row2 = [m.sin(raan)*m.sin(i), -m.cos(raan)*m.sin(i),m.cos(i)]
    return np.array([row0,row1,row2])

def true_anomaly(arr):
    E,e = arr
    return 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
