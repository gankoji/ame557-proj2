import numpy as np
import math

# Functions for classical orbital elements
def classicalElements(r,v,mu):
    # Angular Momentum
    H = np.cross(r,v)
    HMag = np.linalg.norm(H)
    h = H/HMag

    # Energy
    zeta = np.dot(v,v)/2 - mu/np.linalg.norm(r)

    # Semimajor Axis
    a = -mu/(2*zeta)

    # Eccentricity
    a1 = ((np.linalg.norm(v)**2)/mu - 1/np.linalg.norm(r))
    b1 = np.dot(r,v)/mu
    E = a1*r - b1*v
    e = np.linalg.norm(E)

    # Inclination
    I = np.array([1,0,0])
    J = np.array([0,1,0])
    K = np.array([0,0,1])
    inc = (180.0/math.pi)*math.acos(np.dot(K,h))

    # Here we branch for equatorial orbits
    if (inc < 1e-6):
        # Orbit is equatorial
        # Next, check for circular orbits
        if (e < 1e-3):
            # Circular Orbit
            # Longitude of the ascending node and argument of perigee are
            # undefined for circular equatorial orbits
            Omega = 0.0
            omega = 0.0

            # Replace true anomaly with true longitude at epoch, l_0
            l_o = 180.0/math.pi*math.acos(np.dot(I,r)/np.linalg.norm(r))
            # Since the orbit is circular equatorial, we can simply
            # use the angle formed by r with the coordinate system to
            # define the Mean Anomaly
            M = 180.0/math.pi*math.atan2(r[1], r[0]) 
        else:
            # Noncircular Orbit

            # Longitude of the ascending node is
            # undefined for equatorial orbits
            Omega = 0.0
            # Find longitude of perigee (omega_bar)
            omega_bar = 180.0/math.pi*math.acos(np.dot(I,E)/e)
            omega = omega_bar

            # True Anomaly
            f = math.acos(np.dot(E,r)
                          /(e*np.linalg.norm(r)))
            M = 180.0/math.pi*trueAnomalyToMeanAnomaly(f, e)
    else:
        # Orbit is not equatorial
        # Next, check for circular orbits
        if (e < 1e-3):
            # Circular Orbit
            # Line of Nodes
            N = np.cross(K,h)
            n = N/np.linalg.norm(N)

            # Longitude of the ascending node
            Omega = 180.0/math.pi*math.acos(np.dot(I,n))
            # Arugment of perigee is undefined for circular orbits
            omega = 0.0
            # Find arugument of latitude (u) in place of true anomaly
            u = math.acos(np.dot(n,r)/np.linalg.norm(r))
            M = 180.0/math.pi*trueAnomalyToMeanAnomaly(u,e)

        else:
            # Noncircular Orbit
            # Line of Nodes
            N = np.cross(K,h)
            n = N/np.linalg.norm(N)
            #print("Noncircular, nonequatorial")
            # Longitude of the ascending node
            Omega = 180.0/math.pi*math.acos(np.dot(I,n))

            # Argument of perigee
            omega = 180.0/math.pi*math.acos(np.dot(n,E)/e)

            # True Anomaly
            f = math.acos(np.dot(E,r)
                          /(e*np.linalg.norm(r)))
            M = 180.0/math.pi*trueAnomalyToMeanAnomaly(f, e)

    output = [a,e,inc,Omega,omega,M]
    return output

# Conversion functions for mean and true anomaly
def trueAnomalyToMeanAnomaly(f,e):
    # Find initial eccentric anomaly from true anomaly
    inner = math.fabs((1-e)/(1+e))
    
    E = 2*math.atan(inner*math.tan(f/2))
    # Convert initial eccentric anomaly to mean anomaly
    M = E - e*math.sin(E)
    return M

def meanAnomalytoTrueAnomaly(M_f, n):
    # Convert the mean anomaly to eccentric anomaly
    E_new = M_f + e*math.sin(M_f) + (e**2)/2*math.sin(2*M_f)
    err = E_new - e*math.sin(E_new) - M_f
    # Loop to iterate on the transcendtal Kepler's equation until
    # convergence
    while(math.fabs(err)  > 1e-3):
        err = E_new - e*math.sin(E_new) - M_f
        E_last = E_new
        E_deriv = 1 - e*math.cos(E_new)
        E_new = E_last - (err/E_deriv)

    # Converting Eccentric Anomaly to true anomaly
    f_f_rad = 2*math.atan(math.sqrt((1+e)/(1-e))*math.tan(E_new/2))
    f_f = 180.0/math.pi*f_f_rad

    return f_f
