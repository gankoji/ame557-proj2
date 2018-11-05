import numpy as np
import math

from US76 import atmosphere

mu = 3.986e14 # m^3/s^2, Earth's gravitational parameter
Re = 6378137.0 # m, Earth's equatorial radius
J2 = 1.08026267e-3 # Oblateness Coefficient
BC = 1/(12.741621)/(.13071e-3) # Ballistic Coefficient on Tiangong-1
Bstar = 0.00013071
rho_ref = 0.157

__underHeight__ = False
def Oblateness_J2(r,v):
    rMag = np.linalg.norm(r)
    z = r[2]
    

    z_hat = np.array([0,0,1])
    r_hat = r/rMag

    coeff = -3*mu*J2*(Re**2)/(2*(rMag**4))
    #print("Coeff: " + str(coeff))
    rcoeff = (1 - (5*z**2)/(rMag**2))
    #print("Rcoeff: " + str(rcoeff))
    zcoeff = (2*z)/rMag
    #print("Zcoeff: " + str(zcoeff))

    a_d = coeff*(rcoeff*r_hat + zcoeff*z_hat)

    return a_d
    #return np.zeros((3,))

def Drag(r,v):
    global __underHeight__

    rMag = np.linalg.norm(r)
    vMag = np.linalg.norm(v)

    alt = (rMag - Re)/1000.0
    stuff = atmosphere((rMag - Re)/1000.0)
    rho = stuff[0]
    a_d = -(Bstar/rho_ref)*rho*vMag*v
    
    # if(alt < 50) or __underHeight__:
    #     __underHeight__ = True
    #     print("Height: " + str(alt) + " Density: " + str(rho) + " Acceleration: " + str(np.linalg.norm(a_d)))
    #     print(str(1/BC))
    #     print(vMag)
    #     print(r)
    #     print(v)
    #     print(a_d)

    # return np.zeros((3,))
    return a_d
