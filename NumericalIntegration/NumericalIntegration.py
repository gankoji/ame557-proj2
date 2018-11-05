import numpy as np
import math
import matplotlib.pyplot as plt

from OrbitalUtilities import elementsToRV, classicalElements
from Disturbances import Oblateness_J2, Drag

debug = False

if debug:
    #simLength = 8706000
    simLength = 25000
    dt = 0.1
else:
    simLength = 160*24*3600
    dt = 0.5

N= int(math.floor(simLength/dt))
simTime = np.linspace(0,simLength, N)
j = 0

# Constants of the simulation
mu = 3.986e14 # m^3/s^2, Earth's gravitational parameter
Re = 6378137 # m, Earth's equatorial radius
J2 = 1.08026267e-3 # Oblateness Coefficient

r = np.zeros((N,3))
v = np.zeros((N,3))
a = np.zeros((N,3))
aJ2 = np.zeros((N,3))
aD = np.zeros((N,3))
orbEl = np.zeros((N,6))
alt = np.zeros((N))
E = np.zeros((N))

r_c = np.zeros((3,))
v_c = np.zeros((3,))
t_c = 0

d2r = math.pi/180.0

# km, NA, deg, deg, deg, deg
oe0 = np.array([6657391.879, 0.002595, d2r*42.748082,
                d2r*345.325883, d2r*124.412552, d2r*287.718564])

n = math.sqrt(mu/((oe0[0])**3))
print("Mean Motion: " + str(n))
T = 2*math.pi/n
print("Orbital Period: " + str(T))

RV = elementsToRV(oe0, mu)

r0 = RV[0:3]
v0 = RV[3:6]

# These are to help with debugging
# r0 = np.array([-5266895.94243465,  1069106.12828683, -3479540.55705787])
# v0 = np.array([  872.34980757, -6684.83866579, -2791.63486922])

print("Initial Radius: " + str(np.linalg.norm(r0)/1000))
print("Initial Velocity: " + str(np.linalg.norm(v0)/1000))
notCrashed = True
stateFlip = False
for i,x in enumerate(simTime):
    if(i <= 1):
        r[i,:] = r0
        v[i,:] = v0
        rMag = np.linalg.norm(r[0,:])
        a[i,:] = -mu/(rMag**3)*r[0,:]
    elif notCrashed:
        a[i,:] = (-mu/(rMag**3)*r[i-1,:]
                  + Oblateness_J2(r[i-1,:],v[i-1,:])
                  + Drag(r[i-1,:],v[i-1,:]))

        aJ2[i,:] = Oblateness_J2(r[i-1,:],v[i-1,:])
        aD[i,:] = Drag(r[i-1,:],v[i-1,:])
        
        # Use Runge-Kutta 4 to integrate v
        R = r[i-1,:]
        V = v[i-1,:]
        rMag = np.linalg.norm(R)
        
        l1 = dt*(-mu/(rMag**3)*R
                 + Oblateness_J2(R,V)
                 + Drag(R,V))
        k1 = dt*V
        
        R = r[i-1,:] + k1/2
        V = v[i-1,:] + l1/2
        rMag = np.linalg.norm(R)
        l2 = dt*(-mu/(rMag**3)*R
                 + Oblateness_J2(R,V)
                 + Drag(R,V))
        k2 = dt*V
        
        R = r[i-1,:] + (k2/(2.0))
        V = v[i-1,:] + (l2/(2.0))
        rMag = np.linalg.norm(R)
        l3 = dt*(-mu/(rMag**3)*R
                 + Oblateness_J2(R,V)
                 + Drag(R,V))
        k3 = dt*V
        
        R = r[i-1,:] + k3
        V = v[i-1,:] + l3
        rMag = np.linalg.norm(R)
        l4 = dt*(-mu/(rMag**3)*R
                 + Oblateness_J2(R,V)
                 + Drag(R,V))
        k4 = dt*V

        v[i,:] = v[i-1,:] + (1/6.0)*(l1 + 2*l2 + 2*l3 + l4)
        r[i,:] = r[i-1,:] + (1/6.0)*(k1 + 2*k2 + 2*k3 + k4)


            
    else:
        r[i,:] = r[i-1,:]
        v[i,:] = v[i-1,:]
        
    orbEl[i,:] = np.array(classicalElements(r[i,:],v[i,:],mu))

    E[i] = (np.linalg.norm(v[i,:])**2)/(2.0) - (mu/np.linalg.norm(r[i,:]))
    
    alt[i] = np.linalg.norm(r[i,:]) - Re
    if (alt[i] < 0):
        notCrashed = False
        if not stateFlip:
            stateFlip = True
            r_c = r[i,:]
            v_c = v[i,:]
            t_c = simTime[i]
            print("Crashed " + str(simTime[i]/86400) + " days after epoch.")
            print("Position: " + str(r[i,:]))
            print("Velocity: " + str(v[i,:]))
            print("Orbital Elements: " + str(orbEl[i,:]))

plotFigures = True
plotFormat = "svg"

if plotFigures:
    fig1, ax1 = plt.subplots()
    ax1.plot(r[:,0], r[:,1])
    ax1.set(xlabel='R_x (m)', ylabel='R_y (m)',
           title='Planar Orbit Path')
    fig1.savefig('Path',format=plotFormat,
                 bbox_inches='tight', pad_inches=0.2, dpi=4000)

    fig2, ax2 = plt.subplots()
    ax2.plot(simTime, E)
    ax2.set(xlabel='Time (s)', ylabel='Energy (m^2/s^2)',
            title='Orbital Specific Energy')
    fig2.savefig('Energy',format=plotFormat,
                 bbox_inches='tight', pad_inches=0.2, dpi=4000)

    fig6, ax6 = plt.subplots()
    ax6.plot(simTime, alt)
    ax6.set(xlabel='Time (s)', ylabel='Altitude (m)',
           title='Height Above Reference Sphere')
    ax6.legend(['Simulation Elements'])
    fig6.savefig('Altitude',format=plotFormat,
                 bbox_inches='tight', pad_inches=0.2, dpi=4000)

    fig5, ax5 = plt.subplots()
    ax5.plot(simTime, aJ2)
    ax5.set(xlabel='Time(s)', ylabel='Acceleration (m/s^2)',
            title='Acceleration Due to Earth\' Oblateness')
    ax5.legend(['X','Y','Z'])
    fig5.savefig('Oblateness',format=plotFormat,
                 bbox_inches='tight', pad_inches=0.2, dpi=4000)

    fig7,ax7 = plt.subplots()
    ax7.plot(simTime, aD)
    ax7.set(xlabel='Time(s)', ylabel='Acceleration (m/s^2)',
            title='Acceleration Due to Atmospheric Drag')
    ax7.legend(['X','Y','Z'])
    fig7.savefig('Drag',format=plotFormat,
                 bbox_inches='tight', pad_inches=0.2, dpi=4000)
    
    fig16, ax16 = plt.subplots()
    ax16.plot(simTime, orbEl[:,5])
    ax16.set(xlabel='Time (s)', ylabel = 'Argument (deg)',
           title='Argument of Periapsis')
    ax16.legend(['Simulation Elements'])
    fig16.savefig('Periapsis',format=plotFormat,
                  bbox_inches='tight', pad_inches=0.2, dpi=4000)
    
    fig17, ax17 = plt.subplots()
    ax17.plot(simTime, orbEl[:,3])
    ax17.set(xlabel='Time (s)', ylabel='Right Ascension (deg)',
           title='Right Ascension of the Ascending Node')
    ax17.legend(['Simulation Elements'])
    fig17.savefig('Omega',format=plotFormat,
                  bbox_inches='tight', pad_inches=0.2, dpi=4000)
    
    fig18, ax18 = plt.subplots()
    ax18.plot(simTime, orbEl[:,0])
    ax18.set(xlabel='Time (s)', ylabel='Distance (m)',
           title='Semi-Major Axis')
    ax18.legend(['Simulation Elements'])
    fig18.savefig('Semimajor',format=plotFormat,
                  bbox_inches='tight', pad_inches=0.2, dpi=4000)

    fig20, ax20 = plt.subplots()
    ax20.plot(simTime, orbEl[:,1])
    ax20.set(xlabel='Time (s)', ylabel='Eccentricity',
           title='Eccentricity')
    ax20.legend(['Simulation Elements'])
    fig20.savefig('Eccentricity',format=plotFormat,
                  bbox_inches='tight', pad_inches=0.2, dpi=4000)
    
    plt.show()
