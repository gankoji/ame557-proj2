import numpy as np
import math
import matplotlib.pyplot as plt

from OrbitalUtilities import *

debug = False

if debug:
    simLength = 3
    dt = 0.0001
else:
    simLength = 25
    dt = 0.0001

N= int(math.floor(simLength/dt))
simTime = np.linspace(0,simLength, N)
j = 0

# Constants of the simulation
mu = 1

r = np.zeros((N,3,5))
v = np.zeros((N,3,5))
a = np.zeros((N,3,5))
orbEl = np.zeros((N,6,5))
deg2rad = math.pi/(180.0)

# These next few equations are directly from BMW 2.5
# Define the initial orbit in a, e, i, Omega, omega, nu
oe0 = np.array([0.5, 0.1, 20*deg2rad, 10*deg2rad, 10*deg2rad, 0*deg2rad])

# Magnitude of position
p = oe0[0]*(1 - (oe0[1]**2))
rMag = p/(1 + (oe0[1]*math.cos(oe0[5])))

# Position Vector in perifocal (z assumed zero by definition)
rpf = np.zeros((3,))
rpf[0] = rMag*math.cos(oe0[5])
rpf[1] = rMag*math.sin(oe0[5])

# Velocity Vector in perifocal
vpf = np.zeros((3,))
vpf[0] = math.sqrt(mu/p)*(-math.sin(oe0[5]))
vpf[1] = math.sqrt(mu/p)*(oe0[1] + math.cos(oe0[5]))

# Building the DCM from perifocal to equatorial (which is where the
# sim is)
C = np.zeros((3,3))

ci = math.cos(oe0[2])
si = math.sin(oe0[2])
cO = math.cos(oe0[3])
sO = math.sin(oe0[3])
co = math.cos(oe0[4])
so = math.sin(oe0[4])

C[0,0] = cO*co - sO*so*ci
C[0,1] = -cO*so - sO*co*ci
C[0,2] = sO*si

C[1,0] = sO*co + cO*so*ci
C[1,1] = -sO*so + cO*co*ci
C[1,2] = -cO*si

C[2,0] = so*si
C[2,1] = co*si
C[2,2] = ci

# Finally, transform initial velocity/position to equatorial
# coordinates
r0 = np.dot(C, rpf)
v0 = np.dot(C, vpf)

for a_in in [0.000, 0.01]:
    for i,x in enumerate(simTime):
        if(i <= 1):
            r[i,:,j] = r0
            v[i,:,j] = v0
            ad = np.array([0,0,a_in])
            rMag = np.linalg.norm(r[0,:,j])
            a[i,:,j] = -mu/(rMag**3)*r[0,:,j] + ad
        else:
            a[i,:,j] = -mu/(rMag**3)*r[i-1,:,j] + ad

            # Use Runge-Kutta 4 to integrate v
            R = r[i-1,:,j]
            rMag = np.linalg.norm(R)
            l1 = dt*(-mu/(rMag**3)*R + ad)
            k1 = dt*(v[i-1,:,j])
            R = r[i-1,:,j] + k1/2
            rMag = np.linalg.norm(R)
            l2 = dt*(-mu/(rMag**3)*R + ad)
            k2 = dt*(v[i-1,:,j] + l1/2)
            R = r[i-1,:,j] + k2/2
            rMag = np.linalg.norm(R)
            l3 = dt*(-mu/(rMag**3)*R + ad)
            k3 = dt*(v[i-1,:,j] + l2/2)
            R = r[i-1,:,j] + k3
            rMag = np.linalg.norm(R)
            l4 = dt*(-mu/(rMag**3)*R + ad)
            k4 = dt*(v[i-1,:,j] + l3)

            v[i,:,j] = v[i-1,:,j] + (1/6.0)*(l1 + 2*l2 + 2*l3 + l4)
            r[i,:,j] = r[i-1,:,j] + (1/6.0)*(k1 + 2*k2 + 2*k3 + k4)


        #print(classicalElements(r[i,:,j],v[i,:,j],mu))
        orbEl[i,:,j] = np.array(classicalElements(r[i,:,j],v[i,:,j],mu))

    j += 1

# Next, we calculate some averaged quantities for comparison with the
# simulation results
n = math.sqrt(mu/(oe0[0]**3)) # mean motion, rad/s
P = (2*math.pi)/n # period, s
a_D = 0.01
deltaA = 0.0
a_avg = oe0[0] + deltaA

e_dot_bar = ((5.0/4.0)*math.pi*math.cos(oe0[5])
             *math.sqrt(1 - (oe0[1]**2))
             *math.sin(oe0[2])*a_D
             /(n*oe0[0]))

deltaE = -(a_D*math.sin(oe0[2])/((n**2)*oe0[0])
          *(math.pi*math.sqrt(1 - (oe0[1]**2))
            *math.cos(oe0[5])
            - (1 - (oe0[1]**2))/4*math.sin(oe0[5]))
          - e_dot_bar*math.pi/n)

e0 = oe0[1]
ef = e0 + deltaE

e_averaged = np.zeros((simTime.shape[0]))
for i, x in enumerate(e_averaged):
    e_averaged[i] = e0 + simTime[i]*(deltaE/P)
            


plotFigures = True

if plotFigures:
    # fig1, ax1 = plt.subplots()
    # ax1.plot(r[:,0,0], r[:,1,0])
    # ax1.set(xlabel='R_x (m)', ylabel='R_y (m)',
    #        title='Planar Orbit Path - Undisturbed Orbit')
    # fig1.savefig('Path_Undisturbed.png')

    # fig2, ax2 = plt.subplots()
    # ax2.plot(r[:,0,1], r[:,1,1])
    # ax2.set(xlabel='R_x (m)', ylabel='R_y (m)',
    #        title='Planar Orbit Path - a_d = 0.011')
    # fig2.savefig('Path_a1.png')

    # fig6, ax6 = plt.subplots()
    # ax6.plot(simTime,orbEl[:,0:4,0])
    # ax6.set(xlabel='Time (s)', ylabel='Orbital Elements',
    #        title='Classical Orbital Elements - Undisturbed Orbit')
    # ax6.legend(['a','e','i','Omega'])
    # fig6.savefig('Elements_Unidisturbed.png')

    # fig7, ax7 = plt.subplots()
    # ax7.plot(simTime,orbEl[:,0:4,1])
    # ax7.set(xlabel='Time (s)', ylabel='Orbital Elements',
    #        title='Classical Orbital Elements - a_d = 0.01')
    # ax7.legend(['a','e','i','Omega'])
    # fig7.savefig('Elements_a1.png')

    # fig16, ax16 = plt.subplots()
    # ax16.plot(simTime,orbEl[:,5,0])
    # ax16.set(xlabel='Time (s)', ylabel='Orbital Elements',
    #        title='Mean Anomaly - Undisturbed')
    # ax16.legend(['Mean Anomaly'])
    # fig16.savefig('MeanAnomaly_Undisturbed.png')

    # Semi-Major Axis Plots
    # fig17, ax17 = plt.subplots()
    # ax17.plot(simTime, orbEl[:,0,0])
    # ax17.set(xlabel='Time (s)', ylabel='Distance (m)',
    #        title='Semi-Major Axis - Undisturbed')
    # fig17.savefig('semimajor_undisturbed.png')
    
    fig18, ax18 = plt.subplots()
    ax18.plot(simTime, orbEl[:,0,1])
    ax18.plot([simTime[0], simTime[-1]], [a_avg, a_avg])
    ax18.set(xlabel='Time (s)', ylabel='Distance (m)',
           title='Semi-Major Axis - a_d = 0.01')
    ax18.legend(['Simulation Elements','Analytical Average'])
    fig18.savefig('semimajor_a1.png')
    
    # Eccentricity Plots
    # fig19, ax19 = plt.subplots()
    # ax19.plot(simTime, orbEl[:,1,0])
    # ax19.set(xlabel='Time (s)', ylabel='Eccentricity',
    #        title='Eccentricity - Undisturbed')
    # fig19.savefig('Eccentricity_undisturbed.png')    

    fig20, ax20 = plt.subplots()
    ax20.plot(simTime, orbEl[:,1,1])
    ax20.plot(simTime, e_averaged)
    ax20.set(xlabel='Time (s)', ylabel='Eccentricity',
           title='Eccentricity - a_d = 0.01')
    ax20.legend(['Simulation Elements','Analytical Average'])
    fig20.savefig('Eccentricity_a1.png')
    
    plt.show()
