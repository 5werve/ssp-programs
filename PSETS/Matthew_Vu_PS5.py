import numpy as np
from math import *
import matplotlib.pyplot as plt

def GetMagnitude(v):
    return np.sqrt(v[0]**2 + v[1]**2)

def ThreeBodSim(r1, r2, r3, v1, v2, v3, m1, m2, m3):
    delT = 0.001
    animate = 0
    G = 6.674 * 10 ** (-11)

    R1 = []
    R2 = []
    R3 = []

    time = 0

    plt.figure(figsize = (6, 6))
    cntr = 0

    while True:
        a1 = G * ((m2 * (r2 - r1) / GetMagnitude(r2 - r1) ** 3) + (m3 * (r3 - r1) / GetMagnitude(r3 - r1) ** 3))
        a2 = G * ((m1 * (r1 - r2) / GetMagnitude(r1 - r2) ** 3) + (m3 * (r3 - r2) / GetMagnitude(r3 - r2) ** 3))
        a3 = G * ((m1 * (r1 - r3) / GetMagnitude(r1 - r3) ** 3) + (m2 * (r2 - r3) / GetMagnitude(r2 - r3) ** 3))

        r1 = r1 + v1 * delT + 0.5 * a1 * delT ** 2
        r2 = r2 + v2 * delT + 0.5 * a2 * delT ** 2
        r3 = r3 + v3 * delT + 0.5 * a3 * delT ** 2

        na1 = G * ((m2 * (r2 - r1) / GetMagnitude(r2 - r1) ** 3) + (m3 * (r3 - r1) / GetMagnitude(r3 - r1) ** 3))
        na2 = G * ((m1 * (r1 - r2) / GetMagnitude(r1 - r2) ** 3) + (m3 * (r3 - r2) / GetMagnitude(r3 - r2) ** 3))
        na3 = G * ((m1 * (r1 - r3) / GetMagnitude(r1 - r3) ** 3) + (m2 * (r2 - r3) / GetMagnitude(r2 - r3) ** 3))

        v1 = v1 + 0.5 * (a1 + na1) * delT
        v2 = v2 + 0.5 * (a2 + na2) * delT
        v3 = v3 + 0.5 * (a3 + na3) * delT
        
        animate += 1
        
        if (animate == 10):
            
            animate = 0
            R1.append(r1)
            R2.append(r2)
            R3.append(r3)
            
            if (len(R1) > 300):
                R1.pop(0)
            if (len(R2) > 300):
                R2.pop(0)
            if (len(R3) > 300):
                R3.pop(0)
            
            plt.cla()

            plt.plot(r1[0], r1[1], 'ro')
            plt.plot(r2[0], r2[1], 'bo')
            plt.plot(r3[0], r3[1], 'go')

            plt.plot([r[0] for r in R1], [r[1] for r in R1])
            plt.plot([r[0] for r in R2], [r[1] for r in R2])
            plt.plot([r[0] for r in R3], [r[1] for r in R3])
            
            plt.xlabel('x-position')
            plt.ylabel('y-position')
            
            plt.xlim(-6, 6)
            plt.ylim(-6, 6)
            
            plt.pause(delT)

            k1 = 0.5 * m1 * GetMagnitude(v1) ** 2
            k2 = 0.5 * m2 * GetMagnitude(v2) ** 2
            k3 = 0.5 * m3 * GetMagnitude(v3) ** 2

            u1 = (-G * m1 * m2) / GetMagnitude(r2 - r1) + (-G * m1 * m3) / GetMagnitude(r3 - r1)
            u2 = (-G * m2 * m1) / GetMagnitude(r1 - r2) + (-G * m2 * m3) / GetMagnitude(r3 - r2)
            u3 = (-G * m3 * m1) / GetMagnitude(r1 - r3) + (-G * m3 * m2) / GetMagnitude(r2 - r3)

            print('Total energy of system: ', k1 + k2 + k3 + u1 + u2 + u3)
            
        time += delT

ThreeBodSim(np.array([0, 0]), np.array([5, 0]), np.array([5.2, 0]), np.array([0, 0]), np.array([0, 10]), np.array([0, 15.7]), 10**13, 10**11, 10**9) 
