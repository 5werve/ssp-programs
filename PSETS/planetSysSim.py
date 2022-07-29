import numpy as np  
import matplotlib.pyplot as plt

def GetMagnitude(v):
    return np.sqrt(v[0]**2 + v[1]**2)

def PlanetarySystemSim(M, r1, r2, v1, v2):
    rSun = [0, 0]

    delT = 0.01
    G = 6.674 * 10 ** (-11)
    
    time = 0

    plt.figure(figsize = (6, 6))
    cntr = 0

    R1 = []
    R2 = []

    maxX1 = maxX2 = minX1 = minX2 = T1 = T2 = 0
    
    while True:
        a1 = r1 * -G * M / (GetMagnitude(r1) ** 3)
        a2 = r2 * -G * M / (GetMagnitude(r2) ** 3)

        a1Temp = a1
        a2Temp = a2

        r1 = r1 + v1 * delT + 0.5 * a1 * delT**2
        r2 = r2 + v2 * delT + 0.5 * a2 * delT**2

        a1 = r1 * -G * M / (GetMagnitude(r1) ** 3)
        a2 = r2 * -G * M / (GetMagnitude(r2) ** 3)
       

        v1 = v1 + 0.5 * (a1 + a1Temp) * delT
        v2 = v2 + 0.5 * (a2 + a2Temp) * delT

        if r1[0] > maxX1:
            maxX1 = r1[0]
        if r2[0] > maxX2:
            maxX2 = r2[0]
        if r1[0] < minX1:
            minX1 = r1[0]
        if r2[0] < minX2:
            minX2 = r2[0]
        
        if abs(r1[1]) < 0.05 and r1[0] > 0 and time > 2 and T1 == 0:
            T1 = time
            print('set1')

        if abs(r2[1]) < 0.05 and r2[0] > 0 and time > 2 and T2 == 0:
            T2 = time
            print('set2')

        if T1 != 0 and T2 != 0:
            break

        print(time)
        time += delT
        
        R1.append(r1)
        R2.append(r2)
        
        plt.cla()

        plt.plot(0, 0, 'go')
        
        plt.plot(r1[0], r1[1], 'ro')
        plt.plot(r2[0], r2[1], 'bo')

        plt.plot([r[0] for r in R1], [r[1] for r in R1])
        plt.plot([r[0] for r in R2], [r[1] for r in R2])
        
        plt.xlabel('x-position')
        plt.ylabel('y-position')

        plt.pause(delT)

    semiMajA1 = (maxX1 - minX1) / 2
    semiMajA2 = (maxX2 - minX2) / 2

    print('a1', semiMajA1, 'a2', semiMajA2)
    print('T1', T1, 'T2', T2)
    print('a1^3 / T1^2:', semiMajA1**3/T1**2, 'a2^3 / T2^2:', semiMajA2**3/T2**2)


PlanetarySystemSim(10**13, np.array([15, 0]), np.array([10, 0]), np.array([0, 5]), np.array([0, 5]))
