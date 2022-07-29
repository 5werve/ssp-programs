import numpy as np
import matplotlib.pyplot as plt
from math import *

def DrawTrajectoryWithoutAir(v0, angle):
    delT = 0.09
    g = 9.8
    angle = angle*np.pi/180
    
    time = 2 * v0 * np.sin(angle) / g
    times = np.arange(0, time, delT)

    v0x = v0 * np.cos(angle)
    v0y = v0 * np.sin(angle)

    plt.figure(figsize = (6, 6))
    cntr = 0

    for t in times:
        x = v0x * t
        y = v0y * t + (0.5 * -9.8 * t**2)

        plt.plot(x, y, 'ro')

        plt.ylim(0, v0y ** 2 / (2 * g) + 5)
        plt.xlim(0, time * v0x + 5)

        plt.xlabel('x-position')
        plt.ylabel('y-position')

        plt.pause(delT)

    print('time of flight:', time, 'range:', time * v0x)

#DrawTrajectoryWithoutAir(70, 45)

def DrawTrajectoryWithAir(v0, angle):
    delT = 0.05
    g = 9.8
    angle = angle*np.pi/180
    time = 0
    y = x = i = j = 0
    b = 0.005679
    m = 0.4
    timeVac = 0

    v0x = v0 * np.cos(angle)
    v0y = v0 * np.sin(angle)
    vx = v0x
    vy = v0y

    plt.figure(figsize = (6, 6))
    cntr = 0

    while y >= 0 or j >= 0:
        plt.ylim(0, v0y ** 2 / (2 * g) + 5)
        plt.xlim(0, (2 * v0 * np.sin(angle) / g) * v0x + 5)

        if y >= 0:
            ax = (-b * np.sqrt(vx**2 + vy**2) * vx / m)
            ay = (-b * np.sqrt(vx**2 + vy**2) * vy / m) - g

            x = x + vx * delT + 0.5 * ax * delT**2
            y = y + vy * delT + 0.5 * ax * delT**2

            vx = vx + ax * delT
            vy = vy + ay * delT

            plt.plot(x, y, 'ro')

        if j >= 0:
            i = v0x * timeVac
            j = v0y * timeVac + (0.5 * -9.8 * timeVac**2)       
            timeVac += delT
            
            plt.plot(i, j, 'bo')  

        plt.xlabel('x-position')
        plt.ylabel('y-position')

        plt.pause(delT)


    print('Range with air res: ', x)
    print('Range in vacuum: ', i)

    
#DrawTrajectoryWithAir(55, 45)

fruits = np.array([["Apple","Banana","Blueberry","Cherry"], 
["Coconut","Grapefruit","Kumquat","Mango"], 
["Nectarine","Orange","Tangerine","Pomegranate"], 
["Lemon","Raspberry","Strawberry","Tomato"]])

#Extract bottom right element in one command
print(fruits[3, 3])

#Extract the inner 2x2 square in one command
print(fruits[1:3, 1:3])

#Extract the first and third rows in one command
print(fruits[0:3:2, :])

#Extract the inner 2x2 square flipped vertically and horizontally in one command
print(fruits[-2:0:-1, -2:0:-1])

#Swap the first and fourth columns in a few commands
swappedFruits = np.copy(fruits)
fruits[:, 0], fruits[:, 3] = np.copy(fruits)[:, 3], np.copy(fruits)[:, 0]
print(fruits)

#Replace every element in the array above with the string "SLICED!" in one command
fruits[:, :] = 'SLICED!'
print(fruits)

def HMStoDegOrRad(h, m, s, option):
    deg = (h * 15) + (m * (1/60) * 15) + (s * (1/3600) * 15)
    if (option):
        return deg
    else:
        return deg * (math.pi/180)

def FindAngle(s, c):
    return (atan2(s, c) + 2*np.pi) % (2*np.pi)

def RAdecimalToHMS(val):
    h = trunc(val / 15)
    val = (abs(val) / 15) - abs(h)
    m = trunc(val / (1/60))
    val = val - (m * (1/60))
    s = val / (1/3600)
    s = round(s, 0)
    if (s == 60):
        s = 0
        m += 1
    print ("RA:", h, ":", m, ":", s)

def DECdecimalToDMS(dec):
    d = trunc(dec)
    dec = abs(dec) - abs(d)
    m = trunc(dec / (1/60))
    dec -= (m * (1/60))
    s = dec / (1/3600)
    s = round(s, 0)
    if (s == 60):
        s = 0
        m += 1
    print ("Dec:", d, ":", m, ":", s)

def AltAzToRADec(alt, az, lst):
    lst = HMStoDegOrRad(lst[0], lst[1], lst[2], False)
    lat = 40 * np.pi/180
    az = az * np.pi/180
    alt = alt * np.pi/180
    
    dec = np.arcsin(np.sin(lat) * np.sin(alt) + np.cos(lat) * np.cos(alt) * np.cos(az))

    HA = FindAngle((np.sin(2*np.pi - az) * np.cos(alt) / np.cos(dec)), (np.sin(alt) - np.sin(lat) * np.sin(dec)) / (np.cos(lat) * np.cos(dec)))
    RA = lst - HA
    
    RA = RA * 180/np.pi / 15
    if RA < 0:
        RA += 24
    
    dec = dec * 180/np.pi
    
    return RA, dec

RA, dec = AltAzToRADec(233.987, 26.473, 22, 28, 59.1796)
RAdecimalToHMS(RA)
DECdecimalToDMS(dec)
        
        
        
    
