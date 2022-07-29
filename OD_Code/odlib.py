from math import *
import numpy as np

k = 0.0172020989484
cAU = 173.144643267
eps = radians(23.4374)

# Converts hours:min:sec to degrees or radians
def HMStoDegOrRad(h, m, s, option):
    deg = (h * 15) + (m * (1/60) * 15) + (s * (1/3600) * 15)
    if (option):
        return deg
    else:
        return deg * (math.pi/180)

# Converts degrees min' sec" to degrees
def DMStoDeg(d, am, asec):
    return ((float(d)) + (copysign(float(am), float(d)) * (1/60)) + (copysign(float(asec), float(d)) * (1/3600)))

# Finds angle in the correct quadrant (range is positive) given sinx and cosx
def FindAngle(s, c):
    return (atan2(s, c) + 2*np.pi) % (2*np.pi)

# Converts RA in decimal to hours:min:sec
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

# Converts declination in decimal to degrees min' sec"
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

# Gets magnitude of 3D vector
def GetMagnitude(v):
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)

# Converts altitude (deg) and azimuth (deg) of object to RA (hours) and declination (degrees)
def AltAzToRADec(alt, az, lst, lat):
    lst = HMStoDegOrRad(lst[0], lst[1], lst[2], False)
    lat = lat * np.pi/180
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

# Returns 2 3D vectors dotted
def DotVectors(v1, v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
 
# Returns 2 3D vectors crossed  
def CrossVectors(v1, v2):
    return [v1[1]*v2[2] - v1[2]*v2[1], -(v1[0]*v2[2] - v2[0]*v1[2]), v1[0]*v2[1] - v2[0]*v1[1]]

# Returns triple product of v1, v2, v3
def TripleProduct(v1, v2, v3):
    return DotVectors(CrossVectors(v2, v3), v1)

# Returns coordinates of a 3D vector rotated by a rads about the z-axis and b rads about the x-axis
def Rot3DVector(v, a, b):
    rotMatrix = np.array([[np.cos(a), np.sin(a), 0], [-np.sin(a)*np.cos(b), np.cos(b)*np.cos(a), np.sin(b)], [np.sin(a)*np.sin(b), -np.cos(a)*np.sin(b), np.cos(b)]])
    result = np.matmul(rotMatrix, v)
    return result

# Returns the values of a csv file
def InputCSV(path):
    data = np.loadtxt(path, dtype=float, delimiter=',')
    return np.array(data)

# Calulates h of object given its r vector and r. vector
def CalculateAngularMomentum(path):
    ad = InputCSV(path)
    v1 = ad[0:3]
    v2 = ad[3:6]
    conversionRate = 58.13244
    return [round(conversionRate*(v1[1]*v2[2] - v1[2]*v2[1]), 6), round(-conversionRate*(v1[0]*v2[2] - v2[0]*v1[2]), 6), round(conversionRate*(v1[0]*v2[1] - v2[0]*v1[1]), 6)]
    
CalculateAngularMomentum('/home/vumat/Desktop/Python_Code/ang_mom_info.csv')