from math import *
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from matplotlib.ticker import FormatStrFormatter


k = 0.01720209895
cAU = 173.144643267
eps = radians(23.4374)
G = 6.674 * (10 ** -11)
Ep = radians(23.4374)
path = '/home/vumat/Desktop/Python_Code/OD_Input/VuInput.csv'

def HMStoDegOrRad(h, m, s, option):
    deg = (h * 15) + (m * (1/60) * 15) + (s * (1/3600) * 15)
    if (option):
        return deg
    else:
        return deg * (np.pi/180)
        
def DMStoDeg(d, am, asec):
    return ((float(d)) + (copysign(float(am), float(d)) * (1/60)) + (copysign(float(asec), float(d)) * (1/3600)))

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
    return (h, m, s)

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
    return (d, m, s)
    
def GetMagnitude(v):
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)

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
    
def DotVectors(v1, v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
    
def CrossVectors(v1, v2):
    return [v1[1]*v2[2] - v1[2]*v2[1], -(v1[0]*v2[2] - v2[0]*v1[2]), v1[0]*v2[1] - v2[0]*v1[1]]
    
def TripleProduct(v1, v2, v3):
    return DotVectors(CrossVectors(v2, v3), v1)

def NewtonMethod(func, deriv, iguess, error):
    x = iguess
    while 1:
        temp = x
        x = x - (func(x) / deriv(x))
        if abs(temp - x) < error:
            return x

def InputCSV(path):
    data = np.loadtxt(path, dtype=float, delimiter=',')
    return np.array(data)

def CalculateAngularMomentum(v1, v2):
    return [round((v1[1]*v2[2] - v1[2]*v2[1]), 6), round(-(v1[0]*v2[2] - v2[0]*v1[2]), 6), round((v1[0]*v2[1] - v2[0]*v1[1]), 6)]

def PrintOrbitalElements(value, expected, error):
    print('Calculated Value:', value)
    print('Expected Value:', expected)
    print('Error(%):', error)

def CalculateOrbitalElements(path, r, v, t):
    ad = InputCSV(path)
    h = CalculateAngularMomentum(r, v)
    expected = ad[3]
    
    # Calculate orbital elements
    a = 1 / ((2/GetMagnitude(r)) - (GetMagnitude(v) ** 2))
    e = np.sqrt(1 - ((GetMagnitude(h)) ** 2 / a))
    i = np.arccos(h[2] / GetMagnitude(h))
    omega = FindAngle((h[0] / (GetMagnitude(h) * np.sin(i))), (-h[1] / (GetMagnitude(h) * np.sin(i))))
    #print(omega)
    #print((h[0] / (GetMagnitude(h) * np.sin(i))), (-h[1] / (GetMagnitude(h) * np.sin(i))))
    U = FindAngle((r[2] / (GetMagnitude(r) * np.sin(i))), (r[0] * np.cos(omega) + r[1] * np.sin(omega)) / GetMagnitude(r))
    nu = FindAngle((a * (1 - e ** 2) / (GetMagnitude(h) * e) * DotVectors(r, v) / GetMagnitude(r)), (1 / e) * ((a * (1 - e ** 2)) / GetMagnitude(r) - 1))
    w = U - nu
    if w < 0:
        w += 2 * np.pi

    if (1 - e  ** 2) < 0:
        return 0, 0, 0, 0, 0, 0, 0
    E = FindAngle(GetMagnitude(r) * np.sin(nu) / (a * np.sqrt(1 - e ** 2)), (e + np.cos(nu)) / (1 + e * np.cos(nu)))
    M = E - e * np.sin(E)
    T = t - M / (k / a ** 1.5)
    #print(t)
    i = degrees(i)
    omega = degrees(omega)
    w = degrees(w)
    E = degrees(E)
    M = degrees(M)
    
    # Calculate percent errors
    da = abs((a - expected[9]) / expected[9]) * 100
    de = abs((e - expected[0]) / expected[0]) * 100
    di = abs((i - expected[2]) / expected[2]) * 100
    domega = abs((omega - expected[3]) / expected[3]) * 100
    dw = abs((w - expected[4]) / expected[4]) * 100
    dM = abs((M - expected[7]) / expected[7]) * 100
    dT = abs((T - expected[5]) / expected[5]) * 100
    
    print('Semi-major axis(AU):\n')
    PrintOrbitalElements(a, expected[9], da)
    print('\nEccentricity:\n')
    PrintOrbitalElements(e, expected[0], de)
    print('\nInclination(degrees):\n')
    PrintOrbitalElements(i, expected[2], di)
    print('\nOmega(degrees):\n')
    PrintOrbitalElements(omega, expected[3], domega)
    print('\nAngle of perihelion(degrees):\n')
    PrintOrbitalElements(w, expected[4], dw)
    print('\nMean anomaly(degrees):\n')
    PrintOrbitalElements(M, expected[7], dM)
    print('\nTime of perihelion passage(days):\n')
    PrintOrbitalElements(T, expected[5], dT)


    
    return a, e, i, omega, w, M, T

def SEL(taus, Sun2, rhohat2, Ds):
    roots = [0., 0., 0.]
    rhos = [0., 0., 0.]
    
    E = -2 * (DotVectors(rhohat2, Sun2))
    F = GetMagnitude(Sun2) ** 2
    A1 = taus[1] / taus[2]
    B1 = A1 / 6 * (taus[2] ** 2 - taus[1] ** 2)
    A3 = -taus[0] / taus[2]
    B3 = A3 / 6 * (taus[2] ** 2 - taus[0] ** 2)
    
    A = (A1 * Ds[1] - Ds[2] + A3 * Ds[3]) / (-Ds[0])
    B = (B1 * Ds[1] + B3 * Ds[3]) / (-Ds[0])
    
    a = -(A ** 2 + A * E + F)
    b = -(2 * A * B + B * E)
    c = -(B ** 2)
    
    poly = [1, 0, a, 0, 0, b, 0, 0, c]
    roots = [np.real(x) for x in np.roots(poly) if (np.isreal(x) and x > 0)]
    roots.sort()
    
    rhos = [A + B / x ** 3 for x in roots]
    rhos.sort()
    
    return roots, rhos

# Not working yet
def GenerateEphemeris(path):
    ad = InputCSV(path)
    a, e, i, omega, w, M, T, da, de, di, domega, dw, dM, dT = CalculateOrbitalElements(path)
    realRADec = ad[1][0:6]
    rSun = ad[1][6:9]
    
    rRA = HMStoDegOrRad(realRADec[0], realRADec[1], realRADec[2], False)
    rDec = radians(DMStoDeg(realRADec[3], realRADec[4], realRADec[5]))
    inp = input('Enter input as: [year] [month] [day] [hour] [minute] [second] -- ')
    inp = inp.split()
    inp = [int(x) for x in inp]
    
    ts = pd.Timestamp(year = inp[0],  month = inp[1], day = inp[2], hour = inp[3], second = inp[4], tz = 'UTC')
    ts = ts.to_julian_date()
    M = (k * np.sqrt(1 / a ** 3) * (ts - T))
    
    E = NewtonMethod(lambda x: M - (x - e * np.sin(x)), lambda x: e * np.cos(x) - 1, M, 1e-10)
    while abs((E - e * np.sin(E)) - M) > 1e-10:
        E = NewtonMethod(lambda x: M - (x - e * np.sin(x)), lambda x: e * np.cos(x) - 1, M, 1e-10)
    x = a * np.cos(E) - a * e
    y = a * np.sqrt(1 - e ** 2) * np.sin(E)
    z = 0
    v = np.array([x, y, z])

    omega = radians(omega)
    i = radians(i)
    w = radians(w)

    rotMat1 = np.array([[np.cos(omega), -np.sin(omega), 0], [np.sin(omega), np.cos(omega), 0], [0, 0, 1]])   
    rotMat2 = np.array([[1, 0, 0], [0, np.cos(i), -np.sin(i)], [0, np.sin(i), np.cos(i)]])
    rotMat3 = np.array([[np.cos(w), -np.sin(w), 0], [np.sin(w), np.cos(w), 0], [0, 0, 1]])
    rotMat4 = np.array([[1, 0, 0], [0, np.cos(Ep), -np.sin(Ep)], [0, np.sin(Ep), np.cos(Ep)]])
    
    rot1 = np.matmul(rotMat3, v)
    rot2 = np.matmul(rotMat2, rot1)
    v1 = np.matmul(rotMat1, rot2)
    v2 = np.matmul(rotMat4, v1)
    rho = v2 - np.array(rSun)
    rhoHat = rho / GetMagnitude(rho)
    dec = degrees(np.arcsin(rhoHat[2]))
    RA = degrees(FindAngle(rhoHat[1] / np.cos(dec), rhoHat[0] / np.cos(dec)))
    errRA = abs((degrees(rRA) - RA) / degrees(rRA)) * 100
    errDec = abs((degrees(rDec) - dec) / degrees(rDec)) * 100
    
    return RA, dec, errRA, errDec

def ConvertRAsDecs(ad):
    RAs = np.array([HMStoDegOrRad(ad[i][6], ad[i][7], ad[i][8], False) for i in range(len(ad) - 1)])
    Decs = np.array([radians(DMStoDeg(ad[i][9], ad[i][10], ad[i][11])) for i in range(len(ad) - 1)])
    #print(RAs, Decs)
    return RAs, Decs

def CalculateRhoHats(RAs, Decs):
    rhoHats = np.array([np.array([np.cos(Decs[i]) * np.cos(RAs[i]), np.sin(RAs[i]) * np.cos(Decs[i]), np.sin(Decs[i])]) for i in range(3)])
    return RAs, Decs, rhoHats
    
def CalculateDs(path, rhoHats):
    ad = InputCSV(path)
    
    Rsuns = np.array([ad[i][12:] for i in range(len(ad) - 1)])
    
    D0 = DotVectors(rhoHats[0], CrossVectors(rhoHats[1], rhoHats[2]))
    D1s = np.array([DotVectors(CrossVectors(Rsuns[i], rhoHats[1]), rhoHats[2]) for i in range(3)])
    D2s = np.array([DotVectors(CrossVectors(rhoHats[0], Rsuns[i]), rhoHats[2]) for i in range(3)])
    D3s = np.array([DotVectors(rhoHats[0], CrossVectors(rhoHats[1], Rsuns[i])) for i in range(3)])
    return D0, D1s, D2s, D3s

def CalculateFGs(tau, r2, r2dot, flag):
    if flag == 2:
        u = 1 / r2 ** 3
        f = 1 - (u / (2 * r2 ** 3)) * tau ** 2
        g = tau - (u / (6 * r2 ** 3)) * tau ** 3
        return f, g
    
    u = 1 / GetMagnitude(r2) ** 3
    
    z = DotVectors(r2, r2dot) / (GetMagnitude(r2) ** 2)
    q = DotVectors(r2, r2) / (GetMagnitude(r2) ** 2) - u
    
    f = 1 - 0.5 * u * tau ** 2 + 0.5 * u * z * tau ** 3
    g = tau - (1 / 6) * u * tau ** 3
    
    if flag == 4:
        f += (1 / 24) * (3 * u * q - 15 * u * z ** 2 + u ** 2) * tau ** 4
        g += (1 / 4) * u * z * tau ** 4
    return f, g

def Get_RA_Dec_RMS(path):
    table = fits.open(path)[1].data
    n = int(table.shape[0])
    rmsRA = np.sqrt(sum([(table.index_ra[x] - table.field_ra[x]) ** 2 for x in range(n)]) / (n - 1))
    rmsDec = np.sqrt(sum([(table.index_dec[x] - table.field_dec[x]) ** 2 for x in range(n)]) / (n - 1))
    return radians(rmsRA), radians(rmsDec)

def MOG(taus, Sun2, rhoHats, D0, D1s, D2s, D3s, times, ad):
    Ds = [D0, D2s[0], D2s[1], D2s[2]]
    #print(Ds)
    rhohat2 = rhoHats[1]

    #print(taus, Sun2, rhohat2)
    roots, x = SEL(taus, Sun2, rhohat2, Ds)

    timesCopy = np.copy(times)
    times0 = timesCopy[0]
    times1 = timesCopy[1]
    times2 = timesCopy[2]
    
    bestDiff = abs(roots[0] - 1.3)
    for i in range(len(roots)):
        if abs(roots[i] - 1.3) < bestDiff:
            bestDiff = abs(roots[i] - 1.3)
    for i in range(len(roots)):
        if abs(roots[i] - 1.3) == bestDiff:
            bestRoot = roots[i]
    r2Curr = bestRoot
    #print(bestRoot)
    
    f1, g1 = CalculateFGs(taus[0], bestRoot, [], 2)
    f3, g3 = CalculateFGs(taus[1], bestRoot, [], 2)
    
    c1 = g3 / (f1 * g3 - g1 * f3)
    c2 = -1
    c3 = -g1 / (f1 * g3 - g1 * f3)


    rho1 = (c1 * D1s[0] + c2 * D1s[1] + c3 * D1s[2]) / (c1 * D0)
    rho2 = (c1 * D2s[0] + c2 * D2s[1] + c3 * D2s[2]) / (c2 * D0)
    rho3 = (c1 * D3s[0] + c2 * D3s[1] + c3 * D3s[2]) / (c3 * D0)

    r1Vec = rho1 * np.array(rhoHats[0]) - np.array(ad[0][12:])
    r2Vec = rho2 * np.array(rhoHats[1]) - np.array(Sun2)
    r3Vec = rho3 * np.array(rhoHats[2]) - np.array(ad[2][12:])
    
    r2Better = GetMagnitude(r2Vec)
    
    d1 = -f3 / (f1 * g3 - f3 * g1)
    d3 = f1 / (f1 * g3 - f3 * g1)
    
    r2DotVec = d1 * r1Vec + d3 * r3Vec

    timesCopy[0] = times0 - (rho1 / cAU)
    timesCopy[1] = times1 - (rho2 / cAU)
    timesCopy[2] = times2 - (rho3 / cAU)
    taus = [k * (timesCopy[0] - timesCopy[1]), k * (timesCopy[2] - timesCopy[1]), k * (timesCopy[2] - timesCopy[0])]
    f1, g1 = CalculateFGs(taus[0], r2Vec, r2DotVec, 4)
    f3, g3 = CalculateFGs(taus[1], r2Vec, r2DotVec, 4)
    #print(taus[0], r2Vec, r2DotVec)
    
    #print(f1, g1, f3, g3)
    count = 0
    while abs(r2Curr - r2Better) > 1e-50: 
        taus = [k * (timesCopy[0] - timesCopy[1]), k * (timesCopy[2] - timesCopy[1]), k * (timesCopy[2] - timesCopy[0])]
        f1, g1 = CalculateFGs(taus[0], r2Vec, r2DotVec, 4)
        f3, g3 = CalculateFGs(taus[1], r2Vec, r2DotVec, 4)

        #print(f1, g1, f3, g3)


        c1 = g3 / (f1 * g3 - g1 * f3)
        c2 = -1
        c3 = -g1 / (f1 * g3 - g1 * f3)

        rho1 = (c1 * D1s[0] + c2 * D1s[1] + c3 * D1s[2]) / (c1 * D0)
        rho2 = (c1 * D2s[0] + c2 * D2s[1] + c3 * D2s[2]) / (c2 * D0)
        rho3 = (c1 * D3s[0] + c2 * D3s[1] + c3 * D3s[2]) / (c3 * D0)

        r1Vec = rho1 * np.array(rhoHats[0]) - np.array(ad[0][12:])
        r2Vec = rho2 * np.array(rhoHats[1]) - np.array(Sun2)
        r3Vec = rho3 * np.array(rhoHats[2]) - np.array(ad[2][12:])
        r2Curr = r2Better
        r2Better = GetMagnitude(r2Vec)

        d1 = -f3 / (f1 * g3 - f3 * g1)
        d3 = f1 / (f1 * g3 - f3 * g1)
        
        r2DotVec = d1 * r1Vec + d3 * r3Vec

        timesCopy[0] = times0 - (rho1 / cAU)
        timesCopy[1] = times1 - (rho2 / cAU)
        timesCopy[2] = times2 - (rho3 / cAU)

        if count > 100:
            return [0, 0, 0, 0, 0, 0, 0]
        
        count += 1
    
    rotMatToEclip = np.array([[1, 0, 0], [0, np.cos(Ep), np.sin(Ep)], [0, -np.sin(Ep), np.cos(Ep)]])
    r2Vec = np.matmul(rotMatToEclip, r2Vec)
    r2DotVec = np.matmul(rotMatToEclip, r2DotVec)
    #print(r2Vec, r2DotVec)
    
    return CalculateOrbitalElements(path, r2Vec, r2DotVec, timesCopy[1])

def MonteCarlo(path):
    ad = InputCSV(path)
    
    rmsRA1, rmsDec1 = Get_RA_Dec_RMS('OD_Input/corr1.fits')
    rmsRA2, rmsDec2 = Get_RA_Dec_RMS('OD_Input/corr2.fits')
    rmsRA3, rmsDec3 = Get_RA_Dec_RMS('OD_Input/corr3.fits')
    
    rmsRAs = [rmsRA1, rmsRA2, rmsRA3]
    rmsDecs = [rmsDec1, rmsDec2, rmsDec3]
    rRAs, rDecs = ConvertRAsDecs(ad)
    
    N = 1000
    
    aArr = []
    eArr = []
    iArr = []
    omegaArr = []
    wArr = []
    MArr = []
    TArr = []

    times = []
    for i in range(len(ad) - 1):
        ts = pd.Timestamp(year = int(ad[i][0]),  month = int(ad[i][1]), day = int(ad[i][2]), hour = int(ad[i][3]), minute = int(ad[i][4]), second = int(np.trunc(ad[i][5])), tz = 'UTC')
        ts = ts.to_julian_date()
        ts += (ad[i][5] % 1) / 86400
        times.append(ts)

    taus = [k * (times[0] - times[1]), k * (times[2] - times[1]), k * (times[2] - times[0])]

    Sun2 = ad[1][12:]
    x, y, rhoHats = CalculateRhoHats(rRAs, rDecs)
    D0, D1s, D2s, D3s = CalculateDs(path, rhoHats)

    MOG(taus, Sun2, rhoHats, D0, D1s, D2s, D3s, times, ad)
    
    '''
    for i in range(N):
        RAs = [np.random.normal(loc = rRAs[i], scale = rmsRAs[i]) for i in range(len(rmsRAs))]
        Decs = [np.random.normal(loc = rDecs[i], scale = rmsDecs[i]) for i in range(len(rmsDecs))]
        
        x, y, rhoHats = CalculateRhoHats(RAs, Decs)
        D0, D1s, D2s, D3s = CalculateDs(path, rhoHats)

        a, e, i, omega, w, M, T = MOG(taus, Sun2, rhoHats, D0, D1s, D2s, D3s, times, ad)

        if a != 0:
            aArr.append(a)
            eArr.append(e)
            iArr.append(i)
            omegaArr.append(omega)
            wArr.append(w)
            MArr.append(M)
            TArr.append(T)
            N -= 1

    orbElements = [aArr, eArr, iArr, omegaArr, wArr, MArr, TArr]
    aSTD = np.std(aArr) / len(aArr)
    aMean = sum(aArr) / len(aArr)

    orbElementsSTD = [np.std(orbElements[i]) / len(orbElements[i]) for i in range(len(orbElements))]
    orbElementsMean = [sum(orbElements[i]) / len(orbElements[i]) for i in range(len(orbElements))]
    '''
    
MonteCarlo(path)
