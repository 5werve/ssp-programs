import numpy as np
import random as rp
import time
import matplotlib.pyplot as plt

def RandomWalkSim(N):
    rp.seed(time.time())
    
    x = y = 0
    xVals = np.array([0])
    yVals = np.array([0])
    
    for i in range(N):
        delX = rp.random() * 2 - 1
        delY = np.sqrt(1 - delX ** 2)
        if rp.random() > 0.5:
            delY = -delY

        #print(np.sqrt(delX ** 2 + delY ** 2))
        
        x += delX
        y += delY
        
        xVals = np.append(xVals, x)
        yVals = np.append(yVals, y)
        
        #plt.plot(xVals, yVals, 'ro')
        #plt.pause(0.1)

    return np.sqrt(xVals[-1] ** 2 + yVals[-1] ** 2) 

'''
N = 100

Rs = np.array([])   
for i in range(N):
    Rs = np.append(Rs, RandomWalkSim(N))

Rrms = np.sqrt(sum(Rs ** 2) / N)
print('Calculated Rrms:', Rrms)
print('Expected Rrms:', np.sqrt(N))
'''

N = 100

Rs = np.array([])   
for i in range(N):
    Rs = np.append(Rs, RandomWalkSim(N))

Ravg = sum(Rs) / N
print(Ravg)

plt.hist(Rs)
plt.show()

