# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

#START = 4
#END = 5

#bruh
print("==================================================================")
nameOfFile = "dos" + str(5) + ".csv"
df = pd.read_csv(str("C:\\Users\\baron\\Desktop\\" + nameOfFile))
    
print(nameOfFile, ":")
    
#Get Columns
energy = df['energy'][:]
dosUp = df['DOS (up)'][:]
dosDown = df['DOS (down)'][:]
intDosUp = df['integrated DOS (up)'][:]
intDosDown = df['integrated DOS (down)'][:]
eFermi = df['eFermi'][:]


#okay method
numUp = 0
numDwn = 0
num = 0

j = 1
while (j < len(energy)):
    if (energy[j] > eFermi[j]):
        numUp = abs(numUp)
        numDwn = abs(numDwn)
        num = abs(num)
        break
    
    numUp = numUp + dosUp[j] * (energy[j] - energy[j - 1])
    numDwn = numDwn + dosDown[j] * (energy[j] - energy[j - 1])
    num = num + ((dosUp[j] - dosDown[j]) * (energy[j] - energy[j - 1]))

    
    j = j + 1

print(numUp, numDwn, numUp - numDwn, numUp + numDwn)

#better method
numUp = 0
numDwn = 0
num = 0

j = 1
while (j < len(energy)):
    if (energy[j] > eFermi[j]):
        break
    #up:
    fit = np.polyfit([energy[j], energy[j-1]], [dosUp[j], dosUp[j-1]], 1)
    m = fit[0]
    b = fit[1]
    
    area = (m / 2 * (energy[j]**2 - energy[j-1]**2)) + (b * (energy[j] - energy[j-1]))
    numUp = numUp + area
    
    #down:
    fit = np.polyfit([energy[j], energy[j-1]], [dosDown[j], dosDown[j-1]], 1)
    m = fit[0]
    b = fit[1]
    
    area = (m / 2 * (energy[j]**2 - energy[j-1]**2)) + (b * (energy[j] - energy[j-1]))
    numDwn = numDwn + area
    j = j + 1
    
num = numUp - numDwn
print(numUp, numDwn, numUp - numDwn, numUp + numDwn)

    
