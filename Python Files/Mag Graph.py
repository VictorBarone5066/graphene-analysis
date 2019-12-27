# For use with POSCAR files that have already been relaxed (volume relaxation) before being strained

#Imports:
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os

#Definitions:
def Clear (listOfLists = []):
    j = 0
    while (j < len(listOfLists)):
        listOfLists[j].clear()
        j = j + 1
    return

def minPos(list_ = []):
    i = 0
    toReturn = None
    mini = max(list_)
    while (i < len(list_)):
        if (list_[i] < mini):
            mini = list_[i]
            toReturn = i
        i = i + 1
        
    return toReturn
            
#Important Data Starting lattice constant:
LATCNST = 2.46 #angstroms
ENperATM_FOR_PURE = -9.233597199 #eV / atom
NRG_DENSITY_PURE = 3.523713443 #ev / angst^2

half = False
    
#Set up:  
df1 = pd.read_csv('C:\\Users\\baron\\Desktop\\outputM.csv')
df1 = df1.dropna(how='all')
sortedDataFrame = df1.sort_values('name for determining information')  #TODO: figure out why this stuff dosent work
df1 = sortedDataFrame
df1.to_csv('C:\\Users\\baron\\Desktop\\tmp.csv')
df = pd.read_csv('C:\\Users\\baron\\Desktop\\tmp.csv')

#Getting the columns, saving to lists
mag = df['mag']
displacement = df['displacement from original positions'] #displacement from 2.46 angstroms
names = df['name for determining information']

#Declerations
i = 0
x = []
y = []
z = []
lists = []
#Main loop:
try:
    while (i < (df.shape[0] + 1)):
        #As long as the names are the same, continue to add info to the same arrays
        if (names[i] == names[i+1] and displacement[i] != -99):
            x.append(float(mag[i]))
            y.append(float(displacement[i]))
            z.append(names[i])

            
        #If the names are different, we must be about to consider different calculations...so:
        elif (names[i] != names[i+1]) and (names[i] != 'endf'):
           # <codecell> Energy vs. Lattice Constant Graphs
            #append the final values
            if (displacement[i] != -99):
                x.append(float(mag[i]))
                y.append(float(displacement[i]))
                z.append(names[i])
    
            #turn each value into an ordered pair to avoid bad graphs (trust me, dont change this)
            j = 0
            while (j < len(x)):
                lists.append([float(y[j]), float(x[j])])
                j = j+1            
            lists.sort()
            _y_=[y[0] for y in lists]
            _x_=[x[1] for x in lists]
                
            
            # <codecell> Mag Graphs
            print(_y_, _x_)
            #Plot figs
            plt.title('Zig-Zag Strain: Unrelaxed Lattice')
            plt.xlabel(r"Fractional displacement from $\vec{a}$ = 2.46 Å")
            plt.ylabel(r'Magnetization ($\mu_{\mathrm{B}}$)')
            plt.plot(_y_, _x_, color='k', marker='o', linewidth=0) #plot data points
            plt.plot(_y_[:6], _x_[:6], color = 'c')
            plt.plot(_y_[6:], _x_[6:], color = 'c')
            plt.ylim(0,max(_x_) + 0.05)

            plot = plt.gcf()
            plot.savefig('C:\\Users\\baron\\Desktop\\mag.pdf', dpi = 95)            
            
            
            # <codecell> end of loop, cleanup
            
            #Clearing the many lists I have before going into the next iteration       
            Clear([x, y, z, _x_, _y_, lists])
        i = i + 1
except KeyError:
    pass
#Delete temp csv file
os.remove('C:\\Users\\baron\\Desktop\\tmp.csv')












