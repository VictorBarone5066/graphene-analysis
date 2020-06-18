# Creates figures mapping all strain values used for tests.  Much more complicated than necessary,
# since I just copied and pasted most of this code from previous Output Analyzer stuff
#Imports:
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
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
lis = []
lisMin = []
                #cyan               #green          #magenta
for name_ in ["Strain Set A", "Strain Set B"]:  
    #Set up:  
    df1 = pd.read_csv('C:\\Users\\baron\\Desktop\\' + name_ + '\\output.csv')
    df1 = df1.dropna(how='all')
    sortedDataFrame = df1.sort_values('name for determining information')  #TODO: figure out why this stuff dosent work
    df1 = sortedDataFrame
    df1.to_csv('C:\\Users\\baron\\Desktop\\tmp.csv')
    df = pd.read_csv('C:\\Users\\baron\\Desktop\\tmp.csv')
    
    #Getting the columns, saving to lists
    energy = df['energy']
    displacement = df['displacement from original positions'] #displacement from 2.46 angstroms
    names = df['name for determining information']
    cAtoms = df['CAtoms']    
    area = df['area']
    ids = df['id Code']
    xV = df['xV']
    yV = df['yV']
    
    #Also getting extra stuff for the formation energies
    workingDir = df['working directory (where POSCAR is)'] #technically dont need this but useful for checking that output is working correctly
    
    #Declerations
    i = 0
    x = []
    y = []
    z = []
    u = []
    v = []
    w = []
    r = []
    s = []
    lists = []
    
    #Main loop:
    try:
        while (i < (df.shape[0] + 1)):
            #As long as the names are the same, continue to add info to the same arrays
            if (names[i] == names[i+1] and displacement[i] != -99):
                x.append(displacement[i])
                y.append(energy[i])
                z.append(cAtoms[i])
                u.append(area[i])
                v.append(workingDir[i])
                w.append(ids[i])
                r.append(xV[i])
                s.append(yV[i])
                
            #If the names are different, we must be about to consider different calculations...so:
            elif (names[i] != names[i+1]) and (names[i] != 'endf'):
               # <codecell> Energy vs. Lattice Constant Graphs
                #append the final values
                if (displacement[i] != -99):
                    x.append(displacement[i])
                    y.append(energy[i])
                    z.append(cAtoms[i])
                    u.append(area[i])
                    v.append(workingDir[i])
                    w.append(ids[i])
                    r.append(xV[i])
                    s.append(yV[i])
                
                #turn each value into an ordered pair to avoid bad graphs (trust me, dont change this)
                j = 0
                while (j < len(x)):
                    lists.append([x[j], y[j], z[j], u[j], v[j], w[j], r[j], s[j]])
                    j = j+1            
                lists.sort()
                _x_=[x[0] for x in lists]
                _y_=[y[1] for y in lists]
                _z_=[z[2] for z in lists] 
                _u_=[u[3] for u in lists]
                _v_=[v[4] for v in lists]
                _w_=[w[5] for w in lists]
                _r_=[r[6] for r in lists]
                _s_=[s[7] for s in lists]
                    
                # <codecell> Create list of info needed for graphs
                
                #Name Stuff 
                newname_ = names[i].lstrip("1D").lstrip('2D').rstrip('X').rstrip('Y').strip('_').strip("Bulk")    #clear garbage from name
                newname_ = re.sub(r"(?<=\w)([A-Z])", r" \1", newname_) #insert spaces in the file name
                if (newname_[0] == 'i'):
                    newname_ = 'D' + newname_
            
                if (newname_[0] == 'V' and newname_[1] == 'a'):
                    lis.append([newname_, _r_[:], _s_[:]])

                
                # <codecell> end of loop, cleanup
                
                #Clearing the many lists I have before going into the next iteration       
                Clear([x, y, z, u, v, r, s, _x_, _y_, _z_, _u_, _v_, w, _w_, _r_, _s_, lists])
            i = i + 1
    except KeyError:
        pass
#Delete temp csv file
os.remove('C:\\Users\\baron\\Desktop\\tmp.csv')

# <codecell> Make the plot

plt.xlabel("Zig-Zag Supercell Vector (Å)")
plt.ylabel("Armchair Supercell Vector (Å)")

i = 0
colorList = ['c', 'm', 'g']
while (i < len(lis)):
    if (i%3 != 2):
        plt.plot(lis[i][1], lis[i][2], color = colorList[int(np.floor(i / len(colorList)))], marker = '.', linewidth = 0.25)
    i = i + 1

#plt.plot(14.6748131829999, 12.89053689700006, 'kx', markersize = 10) #prerelaxed min (file 4)
#plt.plot(14.7231002294, 12.78251934, 'kx', markersize = 10) #nonprerelaxed min (file 4.5)
#plt.plot() #final vac min (file ?)

plt.legend((r'Set $a$', '_', 'Set $b$')) #<<there has to be a better way of doing this.  Works because 'legend' ignores
                                                        #anything with '_' at the beginning
    
plot = plt.gcf()
plot.savefig('C:\\Users\\baron\\Desktop\\StrainMap.png', dpi = 300)













