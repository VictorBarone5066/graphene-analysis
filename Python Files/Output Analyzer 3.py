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
df1 = pd.read_csv('C:\\Users\\baron\\Desktop\\output.csv')
df1 = df1.dropna(how='all')
sortedDataFrame = df1.sort_values('name for determining information')  #TODO: figure out why this stuff dosent work
df1 = sortedDataFrame
df1.to_csv('C:\\Users\\baron\\Desktop\\tmp.csv')
df = pd.read_csv('C:\\Users\\baron\\Desktop\\tmp.csv')
writer = open("C:\\Users\\baron\\Desktop\\LatticeConstants.csv", "w+")
writer2 = open("C:\\Users\\baron\\Desktop\\FormationEnergies.csv", "w+")
writer.write("Model Name,Min E Lat Const (Thry),Min E Lat Const (file 4),Thry - 2.46")
writer2.write("Model Name,Formation Energy")

#Getting the columns, saving to lists
energy = df['energy']
displacement = df['displacement from original positions'] #displacement from 2.46 angstroms
names = df['name for determining information']
cAtoms = df['CAtoms']    
area = df['area']
ids = df['id Code']
lat = df['lat Const']

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
q = []
gammaInfo = [] #Create a list (of lists) of areas and associated names for use in the gamma value loop.  
latConstantInfo = [] #Create a list (of lists) of ID nums, names, and min energy lattice constants for another graph
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
            q.append(lat[i])
            
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
                q.append(lat[i])
            
            #turn each value into an ordered pair to avoid bad graphs (trust me, dont change this)
            j = 0
            while (j < len(x)):
                lists.append([x[j], y[j], z[j], u[j], v[j], w[j], q[j]])
                j = j+1            
            lists.sort()
            _x_=[x[0] for x in lists]
            _y_=[y[1] for y in lists]
            _z_=[z[2] for z in lists] 
            _u_=[u[3] for u in lists]
            _v_=[v[4] for v in lists]
            _w_=[w[5] for w in lists]
            _q_=[q[6] for q in lists]
            
            if (half):
                #note:  this syntax ( a = b[:] ) creates a copy of b and gives it to a, instead of creating a pointer.  This is important, so do not change it, even though 
                #it looks confusing
                _x_left=_x_[:]
                del _x_left[6:]
                _x_right=_x_[:]
                del _x_right[:6]
                _y_left=_y_[:]
                del _y_left[6:]
                _y_right=_y_[:]
                del _y_right[:6]
            
            #getting number of CC bonds; should all be the same so it dosent (shouldent) matter which index you grab from
            nCAtoms = _z_[4]
            
            #Get the value that corresponds to the starting strain - i.e. file 4
            j = 0
            while (j < len(v)):
                if (int(v[j]) == 4):
                    startingLatConst = q[j]
                j = j + 1
            
            poly_ = np.polyfit(_x_, _y_, 2)
            minX = -poly_[1] / (2 * poly_[0])
            polyMin = (poly_[0] * minX * minX) + (poly_[1] * minX) + (poly_[2])
            print("Name: ", names[i], _v_[minPos(_y_)],"\nPolynomial Minimum: ", polyMin, "\nFile Minimum: ", min(_y_))
            
            #Get _x_ and _y_ into proper units:  _x_ = true strain (not displacement from 2.46) and _y_ = meV diff from min
            _x_ = [float(g) for g in _x_]
            _y_ = [float(1000 * (g - polyMin)) / float(nCAtoms) for g in _y_]
            
            #Get the polynomial fits for the graphs
            polyFit = np.polyfit(_x_, _y_, 2)
            xPoly = np.linspace(_x_[0], _x_[-1], 50)
            xPoly = np.ndarray.tolist(xPoly)
            
            j = 0
            yPoly = []
            while (j < len(xPoly)):
                yPoly.append((xPoly[j] * xPoly[j] * polyFit[0]) + (xPoly[j] * polyFit[1]) + (polyFit[2]))
                j = j + 1
                
            
            # <codecell> Strain Graphs

#-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -  
            if (half):
                #Set left and right components to the correct units; x is true strain and y is meV
                _x_left = [float(g) for g in _x_left]
                _x_right = [float(g) for g in _x_right]
                _y_left = [float(1000 * (g - polyMin)) / float(nCAtoms) for g in _y_left]
                _y_right = [float(1000 * (g - polyMin)) / float(nCAtoms) for g in _y_right]
                
                #Get the polynomial fits
                polyFitL = np.polyfit(_x_left, _y_left, 2)
                polyFitR = np.polyfit(_x_right, _y_right, 2)
                xPolyL = np.linspace(_x_left[0], _x_left[-1], 50)
                xPolyR = np.linspace(_x_right[0], _x_right[-1], 50)
                xPolyL = np.ndarray.tolist(xPolyL)
                xPolyR = np.ndarray.tolist(xPolyR)
                
                j = 0
                yPolyL = []
                while (j < len(xPolyL)):
                    yPolyL.append((xPolyL[j] * xPolyL[j] * polyFitL[0]) + (xPolyL[j] * polyFitL[1]) + (polyFitL[2]))
                    j = j + 1
                j = 0
                yPolyR = []
                while (j < len(xPolyR)):
                    yPolyR.append((xPolyR[j] * xPolyR[j] * polyFitR[0]) + (xPolyR[j] * polyFitR[1]) + (polyFitR[2]))
                    j = j + 1
                    
#-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -                 
            
            #Name Stuff 
            newname_ = names[i].lstrip("1D").lstrip('2D').rstrip('X').rstrip('Y').strip('_').strip("Bulk")    #clear garbage from name
            newname_ = re.sub(r"(?<=\w)([A-Z])", r" \1", newname_) #insert spaces in the file name
            if (newname_[0] == 'i'):
                newname_ = 'D' + newname_
            plot_ = plt.gcf()
            plt.title(newname_, x = 0.5, y = 0.9)
            
            #More name edits
            multi = False
            if (names[i][-1] == "X"):
                plt.xlabel('Zig-Zag Strain')
            elif (names[i][-1] == "Y"):
                plt.xlabel('Armchair Strain')
            else:
                multi = True
                plt.xlabel('MultiDirectional Strain')
            
            #Plot figs
            plt.ylabel('Energy from Min. per C Atom (meV)')
            plt.plot(xPoly, yPoly, color='c') #plot polynomial fit 
            plt.plot(_x_, _y_, color='k', marker='o', linewidth=0) #plot data points
            plt.plot(0, 0, color='#ff00ff', marker='|') #mark zero strain point
            plt.axvline(x = (LATCNST / startingLatConst) - 1, ymax = 0.85, color='#ff00ff', linewidth = 0.25) #Plot line for 2.46 Angst lattice constant
           
            #Special Case for Vac:
            if (((newname_[0] == "V") and (newname_[1] == "a") and (newname_[2] == "c")) and not multi and half):      
                plt.plot(xPolyL, yPolyL, '--', color = 'forestgreen')
                plt.plot(xPolyR, yPolyR, '--', color = 'forestgreen')

            plt.show()
            
            if not os.path.exists('C:\\Users\\baron\\Desktop\\figs\\'):
                os.makedirs('C:\\Users\\baron\\Desktop\\figs\\')
            plot_.savefig('C:\\Users\\baron\\Desktop\\figs\\' + names[i] + 'Strain.jpg', dpi = 95)

            # <codecell> gamma value Prep & Energy vs Area
            
            #To find Gamma Values, need to use the surface area of the minimum energy computational cell.
            #Since the area (should) scale linerally with the applied displacement, can fit all areas to a line:
            # area = m(strain) + b
            # where b is set to the lattice constant before strain.

            #Find the minimum energy area, Get Polynomial for graphing
            polyFitA = np.polyfit(_u_, _y_, 2)
            xPolyA = np.linspace(_u_[0], _u_[-1], 50)
            
            yPolyA = []
            j = 0
            while (j < len(xPolyA)):
                yPolyA.append((xPolyA[j] * xPolyA[j] * polyFitA[0]) + (xPolyA[j] * polyFitA[1]) + (polyFitA[2]))
                j = j + 1
        
            eMinArea = -polyFitA[1] / (2 * polyFitA[0])
            
            #Store areas and names
            gammaInfo.append([names[i], eMinArea])
            
            # <codecell> Formation Energy vs. Area Plots
#            plt.xlabel(r"Cell Area (${\mathring{A}}^2$)")
#            plt.ylabel('Energy from Min. per C Atom (meV)')
#            plt.axvline(x = _u_[minPos(_y_)], color='#ff00ff', ymax = 0.85, linewidth = 0.25)
#            plt.axvline(x = eMinArea, color='#ff00ff', ymax = 0.85)
#            
#            #More name edits
#            if (names[i][-1] == "X"):
#                plt.title(newname_ + ' (Zig-Zag Strain)', x = 0.5, y = 0.9)
#            elif (names[i][-1] == "Y"):
#                plt.title(newname_ + ' (Armchair Strain)', x = 0.5, y = 0.9)
#            else:
#                plt.title(newname_ + ' (MultiDirectional Strain)', x = 0.5, y = 0.9)
#            
#            plt.plot()
#            
#            plt.plot(_u_, _y_, color='k', marker = 'o', linewidth = 0)
#            plt.plot(xPolyA, yPolyA, color = 'c')
#            plt.plot()
#            plt.show()
            
            print("-------------------------------------------------------")
            
            # <codecell> minimum Lattice constant output
            thryMinLatConst = startingLatConst * (1 + minX)
            writer.write("\n" + names[i] + "," + str(thryMinLatConst) + "," + str(startingLatConst) + "," + str(LATCNST - thryMinLatConst))
            
            thisLatInfo = [newname_, _w_[1], thryMinLatConst, _q_[minPos(_y_)]]
            latConstantInfo.append(thisLatInfo)
            
            # <codecell> formation energies
            
            #TODO: Need to find the theoretical energy that we would have at 2.46 Lattice constant?
            
            j = 0
            minE_ = 500
            while (j < len(_x_)):
                if (_y_[j] < minE_):
                    minE_ = _y_[j]
                    minELoc = j
                j = j + 1
            
            defFrmEnrgy = (_y_[minELoc] * nCAtoms / 1000 + polyMin) - (_z_[minELoc] * ENperATM_FOR_PURE)
            writer2.write("\n" + names[i] + "," + str(defFrmEnrgy))
            
            # <codecell> end of loop, cleanup
            
            #Clearing the many lists I have before going into the next iteration       
            Clear([x, y, z, u, v, q, _x_, _y_, _z_, _u_, _v_, w, _w_, _q_, xPoly, yPoly,  lists])
            if (half):
                Clear([_x_left, _x_right, _y_left, _y_right, xPolyL, xPolyR, yPolyL, yPolyR])
        i = i + 1
except KeyError:
    pass
#Delete temp csv file
os.remove('C:\\Users\\baron\\Desktop\\tmp.csv')
writer.close()
writer2.close()

# <codecell> GammaValue Loop

#Calculating Gamma values:  
#need to sort the info so that 1D is first, since g12 relies on g11
df1 = pd.read_csv('C:\\Users\\baron\\Desktop\\output.csv')
df1 = df1.sort_values(by='id Code', ascending=0)   
df1 = df1.dropna(how='all')                                 #need to make a new tmp file since for whatever reason
df1.to_csv("C:\\Users\\baron\\Desktop\\tmp.csv")            #even after sorting, reading in values to the lists
df0 = pd.read_csv('C:\\Users\\baron\\Desktop\\tmp.csv')     #dosent result in the lists being sorted
q = 0
a = []
b = []
c = []
d = []
res = []
writer = open("C:\\Users\\baron\\Desktop\\GammaValues.csv", "w+")
writer.write("Model Name,Gamma(11),Gamma(22),Gamma(12),Young(11),Young(22),Poisson(11),Poisson(22)")


energy = df0['energy']
displacement = df0['displacement from original positions']
names = df0['name for determining information']
idCode = df0["id Code"]

#This is basically the same loop as above, so dont really need comments for the if and elifs
try:
    while (q <= df0.shape[0] + 1):
        if (idCode[q] == idCode[q+1] and (idCode[q] != -1) and (displacement[q] != -99)):
            a.append(displacement[q])
            b.append(energy[q])
            c.append(names[q])
            d.append(idCode[q])
        elif (idCode[q] != idCode[q+1]) and (idCode[q] != -1):
            #append the final values
            if (displacement[q] != -99):
                a.append(displacement[q])
                b.append(energy[q])
                c.append(names[q])
                d.append(idCode[q])
    
            #turn each value into an ordered pair to avoid bad graphs (trust me, dont change this)
            r = 0
            while (r < len(a)):
                res.append([a[r], b[r], c[r], d[r]])
                r = r+1            
            res.sort()
            _a_=[a[0] for a in res]
            _b_=[b[1] for b in res]
            _c_=[c[2] for c in res] #not sure if these last 2 are necessary, but it (in theory) dosent hurt to redefine them 
            _d_=[d[3] for d in res]
            
            #Seperate 1D and 2D stuff into yet more containers
            #syntax: x1DY means, for example, the x axis of the 1-dimensional strain in the Y direction
            x1DX = [] 
            y1DX = []
            x1DY = []
            y1DY = []
            x2D = []
            y2D = []
            j = 0
            if (len(res[0]) == len(res[1]) == len(res[2]) == len(res[3])):
                pass
            else:
                raise Exception("The lengths of 'res []' are not all the same.  They should be the same")
     
            while (j < len(res)):
                if(res[j][2][0] == '1' and (res[j][2][-1] == 'X' or res[j][2][-1] == 'x')): #must be considering 1D x-directional strain
                    x1DX.append(res[j][0]) #append strain
                    y1DX.append(res[j][1]) #append energy
                    #get area
                    k = 0
                    while (k < len(gammaInfo)):
                        if (res[j][2] == gammaInfo[k][0]):
                            eMinArea1DX = gammaInfo[k][1]
                        k = k + 1
                elif(res[j][2][0] == '1' and (res[j][2][-1] == 'Y' or res[j][2][-1] == 'y')): #must be considering 1D y-directional strain
                    x1DY.append(res[j][0]) #append strain
                    y1DY.append(res[j][1]) #append energy
                    #get area
                    k = 0
                    while (k < len(gammaInfo)):
                        if (res[j][2] == gammaInfo[k][0]):
                            eMinArea1DY = gammaInfo[k][1]
                        k = k + 1
                elif(res[j][2][0] == '2'): #must be considering 2D
                    x2D.append(res[j][0]) #append strain
                    y2D.append(res[j][1]) #append energy
                    #get area
                    k = 0
                    while (k < len(gammaInfo)):
                        if (res[j][2] == gammaInfo[k][0]):
                            eMinArea2D = gammaInfo[k][1]
                        k = k + 1
                j = j + 1
                  
            #now have all of the info for a 1D model and its 2D counterpart...
            fit1DX = np.polyfit(x1DX, y1DX, 2)
            fit1DY = np.polyfit(x1DY, y1DY, 2)
            fit2D = np.polyfit(x2D, y2D, 2)
            g11 = (2 * fit1DX[0])/eMinArea1DX
            g22 = (2 * fit1DY[0])/eMinArea1DY
            
            if (g11 == 0 or g22 == 0):
                raise Exception ("either g11 or g22 is zero.  They should both be initialized by now.")
                
            g12 = (fit2D[0] / eMinArea2D) - ((1/2) * (g11 + g22)) 
                
            #Young's moduli and Poisson's ratio
            young1 = g11 * (1.6022E-19) / ((1E-10)*(1E-10)*(3.4E-10))
            young2 = g22 * (1.6022E-19) / ((1E-10)*(1E-10)*(3.4E-10))
            poisson1 = g12 / g22
            poisson2 = g12 / g11
            
            #Write values to file
            writer.write("\n" + str(names[q][2:]) + "," + str(g11) + "," + str(g22) + "," + str(g12) + "," + str(young1) + "," + str(young2) + "," + str(poisson1) + "," + str(poisson2))
            
            #Clear all of the lists
            Clear([a, b, c, d, _a_, _b_, _c_, _d_, res, x1DX, x1DY, x2D, y1DX, y1DY, y2D])
        q = q + 1
except KeyError:
    pass
    
os.remove('C:\\Users\\baron\\Desktop\\tmp.csv')
writer.close()

# <codecell> Min-Energy Lattice Constant Span Graph
#print (latConstantInfo)
i = 0
maxID = -99
while (i < len(latConstantInfo)):
    if (latConstantInfo[i][1] > maxID):
        maxID = latConstantInfo[i][1]
    i = i + 1
    
i = 0
graphX = []
graphY = []
graphZ = []
graphY_ = []
while (i < len(latConstantInfo)):
    graphX.append(latConstantInfo[i][1])
    graphY.append(latConstantInfo[i][2])
    graphZ.append(latConstantInfo[i][0])
    graphY_.append(latConstantInfo[i][3])
    i = i + 1
  
    
plt.show
plot___ = plt.gcf()

plt.plot(graphX, graphY, marker = 'o', linewidth = 0, color = 'k') #plot the lat const as estimated from the polynomial fit
#plt.plot(graphX, graphY_, marker = '+', linewidth = 0, color = '#ff00ff') #plot the lat const from POSCAR file
plt.axhline(y = LATCNST, color='c')
plt.title("Minimum Energy Lattice Constants")        
plt.xticks(graphX, graphZ, rotation = (45))


if not os.path.exists('C:\\Users\\baron\\Desktop\\figs\\'):
    os.makedirs('C:\\Users\\baron\\Desktop\\figs\\')
plot___.savefig('C:\\Users\\baron\\Desktop\\figs\\' + 'Lattice Constant Span.jpg', dpi = 95)


# <codecell> Table of Values:  Gamma Values











