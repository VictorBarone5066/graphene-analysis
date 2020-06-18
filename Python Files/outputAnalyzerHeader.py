#Stuff for Graphene

#Imports:
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import csv
import operator as op
import re

from sklearn.metrics import r2_score

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True

#Constants
UNIT_CELLS_PER_SUPERCELL_X = 6
UNIT_CELLS_PER_SUPERCELL_Y = 3

#Referring to working directory
VAC_X_LEFT_END = 5
VAC_X_RIGHT_BEG = 6
VAC_Y_LEFT_END = 4
VAC_Y_RIGHT_BEG = 5
VAC_XY_LEFT_END = 2
VAC_XY_RIGHT_BEG = 3


# <codecell> FUNCTIONS
def Clear (listOfLists = []):
    j = 0
    while (j < len(listOfLists)):
        listOfLists[j].clear()
        j = j + 1
    return

def MinPos(list_ = []):
    i = 0
    toReturn = None
    mini = max(list_)
    while (i < len(list_)):
        if (list_[i] < mini):
            mini = list_[i]
            toReturn = i
        i = i + 1 
    return toReturn
            
def GetStrainType(string):
    if string[-1] == 'X':
        return 'X'
    elif (string[-1] == 'Y'):
        return 'Y'
    elif (string[0] == '2' and string[1] == 'D'):
        return 'XY'
    else:
        return 'UNDEFINED'  
    
#if order==1: returns xPoly, yPoly, fit params.  if order==2: returns xPoly, yPoly, minX, minY, fit params
def GetPolyEqtn(x, y, order, numPts = 100, findR2 = False, findStdErr = False):
    try:
        fit = np.polyfit(x, y, order, cov=findStdErr)
    except(ValueError):
        fit = [[1,1,1],[[1,1,1,],[1,1,1,],[1,1,1,]]]
    fit_ = np.polyfit(x, y, order)
    #Finds r^2 for a 2nd order polynomial fit
    if(findR2):
        r2 = r2_score(y, [(fit_[0] * x_**2) + (fit_[1] * x_) + (fit_[2]) for x_ in x])

    xPoly = list(np.linspace(min(x), max(x), numPts))
    yPoly = []
    i = 0

    #Finds the standard error (of the coefficient of x^2) for 2nd order polynomial fit
    if(findStdErr):
        cov00 = fit[1][0][0] #the element in the first row and column of the covariance matrix
        stdDev = cov00**(1/2)
    
    if (order == 1):
        while (i < len(xPoly)):
            yPoly.append((fit_[0] * xPoly[i]) + (fit_[1])) #mx + b
            i = i + 1
        return xPoly, yPoly, fit_
    
    if (order == 2):
        while(i < len(xPoly)):
            yPoly.append((fit_[0] * xPoly[i]**2) + (fit_[1] * xPoly[i]) + (fit_[2])) #ax^2 + bx + c
            i = i + 1
        minX = (-fit_[1]) / (2 * fit_[0]) #-b / 2a
        minY = (fit_[0] * minX**2) + (fit_[1] * minX) + (fit_[2])
        if(not findR2):
            return xPoly, yPoly, minX, minY, fit_
        elif(findR2 and not findStdErr):
            return xPoly, yPoly, minX, minY, fit_, r2
        else:
            return xPoly, yPoly, minX, minY, fit_, r2, stdDev
    if (order < 1 or order > 2):
        print('GetPolyEqtn:  Order Undefined')
    
    return -1        

def GetDisplayName(ret, afterMarch = False):
    if(re.search("Pure", ret)):
        return "Pure"
    if(re.search("SingleVac", ret) or re.search("SV", ret)):
        return "SV"
    if(re.search("DiVac", ret) or re.search("DV", ret)):
        if(afterMarch):
            return r"DV$^{b}$"
        else:
            return r"DV$^{a}$"
    if(re.search("555777", ret)):
        return "555777"
    if(re.search("StoneW", ret) or re.search("SW", ret)):
        if(afterMarch):
            return r"STW$^{b}$"
        else:
            return r"STW$^{a}$"
    return "???"
    
        
def FindLATPosition(strain, strainSet, afterMarch, LAT = 2.468217817, initLatConst = True):
    if (strainSet == 'A'):
        initLatConst = LAT
    if (strainSet == 'B'):
        if(afterMarch):         ##the defects made after march have 0 displacement at working dir 2 for set B
            initLatConst = strain.GetLatConstByWD(2) 
        else:
            initLatConst = strain.GetLatConstByWD(4)
    strain.InitializeMinEnergyLatticeConstant()
    minLatConst = strain.minEnergyLatticeConstantPoly
    if(initLatConst):
        return ((LAT / minLatConst) - 1), initLatConst
    else:
        return ((LAT / minLatConst) - 1)
            

def GetAxisName(s):
    ret = '???'
    if (s[-1] == "X"):
        ret = r"Zig-Zag Strain"
    elif (s[-1] == "Y"):
        ret = r"Armchair Strain"
    else:
        ret = r"Isotopic Strain"    
    return ret

def DetermineSkewed(s):
    if(re.search("Skew", s) or re.search("skew", s)):
        return True
    return False

def DetermineAfterMarch(s):
    if(re.search("AfterMarch", s) or re.search("afterMarch", s)):
        return True
    return False

# <codecell> CLASSES
    ##corresponds to a line in the output file
class dataLine:
    directory = None
    energy = None
    workingDir = None
    disp = None
    name = None
    idNum = None    
    area = None
    cAtoms = None
    latConst = None
    xV, yV, zV = None, None, None
    skewed = False
    afterMarch = False
    
    def __init__(self, lis = [None, None, None, None, None, None, None, None, None, None, None, None]):
            self.directory = str(lis[0])
            self.energy = float(lis[1])
            self.workingDir = float(lis[2])
            self.disp = float(lis[3])
            self.name = str(lis[4])
            self.idNum = int(lis[5])    
            self.area = float(lis[6])
            self.cAtoms = int(lis[7])
            self.latConst = float(lis[8])
            self.xV, self.yV, self.zV = float(lis[9]), float(lis[10]), float(lis[11])
            
            if (DetermineSkewed(self.directory)):
                self.skewed = True
                
            if (DetermineAfterMarch(self.directory)):
                self.afterMarch = True
                
    ##corresponds to a strain type (all of the data lines for single vac X strain, for example)
class strainGroup:
    strainGroup = None
    name = None
    idNum = None
    strainType = None
    minEnergyLatticeConstantPoly = None
    minEnergyLatticeConstantFile = None
    minEnergyLatticeConstantBoth = [None, None] #in the [x, y] directions
    minEnergyArea = None
    formationEnergy = None
    vacSplit = False
    skewed = False
    afterMarch = False
    
    
    def __init__(self, lis = [None, None, None, None]):
        self.strainGroup = lis[0]
        self.name = lis[1]
        self.idNum = lis[2]
        self.strainType = lis[3]
        
        if(self.strainType == 'X' or self.strainType == 'XY'):
            for l in self.strainGroup:
                l.latConst = l.xV / UNIT_CELLS_PER_SUPERCELL_X
        if(self.strainType == 'Y'):
            for l in self.strainGroup:
                l.latConst = l.yV / UNIT_CELLS_PER_SUPERCELL_Y / 3**(1/2)
                
        if(self.name[:9] == '1DVacBulk'):
            self.vacSplit = True
            
        if(DetermineSkewed(self.name)):
            self.skewed = True
        
        #python is being dumb.  strainGroup is absolutly NOT an int, and it absolutly IS s subscriptable list
        try:
            if(self.strainGroup[0].afterMarch):
                self.afterMarch = True
        except(TypeError):
            pass
        
    def InitializeMinEnergyLatticeConstant(self):
        #Get min from poly fit and files
        xV, yV, energy = [], [], []
        for l in self.strainGroup:
            xV.append(l.xV)
            yV.append(l.yV)
            energy.append(l.energy)
        if (self.strainType == 'X' or self.strainType == 'XY'):
            x, y, minX, minY, f = GetPolyEqtn(xV, energy, 2)
            self.minEnergyLatticeConstantPoly = minX / UNIT_CELLS_PER_SUPERCELL_X
            self.minEnergyLatticeConstantFile = self.strainGroup[MinPos(energy)].xV / UNIT_CELLS_PER_SUPERCELL_X
        elif (self.strainType == 'Y'):
            x, y, minX, minY, f = GetPolyEqtn(yV, energy, 2)
            self.minEnergyLatticeConstantPoly = minX / UNIT_CELLS_PER_SUPERCELL_Y / 3**(1/2)
            self.minEnergyLatticeConstantFile = self.strainGroup[MinPos(energy)].yV / UNIT_CELLS_PER_SUPERCELL_Y / 3**(1/2)
        else:
            self.minEnergyLatticeConstantPoly = -1
            self.minEnergyLatticeConstantFile = -1
    
    def InitializeMinEnergyArea(self):
        area, energy = [], []
        for l in self.strainGroup:
            area.append(l.area)
            energy.append(l.energy)   
        x, y, minX, minY, f = GetPolyEqtn(area, energy, 2)
        self.minEnergyArea = minX
        
    def InitializeFormationEnergy(self, pureETot, strainSet, pureNAtoms = 72):
        #poly fit to find the theoretical minimum energy
        if(strainSet == 'A'): ##we want to find the formation energy of a defect for pure lattice constant conditions
            for l in self.strainGroup:
                if(float(l.workingDir) == 5.0):
                    self.formationEnergy = l.energy - (pureETot * (l.cAtoms / pureNAtoms))
                    return
        elif(strainSet == 'B'): ##we want to find the minimum energy formation energy
            _, _, _, minEDefect, _ = GetPolyEqtn([l.disp for l in self.strainGroup], [l.energy for l in self.strainGroup], 2)
            self.formationEnergy = minEDefect - (pureETot * (self.strainGroup[0].cAtoms / pureNAtoms))
        else:
            self.formationEnergy = -99
        
    def GetLatConstByWD(self, id_):
        if(self.strainType == 'X' or self.strainType == 'XY'):
            for l in self.strainGroup:
                if (l.workingDir == float(id_)):
                    return l.xV / UNIT_CELLS_PER_SUPERCELL_X
        elif(self.strainType == 'Y'):
            for l in self.strainGroup:
                if (l.workingDir == float(id_)):
                    return l.yV / UNIT_CELLS_PER_SUPERCELL_Y / 3**(1/2)
        else:
            return -1
    
    ##corresponds to a entire defect group (each of the strain groups for single vac X, Y, and XY strains)    
class defectGroup:
    xStrain = None
    yStrain = None
    xyStrain = None
    name = None
    idNum = None
    skewed = False
    afterMarch = False
    
    def __init__(self, lis = [None, None, None, None, None]):
        self.xStrain = lis[0]
        self.yStrain = lis[1]
        self.xyStrain= lis[2]
        self.name = lis[3]
        self.idNum = lis[4]
        
        if (self.xStrain[0].strainGroup[0].skewed or self.yStrain[0].strainGroup[0].skewed):
            self.skewed = True
            
        if (self.xStrain[0].strainGroup[0].afterMarch or self.yStrain[0].strainGroup[0].afterMarch):
            self.afterMarch = True
        
    def InitializeG11(self):
        ##Zig-Zag (X) Direction:
        vec = []
        energies = []
        areas = []
        for l in self.xStrain[0].strainGroup:
            vec.append(l.xV)
            energies.append(l.energy)
            areas.append(l.area)    
        #Fits, minimums, etc
        _, _, vecMin, _, _ = GetPolyEqtn(vec, energies, 2)
        strains_ = [(v / vecMin) - 1 for v in vec]
        xP, yP, strainForeMin, eMin, fit, g11r2, g11StdDev = GetPolyEqtn(strains_, energies, 2, findR2 = True, findStdErr = True)
        aP, eP, areaForeMin, eMin_, fit_ = GetPolyEqtn(areas, energies, 2)
        ##Return: g11, the rSquared for g11, std err for g11    
        return (2 * fit[0] / areaForeMin), (g11r2), (2 * g11StdDev / areaForeMin)   

    def InitializeG22(self):
        ##Armchair (Y) Direction:
        vec = []
        energies = []
        areas = []
        for l in self.yStrain[0].strainGroup:
            vec.append(l.yV)
            energies.append(l.energy)
            areas.append(l.area)    
        #Fits, minimums, etc
        _, _, vecMin, _, _ = GetPolyEqtn(vec, energies, 2)
        strains_ = [(v / vecMin) - 1 for v in vec]
        xP, yP, strainForeMin, eMin, fit, g22r2, g22StdDev = GetPolyEqtn(strains_, energies, 2, findR2 = True, findStdErr = True)
        aP, eP, areaForeMin, eMin_, fit_= GetPolyEqtn(areas, energies, 2)
        ##Return: g22, rSquared for g22, std err for g22    
        return (2 * fit[0] / areaForeMin), (g22r2), (2 * g22StdDev / areaForeMin)
    
    def InitializeG12(self, g11, g22, uncG11, uncG22):        
        ##Isotropic (XY) Directions:
        vec = []
        energies = []
        areas = []
        for l in self.xyStrain[0].strainGroup:
            vec.append(l.xV)
            energies.append(l.energy)
            areas.append(l.area)    
        #Fits, minimums, etc
        _, _, vecMin, _, _ = GetPolyEqtn(vec, energies, 2)
        strains_ = [(v / vecMin) - 1 for v in vec]
        xP, yP, strainForeMin, eMin, fit, g12r2, g12StdDev = GetPolyEqtn(strains_, energies, 2, findR2 = True, findStdErr = True)
        aP, eP, areaForeMin, eMin_, fit_ = GetPolyEqtn(areas, energies, 2)
        ##Return: g12, rSquared for g12, std err for g12   
        ###          sqrt(((partial g12 partial fit)^2 err fit^2) +  (partial g12 partial g11)^2 err g11^2...etc)
        dg12dfit = 1/(2*areaForeMin)
        dg12dg11 = -1/2
        dg12dg22 = -1/2
        g12PropErr = ((dg12dfit * g12StdDev)**2 + (dg12dg11 * uncG11)**2 + (dg12dg22 * uncG22)**2)**(1/2)
        return ((fit[0] / areaForeMin) - (1/2 * (g11 + g22))), (g12r2), g12PropErr

# <codecell> INITIALIZATION    

def GetDataLines(pathToFile_):
    ##Read CSV, put all lines into list
    lines = []
    with open(pathToFile_, 'r') as infile:
        reader = csv.reader(infile, delimiter=",")
        for i, line in enumerate(reader):
            lines.append(list(line))
    
    i = 1
    dataLines = []
    while (i < len(lines) - 1):
        dataLines.append(dataLine(lines[i]))
        i = i + 1
    dataLines.sort(key=op.attrgetter('name'))
    dataLines.append(dataLine([1,2,3,4,'END',5,6,7,8,9,10,11,12]))
  
    return dataLines

def GetStrainGroups(pathToFile_ = None, dataLines_ = None):
    if ((pathToFile_ == None and dataLines_ == None) or (pathToFile_ == None and dataLines_ == None)):
        print("GetStrainGroups:  Bad input")
        return -1
    if(pathToFile_ != None):    
        dataLines = GetDataLines(pathToFile_)
    if (dataLines_ != None):
        dataLines = dataLines_
    ##Go through lines and organize them by strain type.  Put in list of strain types
    strainGroups = []
    lis = []
    i = 1
    while (i < len(dataLines)):
        if (dataLines[i].name == dataLines[i - 1].name):
            lis.append(dataLines[i - 1])
        else:
            lis.append(dataLines[i - 1])
            strainGroups.append(strainGroup([lis, lis[0].name, lis[0].idNum, GetStrainType(lis[0].name)]))
            lis = []
        i = i + 1
    strainGroups.sort(key=op.attrgetter('idNum'))
    strainGroups.append(strainGroup([1,'END',2,3]))
    for a in strainGroups[:-1]:
        a.strainGroup.sort(key=op.attrgetter('disp'))    
    
    return strainGroups
    
def GetDefectGroups(pathToFile__ = None, dataLines__ = None, strainGroups__ = None):
    if(strainGroups__ != None):
        strainGroups = strainGroups__
    elif(dataLines__ != None):
        strainGroups = GetStrainGroups(dataLines_ = dataLines__)
    else:
        strainGroups = GetStrainGroups(pathToFile_ = pathToFile__)
    #Go through the strain types and organize them by defect type.  Put in list of defect types
    defectGroups = []
    lisX, lisY, lisXY = [], [], []
    i = 1
    while (i < len(strainGroups)):
        if (strainGroups[i].idNum == strainGroups[i - 1].idNum):
            if(strainGroups[i - 1].strainType == 'X'):
                lisX.append(strainGroups[i - 1])
            if(strainGroups[i - 1].strainType == 'Y'):
                lisY.append(strainGroups[i - 1])
            if(strainGroups[i - 1].strainType == 'XY'):
                lisXY.append(strainGroups[i - 1])
        else:
            if(strainGroups[i - 1].strainType == 'X'):
                lisX.append(strainGroups[i - 1])
            if(strainGroups[i - 1].strainType == 'Y'):
                lisY.append(strainGroups[i - 1])
            if(strainGroups[i - 1].strainType == 'XY'):
                lisXY.append(strainGroups[i - 1])
            defectGroups.append(defectGroup([lisX, lisY, lisXY, lisX[0].name, lisX[0].idNum]))
            lisX, lisY, lisXY = [], [], []
        i = i + 1
    defectGroups.sort(key=op.attrgetter('idNum'))
       
    return defectGroups  


 
#returns:list of data lines, list of strains, list of defects
def GetOutputInfo(pathToFile):
    dataLines = GetDataLines(pathToFile)
    strainGroups = GetStrainGroups(dataLines_ = dataLines)
    defectGroups = GetDefectGroups(strainGroups__ = strainGroups)
    
    dataLines.pop()
    strainGroups.pop()  

    return dataLines, strainGroups, defectGroups

def GetPureInfo(defect):
    pureEnergies = [[], [], []]
    for d in defect:
        if (GetDisplayName(d.name, d.afterMarch) == 'Pure'):
            nCarbonsPure = d.xStrain[0].strainGroup[0].cAtoms
            for s1, s2, s3 in zip(d.xStrain, d.yStrain, d.xyStrain):
                for l1, l2, l3 in zip(s1.strainGroup, s2.strainGroup, s3.strainGroup):
                    pureEnergies[0].append(l1.energy)
                    pureEnergies[1].append(l2.energy)
                    pureEnergies[2].append(l3.energy)
    pureEnergies_ = []
    for e in pureEnergies:
        for e2 in e:
            pureEnergies_.append(e2)
    return min(pureEnergies_), nCarbonsPure

def SplitVac(pathToFile, pathToOutFile, xLims = [VAC_X_LEFT_END, VAC_X_RIGHT_BEG], yLims = [VAC_Y_LEFT_END, VAC_Y_RIGHT_BEG], xyLims = [VAC_XY_LEFT_END, VAC_XY_RIGHT_BEG]):
    writer = csv.writer(open(pathToOutFile, 'w', newline=''))
    with open(pathToFile, 'r') as infile:
        reader = csv.reader(infile, delimiter=",")
        for i, line in enumerate(reader):
            #Case for X:
            if(line[4] == '1DVacBulkX'):
                end = line[4][-1]
                if (float(line[2]) <= float(xLims[0])):                    
                    line[4] = line[4][:9] + 'Left' + end
                if (float(line[2]) >= float(xLims[1])):
                    line[4] = line[4][:9] + 'Right' + end     
                    line[5] = str(int(line[5]) + 99)
            #Case for Y:
            if(line[4] == '1DVacBulkY'):
                end = line[4][-1]
                if (float(line[2]) <= float(yLims[0])):                    
                    line[4] = line[4][:9] + 'Left' + end
                if (float(line[2]) >= float(yLims[1])):
                    line[4] = line[4][:9] + 'Right' + end     
                    line[5] = str(int(line[5]) + 99)
            #Case for XY:
            if(line[4] == '2DVacBulk'):
                if (float(line[2]) <= float(xyLims[0])):                    
                    line[4] = line[4] + 'Left'
                if (float(line[2]) >= float(xyLims[1])):
                    line[4] = line[4] + 'Right'        
                    line[5] = str(int(line[5]) + 99)
            writer.writerow(line)
            
def IsoVac(pathToFile, pathToOutFile):
    writer = csv.writer(open(pathToOutFile, 'w', newline=''))
    with open(pathToFile, 'r') as infile:
        reader = csv.reader(infile, delimiter=",")
        for i, line in enumerate(reader):    
            if(line[4][2:9] == 'VacBulk' or i == 0):
                writer.writerow(line)

# <codecell> GRAPH FUNCTIONS
         
def EnergyVsStrainGraph(strain, xDat, yDat, xP, yP, eMin, cAtoms, LATpos, outLoc):
    plt.title(GetDisplayName(strain.name, strain.afterMarch))
    plt.xlabel(GetAxisName(strain.name))
    plt.ylabel('Energy from Min. per C Atom (meV)')
    plt.plot(xP, [float(1000 * (e - eMin)) / float(cAtoms) for e in yP], 'c')
    plt.plot(xDat, [float(1000 * (e - eMin)) / float(cAtoms) for e in yDat], 'ko')
    plt.axvline(LATpos, color = '#ff00ff', linewidth = 0.75)
    plt.savefig(outLoc, dpi = 300)
    plt.clf()    
            
def VacGraphs(pathToFile, pathToOut, outType, xLims = [VAC_X_LEFT_END, VAC_X_RIGHT_BEG], yLims = [VAC_Y_LEFT_END, VAC_Y_RIGHT_BEG], xyLims = [VAC_XY_LEFT_END, VAC_XY_RIGHT_BEG]):
    lines, strains, defects = GetOutputInfo(pathToFile)
    
    def GetDisplayName_(ret, afterMarch = False):
        if(re.search("VacBulk", ret)):
            return "SV"
 
    for s in strains:
        #Get displacements and their corresponding energies
        displacements = []
        energies = []
        for l in s.strainGroup:
            nCarbons = l.cAtoms
            displacements.append(l.disp)
            energies.append(l.energy)
        
        #Get initial polynomial information
        xP, yP, zeroStrainDisp, eMin, garbage = GetPolyEqtn(displacements, energies, 2)

        #Find place to split vac up
        if (s.name == '1DVacBulkX'):
            l, r = xLims[0] + 1, xLims[1]
        if (s.name == '1DVacBulkY'):
            l, r = yLims[0] + 1, yLims[1]  
        if (s.name == '2DVacBulk'):
            l, r = xyLims[0] + 1, xyLims[1]
                
        #Define 0 strain to be at the minimum of the energy-displacement curve.  Get energy vs strain
        strains_ = [d - zeroStrainDisp for d in displacements]
        xP, yP, strainForeMin, eMin, fit = GetPolyEqtn(strains_, energies, 2)
        xPL, yPL, strainForeMinL, eMinL, fit = GetPolyEqtn(strains_[:l], energies[:l], 2)
        xPR, yPR, strainForeMinR, eMinR, fit = GetPolyEqtn(strains_[r:], energies[r:], 2)
    
        #Get the position of LAT for true strain
        strainLAT = ((1) - 1) - zeroStrainDisp
        
        #Plot the things.  Save the things
        plt.title(GetDisplayName_(s.name, s.afterMarch))
        plt.xlabel(GetAxisName(s.name))
        plt.ylabel('Energy from Min. per C Atom (meV)')
        plt.plot(xPL, [float(1000 * (e - eMin)) / float(nCarbons) for e in yPL], 'c')
        plt.plot(xPR, [float(1000 * (e - eMin)) / float(nCarbons) for e in yPR], 'c')
        plt.plot(strains_, [float(1000 * (e - eMin)) / float(nCarbons) for e in energies], 'ko')
        plt.axvline(strainLAT, color = '#ff00ff', linewidth = 0.75)
        plt.savefig(pathToOut + str(s.name) + 'Fin' + str(outType), dpi = 300)
        plt.clf()
        
     
def LatConstGraph(pathToIn, pathToOut, outType, defectGroup = None, lims = [2.4, 2.5]):
    if(defectGroup == None):
        lines, strains, defects = GetOutputInfo(pathToIn)
    else:
        defects = defectGroup
    
    #Plot time
    plt.ylim(lims[0], lims[1])
    for d in defects:
        plt.plot(GetDisplayName(d.xStrain[0].name, d.xStrain[0].afterMarch), d.xStrain[0].minEnergyLatticeConstantPoly, 'ko', linewidth=0)##x strain, lat const from x vector
        plt.plot(GetDisplayName(d.yStrain[0].name, d.yStrain[0].afterMarch), d.yStrain[0].minEnergyLatticeConstantPoly, 'c^', linewidth=0) ##y strain, lat const from y vector
    
    plt.xlabel('Defect Type')
    plt.ylabel(r'Lattice Constant ($\mathring{\mathrm{A}})$')
    plt.savefig(str(pathToOut) + str(outType), dpi=300)
    plt.clf()
        
def PolyGridGraph(pathToIn, pathToOut, outType, strainSetA, strainSetB, LAT, strainGroup = None, xLims = [VAC_X_LEFT_END, VAC_X_RIGHT_BEG], yLims = [VAC_Y_LEFT_END, VAC_Y_RIGHT_BEG], xyLims = [VAC_XY_LEFT_END, VAC_XY_RIGHT_BEG]):
    if(strainSetA == strainSetB):
        return
    if(strainGroup == None):
        lines, strains, defects = GetOutputInfo(pathToIn)
    else:
        strains = strainGroup
        
    fig = plt.figure(figsize=(8.5, 11))
    gs = gridspec.GridSpec(7, 3, fig, hspace=0.1)
    
    def GetRowCol(strain):
        row, col = -1, -1
        name = strain.name #Pure, SWS, etc.
        direction = strain.strainType #X, Y, or XY
        am = strain.afterMarch
        
        #Determine the correct row
        if(GetDisplayName(name, am) == "Pure"):
            row = 0
        elif(GetDisplayName(name, am) == "SV" or GetDisplayName(name, am) == "Vac Bulk Right" or GetDisplayName(name, am) == "Vac Bulk Left"):
            row = 1
        elif(GetDisplayName(name, am) == r"DV$^{a}$"):
            row = 2
        elif(GetDisplayName(name, am) == r"DV$^{b}$"):
            row = 3
        elif(GetDisplayName(name, am) == "555777"):
            row = 4
        elif(GetDisplayName(name, am) == r"STW$^{a}$"):
            row = 5
        elif(GetDisplayName(name, am) == r"STW$^{b}$"):
            row = 6

        
        #Determine the correct column
        if(direction == 'X'):
            col = 0
        if(direction == 'Y'):
            col = 1            
        if(direction == 'XY'):
            col = 2      
            
        return row, col

    for s in strains:
        #Get displacements and their corresponding energies
        vec = []
        energies = []
        for l in s.strainGroup:
            if (s.strainType == 'X'):
                vec.append(l.xV)
            if (s.strainType == 'Y'):
                vec.append(l.yV)
            if(s.strainType == 'XY'):
                vec.append(l.xV)
            energies.append(l.energy)
        
        #Get initial polynomial information
        _, _, vecMin, _, _ = GetPolyEqtn(vec, energies, 2)
        #Define 0 strain to be at the minimum of the energy-displacement curve.  Get energy vs strain
        strains_ = [(v / vecMin) - 1 for v in vec]
        xP, yP, strainForeMin, eMin, fit = GetPolyEqtn(strains_, energies, 2)
    
        #Get the position of LAT for true strain
        if (strainSetA):
            strainSet_ = 'A'
        else:
            strainSet_ = 'B'
        strainLAT, initLatConst = FindLATPosition(strain = s, strainSet = strainSet_, afterMarch = s.afterMarch, LAT = LAT)
    
        #Add individual plots to the big boy
        row, col = GetRowCol(s)
        ax = fig.add_subplot(gs[row, col])
        ax.set_facecolor((0,0,0,0))
        
        if(col == 2):
            ax.set_title(GetDisplayName(s.name, s.afterMarch), rotation = 90, x=1.1, y=0.4)       
        if(row == 6):
            ax.set_xlabel(GetAxisName(s.name))
        else:                    
            ax.tick_params('x', labelbottom=False)
        if(row == 3 and col == 0):
            ax.set_ylabel('Energy from Min. per C Atom (meV)', labelpad=10)
           
        if(strainSetA):
            ax.set_xlim(-0.03, 0.03)
        if(strainSetB):
            ax.set_xlim(-0.025, 0.025)
        ##special case for single vacancy, set A    
        if(GetDisplayName(s.name, s.afterMarch) == "SV" and strainSetA): 
            if (s.strainType == 'X'):
                lef, rig = xLims[0] + 1, xLims[1]
            if (s.strainType == 'Y'):
                lef, rig = yLims[0] + 1, yLims[1]  
            if(s.strainType == 'XY'):
                lef, rig = xyLims[0] + 1, xyLims[1]
                
            #Define 0 strain to be at the minimum of the energy-displacement curve.  Get energy vs strain
            strains_ = [d - zeroStrainDisp for d in displacements]
            xP, yP, strainForeMin, eMin, fit = GetPolyEqtn(strains_, energies, 2)
            xPL, yPL, strainForeMinL, eMinL, fit = GetPolyEqtn(strains_[:lef], energies[:lef], 2)
            xPR, yPR, strainForeMinR, eMinR, fit = GetPolyEqtn(strains_[rig:], energies[rig:], 2)
    
            #Get the position of LAT for true strain
            strainLAT = ((1) - 1) - zeroStrainDisp
            ax.plot(xPL, [float(1000 * (e - eMin)) / float(l.cAtoms) for e in yPL], 'c')
            ax.plot(xPR, [float(1000 * (e - eMin)) / float(l.cAtoms) for e in yPR], 'c')
            ax.plot(strains_, [float(1000 * (e - eMin)) / float(l.cAtoms) for e in energies], 'ko')
            ax.axvline(strainLAT, color = '#ff00ff', linewidth = 0.75)
        
        ##Everything but the single vacancy               
        else:
            ax.plot(xP, [float(1000 * (e - eMin)) / float(l.cAtoms) for e in yP], 'c')
            ax.plot(strains_, [float(1000 * (e - eMin)) / float(l.cAtoms) for e in energies], 'ko')
            ax.axvline(strainLAT, color = '#ff00ff', linewidth = 0.75)
    
    fig.savefig(str(pathToOut) + str(outType), dpi=300)
    plt.clf()
    
    
    