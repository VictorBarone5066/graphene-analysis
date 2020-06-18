import outputAnalyzerHeaderChanged as header

import matplotlib.pyplot as plt
import sys

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True

#Constants
##Set A:  non-volume relaxed   ...   Set B:  volume relaxed
strainSetA, strainSetB = True, False
FILE_PATH = "C:\\Users\\baron\\Desktop\\"
outType = ".png"

#X index limits for set A...left, right
xLims = [5, 10]
yLims = [4, 5]
xyLims = [2, 3]

LAT = 2.468217817 #starting lattice constant (angstroms) (for set A)

if (strainSetA == strainSetB):
    print("Choose EITHER strain set A or set B - not both.\n")
    sys.exit()
if (strainSetA):
    strainSet = 'A'
if (strainSetB):
    strainSet = 'B'

FILE_PATH = FILE_PATH + "Strain Set " + str(strainSet) + "\\"

#Need to split vacancy into two seperate groups for strain set A          xLims, yLims, xyLims
if (strainSetA):
    header.SplitVac(FILE_PATH + 'output.csv', FILE_PATH + 'outputSplit.csv', [5, 5.1], [4, 5], [2, 3])
    header.IsoVac(FILE_PATH + 'output.csv', FILE_PATH + 'vacOnly.csv')
    header.VacGraphs(FILE_PATH + 'vacOnly.csv', FILE_PATH, outType, xLims, yLims, xyLims)

#Open csvs for writing
latConstWriter = open(FILE_PATH + 'LatticeConstantsSet' + strainSet + '.csv', "w+")
formationEnergyWriter = open(FILE_PATH + 'FormationEnergiesSet' + strainSet + '.csv', "w+")
gammaValueWriter = open(FILE_PATH + 'GammaValueSet' + strainSet + '.csv', "w+")

latConstWriter.write("Model Name,Min E Lat Const (Polynomial),Min E Lat Const (from File),|Polynomial - File|,Init Lat Const,|Polynomial - Init|,|File - Init|,All in Angstroms")
formationEnergyWriter.write("Model Name,Formation Energy (eV)")
gammaValueWriter.write("Model Name,g11 (ev/A^2),g22(eV/A^2),g12(ev/A^2),y11(Pa),y22(Pa),p11,p22,y11/yp11,y22/yp22,p11/pp11,p22/pp22,rSqrdg11,rSqrdg22,rSqrdg12,stdErrg11,stdErrg22,stdErrg12")

# <codecell> MAIN
if (strainSetA):
    lines, strains, defects = header.GetOutputInfo(FILE_PATH + 'outputSplit.csv')
else:
    lines, strains, defects = header.GetOutputInfo(FILE_PATH + 'output.csv')

#Initialization:  Need to find chemical potential of carbon for the formation energies
##Corresponds to fitting energy / atom for a range of areas for pure graphene
pureETot, nCarbonsPure = header.GetPureInfo(defects)

#Individual Strain Loops:  
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
    _, _, vecMin, _, _ = header.GetPolyEqtn(vec, energies, 2)
    #Define 0 strain to be at the minimum of the energy-displacement curve.  Get energy vs strain
    strains_ = [(v / vecMin) - 1 for v in vec]
    xP, yP, strainForeMin, eMin, fit = header.GetPolyEqtn(strains_, energies, 2)
    print(s.name, vecMin)
    if(s.name == "1DVacBulkLeftX"):
        print(strains_)
    #Get the position of LAT for true strain
    strainLAT, initLatConst = header.FindLATPosition(strain = s, strainSet = strainSet, afterMarch = s.afterMarch, LAT = LAT)

    #Plot individual energy v strains
    header.EnergyVsStrainGraph(s, strains_, energies, xP, yP, eMin, l.cAtoms, strainLAT, str(FILE_PATH) + str(s.name) + str(outType))

    #Print the Minimum energy lattice constants
    s.InitializeMinEnergyLatticeConstant()
    file, poly = s.minEnergyLatticeConstantFile, s.minEnergyLatticeConstantPoly
    latConstWriter.write('\n' + str(s.name) + ',' + str(poly) + ',' + str(file) + ',' + str(abs(file - poly)) + ',' + str(initLatConst) + ',' + str(abs(initLatConst - poly)) + ',' + str(abs(initLatConst - file)))
    
    #Formation Energies:  Ef = Ed - N*mu; Ed = eTot of defected, N = num atoms in defected, mu = chem potential (energy / atom for pure model of the same size as the defected one)
    s.InitializeFormationEnergy(pureETot, strainSet, pureNAtoms = nCarbonsPure)
    formationEnergyWriter.write('\n' + str(s.name) + ',' + str(s.formationEnergy))

latConstWriter.close()
formationEnergyWriter.close()

#Gamma Value, Youngs Moduli, Poisson Ratio Loop:
for d in defects:

    #Get all gamma values
    g11, g11r2, g11StdDev = d.InitializeG11()
    g22, g22r2, g22StdDev = d.InitializeG22()
    g12, g12r2, g12StdDev = d.InitializeG12(g11 = g11, g22 = g22, uncG11 = g11StdDev, uncG22 = g22StdDev)

    #Youngs Moduli:  Y = gamma (eV to J) / (Angstroms to Meters * Interplanar spacing for graphene)
    young11 = g11 * (1.6022E-19) / ((1E-10)*(1E-10)*(3.4E-10))
    young22 = g22 * (1.6022E-19) / ((1E-10)*(1E-10)*(3.4E-10))
    #Poissons Ratio: ratio of gamma values
    poisson11 = g12 / g22
    poisson22 = g12 / g11

    #Ratios
    if(d.name[2:6] == 'Pure'):
        young11Pure = young11
        young22Pure = young22
        poisson11Pure = poisson11
        poisson22Pure = poisson22

    #Write Values:
    gammaValueWriter.write('\n'+ header.GetDisplayName(d.name, d.afterMarch) + ',' + str(g11) + ',' + str(g22) + ',' + str(g12) + ',' + str(young11) + ',' + str(young22) + ',' + str(poisson11) + ',' + str(poisson22) + ',' + str(young11/young11Pure) + ',' + str(young22/young22Pure) + ',' + str(poisson11/poisson11Pure) + ',' + str(poisson22/poisson22Pure) + ',' + str(g11r2) + ',' + str(g22r2) + ',' + str(g12r2)+ ',' + str(g11StdDev) + ',' + str(g22StdDev) + ',' + str(g12StdDev))
         
    
gammaValueWriter.close()
header.LatConstGraph(FILE_PATH + 'output.csv', FILE_PATH + 'LatConstGraph' + strainSet, outType, defectGroup = defects, lims = [2.4, 2.50])
header.PolyGridGraph(FILE_PATH + 'output.csv', FILE_PATH + 'PolyGrid' + strainSet, outType, strainSetA, strainSetB, LAT, xLims = xLims, yLims = yLims, xyLims = xyLims)

