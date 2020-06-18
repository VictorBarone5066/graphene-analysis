Python files specifically for use with output from shell scripts.  

'output analyzier' reads a csv file and generates elastic data and graphs for graphene defects.  Supports 'A' and 'B' type files - 'A' is for non-defect lattice constants and -B- is for volume relaxed lattice constants.  Most of the main code depends on the header file
INPUT:	output.csv, from shell scripts.
OUTPUT:	FormationEnergies.csv, GammaValues.csv, LatticeConstants.csv, figs 

Temp graphs are for analysis of MD data

'Mag Graph.py' and 'Strain Map Graphs.py' simply create extra graphs - for magnetization and strain amounts, specifically.  Mag graph reads in an edited version of 'output.csv', which can be generated from the corresponding shell script.  Strain map graphs uses the same input file as the output analyzers.  

'nElect.py' counts the number of electrons in a system using DOSCAR - can be used for consistency checks, or calculating magnetization, as long as DOSCAR reports spin-polarized values.  
