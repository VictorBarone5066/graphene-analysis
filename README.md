# graphene-analysis
Python files for calculations related to elastic and magnetic properties of graphene.  To be used with output files created from shell scripts available on github.

Python files specifically for use with output from shell scripts.  

'output analyzier' (2 and 3) read data and generate a bunch of information on elastic properties of graphene, and generate graphs of the data points.  2 is different than 3 in the following ways:
	-2 assumes that the initial starting supercell size is 14.76 x 12.78 angstroms.  3 allows for variable supercell sizes, but requires extra data in the read in file to compute cell size.  This means that one would have to use the later versions of 'writeFromWork.sh' to generate the necessary columns.
	-The vertical line on 2's graphs represent a calculated zero strain point, where the the vertical lines on 3's graphs represent 2.46 angstrom (zig-zag) lattice constant.
INPUT:	output.csv, from shell scripts.  poscars/, a folder of poscars, also from shell scripts.
OUTPUT:	FormationEnergies.csv, GammaValues.csv, LatticeConstants.csv, figs/.  Elastic moduli can be found in the gamma values file, and graphs can be found in the figs/ folder.  


'Mag Graph.py' and 'Strain Map Graphs.py' simply create extra graphs - for magnetization and strain amounts, specifically.  Mag graph reads in an edited version of 'output.csv', which can be generated from the corresponding shell script.  Strain map graphs uses the same input file as the output analyzers.  

'nElect.py' counts the number of electrons in a system using DOSCAR - can be used for consistency checks, or calculating magnetization, as long as DOSCAR reports spin-polarized values.  
