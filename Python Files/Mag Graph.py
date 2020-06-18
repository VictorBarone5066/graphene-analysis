# Magnetization graph

#Imports:
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
    
FILE_PATH = 'C:\\Users\\baron\\Desktop\\School\\Research\\Fa 2018 - Spring 2020\\Vac-Defect Specific\\Magnetization\\'

def polyEval(x, coeffs):
    order = len(coeffs) - 1
    
    i = 0
    s = 0
    while(order >= 0):
        s = s + coeffs[i]*(x**order)
        i = i + 1
        order = order - 1
  
    return s

def poly(lisX, lisY, degEdit = 0):
    x, y = [], []
    fit = np.polyfit(lisX, lisY, len(lisX) + degEdit)

    i = lisX[0]
    while(i <= lisX[-1] + ((lisX[-1]-lisX[0])/100)):
        x.append(i)
        y.append(polyEval(x[-1], fit))
        i = i + (lisX[-1]-lisX[0])/100
    
    return x, y
    

#Set up:  
df = pd.read_csv(FILE_PATH + 'outputM.csv')
df.sort_values(by = ['displacement from original positions'], inplace = True)

#Getting the columns, saving to lists
mag = df['mag']
displacement = df['xV'] #displacement from starting 2.468217 angstroms
names = df['name for determining information']
energy = df['energy']

#Declerations
left, right = 6, 6


#Plots for strain vs mag
x = list(displacement)[:]
x = [(v/14.686191521410283) - 1 for v in x]
e = list(energy)[:]
y = list(mag)[:]

#Plots for dM/de vs strain
i = 1
dx = []
dy = []
while(i < len(x)-1): ##previously while(i < len(x)):
    if (True): ##do forward diff ##previously if(i == 0):
        dx.append(x[i])
        dy.append((y[i+1] - y[i]) / (x[i+1] - x[i]))
#    elif (i == len(x) - 1): ##do backwards diff
#        dy.append((y[i] - y[i-1]) / (x[i] - x[i-1]))
#    else: #Average forwards and backwards diff
#        forw = (y[i+1] - y[i]) / (x[i+1] - x[i])
#        back = (y[i] - y[i-1]) / (x[i] - x[i-1])
#        dy.append((forw + back) / 2)
    i = i + 1          
      
# <codecell> Mag Graphs
#Plot figs
    
##Full Image------------------------------------------------------------------
#plt.title('Full')
plt.xlabel(r"Zig-Zag Strain")
plt.ylabel(r'Net Magnetic Moment ($\mu_{\mathrm{B}}$)')
plt.ylim(1.5, 2)

leftX, leftY = poly(x[:left], y[:left])
rightX, rightY = poly(x[right:], y[right:], -5)
plt.plot(leftX, leftY, color = 'c')
plt.plot(rightX, rightY, color = 'c')
plt.plot(x, y, 'ko')

plot = plt.gcf()
plot.savefig(FILE_PATH + 'mag.png', dpi = 300)            
plt.show()
plt.clf()

#Derivitive plots----------------------------------------------------------
#plt.title('Derivitive (forward difference)')
plt.xlabel(r"Zig-Zag Strain")
plt.ylabel(r'$\partial M/\partial\varepsilon$ ($\mu_{\mathrm{B}}$)')

leftX, leftY = poly(dx[:left-1], dy[:left-1])
rightX, rightY = poly(dx[right-1:], dy[right-1:], -5)
plt.plot(leftX, leftY, color = 'c')
plt.plot(rightX, rightY, color = 'c')
plt.plot(dx, dy, 'ko')
plot = plt.gcf()
plot.savefig(FILE_PATH + 'magDerivitive.png', dpi = 300) 
plt.show()
plt.clf()

##Zoomed image--------------------------------------------------------------
centX, centY = 0.0087, 1.856
offX, offY = 0.006, 0.07
#plt.title('Zoomed')
leftX, leftY = poly(x[:left], y[:left])
rightX, rightY = poly(x[right:], y[right:], -2)
plt.xlabel(r"Zig-Zag Strain")
plt.ylabel(r'Net Magnetic Moment ($\mu_{\mathrm{B}}$)')
plt.plot(leftX, leftY, color = 'c')
plt.plot(rightX, rightY, color = 'c')
plt.plot(x, y, color='k', marker='o', linewidth=0) #plot data points
plt.xlim(centX - offX, centX + offX)
plt.ylim(centY - offY, centY + offY)

plot = plt.gcf()
plot.savefig(FILE_PATH + 'magZoom.png', dpi = 300)            
plt.show()
plt.clf()

            

#plt.axvline(x=(-0.00125 / 2), ymin = 0.02, ymax = 0.98, color = 'g', linestyle = '--', linewidth = 0.75)
#plt.text(x = -0.022, y = 1.8, s = (r'$\Delta M \approx {:.1f} {}$').format(float(y[7]) - float(y[8]), r'\mu_{\mathrm{B}}'))
#plt.axvline(x = -0.02, ymin = .53, ymax = .57, color ='g', linestyle = '--', linewidth = 0.75)
#plt.axvline(x = -0.02, ymin = .665, ymax = .705, color ='g', linestyle = '--', linewidth = 0.75)
#plt.axhline(y=y[7], xmin = 0.10, xmax = .465, color ='g', linestyle = '--', linewidth = 0.75)
#plt.axhline(y=y[8], xmin = 0.10, xmax = .485, color ='g', linestyle = '--', linewidth = 0.75)








