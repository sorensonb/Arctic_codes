#!/usr/bin/env python

"""
NAME:

PURPOSE:

MODIFICATIONS:
  Blake Sorenson  - 170628: Written 

"""

# Import functions
import numpy as np
import sys
import matplotlib.pyplot as plt
from readsounding import *
from compare_cn2 import *

#if(len(sys.argv)!=2):
#    print "SYNTAX: ./static_stability.py Graw_RTS_table"
#    sys.exit()

graw = readsounding(sys.argv[1])


#graw = smooth(graw)
graw = cn2_calc(graw)
d0dz = np.zeros(len(graw['CN2']))
Rb = np.zeros(len(graw['CN2']))


# Calculate theta
#thetaO = graw['PotTe']+273.15
tempK = graw['TEMP']+273.15
thetaCalc = tempK*((1000./graw['PRESS'])**0.286)

counter=0
for i in range(1,len(thetaCalc)-1):
    d0dz[counter] = (thetaCalc[i+1]-thetaCalc[i-1])/(graw['ALT'][i+1]-graw['ALT'][i-1]) 
    counter+=1

# Calculate Bulk Richardson Number
#alt_km = (graw['ALT']-graw['ALT'][0])/1000.
#theta = graw['TEMP']+9.8*alt_km
#tK = graw['TEMP']+273.15
for i in range(0,len(graw['TEMP'])-1):
    Rb[i] = ((9.8/(np.average([tempK[i+1],tempK[i]])))*(thetaCalc[i+1]-thetaCalc[i])*(graw['ALT'][i+1]-graw['ALT'][i]))/\
            ((graw['UWIND'][i+1]-graw['UWIND'][i])**2. + (graw['VWIND'][i+1]-graw['VWIND'][i])**2.) 
logCn2 = np.log10(graw['CN2'])

# Plot scatter plot of d0dz vs Cn2
fig1, ax = plt.subplots()
fit = np.polyfit(d0dz,logCn2,deg=1)
plt.scatter(d0dz,logCn2)

plt.ylabel('log$_{10}$ $C_{n}^{2}$')
plt.xlabel('d$\Theta$/dz')

# Plot time series of d0dz and  Cn2
fig2 = plt.figure()
ax1 = plt.gca() 
ax1.plot(np.arange(0,len(logCn2)),logCn2,color='blue')
ax1.spines['left'].set_color('blue')
ax1.set_ylabel('Cn2',color='blue')
ax2 = ax1.twinx()
ax2.plot(np.arange(0,len(d0dz)),d0dz,color='black')
ax2.spines['left'].set_color('black')
ax2.set_ylabel('d0dz',color='black')
plt.show()
sys.exit()

# Get rid of infinite values in Rb before proceeding
logCn2 = np.delete(logCn2,np.where(np.isinf(Rb)))
Rb     = np.delete(Rb,np.where(np.isinf(Rb)))
# Plot scatter plot of Rb vs Cn2
fig3, ax3 = plt.subplots()
fit = np.polyfit(Rb,logCn2,deg=1)
plt.scatter(Rb,logCn2)

# Plot r2 line
ax3.plot(Rb,fit[0]*Rb+fit[1],color='red')
plt.ylabel('log$_{10}$ $C_{n}^{2}$')
plt.xlabel('Bulk Richardson Number')

# Plot time series of d0dz and  Cn2
fig4 = plt.figure()
ax4 = plt.gca() 
ax4.plot(np.arange(0,len(logCn2)),logCn2,color='blue')
ax4.spines['left'].set_color('blue')
ax4.set_ylabel('Cn2',color='blue')
ax5 = ax4.twinx()
ax5.plot(np.arange(0,len(Rb)),Rb,color='black')
ax5.spines['left'].set_color('black')
ax5.set_ylabel('Bulk Richardson Number',color='black')

plt.show()
