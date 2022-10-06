#!/usr/bin/env python

"""
NAME:
  ascent_plotter.py

PURPOSE:
  Calculate the refractive index structure parameter for a maximum of 

MODIFICATIONS:
  Blake Sorenson  - 170628: Written 

"""

# Import functions
import numpy as np
import sys
import matplotlib.pyplot as plt
from adpaa import ADPAA

if(len(sys.argv)!=2):
    print("SYNTAX: ./thermo_test.py thermo_file")
    sys.exit()

thermo = ADPAA()
thermo.ReadFile(sys.argv[1])
# Mask any missing data in the thermosonde data
thermo.mask_MVC()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(thermo.data['Time'], thermo.data['Alt'])
plt.show()


#start_time = 86332.0
#stop_time  = 86606.0
#
#ascent_start = np.where(thermo.data['Time']==start_time)[0][0]
#ascent_stop  = np.where(thermo.data['Time']==stop_time)[0][0]
#
#thermo.data['Alt']-=thermo.data['Alt'][ascent_start]
#
#fig1 = plt.figure()
#plt.xlabel('Time [UTC]')
#ax1 = plt.gca()
#ax1.plot(thermo.data['Time'][ascent_start:ascent_stop+1],thermo.data['Alt'][ascent_start:ascent_stop+1],color='black')
#ax1.set_ylabel('Altitude AGL [m]',color='black')
#
#ax2 = ax1.twinx()
#ax2.set_ylim(0.001,0.035)
#ax2.plot(thermo.data['Time'][ascent_start:ascent_stop+1],thermo.data['TempDiff'][ascent_start:ascent_stop+1],color='blue')
#ax2.set_ylabel('Temperature Difference [degC]',color='blue')
#ax2.spines['right'].set_color('blue')
#
##plt.title('Tethered Test Ascent')
##fig1.savefig('17_11_17_22_51_30_Graw_ascent_altTemp.png',dpi=300)
#plt.show()
#
