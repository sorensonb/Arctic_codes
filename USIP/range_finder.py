#!/usr/bin/env python

"""
NAME:
  plot_cn2.py

PURPOSE:
  Calculate the refractive index structure parameter for a maximum of 

MODIFICATIONS:
  Blake Sorenson  - 170628: Written 

"""

# Import functions
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter
from readsounding import *
from compare_cn2 import *

if(len(sys.argv)!=1):
    print "SYNTAX: ./range_finder.py" 
    sys.exit()


p    = 100000.
t    = 288.15 
diff = np.arange(0.005,1.0001,0.005)
ct2 = diff**2.
x =[(((76.*10.**-8.)*(p/(t**2)))**2)*t2 for t2 in ct2]

fig1 = plt.figure()
plt.plot(diff,np.log10(x))
##plt.yscale('log')
#plt.ylim(-17,-12)
plt.title('$C_{n}^{2}$ with $P$=1000 hPa and $T$=288.15 K')
plt.ylabel('log$_{10}$ $C_n^2$ [$m^{2/3}$]')
plt.xlabel('$\Delta$T [K]')
fig1.savefig('cn2_temp_relation.png',dpi=300)
print "Saved image: cn2_temp_relation.png"
