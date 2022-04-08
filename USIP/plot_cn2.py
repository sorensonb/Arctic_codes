#!/usr/bin/env python

"""
  NAME:
    plot_cn2.py
  
  PURPOSE:
    Plot any number of refractive index structure parameter profiles on a single
    graph. 
  
  MODIFICATIONS:
    Blake Sorenson  - 170714: Written 
  
"""
def syntax():
    print "SYNTAX: ./plot_cn2_single.py <soundings> [raw] [old / model] [save]"
    print "         raw       - do no smoothing"
    print "         old/model - use the data directly above and below"
    print "                     each level to calculate Cn2"
    print "         save      - saves the image"

# Import functions
import numpy
import sys
import matplotlib.pyplot as plt
from readsounding import *
from compare_cn2 import *

if(len(sys.argv)<2):
    syntax()
    sys.exit()

argmax=9
# Check if the data should be plotted without any smoothing
raw=False
old=False
save=False
# Assume that if the user wants to save an image, there will only be one
# sounding passed
savename = ''
method="thermo"
for arg in sys.argv:    
    if((arg=='-h') | (arg=='help')):
        syntax()
        sys.exit()
    if(arg=="raw"):
        raw=True
    elif((arg=="old") | (arg=="model")): 
        old=True
        method="model"
    elif(arg=='save'):
        save=True
# If "raw" was one of the arguments, remove it from the argument list.
# Otherwise, it will be passed as an argument to readsounding, which won't
# work.
if(raw==True):
    sys.argv = numpy.delete(sys.argv,-1)
if(old==True):
    sys.argv = numpy.delete(sys.argv,-1)
if(save==True):
    sys.argv = numpy.delete(sys.argv,-1)
colors=['orange','black','red','cyan','green','olive','blue','purple','yellow']
# Use plot_cn2() to make a plot with a logarithmic x axis
plot_cn2()
# Loop over the arguments and plot the cn2 profile for each sounding
for i in range(1, len(sys.argv)):
    # Use readsounding to read the current sounding file into a dictionary
    a = readsounding(sys.argv[i])
    # Convert wind speeds from knots to meters/second
    #a['UWIND'] = kts2ms(a['UWIND'])
    #a['VWIND'] = kts2ms(a['VWIND'])
    # Apply an 11-point smoother to the T, u, and v data in the sounding
    # if it is over 900 elements long 
    if((len(a['TEMP'])>900) & (raw==False)):
        a = smooth(a) 
        a = cn2_calc(a,method=method)
    else:
        a = cn2_calc(a,method=method)
    # Plot cn2 with height
    #olabel=a['TITLE'].split(" ")[0]+" "+a['TITLE'].split(" ")[1]+" "+\
    #       a['TITLE'].split(" ")[-1]
    olabel=a['LABEL']
    savename = a['NAME']+'_CN2.png'
    plt.plot(np.log10(a['CN2']),(a['ALT']/1000.),color=colors[i%argmax],label=olabel)
plt.legend(loc='upper right', fontsize=10)
plt.ylim(0,12)
if save is True:
    plt.savefig(savename,dpi=300)
    print "Saved image: "+savename
else:
    plt.show()
