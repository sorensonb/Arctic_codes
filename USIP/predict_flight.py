#!/usr/bin/env python

"""
  NAME:
    predict_landing.py

  PURPOSE:
    Calculate the predicted landing point of a radiosonde package. 

  CALLS:
    - Python modules np, sys, re, matplotlib.pyplot, matplotlib.gridspec,
      mpl_toolkits.basemap, scipy.interpolate
      (Requires scipy version 0.17.0 or greater)
    - Custom module readsounding

  SYNTAX:
    /nas/home/bsorenson/prg/thermosonde1617/predict_landing.py <GRAW RTS file> <whole or end>
    
    This program needs a Graw RTS table (/nas/Radiosonde/20170503/170503_204757_GFK_GRAW.txt) 
    
  MODIFICATIONS:
    Blake Sorenson  - 170616: Written

"""

import numpy as np
import sys
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import commands
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
from readsounding import readsounding

if(len(sys.argv)!=5):
    print "SYNTAX: ./predict_flight.py sounding_file starting_lat starting_lon whole/final > <outputfilename>.txt"
    print "         Whole: Predicts the entire balloon path starting with the 100th point"
    print "                and using the u and v winds to calculate the path"
    print "         Final: Predicts the final few minutes of the balloon path starting "
    print "                at the last data point received from the radiosonde."
    sys.exit(1)

a = readsounding(sys.argv[1])
inlat = float(sys.argv[2])
inlon = float(sys.argv[3])
wf = sys.argv[4]

"""
# Extrapolation

# To predict landing location, extrapolate altitude to surface. Assume a
# constant rate of descent 
tx = TIME
ty = ALT
tlat = LAT
tlon = LON
altf = interpolate.interp1d(tx,ty,fill_value='extrapolate')
i = int(tx[1])

# Append the extrapolated altitudes and times to arrays
while(altf(i)>ALT[0]):
    fx = np.append(fx[:],float(i))
    fy = np.append(fy[:], altf(i))
    i+=1 
#latf = interpolate.interp1d(tx,tlat,fill_value='extrapolate')
#lonf = interpolate.interp1d(tx,tlon,fill_value='extrapolate')
#i = int(tx[0])
"""
# Generate basemap
fig5 = plt.figure(figsize=(10,8))
m = Basemap(llcrnrlon=-98., llcrnrlat=45., urcrnrlon=-90., urcrnrlat=50.,\
           resolution='f',projection='cass',lon_0=-95.,lat_0=45.)
#m.etopo()
#m.shadedrelief()
m.drawcoastlines()
m.drawcountries()
m.drawstates()
#m.drawcounties()
# Convert balloon latitude and longitude to map projection coordinates
# (in meters)
plotx, ploty = m(inlon,inlat)

# Predict balloon path from point 100
# Use readsounding as a shortcut to get the U and V wind components 
a['UWIND']*=0.514444
a['VWIND']*=0.514444

count=0
if(wf.lower()=="whole"):
    # Predict the entire balloon path starting at the 100th point
    #closealt = a['ALT'][:].flat[np.abs(a['ALT'][:]-fy[i]).argmin()]
    # Proof of concept: predict the entire flight path
    startx = plotx
    starty = ploty
    newx = np.array([])
    newy = np.array([])

    # Interpolate the U and V winds, along with altitudes
    # Assume 5 m/s ascent and descent rate
    interp_alts = np.arange(a['ALT'][0],a['ALT'][-1],5.0) 
    interp_u = interpolate.interp1d(a['ALT'],a['UWIND'])
    interp_v = interpolate.interp1d(a['ALT'],a['VWIND'])
    new_u = interp_u(interp_alts)
    new_v = interp_v(interp_alts)
    
    for i in range(1, new_u):
        startx+=new_u[i]#*1sec
        starty+=new_v[i]#*1sec
        print startx, starty
        newx = np.append(newx[:],startx) 
        newy = np.append(newy[:],starty) 
    """
    for i in range(len(a['UWIND']),len(fy)):
        closealt = a['ALT'][:].flat[np.abs(a['ALT'][:]-fy[i]).argmin()]
        caindex = np.where(a['ALT']==closealt)     
        startx+=a['UWIND'][caindex]
        starty+=a['VWIND'][caindex]
        newx[count] = startx 
        newy[count] = starty 
        count+=1
    """
else:
    # Predict the unknown portion of the balloon path starting at the last
    # known point. To do this, for each extrapolated altitude, a comparible
    # altitude will be found from the ascent. The u and v wind components
    # from that altitude on the ascent will be used to find the next 
    # balloon position. Assumes a similar wind pattern on the descent 
    startx = plotx[-1]
    starty = ploty[-1]
   
    newx = np.zeros(len(fy)-len(LAT))
    newy = np.zeros(len(fy)-len(LAT))
    
    for i in range(len(LAT),len(fy)):
        # Find the altitude from the ascent profile that is closest to 
        # the current extrapolated descent altitude
        closealt = a['ALT'][:].flat[np.abs(a['ALT'][:]-fy[i]).argmin()]
        # Find the index of that ascent altitude 
        caindex = np.where(a['ALT']==closealt)     
        # Use the index of the ascent altitude to find the u and v
        # components at that altitude.
    
        ###########################################################################
        # 
        # Since it is assumed that the winds on the descent are similar
        # to the winds on the ascent and the movement of the balloon is only due
        # to wind, the next x and y coordinates can be found by adding the current 
        # u and v wind components to the current x and y positions, respectively.  
        #
        # To find the next balloon position, these equations are used:
        #   new_x_position = old_x_position + (current_u_speed*1s)
        #   new_y_position = old_y_position + (current_v_speed*1s)
        #
        # Because the wind speeds are in meters per second and the projection
        # coordinates are in meters, the new x and y coordinates can be found
        # by first multiplying the u and v winds by 1 second to find the change
        # in the x and y directions caused by the u and v winds (the 1 Hz data from
        # the radiosonde allows this). Adding the change in x and y coordinates
        # to the current x and y coordinates gives the new x and y coordinates.
        #
        ###########################################################################
        startx+=a['UWIND'][caindex]
        starty+=a['VWIND'][caindex]
        newx[count] = startx 
        newy[count] = starty 
        count+=1
    
# Convert the predicted x and y projection coordinates to lat/lon
printlon, printlat = m(newx, newy, inverse=True)
print "Latitude,Longitude"
for x, y, in zip(printlat, printlon):
    print str(x)+","+str(y)

#m.drawrivers()
#plotfx, plotfy = m(flon,flat)
#m.plot(plotfx, plotfy, color='red')
m.plot(plotx, ploty)
m.plot(newx,newy)
fig5.show()
plt.show()
"""
x = np.arange(0,10)
y = np.exp(-x/3.0)
f = interpolate.interp1d(x, y, fill_value='extrapolate')

print f(9)
print f(11)
"""
