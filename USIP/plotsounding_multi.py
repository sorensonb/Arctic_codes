#!/usr/bin/env python

"""
  NAME:
    plotsounding_multi.py

  PURPOSE:
    Plot up to six soundings on a single skew-T diagram

  CALLS:
    - Built-in Python modules np, sys
    - matplotlib.pyplot, matplotlib.gridspec
    - Metpy:
      - metpy.plots:
        - SkewT, Hodograph
      - metpy.units
        - units
    - Custom module readsounding

  PYTHON VERSION:
    2.7.5

  MODIFICATIONS:
    Blake Sorenson  - 170601: Written

  USAGE:
    To run the script, use this syntax:
    >>>./plotsounding_multi.py <up to seven soundings>

  NOTES:
    This script makes use of readsounding to read the soundings for plotting.
    Thus, this program is compatible with COAMPS forecast soundings,
    GRAW ProfileDataTables, NWS soundings (UWYO upper-air database),
    MCR Snow White soundings (CAPE2015), and any model sounding from
    rucsoundings.noaa.gov (RAP.Op40, NAM, GFS, FIM, etc.)
    
    For more documentation on readsounding, please read the comments on 
    /nas/home/bsorenson/prg/thermosonde1617/readsounding.py

    This script also makes use of the metpy module for plotting the 
    soundings. As of the last modification, metpy was installed on
    mammatus, aurora, and haboob. Thus, this script must be run on
    these three computers.

"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from metpy.plots import SkewT 
from metpy.units import units
from readsounding import readsounding

argmax=8
if(len(sys.argv)==1):
    print "SYNTAX: ./plotsounding_multi.py any_number_of_soundings"
    sys.exit(1)

# colors is used to color each sounding and its associated wind barbs
colors=['orange','blue','green','red','cyan','purple','olive','black']

plt.rcParams['figure.figsize'] = (9, 9)
fig1 = plt.figure()

# Determine if the skew-T needs to be saved
save=False
try:
    sys.argv.index('save')
    sys.argv = np.delete(sys.argv,sys.argv.index('save'))
    save=True 
except ValueError:
    # do nothing
    x=1  
"""
if('save' in sys.argv):
    save=True
    filename = sys.argv[0]['NAME']+'_soundings.png'
    sys.argv = np.delete(sys.argv,sys.argv.index('save'))
"""
# Create skew-T plot
skew = SkewT(fig1, rotation=45)

for i in range(1, len(sys.argv)):
    # Use readsounding to read in sounding
    temp=readsounding(sys.argv[i])
    # Generate a label for the current sounding using the sounding's title
    # The name consists of the type and time/forecast of the sounding
    t1label=temp['LABEL'].split()
    tlabel=t1label[0]+' '+t1label[1]
    # Add units to the data
    temp['PRESS']=temp['PRESS']*units.mbar
    temp['TEMP'] = temp['TEMP'] * units.degC
    temp['DP'] = temp['DP'] * units.degC
    temp['UWIND'] = temp['UWIND'] * units.knots
    temp['VWIND'] = temp['VWIND'] * units.knots
    # Plot the temperature and dewpoint profiles
    skew.plot(temp['PRESS'], temp['TEMP'], colors[i%argmax], label=tlabel)
    skew.plot(temp['PRESS'], temp['DP'], colors[i%argmax])
    # Generate wind barbs at every 25 millibars
    plotP = []
    my_interval=np.arange(100, 1000, 25) * units('mbar')
    for center in my_interval:
        index=(np.abs(temp['PRESS']-center)).argmin()
        if index not in plotP:
            plotP.append(index) 
    # Plot the wind barbs using the same color as the temp/dp profiles
    skew.plot_barbs(temp['PRESS'][plotP], temp['UWIND'][plotP],                \
                    temp['VWIND'][plotP], color=colors[i%argmax])

    if((i==1) and (save is True)):
        filename = temp['NAME']+'_soundings.png'

# Add adiabats and mixing ratio lines and set axis limits
skew.ax.set_ylim(1050, 50)
skew.ax.set_xlim(-50,50)
skew.plot_dry_adiabats(np.arange(-50,151,10)*units.degC)
skew.plot_moist_adiabats(np.arange(-50,41,10)*units.degC)
skew.plot_mixing_lines()
plt.title('2019/05/04 03:00 UTC',fontsize=14)
plt.xlabel('Temperature [Degrees C]',fontsize=12)
plt.ylabel('Pressure [mb]',fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

#plt.title(sounding1['TITLE']+" vs "+sounding2['TITLE'])
plt.legend(loc='upper left', fontsize=12)
if save is True:
    plt.savefig(filename,dpi=300)
    print "Saved image: "+filename
else: 
    plt.show()
