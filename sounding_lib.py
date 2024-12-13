"""
  NAME:

  PURPOSE:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>         - 2022/04/09:
        Written
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from metpy.plots import SkewT 
from metpy.units import units
import matplotlib.gridspec as gridspec
from datetime import datetime

# User-defined modules
from readsounding import readsounding
from python_lib import *

# colors is used to color each sounding and its associated wind barbs
colors=['tab:cyan','tab:red','tab:red','cyan','purple','olive','black']

lwidth = 1.0

def plot_sounding_figure(in_data, fig = None, skew = None, save = False, \
        ptitle = None, color = None):

    dt_date_str = datetime.strptime(in_data[0].strip().split('/')[-1][:13],'%y%m%d_%H%M%S')

    # Set up figure, if no axis passed in
    # -----------------------------------
    in_ax = True
    ##!#if(skew is None):
    ##!#    # Create skew-T plot
    ##!#    fig = plt.figure()
    ##!#    skew = SkewT(fig, rotation=45)
    ##!#else:
    ##!#    # Create skew-T plot
    ##!#    print('here')
    skew = SkewT(fig, rotation=45, subplot = skew)
   
    if(isinstance(in_data, str)):
        in_data = [in_data]
 
    for ii, ifile in enumerate(in_data):
        # Use readsounding to read in sounding
        data=readsounding(ifile)

        plot_sounding(data, skew, idx = ii, color = color)
    
    # Add adiabats and mixing ratio lines and set axis limits
    skew.ax.set_ylim(1050, 50)
    skew.ax.set_xlim(-50,50)
    skew.plot_dry_adiabats(np.arange(-50,200,10)*units.degC, linewidth = lwidth - 0.2)
    skew.plot_moist_adiabats(np.arange(-50,41,10)*units.degC, linewidth = lwidth - 0.2)
    skew.plot_mixing_lines(linewidth = lwidth - 0.2)
    if(ptitle is None):
        skew.ax.set_title(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    else:
        skew.ax.set_title(ptitle, fontsize = 10)
        skew.ax.set_title('  a)', loc = 'left', weight = 'bold')
    #skew.ax.text(-158, 46, 'a)', weight = 'bold', fontsize = 12)
    #plt.title('2019/05/04 03:00 UTC',fontsize=14)
    skew.ax.set_xlabel('Temperature [Degrees C]', fontsize = 10)
    #plt.xlabel('Temperature [Degrees C]',fontsize=12)
    skew.ax.set_ylabel('Pressure [mb]', fontsize = 10)
    #plt.ylabel('Pressure [mb]',fontsize=12)
    skew.ax.tick_params(axis = 'both', labelsize=8)
    
    #plt.title(sounding1['TITLE']+" vs "+sounding2['TITLE'])
    skew.ax.legend(loc='upper left', prop={'size': 8}, framealpha = 1)
    #plt.legend(loc='upper left', fontsize=12)
    if(not in_ax):
        if(save):
            filename = data['NAME']+'_soundings.png'
            fig.savefig(filename,dpi=300)
            print("Saved image: "+filename)
        else: 
            plt.show()

def plot_sounding(data, skew, idx = 0, color = None):
    # Generate a label for the current sounding using the sounding's title
    # The name consists of the type and time/forecast of the sounding
    t1label=data['LABEL'].split()
    tlabel=t1label[0]
    if(tlabel == 'GRAW'):
        tlabel = 'RAOB'
    #tlabel=t1label[0]+' '+t1label[1]
    # Add units to the data
    data['PRESS']=data['PRESS']*units.mbar
    data['TEMP'] = data['TEMP'] * units.degC
    data['DP'] = data['DP'] * units.degC
    data['UWIND'] = data['UWIND'] * units.knots
    data['VWIND'] = data['VWIND'] * units.knots
    # Plot the dataerature and dewpoint profiles
    if(color is None):
        color = colors[idx]
    skew.plot(data['PRESS'], data['TEMP'], color, label=tlabel, linewidth = lwidth + 0.2)
    skew.plot(data['PRESS'], data['DP'], color, linewidth = lwidth + 0.2)
    # Generate wind barbs at every 25 millibars
    plotP = []
    my_interval=np.arange(100, 1000, 75) * units('mbar')
    #my_interval=np.arange(100, 1000, 25) * units('mbar')
    for center in my_interval:
        index=(np.abs(data['PRESS']-center)).argmin()
        if index not in plotP:
            plotP.append(index) 
    # Plot the wind barbs using the same color as the data/dp profiles
    skew.plot_barbs(data['PRESS'][plotP], data['UWIND'][plotP],                \
                    data['VWIND'][plotP], color = colors[idx], linewidth = lwidth - 0.2)
