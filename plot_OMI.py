#!/usr/bin/env python
"""
  NAME:
    plot_OMI.py

  PURPOSE:
    Plot trends in OMI-measured aerosol index across the northern hemisphere.
    The average AI values are read from files produced by 
    Average_AerosolIndexCalculator.pro (included in ai_trend_codes_YYMMDD.tar.gz)

  PYTHON VERSION:
    2.6.6

  MODULES:
    - Custom module AILib
    - Matplotlib
    - mpl_toolkits.basemap
    - datetime
    - sys
    - numpy
    
"""

import numpy as np
import sys
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.colors as color
##import matplotlib.colors as colors
from matplotlib.colors import rgb2hex,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
from scipy import stats
from OMILib import readOMI,plotOMI,plotOMI_Climo,plotOMI_MK


# Check arguments
if(len(sys.argv)<4):
    print("SYNTAX: ./plot_OMI.py file_type start_date end_date [key] [trend/pval/climo] [ttype=<standard,thiel-sen>] [season=<season_option>] [save]")
    print("         ttype : NXAR  = no xtrack all rows")
    print("                 XAR   = xtrack all rows")
    print("                 XR123 = xtrack rows 1-23")
    print("         key   : (optional) :  formatted like LATxLON")
    print("              Example: 48x-97")
    print("         ttype : (optional) :  standard (default),thiel-sen")
    print("         season: (optional) :  spring,summer,autumn,winter")
    sys.exit()

ftype = sys.argv[1]
start_date = sys.argv[2]
end_date = sys.argv[3]
save_arg=False
trend = True
pval  = False
climo = False
ttype = 'standard'
#ftype = 'XR123'
season_arg=''
for arg in sys.argv:
    split_arg = arg.split('=')
    if(arg=='climo'):
        trend = False
        climo = True
    elif(arg=='pval'):
        trend = False
        pval = True
    elif(arg=='save'):
        save_arg=True
    elif(split_arg[0]=='season'):
        season_arg=split_arg[1]
    elif(split_arg[0]=='ttype'):
        ttype=split_arg[1]
file_path = '/home/bsorenson/HighLatitudeStudy/OMI/30to90/'
# NXAR = No Xtrack, All Rows
if(ftype=='NXAR'):
    input_file = file_path+'no_xtrack_all_rows/omi_monthly_average_AI_noxtrack_allrows.txt.gz'
# XAR = Xtrack, All Rows
elif(ftype=='XAR'):
    input_file = file_path+'xtrack_all_rows/omi_monthly_average_AI_xtrack_allrows.txt.gz'
# XR123 = XTrack, Rows 1-23
elif(ftype=='XR123'):
    input_file = file_path+'omi_monthly_average_AI_xtrack_row023.txt.gz'
else:
    print('ERROR: Invalid file type. See syntax')
    sys.exit()

global AIDict

AIDict = readOMI(input_file,int(start_date),int(end_date))
print('Data read successful')

if(trend is True):
    plotOMI(AIDict,save=save_arg,trend_type=ttype,file_type=ftype,season=season_arg)
elif(pval is True):
    plotOMI_MK(AIDict,save=save_arg,file_type='XR123',season='')
else:
    plotOMI_Climo(AIDict,save=save_arg,trend_type=ttype,file_type=ftype,season=season_arg)


