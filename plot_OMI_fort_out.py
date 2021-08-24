#!/usr/bin/env python
"""
  NAME:
    plot_OMI_fort_out.py

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2021/04/01:
      Written (modification of AerosolIndexPlotting_Procedure_August112014t0046.pro by rcontreras) 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd

if(len(sys.argv)<3):
    print("SYNTAX: python plot_OMI_fort_out.py out_file min_lat")
    print("        name is of format: omi_VERSION_counts_YYYY_YYYY_ZZZ.txt")
    #print("        name is of format: omi_VERSION_counts_YYYY_YYYY_ZZZ.txt")
    print("           VERSION = VSJ2")
    print("           ZZZ  = 100 * AI thresh (i.e., 060)")
    print("        min_lat is of format: 60,65,70,75,80,85")
    sys.exit()

file_name = sys.argv[1]
min_lat   = sys.argv[2]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extract date information from the file name
name_split = file_name.strip().split('/')[-1].split('_')
dtype      = name_split[1]
vtype      = name_split[2]  # counts or areas?
start_year = name_split[3]
end_year   = name_split[4]
ai_thresh  = float(int(name_split[5].split('.')[0])/100)

#if((dtype != 'vsj22') & (dtype != 'vjz2112')):
#    time_fmt = "%Y%m%d"
#else:
#    time_fmt = "%Y%m%d%H"

if(vtype == 'areas'):
    data_str = 'Area'
    divider = 1e5
    axis_label = 'Area of high AI [10$^{5}$ km$^{2}$]'
elif(vtype == 'counts'):
    data_str = 'Cnt'
    divider = 1
    axis_label = 'Count of quarter-degree lat/lon grids of high AI'
else:
    print("INVALID VARIABLE TYPE. Must be areas or counts")
    sys.exit()

print(ai_thresh)

in_data = pd.read_csv(file_name, delim_whitespace=True)
dates  = in_data['Date'].values
if(len(str(dates[0])) == 10):
    time_fmt = "%Y%m%d%H"
    day_avgs = False
elif(len(str(dates[0])) == 8):
    time_fmt = "%Y%m%d"
    day_avgs = True
else:
    print("INVALID DATE FORMAT")
    sys.exit()
dt_dates = [datetime.strptime(str(tmpdate),time_fmt) \
    for tmpdate in dates]
count65  = in_data[data_str + min_lat].values / divider
#count70  = in_data['Cnt70'].values
#count75  = in_data['Cnt75'].values
#count80  = in_data['Cnt80'].values
#count85  = in_data['Cnt85'].values
x_range = np.arange(len(dates))

# Calculate daily totals
#daily_counts_65 = [np.average(tmparr) for tmparr in \
if(not day_avgs):
    daily_counts_65 = [np.sum(tmparr) for tmparr in \
        np.array_split(count65,len(count65)/4)]
##daily_counts_70 = [np.average(tmparr) for tmparr in \
##    np.array_split(count70,len(count70)/4)]
##daily_counts_75 = [np.average(tmparr) for tmparr in \
##    np.array_split(count75,len(count75)/4)]
##daily_counts_80 = [np.average(tmparr) for tmparr in \
##    np.array_split(count80,len(count80)/4)]
##daily_counts_85 = [np.average(tmparr) for tmparr in \
##    np.array_split(count85,len(count85)/4)]
    daily_dt_dates = dt_dates[::4]
    daily_dates = dates[::4]/100
    daily_xrange = x_range[::4]
else:
    daily_counts_65 = count65
    daily_dates = dates
    daily_dt_dates = dt_dates

## Split by synoptic time
#count_00 = count[0::4]
#count_06 = count[1::4]
#count_12 = count[2::4]
#count_18 = count[3::4]
#xrange_00 = x_range[0::4]
#xrange_06 = x_range[1::4]
#xrange_12 = x_range[2::4]
#xrange_18 = x_range[3::4]

fig1, ax = plt.subplots()
#if(dtype != 'vsj22'):
#    ax.plot(dt_dates,count65,label='synoptic')
ax.plot(daily_dt_dates,daily_counts_65,label='daily')
#ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
ax.legend()
ax.set_ylabel(axis_label)
ax.set_title('AI ' + vtype.title() + ' ' + dtype + ': Threshold of '+str(ai_thresh)+\
    '\nNorth of '+min_lat+'$^{o}$')
ax.grid()

save = False
if(save == True):
    outname = "ai_" + vtype + "_total_" + dtype + ".png"
    plt.savefig(outname,dpi=300)
    print("Saved image",outname) 
plt.show()

#fig2 = plt.figure()
#plt.plot(xrange_00,count_00,label='00Z')
#plt.plot(xrange_06,count_06,label='06Z')
#plt.plot(xrange_12,count_12,label='12Z')
#plt.plot(xrange_18,count_18,label='18Z')
#plt.legend()
#plt.show()
