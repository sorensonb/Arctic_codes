#!/usr/bin/env python
"""
  NAME:
    plot_single_OMI.py

  PURPOSE:
    Plot data from a raw OMI HDF5 file using cartopy. See the 'path_dict'
    dictionary keys for accepted variables.

    NOTE: If running on a computer that is not Blake Sorenson's office computer,
    change the 'base_path' variable to point to the directory containing the
    HDF5 OMI files.

  SYNTAX:
    $ python plot_single_OMI.py date variable row_max

       - date: for a single-swath image, use a format of YYYYMMDDHH
               for a daily average, use a format of YYYYMMDD
       - variable: variable name from the HDF5 file. This code is only
               set up to handle certain variables right now, so to see
               which variables are currently allowed, see either the syntax
               message or the keys to any of the dictionaries found below
       - row_max:  the maximum row to plot data through. For example, if you
               only want to plot rows 1 to 23, use a row_max value of 23.
               Otherwise, use a value of 60 to plot all rows.

  EXAMPLE:
    $ python plot_single_OMI.py 201807260741 UVAerosolIndex 60 

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
    print("        name is of format: omi_counts_YYYY_YYYY.txt")
    print("        min_lat is of format: 65,70,75,80,85")
    sys.exit()

file_name = sys.argv[1]
min_lat   = sys.argv[2]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extract date information from the file name
name_split = file_name.strip().split('_')
start_year = name_split[2]
end_year   = name_split[3]
ai_thresh  = float(int(name_split[4].split('.')[0])/100)
print(ai_thresh)

in_data = pd.read_csv(file_name, delim_whitespace=True)
dates  = in_data['Date'].values
dt_dates = [datetime.strptime(str(tmpdate),"%Y%m%d%H") \
    for tmpdate in dates]
count65  = in_data['Cnt'+min_lat].values
#count70  = in_data['Cnt70'].values
#count75  = in_data['Cnt75'].values
#count80  = in_data['Cnt80'].values
#count85  = in_data['Cnt85'].values
x_range = np.arange(len(dates))

# Calculate daily averages
daily_counts_65 = [np.average(tmparr) for tmparr in \
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
ax.plot(dt_dates,count65,label='synoptic')
ax.plot(daily_dt_dates,daily_counts_65,label='daily')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.legend()
ax.set_ylabel('Counts')
ax.set_title('AI Counts: Threshold of '+str(ai_thresh)+\
    '\nNorth of '+min_lat+'$^{o}$')
ax.grid()
plt.show()

#fig2 = plt.figure()
#plt.plot(xrange_00,count_00,label='00Z')
#plt.plot(xrange_06,count_06,label='06Z')
#plt.plot(xrange_12,count_12,label='12Z')
#plt.plot(xrange_18,count_18,label='18Z')
#plt.legend()
#plt.show()
