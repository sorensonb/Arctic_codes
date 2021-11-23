#!/usr/bin/env python
"""
  NAME:
    plot_OMI_fort_out.py

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2021/04/01:
      Written 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd

if(len(sys.argv)<3):
    print("SYNTAX: python plot_OMI_fort_out.py out_file [optional second file] min_lat")
    print("        name is of format: omi_VERSION_counts_YYYY_YYYY_ZZZ.txt")
    #print("        name is of format: omi_VERSION_counts_YYYY_YYYY_ZZZ.txt")
    print("           VERSION = VSJ2")
    print("           ZZZ  = 100 * AI thresh (i.e., 060)")
    print("        min_lat is of format: 60,65,70,75,80,85")
    sys.exit()

second_file = False
file_name = sys.argv[1]
min_lat   = sys.argv[2]
if(len(min_lat.strip().split('.')) != 1):
    # Second file provided
    second_file = True
    file_name2 = sys.argv[2]
    min_lat = sys.argv[3]

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

if second_file:
    name_split  = file_name2.strip().split('/')[-1].split('_')
    dtype2      = name_split[1]
    vtype2      = name_split[2]  # counts or areas?
    start_year2 = name_split[3]
    end_year2   = name_split[4]
    

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

if second_file:
    in_data = pd.read_csv(file_name2, delim_whitespace=True)
    dates2  = in_data['Date'].values
    if(len(str(dates2[0])) == 10):
        time_fmt = "%Y%m%d%H"
        day_avgs = False
    elif(len(str(dates2[0])) == 8):
        time_fmt = "%Y%m%d"
        day_avgs = True
    else:
        print("INVALID DATE FORMAT")
        sys.exit()
    dt_dates_2 = [datetime.strptime(str(tmpdate),time_fmt) \
        for tmpdate in dates2]
    count65_2  = in_data[data_str + min_lat].values / divider
    x_range_2 = np.arange(len(dates2))
    

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
    # Just use the daily totals from the file
    daily_counts_65 = count65
    daily_dates = dates
    daily_dt_dates = dt_dates

if second_file:
    # Calculate daily totals
    #daily_counts_65 = [np.average(tmparr) for tmparr in \
    if(not day_avgs):
        daily_counts_65_2 = [np.sum(tmparr) for tmparr in \
            np.array_split(count65_2,len(count65_2)/4)]
    ##daily_counts_70 = [np.average(tmparr) for tmparr in \
    ##    np.array_split(count70,len(count70)/4)]
    ##daily_counts_75 = [np.average(tmparr) for tmparr in \
    ##    np.array_split(count75,len(count75)/4)]
    ##daily_counts_80 = [np.average(tmparr) for tmparr in \
    ##    np.array_split(count80,len(count80)/4)]
    ##daily_counts_85 = [np.average(tmparr) for tmparr in \
    ##    np.array_split(count85,len(count85)/4)]
        daily_dt_dates_2 = dt_dates_2[::4]
        daily_dates_2 = dates2[::4]/100
        daily_xrange_2 = x_range_2[::4]
    else:
        # Just use the daily totals from the file
        daily_counts_65_2 = count65_2
        daily_dates_2 = dates2
        daily_dt_dates_2 = dt_dates_2


## Split by synoptic time
#count_00 = count[0::4]
#count_06 = count[1::4]
#count_12 = count[2::4]
#count_18 = count[3::4]
#xrange_00 = x_range[0::4]
#xrange_06 = x_range[1::4]
#xrange_12 = x_range[2::4]
#xrange_18 = x_range[3::4]

# 


# Add up the areas for each year
daily_dt_dates = np.array(daily_dt_dates)
daily_counts_65 = np.array(daily_counts_65)

years = np.arange(int(start_year),int(end_year) + 1)
yearly_totals = np.zeros(len(years))

# Create a 2d array to hold all the area data from each day
day_values = np.full((int(datetime(year=2008,month=9,day=30).strftime('%j')) \
    - int(datetime(year=2008,month=4,day=1).strftime('%j'))+1, len(years)+1),-9.)

# Insert each day's area data into the array
for xx in range(len(daily_dt_dates)):
    indices = int(daily_dt_dates[xx].strftime('%j')) - \
        int(datetime(year=daily_dt_dates[xx].year,month=4,\
        day=1).strftime('%j')), daily_dt_dates[xx].year - \
        daily_dt_dates[0].year 
    day_values[indices] = daily_counts_65[xx]

# Mask any -9s, which are caused by days without data
mask_values = np.ma.masked_where(day_values == -9,day_values)
# Calculate the mean and standard deviation of the areas for each
# day of the year.
mask_day_avgs = np.nanmean(mask_values,axis=1)
mask_day_stds = np.nanstd(mask_values,axis=1)

mean_total_std = np.nanstd(mask_values)
#mean_std = np.nanmean(mask_day_stds)

# Calculate the
# This adds and subtracts the standard deviation for each day
# ---------------------------------------------------------------------
lower_range = mask_day_avgs - mask_day_stds * 0.75
upper_range = mask_day_avgs + mask_day_stds * 0.75

# This adds and subtracts the total standard deviation of all the areas
# ---------------------------------------------------------------------
#lower_range = mask_day_avgs - mean_total_std * 0.75
#upper_range = mask_day_avgs + mean_total_std * 0.75

# Use a 1-sigma check to mask any daily_counts_65 values that are
# outside of the average
new_daily_counts = np.full((len(daily_dt_dates)),-9.)
for xx in range(len(daily_dt_dates)):
    indices = int(daily_dt_dates[xx].strftime('%j')) - \
        int(datetime(year=daily_dt_dates[xx].year,month=4,\
        day=1).strftime('%j')), daily_dt_dates[xx].year - \
        daily_dt_dates[0].year 
    if((daily_counts_65[xx] - mask_day_avgs[indices[0]]) > 1. * \
    #if((daily_counts_65[xx] - mask_day_avgs[indices[0]]) < 0.5 * \
        mask_day_stds[indices[0]]):
        #daily_counts_65[xx] = -9
        new_daily_counts[xx] = daily_counts_65[xx]

daily_counts_65 = np.ma.masked_where(daily_counts_65 == -9,daily_counts_65)

fig0 = plt.figure()
for ii, year in enumerate(years):
    yearly_totals[ii] = np.sum(daily_counts_65[np.where( \
        (daily_dt_dates >= datetime(year,4,1)) & \
        (daily_dt_dates <= datetime(year,9,30)))])
    print(len(daily_dt_dates[np.where( \
        (daily_dt_dates >= datetime(year,4,1)) & \
        (daily_dt_dates <= datetime(year,9,30)))]))
    plt.plot(np.arange(len(daily_dt_dates[np.where( \
        (daily_dt_dates >= datetime(year,4,1)) & \
        (daily_dt_dates <= datetime(year,9,30)))])),\
        daily_counts_65[np.where( \
        (daily_dt_dates >= datetime(year,4,1)) & \
        (daily_dt_dates <= datetime(year,9,30)))],label=str(year))
plt.plot(mask_day_avgs,label='avg',color='black')
plt.plot(upper_range,label='+avg σ',linestyle='--',color='black')
plt.plot(lower_range,label='-avg σ',linestyle='--',color='black')
#daily_dt_dates[np.where( (daily_dt_dates >= datetime(2018,4,1)) & (daily_dt_dates <= datetime(2018,9,30)))
plt.legend()

fig2 = plt.figure()
plt.plot(years,yearly_totals)

fig1, ax = plt.subplots()
#if(dtype != 'vsj22'):
#    ax.plot(dt_dates,count65,label='synoptic')
ax.plot(daily_dt_dates,daily_counts_65,label='daily '+dtype)
if second_file:
    ax2 = ax.twinx()
    ax2.plot(daily_dt_dates_2,daily_counts_65_2,label='daily '+dtype2,color='tab:orange')
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
