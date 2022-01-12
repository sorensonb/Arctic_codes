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
import matplotlib.cm as cm
import matplotlib.colors as mc
from datetime import datetime
import pandas as pd
from scipy.signal import argrelextrema, find_peaks

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
thresh     = name_split[5].split('.')[0]
ai_thresh  = float(int(thresh)/100)

if second_file:
    name_split  = file_name2.strip().split('/')[-1].split('_')
    dtype2      = name_split[1]
    vtype2      = name_split[2]  # counts or areas?
    start_year2 = name_split[3]
    end_year2   = name_split[4]
    
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

if(int(min_lat) == 75):
    interval = 1
    h_val = 1
elif(int(min_lat) > 75):
    interval = 0.2
    h_val = 0.1
else:
    interval = 2
    h_val = 1

print(ai_thresh)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Read area data from the input file
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
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

# Convert the areas from square kilometers to 1e5 square kilometers
count65  = in_data[data_str + min_lat].values / divider
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
        daily_dt_dates_2 = dt_dates_2[::4]
        daily_dates_2 = dates2[::4]/100
        daily_xrange_2 = x_range_2[::4]
    else:
        # Just use the daily totals from the file
        daily_counts_65_2 = count65_2
        daily_dates_2 = dates2
        daily_dt_dates_2 = dt_dates_2

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

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Plot the individual years of area time series
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

fig0 = plt.figure(figsize = (10,4))
ax0 = fig0.add_subplot(1,1,1)
plot_c = cm.turbo((years-np.min(years))/(np.max(years)-np.min(years)))

# Set up an array to hold the event counts
# ----------------------------------------
d_val = 4
#h_val = 1
if(interval == 0.2):
    max_val = np.round(np.max(daily_counts_65), 1)+interval
else:
    if(int(min_lat) == 75):
        max_val = int(np.max(daily_counts_65))+2*interval
    else:
        max_val = int(np.max(daily_counts_65))+interval

event_sizes = np.arange(h_val, max_val , interval)
events_yearly = np.zeros((years.shape[0], event_sizes.shape[0]))

for ii, year in enumerate(years):
    yearly_totals[ii] = np.sum(daily_counts_65[np.where( \
        (daily_dt_dates >= datetime(year,4,1)) & \
        (daily_dt_dates <= datetime(year,9,30)))])

    # Test the local maxima
    # ---------------------
    for jj in range(len(event_sizes) - 1):
        # Find the peaks for this size
        # ----------------------------
        peaks, _ = find_peaks(daily_counts_65[np.where( \
            (daily_dt_dates >= datetime(year,4,1)) & \
            (daily_dt_dates <= datetime(year,9,30)))], \
            height = [event_sizes[jj], event_sizes[jj+1]], distance = d_val)
        events_yearly[ii,jj] = len(peaks)        
 
    # Determine the counts of peaks greater than each height range
    # ------------------------------------------------------------

    # Calculate the days since April 1
    # --------------------------------
    loop_dates = daily_dt_dates[np.where( \
        (daily_dt_dates >= datetime(year,4,1)) & \
        (daily_dt_dates <= datetime(year,9,30)))]

    plot_dates = datetime(year = 2000, month = 4, day = 1) + \
        (loop_dates - datetime(year = year, month = 4, day = 1))

    ax0.plot(plot_dates,\
        daily_counts_65[np.where( \
        (daily_dt_dates >= datetime(year,4,1)) & \
        (daily_dt_dates <= datetime(year,9,30)))], c = plot_c[ii])

ax0.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
norm = mc.Normalize(vmin=np.min(years), vmax=np.max(years))

cmap = cm.turbo
cbar = fig0.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax0, orientation='vertical', label='Year')

#ax0.plot(plot_dates, mask_day_avgs,label='avg',color='black')
#ax0.plot(plot_dates, upper_range,label='+avg σ',linestyle='--',color='black')
ax0.set_ylabel(axis_label,weight='bold',fontsize=12)
ax0.set_title('AI ' + vtype + ': Threshold of '+str(ai_thresh)+\
    '\nNorth of '+min_lat+'$^{o}$', fontsize = 14, weight = 'bold')
ax0.tick_params(axis='both', labelsize = 12)
ax0.grid()
fig0.tight_layout()
save = True 
if(save == True):
    outname = "ai_" + vtype + "_indiv_yearly_" + dtype + "_thresh"+thresh+"_minlat"+min_lat+".png"
    fig0.savefig(outname,dpi=300)
    print("Saved image",outname) 

plt.legend()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Plot the total yearly areas on a graph
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
fig2 = plt.figure()
ax4 = fig2.add_subplot(1,1,1)
ax4.plot(years,yearly_totals)

if(save == True):
    outname = "ai_" + vtype + "_total_" + dtype + "_thresh"+thresh+"_minlat"+min_lat+".png"
    fig2.savefig(outname,dpi=300)
    print("Saved image",outname) 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Plot the total time series of daily areas with peaks
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
fig1 = plt.figure(figsize = (13,5)) 
ax1 = fig1.add_subplot(1,2,1)
ax3 = fig1.add_subplot(1,2,2)
ax1.plot(daily_dt_dates,daily_counts_65,label='daily '+dtype)

# Test peaks here
# ---------------
peaks, _ = find_peaks(daily_counts_65, height = h_val, distance = d_val)
ax1.plot(daily_dt_dates[peaks], daily_counts_65[peaks], 'x', color='black')

if second_file:
    ax2 = ax.twinx()
    ax2.plot(daily_dt_dates_2,daily_counts_65_2,label='daily '+dtype2,color='tab:orange')
ax1.legend()
ax1.set_ylabel(axis_label)
ax1.set_title('AI ' + vtype.title() + ': Threshold of '+str(ai_thresh)+\
    '\nNorth of '+min_lat+'$^{o}$')
ax1.grid()
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Plot the yearly event size box graph
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#fig3 = plt.figure()
#ax3  = fig3.add_subplot(1,1,1)
for ii in range(len(event_sizes)-1):
    ax3.bar(years, events_yearly[:,ii], \
        bottom = np.sum(events_yearly[:,0:ii],axis=1), \
        label=str(np.round(event_sizes[ii],1)) + ' - ' + str(np.round(event_sizes[ii+1],1)) + ' * 10$^{5}$ km$^{2}$')
ax3.legend()
ax3.set_ylabel('# of AI peaks in each size range')
ax3.set_title('AI ' + vtype.title() + ' : Threshold of '+str(ai_thresh)+\
    '\nNorth of '+min_lat+'$^{o}$')
#save = False
#if(save):
#    outname = 'ai_bar_minhgt'+str(int(event_sizes[0])) + '_dist'+str(d_val)+'_thresh'+thresh+"_minlat"+min_lat+'.png'
#    fig3.savefig(outname,dpi=300)
#    print("Saved image",outname)

if(save):
    outname = 'ai_daily_areas_peaks_thresh'+thresh+"_minlat"+min_lat+'.png'
    fig1.savefig(outname, dpi = 300)
    print("Saved image",outname)


plt.show()

#fig2 = plt.figure()
#plt.plot(xrange_00,count_00,label='00Z')
#plt.plot(xrange_06,count_06,label='06Z')
#plt.plot(xrange_12,count_12,label='12Z')
#plt.plot(xrange_18,count_18,label='18Z')
#plt.legend()
#plt.show()
