#!/usr/bin/env python
"""


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import sys
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from MODISLib import plot_limits_dict

# NOTE: sus index of 3:15 pm central time = 379
#       fyg index of 3:55 pm central time = 114

# SUS coordinates:
#       -90.6558, 38.6572
# FYG coordinates:
#       -90.9938, 38.5876

if(len(sys.argv) != 2):
    print("SYNTAX: ./process_asos.py asos_file")
    sys.exit()

infile = sys.argv[1]
df = pd.read_csv(infile)

# Make some figures
fig1,ax = plt.subplots(figsize=(12,5))
#fig1.set_size_inches()

mapcrs = ccrs.LambertConformal()
datacrs = ccrs.PlateCarree()
fig2 = plt.figure()
ax2 = fig2.add_subplot(projection = mapcrs)

#fig3,ax3 = plt.subplots()
#fig3.set_size_inches(8,5)

colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple',\
    'tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']

station_names = sorted(set(df['station'].values))
for ii, station in enumerate(station_names):
    time_stn = df['valid'][df['station'] == station]
    tmpc_stn = pd.to_numeric(df['tmpc'][df['station'] == station], errors='coerce').values
    drct_stn = pd.to_numeric(df['drct'][df['station'] == station], errors='coerce').values
    pwxc_stn = df['wxcodes'][df['station'] == station].values
    lat_stn = df['lat'][df['station'] == station].values[0]
    lon_stn = df['lon'][df['station'] == station].values[0]

    # Remove masked data
    # ------------------
    indices = np.where(~np.isnan(tmpc_stn))
    time_stn = time_stn.values[indices] 
    tmpc_stn = np.ma.masked_invalid(tmpc_stn).compressed()
    drct_stn = np.ma.masked_invalid(drct_stn).compressed()
    pwxc_stn = pwxc_stn[indices]

    # Plot the temp data differently for when it reports smoke/haze
    # -------------------------------------------------------------
    tmpc_stn_hz = np.copy(tmpc_stn)
    tmpc_stn_nohz = np.copy(tmpc_stn)
    tmpc_stn_hz   = np.ma.masked_where((pwxc_stn != 'HZ') & \
        (pwxc_stn != 'FU'), tmpc_stn_hz)
    tmpc_stn_nohz = np.ma.masked_where((pwxc_stn == 'HZ') | \
        (pwxc_stn == 'FU'), tmpc_stn_nohz)

    ### Plot the drct data differently for when it reports smoke/haze
    ### -------------------------------------------------------------
    ##drct_stn_hz   = np.copy(drct_stn)
    ##drct_stn_nohz = np.copy(drct_stn)
    ##print(pwxc_stn, drct_stn_hz)
    ##drct_stn_hz   = np.ma.masked_where(pwxc_stn != 'HZ', drct_stn_hz)
    ##drct_stn_nohz = np.ma.masked_where(pwxc_stn == 'HZ', drct_stn_nohz)

    # Convert times to datetime objects
    # In the file, the times are in UTC, so convert to CDT (-5 hrs)
    dtime_stn = [datetime.strptime(ttime,"%Y-%m-%d %H:%M") - timedelta(hours=5) \
    #dtime_sus = [datetime.strptime(ttime,"%Y-%m-%d %H:%M") - timedelta(hours=5) \
            for ttime in time_stn]

    ax.plot(dtime_stn,tmpc_stn,label=station, color=colors[ii])
    #ax.plot(dtime_stn,tmpc_stn_nohz,label=station, color=colors[ii])
    #ax.plot(dtime_stn,tmpc_stn_hz,'--', label=station, color=colors[ii])
    print(station, lat_stn, lon_stn)
    ax2.text(lon_stn, lat_stn, station, transform=datacrs, color=colors[ii])

    #ax3.plot(dtime_stn,drct_stn_nohz,label=station, color=colors[ii])
    #ax3.plot(dtime_stn,drct_stn_hz,'--', label=station, color=colors[ii])

# Convert the file time to a plot_limits_dict format
event_date = datetime.strptime(infile.split('/')[-1].split('_')[-1][:8], "%Y%m%d")
grabber_date = event_date.strftime('%Y-%m-%d')
first_time = list(plot_limits_dict[grabber_date].keys())[0]

# Pull the event time from the plot_limits_dict
event_dtime = event_date + timedelta(hours = int(first_time[:2]), \
    minutes = int(first_time[2:4]))

# Plot a vertical line at the time of the crash
# ---------------------------------------------
ax.axvline(event_dtime,color='tab:red',lw=2,alpha=0.75)
#ax3.axvline(event_dtime,color='tab:red',lw=2,alpha=0.75)

ax.set_xlabel('Time [UTC]')
ax.set_ylabel('Temperature [degC]')
ax.legend()
 
ax2.set_extent([np.min(df['lon']) - 1.0, np.max(df['lon']) + 1.0,\
     np.min(df['lat']) - 1.0, np.max(df['lat']) + 1.0],datacrs)
ax2.coastlines()
ax2.add_feature(cfeature.STATES)
ax2.add_feature(cfeature.LAKES)
ax2.add_feature(cfeature.RIVERS)
ax2.set_title(event_dtime.strftime('%Y-%m-%d %H:%M'))

 
series_name = "asos_time_series_"+event_dtime.strftime('%Y%m%d%H%M')+'.png'
map_name    = "asos_site_map_"+event_dtime.strftime('%Y%m%d%H%M')+'.png'
fig1.savefig(series_name,dpi=300)
print("Saved image",series_name)
fig2.savefig(map_name,dpi=300)
print("Saved image",map_name)
plt.show() 

