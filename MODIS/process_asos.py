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
fig1,ax = plt.subplots()
fig1.set_size_inches(8,5)

mapcrs = ccrs.LambertConformal()
datacrs = ccrs.PlateCarree()
fig2 = plt.figure()
ax2 = fig2.add_subplot(projection = mapcrs)

colors = ['tab:blue','tab:orange','tab:olive','tab:purple','tab:red']

station_names = set(df['station'].values)
for ii, station in enumerate(station_names):
    time_stn = df['valid'][df['station'] == station]
    tmpc_stn = pd.to_numeric(df['tmpc'][df['station'] == station], errors='coerce').values
    pwxc_stn = df['wxcodes'][df['station'] == station].values
    lat_stn = df['lat'][df['station'] == station].values[0]
    lon_stn = df['lon'][df['station'] == station].values[0]

    # Plot the data differently for when it reports smoke/haze
    print(pwxc_stn.shape, tmpc_stn.shape, time_stn.shape)
    tmpc_stn_hz = np.copy(tmpc_stn)
    tmpc_stn_nohz = np.copy(tmpc_stn)
    print(pwxc_stn.shape, tmpc_stn_hz.shape, tmpc_stn_nohz.shape, time_stn.shape)
    tmpc_stn_hz   = np.ma.masked_where(tmpc_stn[np.where(pwxc_stn == 'HZ')], tmpc_stn_hz)
    tmpc_stn_nohz = np.ma.masked_where(tmpc_stn[np.where(pwxc_stn != 'HZ')], tmpc_stn_nohz)

    # Convert times to datetime objects
    # In the file, the times are in UTC, so convert to CDT (-5 hrs)
    dtime_stn = [datetime.strptime(ttime,"%Y-%m-%d %H:%M") - timedelta(hours=5) \
    #dtime_sus = [datetime.strptime(ttime,"%Y-%m-%d %H:%M") - timedelta(hours=5) \
            for ttime in time_stn]

    ax.plot(dtime_stn,tmpc_stn_nohz,label=station, color=colors[ii])
    ax.plot(dtime_stn,tmpc_stn_hz,label=station,linestlye='--', color=colors[ii])
    print(station, lat_stn, lon_stn)
    ax2.text(lon_stn, lat_stn, station, transform=datacrs)

ax.set_xlabel('Time [UTC]')
ax.set_ylabel('Temperature [degC]')
ax.legend()
 
ax2.set_extent([np.min(df['lon']) - 1.0, np.max(df['lon']) + 1.0,\
     np.min(df['lat']) - 1.0, np.max(df['lat']) + 1.0],datacrs)
ax2.coastlines()
ax2.add_feature(cfeature.STATES)
 
plt.show() 

# Plot the station locations

sys.exit()
# Grab the precip data for each station
time_sus  = df['valid'][df['station'] == 'SUS'].values
prcp_sus  = df['p01i'][df['station'] == 'SUS'].values
vsby_sus  = df['vsby'][df['station'] == 'SUS'].values
codes_sus = df['wxcodes'][df['station'] == 'SUS'].values

time_fyg  = df['valid'][df['station'] == 'FYG'].values
prcp_fyg  = df['p01i'][df['station'] == 'FYG'].values
vsby_fyg  = df['vsby'][df['station'] == 'FYG'].values
codes_fyg = df['wxcodes'][df['station'] == 'FYG'].values

# Mask any missing precip values
prcp_sus[(prcp_sus == 'M') | (prcp_sus == 'T')] = '0.000'
prcp_fyg[(prcp_fyg == 'M') | (prcp_fyg == 'T')] = '0.000'

prcp_sus = np.asfarray(prcp_sus,float)
prcp_fyg = np.asfarray(prcp_fyg,float)

# Convert times to datetime objects
# In the file, the times are in UTC, so convert to CDT (-5 hrs)
dtime_sus = [datetime.strptime(ttime,"%Y-%m-%d %H:%M") - timedelta(hours=5) \
#dtime_sus = [datetime.strptime(ttime,"%Y-%m-%d %H:%M") - timedelta(hours=5) \
        for ttime in time_sus]
dtime_fyg = [datetime.strptime(ttime,"%Y-%m-%d %H:%M") - timedelta(hours=5) \
#dtime_fyg = [datetime.strptime(ttime,"%Y-%m-%d %H:%M") - timedelta(hours=5) \
        for ttime in time_fyg]


##!## Loop over arrays and calculate 5-minute precip values
##!## from the 1-hour total precip values
##!## For SUS, the precip cycles at the 55 minute obs
##!#prcp_5min_sus = np.zeros(len(prcp_sus))
##!#prcp_5min_fyg = np.zeros(len(prcp_fyg))
##!#for pi in range(1,len(prcp_sus)):
##!#    if(prcp_sus[pi] < prcp_sus[pi-1]):
##!#    #if(dtime_sus[pi].minute == 55):
##!#        prcp_5min_sus[pi] = prcp_sus[pi]
##!#    else:
##!#        prcp_5min_sus[pi] = prcp_sus[pi] - prcp_sus[pi-1]
##!## For FYG, the precip cycles at the 15 minute obs
##!#for pi in range(1,len(prcp_fyg)):
##!#    if(prcp_fyg[pi] < prcp_fyg[pi-1]):
##!#    #if(dtime_fyg[pi].minute == 15):
##!#        prcp_5min_fyg[pi] = prcp_fyg[pi]
##!#    else:
##!#        prcp_5min_fyg[pi] = prcp_fyg[pi] - prcp_fyg[pi-1]


##!## Calculate cumulative sum of prcip
##!#prcp_csum_sus = np.cumsum(prcp_5min_sus)
##!#prcp_csum_fyg = np.cumsum(prcp_5min_fyg)
##!#
##!## Set up a line for the time of the incident
##!## Use local time of 3:15 pm CDT
##!#event_dtime = datetime.strptime("2017-04-29 15:15","%Y-%m-%d %H:%M") 

# Plot the data
# -------------
ax.plot(dtime_sus,prcp_csum_sus,label='[SUS] St.LOUIS/SPIRIT')
ax.plot(dtime_fyg,prcp_csum_fyg,label='[FYG] Washington')
#plt.plot(dtime_sus,prcp_5min_sus,label='SUS')
#plt.plot(dtime_fyg,prcp_5min_fyg,label='FYG')
# Plot a vertical line at the time of the crash
# ---------------------------------------------
ax.axvline(event_dtime,color='tab:red',lw=2,alpha=0.75,label='Crash time')
# Format the x axis
# -----------------
date_fmt = mdates.DateFormatter("%H:%M %b %d")
ax.xaxis.set_major_formatter(date_fmt)
fig1.autofmt_xdate()
# Format other plot characteristics
# ---------------------------------
ax.grid()
ax.legend()
ax.set_title("Regional Airport Precipitation Measurements")
ax.set_ylabel('Precipitation (inches)')
plt.savefig('precip_amounts.png',dpi=300)
plt.show()
