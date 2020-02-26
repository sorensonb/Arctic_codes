#!/usr/bin/env python

# Copyright (c) 2013-2015 Siphon Contributors.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
"""
================
NCSS Time Series
================

Use Siphon to query the NetCDF Subset Service for a timeseries.
"""
from datetime import datetime, timedelta

import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import num2date
import numpy as np
from siphon.catalog import TDSCatalog
from metpy.calc import dewpoint_from_relative_humidity,wind_direction,wind_speed
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo
from metpy.units import units

if(len(sys.argv) != 4):
    print("SYNTAX: python model_meteogram.py lat lon model")
    print("        GFK: 47.9253, -97.0329")
    print("        model: GFS, GEFS, NAM, RAP, HRRR")
    sys.exit()

# Make meteogram plot
class Meteogram(object):
    """ Plot a time series of meteorological data from a particular station as a
    meteogram with standard variables to visualize, including thermodynamic,
    kinematic, and pressure. The functions below control the plotting of each
    variable.
    TO DO: Make the subplot creation dynamic so the number of rows is not
    static as it is currently. """

    def __init__(self, fig, dates, probeid, time=None, axis=0):
        """
        Required input:
            fig: figure object
            dates: array of dates corresponding to the data
            probeid: ID of the station
        Optional Input:
            time: Time the data is to be plotted
            axis: number that controls the new axis to be plotted (FOR FUTURE)
        """
        if not time:
            time = datetime.utcnow()
        self.start = dates[0]
        self.fig = fig
        self.end = dates[-1]
        self.axis_num = 0
        self.dates = mpl.dates.date2num(dates)
        self.time = time.strftime('%Y-%m-%d %H:%M UTC')
        self.title = 'Latest Ob Time: {0}\nProbe ID: {1}'.format(self.time, probeid)

    def plot_winds(self, ws, wd, model, plot_range=None):
    #def plot_winds(self, ws, wd, wsmax, plot_range=None):
        """
        Required input:
            ws: Wind speeds (knots)
            wd: Wind direction (degrees)
            wsmax: Wind gust (knots)
        Optional Input:
            plot_range: Data range for making figure (list of (min,max,step))
        """
        # PLOT WIND SPEED AND WIND DIRECTION
        self.ax1 = fig.add_subplot(4, 1, 1)
        ln1 = self.ax1.plot(self.dates, ws, label='Wind Speed')
        #self.ax1.fill_between(self.dates, ws, 0)
        self.ax1.set_xlim(self.start, self.end)
        if not plot_range:
            plot_range = [0, 20, 1]
        self.ax1.set_ylabel('Wind Speed (knots)', multialignment='center')
        #self.ax1.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        self.ax1.grid(b=True, which='major', axis='y', color='k', linestyle='--',
                      linewidth=0.5)
        #ln2 = self.ax1.plot(self.dates, wsmax, '.r', label='3-sec Wind Speed Max')

        ax7 = self.ax1.twinx()
        ln3 = ax7.plot(self.dates, wd, '.k', linewidth=0.5, label='Wind Direction')
        ax7.set_ylabel('Wind\nDirection\n(degrees)', multialignment='center')
        ax7.set_ylim(0, 360)
        ax7.set_yticks(np.arange(45, 405, 90), ['NE', 'SE', 'SW', 'NW'])
        lns = ln1 + ln3
        #lns = ln1 + ln2 + ln3
        labs = [l.get_label() for l in lns]
        ax7.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H UTC'))
        if(model!='GEFS'):
            ax7.legend(lns, labs, loc='upper center',
                       bbox_to_anchor=(0.5, 1.2), ncol=3, prop={'size': 12})

    def plot_thermo(self, t, td, model, plot_range=None):
        """
        Required input:
            T: Temperature (deg F)
            TD: Dewpoint (deg F)
        Optional Input:
            plot_range: Data range for making figure (list of (min,max,step))
        """
        # PLOT TEMPERATURE AND DEWPOINT
        #if not plot_range:
        #    plot_range = [10, 90, 2]
        self.ax2 = fig.add_subplot(4, 1, 2, sharex=self.ax1)
        ln4 = self.ax2.plot(self.dates, t, 'r-', label='Temperature')
        #self.ax2.fill_between(self.dates, t, td, color='r')

        self.ax2.set_ylabel('Temperature\n(F)', multialignment='center')
        self.ax2.grid(b=True, which='major', axis='y', color='k', linestyle='--',
                      linewidth=0.5)
        #self.ax2.set_ylim(plot_range[0], plot_range[1], plot_range[2])

        ln5 = self.ax2.plot(self.dates, td, 'g-', label='Dewpoint')
        #self.ax2.fill_between(self.dates, td, self.ax2.get_ylim()[0], color='g')

        ax_twin = self.ax2.twinx()
        #ax_twin.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        lns = ln4 + ln5
        labs = [l.get_label() for l in lns]
        ax_twin.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H UTC'))

        if(model!='GEFS'):
            self.ax2.legend(lns, labs, loc='upper center',
                            bbox_to_anchor=(0.5, 1.2), ncol=2, prop={'size': 12})

    def plot_rh(self, rh, model, plot_range=None):
        """
        Required input:
            RH: Relative humidity (%)
        Optional Input:
            plot_range: Data range for making figure (list of (min,max,step))
        """
        # PLOT RELATIVE HUMIDITY
        if not plot_range:
            plot_range = [0, 100, 4]
        self.ax3 = fig.add_subplot(4, 1, 3, sharex=self.ax1)
        self.ax3.plot(self.dates, rh, 'g-', label='Relative Humidity')
        if(model!='GEFS'):
            self.ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22), prop={'size': 12})
        self.ax3.grid(b=True, which='major', axis='y', color='k', linestyle='--',
                      linewidth=0.5)
        #self.ax3.set_ylim(plot_range[0], plot_range[1], plot_range[2])

        #self.ax3.fill_between(self.dates, rh, self.ax3.get_ylim()[0], color='g')
        self.ax3.set_ylabel('Relative Humidity\n(%)', multialignment='center')
        self.ax3.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H UTC'))
        axtwin = self.ax3.twinx()
        #axtwin.set_ylim(plot_range[0], plot_range[1], plot_range[2])

    def plot_pressure(self, p, model, plot_range=None):
        """
        Required input:
            P: Mean Sea Level Pressure (hPa)
        Optional Input:
            plot_range: Data range for making figure (list of (min,max,step))
        """
        # PLOT PRESSURE
        #if not plot_range:
        #    plot_range = [970, 1030, 2]
        self.ax4 = fig.add_subplot(4, 1, 4, sharex=self.ax1)
        self.ax4.plot(self.dates, p, 'm', label='Mean Sea Level Pressure')
        self.ax4.set_ylabel('Mean Sea\nLevel Pressure\n(mb)', multialignment='center')
        #self.ax4.set_ylim(plot_range[0], plot_range[1], plot_range[2])

        axtwin = self.ax4.twinx()
        #axtwin.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        #axtwin.fill_between(self.dates, p, axtwin.get_ylim()[0], color='m')
        axtwin.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H UTC'))

        if(model!='GEFS'):
            self.ax4.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), prop={'size': 12})
        self.ax4.grid(b=True, which='major', axis='y', color='k', linestyle='--',
                      linewidth=0.5)
        # OTHER OPTIONAL AXES TO PLOT
        # plot_irradiance
        # plot_precipitation
    def plot_precipitation(self, precip, model, plot_range=None):
        """
        Required input:
            precip: Precipitation
        Optional Input:
            plot_range: Data range for making figure (list of (min,max,step))
        """
        # PLOT RELATIVE HUMIDITY
        if not plot_range:
            plot_range = [0, 100, 4]
        self.ax3 = fig.add_subplot(4, 1, 3, sharex=self.ax1)
        self.ax3.plot(self.dates, precip, 'g-', label='Precipitation')
        if(model!='GEFS'):
            self.ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22), prop={'size': 12})
        self.ax3.grid(b=True, which='major', axis='y', color='k', linestyle='--',
                      linewidth=0.5)
        #self.ax3.set_ylim(plot_range[0], plot_range[1], plot_range[2])

        #self.ax3.fill_between(self.dates, rh, self.ax3.get_ylim()[0], color='g')
        self.ax3.set_ylabel('Precipitation Accumulation\n(kg m-2 s-1)', multialignment='center')
        self.ax3.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H UTC'))
        axtwin = self.ax3.twinx()
        #axtwin.set_ylim(plot_range[0], plot_range[1], plot_range[2])


# Create dictionary to hold all 
var_names = {
    'GFS': {
        'precip': 'Total_precipitation_surface_Mixed_intervals_Accumulation',
        'temp':   'Temperature_isobaric',
        'pmsl':   'Pressure_reduced_to_MSL_msl',
        'rh':     'Relative_humidity_isobaric',
        'uwind':  'u-component_of_wind_isobaric',
        'vwind':  'v-component_of_wind_isobaric',
     },
    'GEFS': {
        'precip': 'Total_precipitation_surface_6_Hour_Accumulation_ens',
        'temp':   'Temperature_isobaric_ens',
        'pmsl':   'Pressure_reduced_to_MSL_msl_ens',
        'rh':     'Relative_humidity_isobaric_ens',
        'uwind':  'u-component_of_wind_isobaric_ens',
        'vwind':  'v-component_of_wind_isobaric_ens',
     },
    'NAM': {
        'precip': 'Total_precipitation_surface_3_Hour_Accumulation',
        'temp':   'Temperature_isobaric',
        'pmsl':   'Pressure_reduced_to_MSL_msl',
        'rh':     'Relative_humidity_isobaric',
        'uwind':  'u-component_of_wind_isobaric',
        'vwind':  'v-component_of_wind_isobaric',
     },
    'RAP': {
        'precip': 'Total_precipitation_surface_1_Hour_Accumulation',
        'temp':   'Temperature_isobaric',
        'pmsl':   'Pressure_surface',
        'rh':     'Relative_humidity_isobaric',
        'uwind':  'u-component_of_wind_isobaric',
        'vwind':  'v-component_of_wind_isobaric',
     }
##    'HRRR': {
##        'precip': 'Total_precipitation_surface_1_Hour_Accumulation',
##        'temp':   'Temperature_isobaric',
##        'pmsl':   'Pressure_reduced_to_MSL_msl',
##        'rh':     'Relative_humidity_isobaric',
##        'uwind':  'u-component_of_wind_isobaric',
##        'vwind':  'v-component_of_wind_isobaric',
##     }
}

in_lat = float(sys.argv[1])
in_lon = float(sys.argv[2])
model = sys.argv[3]
if(model=='GFS'):
    data_link = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/' + \
                      'Global_0p5deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p5deg/Best'
elif(model=='GEFS'):
    data_link = 'https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GEFS/Global_1p0deg_Ensemble/' + \
                        'members/catalog.html?dataset=grib/NCEP/GEFS/Global_1p0deg_Ensemble/members/Best'
elif(model=='NAM'):
    data_link = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/NAM/' + \
                        'CONUS_12km/catalog.html?dataset=grib/NCEP/NAM/CONUS_12km/Best'
elif(model=='RAP'):
    data_link = 'https://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/' + \
                        'catalog.html?dataset=grib/NCEP/RAP/CONUS_13km/Best'
elif(model=='HRRR'):
    data_link = 'https://thredds.ucar.edu/thredds/catalog/grib/NCEP/HRRR/' + \
                        'CONUS_2p5km/catalog.html?dataset=grib/NCEP/HRRR/CONUS_2p5km/Best'
    

###########################################
# First we construct a TDSCatalog instance pointing to our dataset of interest, in
# this case TDS' "Best" virtual dataset for the GFS global 0.5 degree collection of
# GRIB files. We see this catalog contains a single dataset.
best_gfs = TDSCatalog(data_link)
print(best_gfs.datasets)

# Grab the NAM data

###########################################
# We pull out this dataset and get the NCSS access point
best_ds = best_gfs.datasets[0]
ncss = best_ds.subset()

###########################################
# We can then use the `ncss` object to create a new query object, which
# facilitates asking for data from the server.
query = ncss.query()

###########################################
# We construct a query asking for data corresponding to latitude 40N and longitude 105W,
# for the next 7 days. We also ask for NetCDF version 4 data, for the variable
# 'Temperature_isobaric', at the vertical level of 100000 Pa (approximately surface).
# This request will return all times in the range for a single point. Note the string
# representation of the query is a properly encoded query string.
now = datetime.utcnow()
query.lonlat_point(in_lon, in_lat).vertical_level(100000).time_range(now, now + timedelta(days=7))
query.variables(var_names[model]['temp'],
                var_names[model]['pmsl'],
                var_names[model]['rh'],
                var_names[model]['uwind'],
                var_names[model]['vwind'],
                var_names[model]['precip']
                ).accept('netcdf')

###########################################
# We now request data from the server using this query. The `NCSS` class handles parsing
# this NetCDF data (using the `netCDF4` module). If we print out the variable names, we
# see our requested variables, as well as a few others (more metadata information)
data = ncss.get_data(query)
list(data.variables)

###########################################
# We'll pull out the temperature  and time variables.
temp = data.variables[var_names[model]['temp']]
pmsl = data.variables[var_names[model]['pmsl']]
rh   = data.variables[var_names[model]['rh'] ]
uwind= data.variables[var_names[model]['uwind']]
vwind= data.variables[var_names[model]['vwind']]
prcip= data.variables[var_names[model]['precip']]
time = data.variables['time']

###########################################
# The time values are in hours relative to the start of the entire model collection.
# Fortunately, the `netCDF4` module has a helper function to convert these numbers into
# Python `datetime` objects. We can see the first 5 element output by the function look
# reasonable.
time_vals = num2date(time[:].squeeze(), time.units)
print(time_vals[:5])

##############################################
#### Now we can plot these up using matplotlib, which has ready-made support for `datetime`
#### objects.
###fig, ax = plt.subplots(1, 1, figsize=(9, 8))
###ax.plot(time_vals, temp[:].squeeze(), 'r', linewidth=2)
###ax.plot(time_vals, dp[:].squeeze(), 'g', linewidth=2)
###ax.set_ylabel('{} ({})'.format(temp.standard_name, temp.units))
###ax.set_xlabel('Forecast Time (UTC)')
###ax.grid(True)
###plt.show()

# Calculate wind speed and direction from UV wind
ws = wind_speed(np.array(uwind) * units('m/s'),np.array(vwind) * units('m/s'))
wd = wind_direction(np.array(uwind) * units('m/s'),np.array(vwind) * units('m/s'))

# ID For Plotting on Meteogram
probe_id = '0102A'

ddata = {'wind_speed': (np.array(ws).squeeze() * units('m/s')).to(units('knots')),
        #'wind_speed_max': (np.array(wsmax) * units('m/s')).to(units('knots')),
        'wind_direction': np.array(wd).squeeze() * units('degrees'),
        'dewpoint': dewpoint_from_relative_humidity((np.array(temp).squeeze() * units.K),
                                                    np.array(rh).squeeze() / 100.).to(units('degF')),
        'air_temperature': (np.array(temp).squeeze() * units('K')).to(units('degF')),
        'mean_slp': (np.array(pmsl).squeeze() * units('Pa')).to(units('hPa')),
        #'mean_slp': calc_mslp(np.array(temp), np.array(pres), hgt_example) * units('hPa'),
        'relative_humidity': np.array(rh).squeeze(), 
        'precipitation': np.array(prcip).squeeze(), 
        'times': np.array(time_vals)}

fig = plt.figure(figsize=(20, 16))
add_metpy_logo(fig, 250, 180)
meteogram = Meteogram(fig, ddata['times'], probe_id)
meteogram.plot_winds(ddata['wind_speed'], ddata['wind_direction'],model)
#meteogram.plot_winds(data['wind_speed'], data['wind_direction'], data['wind_speed_max'])
meteogram.plot_thermo(ddata['air_temperature'], ddata['dewpoint'],model)
meteogram.plot_precipitation(ddata['precipitation'],model)
#meteogram.plot_rh(ddata['relative_humidity'])
meteogram.plot_pressure(ddata['mean_slp'],model)
fig.subplots_adjust(hspace=0.5)
plt.suptitle(model)
plt.show()




