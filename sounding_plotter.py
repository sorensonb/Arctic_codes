#!/usr/bin/env python
"""
  NAME:
    sounding_plotter.py

  PURPOSE:
    Plots the current sounding from a specified station using siphon and metpy

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>  - 2020/02/19
      Written (based on Wyoming_request.py and Skew-T_Layout.py)

"""

# Copyright (c) 2017 Siphon Contributors.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
"""
Wyoming Upper Air Data Request
==============================

This example shows how to use siphon's `simplewebswervice` support to create a query to
the Wyoming upper air archive.
"""

import sys
import matplotlib.pyplot as plt
from datetime import datetime
from metpy.units import units
from siphon.simplewebservice.wyoming import WyomingUpperAir

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pandas as pd

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, Hodograph, SkewT
from metpy.units import units

if(len(sys.argv)!=2):
    print("SYNTAX: python sounding_plotter.py stationID")
    print("        stationID: station identifier (ex., BIS, ABR")
    sys.exit()

station = sys.argv[1]

####################################################
# Create a datetime object for the sounding and string of the station identifier.
#date = datetime(2017, 9, 10, 12)
#station = 'BIS'
date = datetime.utcnow().replace(microsecond=0,second=0,minute=0)
if(date.hour>=12):
    print("retrieving for ",date)
    date = datetime.utcnow().replace(microsecond=0,second=0,minute=0,hour=12)
else:
    print("retrieving for ",date)
    date = datetime.utcnow().replace(microsecond=0,second=0,minute=0,hour=0)

####################################################
# Make the request (a pandas dataframe is returned).
df = WyomingUpperAir.request_data(date, station)

####################################################
# Inspect data columns in the dataframe.
print(df.columns)

####################################################
# Pull out a specific column of data.
print(df['pressure'])

####################################################
# Units are stored in a dictionary with the variable name as the key in the `units` attribute
# of the dataframe.
print(df.units)

####################################################
print(df.units['pressure'])

####################################################
# Units can then be attached to the values from the dataframe.
pressure = df['pressure'].values * units(df.units['pressure'])
temperature = df['temperature'].values * units(df.units['temperature'])
dewpoint = df['dewpoint'].values * units(df.units['dewpoint'])
u_wind = df['u_wind'].values * units(df.units['u_wind'])
v_wind = df['v_wind'].values * units(df.units['v_wind'])


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# This code below is taken from Skew-T_Layout.py from the metpy examples
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Create a new figure. The dimensions here give a good aspect ratio
fig = plt.figure(figsize=(11, 9))
add_metpy_logo(fig, 750, 80, size='large')

# Grid for plots
gs = gridspec.GridSpec(3, 3)
skew = SkewT(fig, rotation=45, subplot=gs[:, :2])

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(pressure, temperature, 'r')
skew.plot(pressure, dewpoint, 'g')
skew.plot_barbs(pressure, u_wind, v_wind)
skew.ax.set_ylim(1000, 100)

# Add the relevant special lines
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()

# Good bounds for aspect ratio
skew.ax.set_xlim(-40, 60)

# Create a hodograph
ax = fig.add_subplot(gs[0, -1])
h = Hodograph(ax, component_range=60.)
h.add_grid(increment=20)
h.plot(u_wind, v_wind)

# Add a title
skew.ax.set_title(station + ' '+ date.strftime("%d %b %Y %HZ"))

# Show the plot
plt.show()
