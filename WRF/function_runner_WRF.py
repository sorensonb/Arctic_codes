#!/usr/bin/env python

"""


"""

from WRFLib import *

infile = 'wrfout_d01_2015-06-28_00:00:00_nikkiupdate'
#infile = 'wrfout_d01_2015-06-28_00:00:00_mycontrol'
#data = Dataset(infile)

# Plot the 500 mb wind field
# --------------------------
#plot_500_hgt_wnd(data, 45)

# Plot a surface parameter
# ------------------------
#plot_sfc_field(data, 'T2', 43)
##plot_sfc_field(data, 'PSFC', 45)
#plot_sfc_field(data, 'BC1', 43)
#plot_sfc_field(data, 'SWDOWN', 45)

# Plots 4 panels (2-m temperature, surface downward SWF, 
# total column integrated smoke, cloud fraction) for two simulations:
# a control (clear) and smoky simulation.
#plot_WRF_combined_output(data, 43, out_add = 'clear', save = False)

plot_WRF_combined_output_2file('wrfout_d01_2015-06-28_00:00:00_mycontrol', \
    'wrfout_d01_2015-06-28_00:00:00_nikkiupdate', 43, \
    save = True)

#data.close()
sys.exit()

# Plot a sounding at a given town name, formatted City, State
# -----------------------------------------------------------

# Breckenridge?
lat = 45.1304
lon = -96.2314

# About St. Cloud
#lat = 45.0407
#lon = -94.9825

# Dense smoke. W of GF
#lat = 47.7269
#lon = -98.1803

#skew = SkewT(rotation=45)
#plot_WRF_sounding(data, 43, skew = skew, lat = lat, lon = lon, linestyle = ':')
#
#data.close()
#
#infile = 'wrfout_d01_2015-06-28_00:00:00_nikkiupdate'
#infile = 'wrfout_d01_2015-06-28_00:00:00_mycontrol'
#data = Dataset(infile)

#plot_WRF_sounding(data, 43, skew = skew, lat = lat, lon = lon, linestyle = '--')

#data.close()

plot_WRF_sounding_compare_data('wrfout_d01_2015-06-28_00:00:00_nikkiupdate', \
                               'wrfout_d01_2015-06-28_00:00:00_mycontrol', \
                                lat, lon, save = False)

##!#infile = 'wrfout_d01_2015-06-28_00:00:00_nikkiupdate'
##!##infile = 'wrfout_d01_2015-06-28_00:00:00_mycontrol'
##!#data = Dataset(infile)
##!#
##!#skew = SkewT(rotation=45)
##!#skew.plot(pres_prof,t_prof,'r', linestyle = ':')
##!#skew.plot(pres_prof,td_prof,'g', linestyle = ':')
##!#skew.plot_dry_adiabats()
##!#skew.plot_moist_adiabats()
##!#skew.plot_mixing_lines()
##!#skew.ax.set_ylim(1000,100)
##!#skew.ax.set_xlim(-40,50)
##!##plot_WRF_sounding(in_data,time_index, skew = None, city_name = None, \
##!##        lat = None, lon = None)

#data.close()

plt.show()
