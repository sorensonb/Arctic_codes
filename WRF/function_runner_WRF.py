#!/usr/bin/env python

"""


"""

from WRFLib import *

#infile = 'wrfout_d01_2015-06-28_00:00:00_nikkiupdate'
infile = 'wrfout_d01_2015-06-28_00:00:00_mycontrol'
data = Dataset(infile)

# Plot the 500 mb wind field
# --------------------------
#plot_500_hgt_wnd(data, 45)

# Plot a surface parameter
# ------------------------
#plot_sfc_field(data, 'T2', 43)
##plot_sfc_field(data, 'PSFC', 45)
#plot_sfc_field(data, 'BC1', 43)
#plot_sfc_field(data, 'SWDOWN', 45)

plot_WRF_combined_output(data, 43, out_add = 'clear', save = True)

# Plot a sounding at a given town name, formatted City, State
# -----------------------------------------------------------
#plot_WRF_sounding(data, 41, 'Fargo, North Dakota')

data.close()
