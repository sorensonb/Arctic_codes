#!/usr/bin/env python

"""


"""

from MSPAM_Lib import *

# Source the source file
# ----------------------
os.system('source source_this')

#for ii in range(1000, 2100, 100):
#    print(ii)
#    run_MSPAM_yearly_yields(year = None, county = None, \
#        crop_type = 11, solar_units = ii, \
#        plant_month = 4, plant_day = 10, \
#        harvest_month = 9, harvest_day = 30, \
#        crop_density = 50,
#        output_unit = 'BPA')

#sys.exit()
#county = 'GRAND FORKS'
#
## Calculate M-SPAM yields for a single county, but for ranges of
## heat units and crop densities. Output is in an HDF5 file
## --------------------------------------------------------------
#run_MSPAM_yearly_auto_singlecounty(county, crop_type = 11, \
#    min_heat_unit = 1000, max_heat_unit = 2000, delta_heat_unit = 100, 
#    min_crop_density = 30, max_crop_density = 80, delta_crop_density = 10,
#    plant_month = 4, plant_day = 10, \
#    harvest_month = 9, harvest_day = 30, \
#    output_unit = 'BPA', \
#    REPROCESS = False)
#
#sys.exit()

# Read the HDF5 single-county results
# -----------------------------------
#infile = 'mspam_out_ctyGRAND_FORKS.hdf5'
infile = 'mspam_out_ctyGRAND_FORKS.hdf5'
mspam_gfk_dict = calc_model_errors_singlecounty(infile)

# Plot the yearly-averaged percent errors as a function of both the
# heat units and crop density for the given county
# -----------------------------------------------------------------
plot_MSPAM_errors_yearavg_singlecounty_2D(mspam_gfk_dict, save = True)




sys.exit()

# Load the observed yields, modeled yields, and percent errors for
# all counties in North Dakota, with calculations varying only
# the heat units.
# ------------------------------------------------------------
mspam_dict = calc_model_errors_allcounties(crop_density = 50)

# Plot a time series of wheat yields and M-SPAM modeled yields
# for a given county with a given heat unit
# ------------------------------------------------------------
plot_yield_time_series(mspam_dict, 'WALSH', 1700, save = True)

# Plot yearly-averaged M-SPAM errors as a function of heat units for
# a single county
# ------------------------------------------------------------------
#plot_MSPAM_errors_yearavg_single(mspam_dict, 'GRAND FORKS', save = False)

# Plot yearly-averaged M-SPAM errors as a function of heat units for 
# four counties
# ------------------------------------------------------------------
#plot_MSPAM_errors_yearavg_fourpanel(mspam_dict, 'GRAND FORKS', 'WALSH', \
#    'CASS', 'TRAILL', save = False)

# Plot yearly-averaged M-SPAM errors for all counties and heat units
# in a 2-d colormesh
# ------------------------------------------------------------------
#plot_MSPAM_errors_yearavg_2D(mspam_dict, save = False)
sys.exit()

# Code to automate the running of the M-SPAM model for a range of heat units
# --------------------------------------------------------------------------
for ii in range(1000, 2100, 100):
    print(ii)
    run_MSPAM_yearly_yields(year = None, county = None, \
        crop_type = 11, solar_units = ii, \
        plant_month = 4, plant_day = 10, \
        harvest_month = 9, harvest_day = 30, \
        output_unit = 'BPA')
