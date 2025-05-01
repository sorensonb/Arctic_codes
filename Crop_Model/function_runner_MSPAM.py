#!/usr/bin/env python

"""


"""

from MSPAM_Lib import *

# Source the source file
# ----------------------
os.system('source source_this')

# Load the observed yields, modeled yields, and percent errors
# ------------------------------------------------------------
mspam_dict = calc_model_errors()

# Plot yearly-averaged M-SPAM errors as a function of heat units for
# a single county
# ------------------------------------------------------------------
#plot_MSPAM_errors_yearavg_single(mspam_dict, 'GRAND FORKS', save = False)

# Plot yearly-averaged M-SPAM errors as a function of heat units for 
# four counties
# ------------------------------------------------------------------
plot_MSPAM_errors_yearavg_fourpanel(mspam_dict, 'GRAND FORKS', 'WALSH', \
    'BOWMAN', 'GRANT', save = True)
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
