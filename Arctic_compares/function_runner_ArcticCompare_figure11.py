#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *
import random
#from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#                      BEGIN FORCING VERSION NN STUFF
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#                        PLOTTING/ANALYSIS FUNCTIONS
#
#  select_idxs(test_dict, ai_min, ice_min, ice_max:
#           This function identifies the indices that give the NN data
#           that meet the AI, ICE, COD, and SZA requirements
#
#  combine_NN_data(sim_name):
#           This function selects all the neural network output files 
#           for the given simulation , while also applying some QC criteria,
#           and returns a dictionary containing the observed and NN-estimated
#           SWF, in addition to the binning variables.
#
#  calc_NN_force_slope_intcpt(test_dict, ice_bin_edges, ...):
#           This function calculates the slope and intercepts for the 
#           forcing vs AI values. NOTE: here, forcing is calculated as NN - OBS.
#           
#           NOTE: This is also where the error stuff can be done by 
#                 switching the slope found by Theil-Sen (res[0])
#                 to either the lower bound slope (res[2]) or 
#                 upper bound slope (res[3])
#
#  test_calculate_type_forcing_v4:
#           Calculate the daily observation-estimated forcing
#           given the NN output-derived forcing slopes 
#           (and intercepts, if desired)
#
#  calculate_type_forcing_v4_monthly: calculates the single-day forcing 
#            estimates and averages them into monthly averages.
#
#
# calc_forcing_slopes_v4_all_months_arctic_avg_manyrefice(all_month_vals):
#           Calculates the trends (original and modified) for either
#           the refice or refcld simulations. The results are returned
#           in a dictionary.
#
#
# calc_forcing_slopes_v4_all_months_arctic_avg_uncert(all_month_files):
#           Same as above, but for the specified uncertainty runs.
#
#  plot_NN_architecture(save = False):
#           Plots a graphic showing the NN architecture
#
#  plot_NN_scatter_multiCOD(test_dict, cod_bin_edges,...):
#           This function plots the NN output for given COD bin edges
#           and for the specified AI, ICE, and SZA min/max values. For
#           each COD bin, the selected, scattered forcing 
#           (calculated as NN - OBS) data are plotted as a function of AI.
#           There is one plot for each COD bin.
#
#  plot_NN_bin_slopes(slope_dict):
#           Plot the binned NN/AI slopes for the 4 surface types as well
#           as the bin counts
#
#
#  plot_NN_bin_slopes_codmean(slope_dict, bin_dict, min_ob = 50, save = False):
#           Plots the NN/AI slopes for the 4 surface types averaged ALONG
#           the COD bins. Thus, the meaned forcing slopes are functions
#           of the COD.
#    
#  plot_NN_forcing_daily(...):
#           This function calculates the version 4 (NN) estimated forcing
#           for a single day based on the forcing slopes (and bins)
#           as well as the daily OMI AI, MODIS COD, NSIDC ICE, and
#           the calculated solar zenith angle for the day. 
#
#  plot_type_forcing_v4_all_months_arctic_avg(all_month_files, \)
#           Calculate Arctic-wide average for all of the provided monthly forcings.
#           Meant to be a plot of the Monte Carlo-esque uncertainty analysis.
#           The first file in the list is the control
#
# plot_type_forcing_v3_all_months_arctic_avg_combined(all_month_values, \)
#           Plots three Arctic-averaged monthly results in one image:
#           The first row is the 65 - 87, the second row is 65 - 75, and
#           the third row is 75 - 87.
#
# plot_compare_NN_output(sys.argv[1], save = False):
#           Makes a plot of Observed SWF, NN SWF, and Regression
#           Designed for showing the validation of the NN under
#           aerosol-free conditions
# 
# plot_compare_NN_output_v2(sys.argv[1], auto_zoom = False, save = False):
#           Plots OMI UVAI, MODIS true color, CERES obs, NN output,
#           and forcing (NN - obs) for either the whole Arctic
#           (auto_zoom = False) or just the plume area (auto_zoom = True)
# 
# plot_compare_NN_output_noaer(sys.argv[1], save = False):
#           Makes a plot of Observed SWF, NN SWF, and Regression
#           Designed for showing the validation of the NN under
#           aerosol-free conditions
# 
# plot_compare_NN_output_overlay(sys.argv[1], auto_zoom = True, save = False)
#           Plots zoomed in images of OMI AI, MODIS true color, and 
#           forcing, with hatched regions showing AI greater than xxx.
#           Function written to determine if there was significant drift
#           between the Aura and Aqua overpasses when viewing a smoke plume.
#           
# compare_nn_version_output(date_str, skip_version = skip_version, save = False)
#           This function reads in all NN output files for a given date_str,
#           calculates the mean and standard deviations of both the nn SWF
#           and the forcing (NN - obs), and plots the CERES obs, 
#           NN mean, and standard deviation of the NN values at each
#           swath point.
# 
#                        L2/L3 ERROR FUNCTIONS
#
# plot_L2_validate_regress_all(sim_name, slope_dict_lin, bin_dict,...
#           This function grabs all the NN output files for the 
#           given simualtion name (i.e. noland50), calculates L2-style
#           (NN - CERES) and L3-style (regressions & AI) forcings for
#           each aerosol-containing point in each swath and combines
#           them into arrays. These are also plotted using...
#
# plot_scatter_hist_L2_L3_errors(direct_forcings, calc_forcings,...
#           This function takes the direct (L2-style) and calculated 
#           (L3-style) forcings from plot_L2_validate_regress_all and
#           generates a scatter plot of the two, with error bars, and
#           a histogram of all the errors. 
#
#plot_NN_forcing_daily_L2L3_errors(date_str, daily_VSJ4, \
#    OMI_daily_VSJ4, slope_dict_lin, bin_dict, L2L3_err_mean, L2L3_err_std, \
#           This funtion calculates daily forcing values originally, and then
#           calculates (num_calcs) additional daily forcing values in which
#           each grid point with a forcing value calculated has L2L3 error
#           added that fits in a Gaussian distribution with the mean and 
#           standard deviation given above.
#
#                   STEPS FOR NEW DAILY FORCING ESTIMATION
#
#
# 1. Grab all the NN data (combine_NN_data)
# 2. Calculate forcing slopes and intercepts on the COD, SZA, and ICE bins
# 3. For each day, loop over each grid box
# 4. Use the COD, SZA, and ICE values to select the matching binned forcing 
#    slope and intercept. NOTE: This is where the "refice" and "refcld"
#    simulations can be easily redone.
# 5. Use the regression slope/intercept values to compute the forcing
#    given the daily AI at that grid box.
#
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# check_omi_ceres_time_offset(dir_list, num_points_per_file = 10):
#           This function grabs the OMI and CERES HDF5 files for given DTGs,
#           grabs 10 (or other specified number) of OMI points, finds the closest
#           CERES grid point to that OMI point, and figures out how many minutes
#           are between the OMI and CERES observations at that point.
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Load in all the NN-estimated SWF values (and other values) for
# simulation 'noland48'
# ----------------------------------------------------------------
#test_dict = combine_NN_data('noland48')
#test_dict = combine_NN_data('noland50')
#print("AS OF 2024/08/27, USING NOLAND72")
#print("AS OF 2024/09/08, USING NOLAND50 AGAIN")
#sim_name = 'noland48'
#sim_name = 'noland50'
#sim_name = 'noland72'
#sim_name = 'noland73'
#sim_name = 'noland74'
#sim_name = 'noland103'
sim_name = 'noland105'
#sim_name = 'noland75'
print("AS OF 2024/09/09, USING ", sim_name)
test_dict = combine_NN_data(sim_name)


#calc_pcnt_aerosol_over_type_v2(sim_name, 1.5, minlat = 70., dtype = 'PERT', ax = None)
#sys.exit()

# Identify the dictionary array indices that give the NN data
# associated with:
# - AI >= 1
# - 0 <= NSIDC ICE < 20
# - 2 <= MODIS COD < 4
# - 50 <= OMI SZA < 55
# -----------------------------------------------------------
#good_idxs = select_idxs(test_dict, 1, \
#    0, 20, 2, 4, 50, 55)

# Define COD, SZA, and NSIDC ICE bin ranges
# -----------------------------------------
cod_bin_edges = np.array([0,0.5,2,4,8,12,20,30,50])
cod_bin_means = (cod_bin_edges[1:] + cod_bin_edges[:-1]) / 2
sza_bin_edges = np.arange(40, 85, 5)
sza_bin_means = (sza_bin_edges[1:] + sza_bin_edges[:-1]) / 2
#ice_bin_edges = np.array([0, 20, 80, 100.2, 255])
ice_bin_edges = np.array([0, 20, 40, 60, 80, 100.2, 255])
ice_bin_means = (ice_bin_edges[1:] + ice_bin_edges[:-1]) / 2

bin_dict = {
    'cod_bin_edges': cod_bin_edges, \
    'cod_bin_means': cod_bin_means, \
    'sza_bin_edges': sza_bin_edges, \
    'sza_bin_means': sza_bin_means, \
    'ice_bin_edges': ice_bin_edges, \
    'ice_bin_means': ice_bin_means, \
}



# Calculate the regression slopes and intercepts from the NN data
# ---------------------------------------------------------------
min_ob = 50
ai_min_forslopes = 0.0
#slope_dict = calc_NN_force_slope_intcpt(test_dict, ice_bin_edges, \
#        sza_bin_edges, cod_bin_edges, ai_min = 0, min_ob = min_ob, \
#        trend_type = 'theil-sen')
slope_dict_lin = calc_NN_force_slope_intcpt(test_dict, ice_bin_edges, \
        sza_bin_edges, cod_bin_edges, ai_min = ai_min_forslopes, min_ob = min_ob, \
        trend_type = 'linregress')

# Plot the scattered NN output as function of AI for the bins given above
# -----------------------------------------------------------------------
#plot_NN_scatter_multiCOD(test_dict, cod_bin_edges, 0, 101, 255, 50, 55, save = False)
#sys.exit()

sza_min = 40
sza_max = 90
#plot_NN_scatter_combined_alltypes(test_dict, bin_dict, \
#    ai_min_forslopes, sza_min, sza_max, trend_type = 'linregress', \
#    show_specific_cod = 0, min_ai_for_stats = 2.0, \
#    plot_bounds = False, save = False)
#sys.exit()

# See what happens with differnet SZA bins. Is it necessary to bin
# by the SZA?
# -----------------------------------------------------------------
sfc_idx = 0
maxerr = 1.5

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Calculate L2- and L3-style forcing estimates for all high AI colocation 
# pixels
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#####sim_name = 'noland72'
#infile = 'validate_values_' + sim_name + '.hdf5'
infile = 'validate_values_' + sim_name + '_useintcptTrue.hdf5'
in_data = h5py.File(infile)
direct_forcings = in_data['direct_forcings'][:]
calc_forcings   = in_data['calc_forcings'][:]
in_data.close()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# TEST DAILY FORCING ESTIMATION 
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Read in the daily NSIDC, MODIS, and OMI data

begin_date = '200504'
end_date   = '202009'
season     = 'sunlight'
minlat = 65
maxlat = 87
shawn_file = home_dir + '/Research/OMI/omi_shawn_daily_2005_2020_v2.hdf5'
#jz_file    = home_dir + '/Research/OMI/omi_VJZ211_daily_2005_2020.hdf5'

# Is this the "clean-sky" background version?
OMI_daily_VSJ4  = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -20.00, max_AI = 0.7, minlat = minlat, maxlat = maxlat)
# Is this the OMI AI monthly trend version?
#OMI_daily_VSJ4  = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -0.10, max_AI = 20.0, minlat = minlat, maxlat = maxlat)

#OMI_daily_VJZ211 = calcOMI_MonthAvg_FromDaily(jz_file, min_AI = -0.10, minlat = minlat, maxlat = maxlat)
daily_VSJ4 = readOMI_daily_HDF5(shawn_file, minlat = minlat, maxlat = maxlat)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Test the L2L3_err stuff for a whole study period
#
# NOTE: RUNNING WITH filter_bad_vals SET TO False
# 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
ai_thresh = 0.7
maxerr = 1.5

##!## Control daily values
##!## --------------------
##!##daily_filename = 'arctic_daily_est_forcing_v1.hdf5'        
##!##daily_filename = 'arctic_daily_est_forcing_numsfcbins6.hdf5' # noland72
##!##daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v1.hdf5' # noland50
##!##daily_filename = 'arctic_daily_est_forcing_numsfcbins4.hdf5' # noland50
##!##daily_filename = 'arctic_daily_est_forcing_numsfcbins4_v1.hdf5' # noland72
##!#daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v2.hdf5' # noland74
##!##daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v3.hdf5' # noland103
##!#
##!## Daily values with ref_cld
##!## -------------------------
##!#refcld_filename = 'arctic_daily_est_forcing_numsfcbins6_refcld2005.hdf5' # noland74 , new error (doesn't matter)
##!#
##!## Daily values with ref_ice
##!## -------------------------
##!#refice_filename = 'arctic_daily_est_forcing_numsfcbins6_refice2005.hdf5' # noland74 , new error (doesn't matter)
##!#
##!## Daily values with ice modifiations
##!## ----------------------------------
##!##ice_filename = 'arctic_daily_est_forcing_iceerr_v1.hdf5'
##!##ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr.hdf5' # noland72
##!##ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v1.hdf5' # noland50
##!##ice_filename = 'arctic_daily_est_forcing_numsfcbins4_iceerr.hdf5' # noland50
##!##ice_filename = 'arctic_daily_est_forcing_numsfcbins4_iceerr_v1.hdf5' # noland72
##!#ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v2.hdf5' # noland74
##!##ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v3.hdf5' # noland103
##!#
##!## Daily values with COD modifiations
##!## ----------------------------------
##!##cod_filename = 'arctic_daily_est_forcing_coderr.hdf5'
##!##cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr.hdf5' # std = 5, noland72
##!##cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_v1.hdf5' # std = 5, noland50
##!##cod_filename = 'arctic_daily_est_forcing_numsfcbins4_coderr.hdf5' # std = 5, noland50
##!##cod_filename = 'arctic_daily_est_forcing_numsfcbins4_coderr_v1.hdf5' # std = 5, noland72
##!#cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_v2.hdf5' # std = 5, noland74
##!##cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_v3.hdf5' # std = 5, noland103


# HERERE:
#run_type = 'final'
#run_type = 'newerr'
#run_type = 'noback1'
run_type = 'noback2' # noback2: same as noback, but with the sign of the NN - CERES error flipped


# Control daily values
# --------------------
#daily_filename = 'arctic_daily_est_forcing_v1.hdf5'        
if(sim_name == 'noland74'):
    if(run_type == 'final'):
        daily_filename = 'arctic_daily_est_forcing_numsfcbins6_final.hdf5' # noland74, redone for validation
    else:
        daily_filename = 'arctic_daily_est_forcing_numsfcbins6.hdf5' # noland72
        daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v1.hdf5' # noland50
        daily_filename = 'arctic_daily_est_forcing_numsfcbins4.hdf5' # noland50
        daily_filename = 'arctic_daily_est_forcing_numsfcbins4_v1.hdf5' # noland72
        daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v2.hdf5' # noland74


    # Daily values with ref_cld
    # -------------------------
    refcld_filename = 'arctic_daily_est_forcing_numsfcbins6_refcld2005.hdf5' # noland74 , new error (doesn't matter)
    
    # Daily values with ref_ice
    # -------------------------
    refice_filename = 'arctic_daily_est_forcing_numsfcbins6_refice2005.hdf5' # noland74 , new error (doesn't matter)

    if(run_type == 'final'):
        ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_final.hdf5' # noland74, redone for validation
    else:
        # Daily values with ice modifiations
        # ----------------------------------
        #ice_filename = 'arctic_daily_est_forcing_iceerr_v1.hdf5'
        ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr.hdf5' # noland72
        ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v1.hdf5' # noland50
        ice_filename = 'arctic_daily_est_forcing_numsfcbins4_iceerr.hdf5' # noland50
        ice_filename = 'arctic_daily_est_forcing_numsfcbins4_iceerr_v1.hdf5' # noland72
        ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v2.hdf5' # noland74

    if(run_type == 'final'):
        cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_final.hdf5' # noland74, redone for validation
   
    else: 
        # Daily values with COD modifiations
        # ----------------------------------
        #cod_filename = 'arctic_daily_est_forcing_coderr.hdf5'
        cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr.hdf5' # std = 5, noland72
        cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_v1.hdf5' # std = 5, noland50
        cod_filename = 'arctic_daily_est_forcing_numsfcbins4_coderr.hdf5' # std = 5, noland50
        cod_filename = 'arctic_daily_est_forcing_numsfcbins4_coderr_v1.hdf5' # std = 5, noland72
        cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_v2.hdf5' # std = 5, noland74


elif(sim_name == 'noland103'):
    daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v3.hdf5' # noland103

    ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v3.hdf5' # noland103

    cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_v3.hdf5' # std = 5, noland103

elif(sim_name == 'noland105'):
    if(run_type == 'noback2'):
        daily_filename = 'arctic_daily_est_forcing_numsfcbins6_noland105_noback1.hdf5' # noland103 run with no background 

        ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_noland105_noback1.hdf5' # noland103

        cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_noland105_noback1.hdf5' # std = 5, noland103
    elif(run_type == 'noback1'):
        daily_filename = 'arctic_daily_est_forcing_numsfcbins6_noland105_noback1.hdf5' # noland103 run with no background 

        ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_noland105_noback1.hdf5' # noland103

        cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_noland105_noback1.hdf5' # std = 5, noland103
    elif(run_type == 'newerr'):
        daily_filename = 'arctic_daily_est_forcing_numsfcbins6_noland105_errtest.hdf5' # noland103 run with OMI_daily_data used
                                                                                       # with min_AI = -0.1, and max_AI = 20

        ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_noland105_errtest.hdf5' # noland103

        cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_noland105_errtest.hdf5' # std = 5, noland103
    else:
        daily_filename = 'arctic_daily_est_forcing_numsfcbins6_noland105.hdf5' # noland103

        ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_noland105.hdf5' # noland103

        cod_filename = 'arctic_daily_est_forcing_numsfcbins6_coderr_noland105.hdf5' # std = 5, noland103

else:
    print("INVALID SIM NAME")
    sys.exit()










#daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)
#ice_dict   = read_daily_month_force_L2L3_error_from_HDF5(ice_filename)
#cod_dict   = read_daily_month_force_L2L3_error_from_HDF5(cod_filename)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Functions for 4-panel plot with all analyzed errors
# - NN - obs
# - L2 - L3
# - Ice errors
# - COD errors
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#num_bins = 200
#plot_error_components_combined(direct_forcings, calc_forcings, sim_name, \
#    daily_filename, ice_filename, cod_filename, num_bins, \
#    astrofit = True, log_scale = True, save = False)
#sys.exit()


# Using 6 sfc bins, 
# ice errors: 
#   Applied ice errors: mean = 0, std dev = 15
#   Resulting forcing errors: mean = 0.011 W/m2, std dev = 3.11 W/m2
# cod errors: 
#   Applied ice errors: mean = 0, std dev = 5
#   Resulting forcing errors: mean = 0.848 W/m2, std dev = 15.9 W/m2
# L2L3 errors:
##    noland50
#       mean error (astrofit) = -2.57
#       stdv error (astrofit) = 22.2
#
# Combined error = sqrt(

# = = = = = =
#
# 2024/09/08
#
# USING 6 SURFACE BINS AND NOLAND50,
#
# Impacts of L3 ice errors on L3 forcing, from histogram
#   Applied ice errors: mean = 0, std dev = 15
#   Resulting forcing errors: mean = 0.031 W/m2, std dev = 3.24 W/m2
# Impacts of L3 COD errors on L3 forcing, from histogram
#   Applied COD errors: mean = 0, std dev = 5
#   Resulting forcing errors: mean = 1.012 W/m2, std dev = 14.374 W/m2
# L2L3 error numbers, from histogram
#   Mean L2L3 error, with astrofit = -2.57 
#   STDV L2L3 error, with astrofit = 22.2 
#
# Mean combined error = ice_error_mean + COD_error_mean + L2L3_error_mean = -1.527
# STD  combined error = sqrt(ice_error_std**2 + COD_error_std**2 + L2L3_error**2) = 26.64
# = = = = = =

# = = = = = =
#
# 2024/09/09
#
# USING 6 SURFACE BINS AND NOLAND74,
# 
# NN vs Obs errors, from histogram
#   Mean NN-ob error, with astrofit = 0.4
#   StDv NN-ob error, with astrofit = 17.3
# L2L3 error numbers, from histogram
#   Mean L2L3 error, with astrofit = -3.74 
#   STDV L2L3 error, with astrofit = 22.2 
# Impacts of L3 ice errors on L3 forcing, from histogram
#   Applied ice errors: mean = 0, std dev = 15
#   Resulting forcing errors: mean = 0.002 W/m2, std dev = 3.16 W/m2
# Impacts of L3 COD errors on L3 forcing, from histogram
#   Applied COD errors: mean = 0, std dev = 5
#   Resulting forcing errors: mean = 0.638 W/m2, std dev = 14.51 W/m2
#
# Mean combined error = ice_error_mean + COD_error_mean + L2L3_error_mean = -3.1
# STD  combined error = sqrt(ice_error_std**2 + COD_error_std**2 + L2L3_error**2) = 26.70
#
# With NN-ob errors,
# Mean combined error = nn_ob_error_mean + ice_error_mean + COD_error_mean + L2L3_error_mean = -2.7
# STD  combined error = sqrt(nn_ob_error_mean**2 + ice_error_std**2 + COD_error_std**2 + L2L3_error**2) = 31.82
#
# = = = = = =


# = = = = =
#
# TEST BULK ERROR ADDITIONS
#
# = = = = =

# NEW ERROR VALUES AS OF 2024/10/30. These use noland103 and account for
# use_intercepts in L2_L3 errors. The "...count300_noland103.hdf5" and 
# "...count300_noland103_v1.hdf5" files were generated with these new error 
# values.
# -------------------------------------------------------------------------
if(sim_name == 'noland105'):
    if(run_type == 'noback2'):
        total_err_mean = 0.8
        total_err_std = 31.61
    elif(run_type == 'noback1'):
        total_err_mean = -2.0
        total_err_std = 31.61
    else:
        total_err_mean = -2.2
        total_err_std = 29.76

elif(sim_name == 'noland103'):
    total_err_mean = 0.1
    total_err_std = 30.38

elif(sim_name == 'noland74'):
    if(run_type == 'final'):
        # VALIDATION (_final) ERROR VALUES AS OF 2024/11/04, accounting for use_intercepts in the L2_L3 
        # errors. The "...count300_final.hdf5" and "...count300_final_v1.hdf5"
        # files were generated with these new error values. These also are for noland74
        # -----------------------------------------------------------------------------
        total_err_mean = -0.1
        total_err_std = 29.13

    else:
        # NEW ERROR VALUES AS OF 2024/10/29, accounting for use_intercepts in the L2_L3 
        # errors. The "...count300_newerror.hdf5" and "...count300_newerror_v1.hdf5"
        # files were generated with these new error values. These also are for noland74
        # -----------------------------------------------------------------------------
        total_err_mean = 0.1
        total_err_std = 30.93
        
        # OLD ERROR VALUES: for noland74
        #total_err_mean = -2.7
        #total_err_std = 31.82

num_sims = 300
#num_sims = 300
minlat = 65.5
maxlat = 87.5
sim_values = None

daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)

# = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = =

read_force_sim_vals = False
save_force_vals = False
read_trend_sim_vals = True 
save_trend_vals = False

calc_region_avg_force_vals = False
calc_force_trends_from_file = False

if(read_trend_sim_vals):
   
    # NOTE: Only using 2 files here. Can change ":2" to allow it to read more files 
    # -----------------------------------------------------------------------------
    if(sim_name == 'noland105'):
        if(run_type == 'noback2'):
            file_start = 'arctic_monthly_force_trends_count300_noland105_noback2'  # noland105, but forcing calced without background
        elif(run_type == 'noback1'):
            file_start = 'arctic_monthly_force_trends_count300_noland105_noback1'  # noland105, but forcing calced without background
        elif(run_type == 'newerr'):
            file_start = 'arctic_monthly_force_trends_count300_noland105_errtest'  # noland105, but with OMI_daily_data used
                                                                                       # with min_AI = -0.1, and max_AI = 20
        else:
            file_start = 'arctic_monthly_force_trends_count300_noland105'  # noland103
    elif(sim_name == 'noland103'):
        file_start = 'arctic_monthly_force_trends_count300_noland103'  # noland103, correct errors
    elif(sim_name == 'noland74'):
        if(run_type == 'final'):
            file_start = 'arctic_monthly_force_trends_count300_final'  # noland74, final check
        else:
            file_start = 'arctic_monthly_force_trends_count300_newerror'  # noland74, correct errors
            #file_start = 'arctic_monthly_force_trends_count'
    force_files = glob(file_start + '*.hdf5')[:2]

    # Figure out how many simulations are in all the files
    # ----------------------------------------------------
    total_sims = 0
    for ffile in force_files:
        print(ffile)
        local_numsim = ffile.strip().split('/')[-1].split('_')[4].split('_')[0].split('.')[0][5:]   
        total_sims += int(local_numsim)

    print("Total number of sims:", total_sims)

    # Create an array to hold the values (Assume 96 months)
    # -----------------------------------------------------
    forcing_trends = np.full( (total_sims, 6, daily_dict['latitude'].shape[0], daily_dict['longitude'].shape[0]), np.nan)

    # Insert the values into the array
    # --------------------------------
    beg_idx = 0
    end_idx = 0
    for ii in range(len(force_files)):
        local_numsim = int(force_files[ii].strip().split('/')[-1].split('_')[4].split('_')[0].split('.')[0][5:])
        data = h5py.File(force_files[ii])
        end_idx = beg_idx + local_numsim
        forcing_trends[beg_idx:end_idx,:,:,:] = data['force_trends'][:,:,:,:]
        data.close()
        print(beg_idx, end_idx)
        beg_idx = end_idx
    
else:
    if(read_force_sim_vals):
        
        # Grab the filenames
        # NOTE: Only using 2 files here. Can change ":2" to allow it to read more files 
        # -----------------------------------------------------------------------------
        if(sim_name == 'noland105'):
            if(run_type == 'noback2'):
                file_start = 'arctic_monthly_force_values_count300_noland105_noback2'
            elif(run_type == 'noback1'):
                file_start = 'arctic_monthly_force_values_count300_noland105_noback1'
            elif(run_type == 'newerr'):
                file_start = 'arctic_monthly_force_values_count300_noland105_errtest'
            else:
                file_start = 'arctic_monthly_force_values_count300_noland105'
        elif(sim_name == 'noland103'):
            file_start = 'arctic_monthly_force_values_count300_noland103'
        elif(sim_name == 'noland74'):
            if(run_type == 'final'):
                file_start = 'arctic_monthly_force_values_count300_final' # noland74, final check
            else:
                file_start = 'arctic_monthly_force_values_count300_newerror' # noland74, correct errors
            file_start = 'arctic_monthly_force_values_count'

        force_files = glob(file_start +'*.hdf5')[:2]
    
        # Figure out how many simulations are in all the files
        # ----------------------------------------------------
        total_sims = 0
        for ffile in force_files:
            print(ffile)
            local_numsim = ffile.strip().split('/')[-1].split('_')[4].split('_')[0].split('.')[0][5:]   
            total_sims += int(local_numsim)
    
        print("Total number of sims:", total_sims)
   
        if(calc_region_avg_force_vals): 
            # Create an array to hold the region-averaged values 
            # (Assume 96 months & 3 regions (65 - 90, 65 - 75, and 75 - 90)
            # -------------------------------------------------------------
            sim_values = np.full( (total_sims, 3, 96), np.nan)
        else:
            # Create an array to hold the values (Assume 96 months)
            # -----------------------------------------------------
            sim_values = np.full( (total_sims, 96, daily_dict['latitude'].shape[0], daily_dict['longitude'].shape[0]), np.nan)
    
        # Insert the values into the array
        # --------------------------------
        beg_idx = 0
        end_idx = 0
        for ii in range(len(force_files)):
            local_numsim = int(force_files[ii].strip().split('/')[-1].split('_')[4].split('_')[0].split('.')[0][5:])
            data = h5py.File(force_files[ii])
            end_idx = beg_idx + local_numsim
            
            if(calc_region_avg_force_vals):
                lat_mins = [65.5, 65.5, 75.5]
                lat_maxs = [89.5, 75.5, 89.5]
                for jj in range(len(lat_mins)):
                    lat_idxs = np.where( (data['latitude'][:] >= lat_mins[jj]) & (data['latitude'][:] < lat_maxs[jj])) 
                    lat_beg_idx = lat_idxs[0][0]
                    lat_end_idx = lat_idxs[0][-1] + 1
                    sim_values[beg_idx:end_idx,jj,:] = \
                        np.nanmean(data['monthly_force_vals'][:,:,lat_beg_idx:lat_end_idx,:], axis = (2,3))
            else:
                sim_values[beg_idx:end_idx,:,:,:] = data['monthly_force_vals'][:,:,:,:]

            data.close()
            print(beg_idx, end_idx)
            beg_idx = end_idx
        
    
    else:
        #for nn in range(2):
            
        return_std = False
        if(return_std):
            sim_values, std_values = calc_force_vals_bulk(daily_dict, total_err_mean, total_err_std, \
                minlat = minlat, maxlat = maxlat, num_sims = num_sims, return_std = True)
        else:
            sim_values = calc_force_vals_bulk(daily_dict, total_err_mean, total_err_std, \
                minlat = minlat, maxlat = maxlat, num_sims = num_sims, return_std = False)
            #sim_values_2 = calc_force_vals_bulk(daily_dict, total_err_mean, total_err_std, \
            #    minlat = minlat, maxlat = maxlat, num_sims = num_sims, return_std = False)
    
        if(save_force_vals):
            # Save error monthly forcing values to an output file
            # --------------------------------------------------- 
            if(sim_name == 'noland74'):
                if(run_type == 'final'):
                    name_add = '_final'
                else:
                    name_add = '_newerr'
            elif(sim_name == 'noland103'):
                name_add = '_noland103'

            elif(sim_name == 'noland105'):
                if(run_type == 'noback2'):
                    name_add = '_noland105_noback2'
                elif(run_type == 'noback1'):
                    name_add = '_noland105_noback1'
                elif(run_type == 'newerr'):
                    name_add = '_noland105_errtest'
                else:
                    name_add = '_noland105'

            write_monthly_force_vals_sims_to_HDF5(daily_dict, sim_values, \
                total_err_mean, total_err_std, save_path = './', name_add = name_add)
                #total_err_mean, total_err_std, save_path = './', name_add = '_final') #noland74, final check
                #total_err_mean, total_err_std, save_path = './', name_add = '_newerror') #noland74, correct errors
                #total_err_mean, total_err_std, save_path = './', name_add = '_noland103')
    
    
        # Test calculating trends across all calculations and across the whole grid
        # -------------------------------------------------------------------------
        forcing_trends = np.full( (sim_values.shape[0], 6, sim_values.shape[2], sim_values.shape[3]), np.nan)
        #forcing_pvals  = np.full( (sim_values.shape[0], 6, sim_values.shape[2], sim_values.shape[3]), np.nan)
        #forcing_uncert = np.full( (sim_values.shape[0], 6, sim_values.shape[2], sim_values.shape[3]), np.nan)
        
        for ii in range(6):
            for jj in range(sim_values.shape[0]):
                print(ii,jj)
                forcing_trends[jj,ii,:,:], _, \
                    _ = calc_forcing_grid_trend(\
                    sim_values[jj,ii::6,:,:], 'standard')
                #forcing_trends[jj,ii,:,:], forcing_pvals[jj,ii,:,:], \
                #    forcing_uncert[jj,ii,:,:] = calc_forcing_grid_trend(\
                #    sim_values[jj,ii::6,:,:], 'standard')
        
        
        if(save_force_vals): 
            # Save error trends to an output file
            # -----------------------------------
            if(sim_name == 'noland74'):
                if(run_type == 'final'):
                    name_add = '_final'
                else:
                    name_add = '_newerr'
            elif(sim_name == 'noland103'):
                name_add = '_noland103'

            elif(sim_name == 'noland105'):
                if(run_type == 'noback2'):
                    name_add = '_noland105_noback2'
                elif(run_type == 'noback1'):
                    name_add = '_noland105_noback1'
                elif(run_type == 'newerr'):
                    name_add = '_noland105_errtest'
                else:
                    name_add = '_noland105'

            write_monthly_force_trend_sims_to_HDF5(daily_dict, forcing_trends, \
                total_err_mean, total_err_std, save_path = './', name_add = name_add)
                #total_err_mean, total_err_std, save_path = './', name_add = '_final') #noland74, final check
                #total_err_mean, total_err_std, save_path = './', name_add = '_newerror') #noland74, correct errors
                #total_err_mean, total_err_std, save_path = './', name_add = '_noland103')
        
        #del(sim_values)
        #del(forcing_trends) 

# Plot the distribution of trend estimates at a lat/lon idx and month
# -------------------------------------------------------------------
test_error_dist(daily_dict, forcing_trends, 4, 1, 287, 20, \
    sim_name = sim_name, run_type = run_type, save = True)
sys.exit()
#test_error_dist(daily_dict, OMI_data, forcing_trends, 3, 10, 340, 30)
OMI_data = calcOMI_MonthAvg_FromDaily(shawn_file, \
    min_AI = 0, max_AI = 20, minlat = 65.5, maxlat = 90.5)
#plot_many_trend_test(daily_dict, OMI_data, forcing_trends, 3, 10, 340, \
plot_many_trend_test(daily_dict, OMI_data, forcing_trends, 4, 2, 285, \
    30, conf_level = 90, save = False)
sys.exit()

