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
sim_name = 'noland74'
#sim_name = 'noland103'
#sim_name = 'noland104'
#sim_name = 'noland75'
print("AS OF 2024/09/09, USING ", sim_name)
test_dict = combine_NN_data(sim_name)

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

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Calculate L2- and L3-style forcing estimates for all high AI colocation 
# pixels
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
infile = 'validate_values_' + sim_name + '_useintcptTrue.hdf5'
in_data = h5py.File(infile)
direct_forcings = in_data['direct_forcings'][:]
calc_forcings   = in_data['calc_forcings'][:]

#direct_forcings, calc_forcings = \
#    plot_L2_validate_regress_all(sim_name, slope_dict_lin, bin_dict, \
#    ai_thresh = 0.7, mod_slopes = None, mod_intercepts = None, \
#    mod_cod = None, mod_ice = None, use_intercept = False, \
#    min_cod = None, max_cod = None,\
#    min_sza = None, max_sza = None, \
#    min_ice = None, max_ice = None, \
#    save = False, return_values = True, save_values = False)
###write_L2_L3_validation_values(direct_forcings, calc_forcings, sim_name, ai_min, bin_dict)
#plot_scatter_hist_L2_L3_errors(direct_forcings, calc_forcings, \
#    sim_name, num_bins = 200, delta_calc = 20, astrofit = True, \
#    screen_outliers = False, save = False)
#plot_hist_L2_L3_errors(direct_forcings, calc_forcings, sim_name, \
#    ax = None, num_bins = 200, astrofit = False, save = False)
#sys.exit()


## From np.mean and np.std on the errors
## -------------------------------------
#errors = calc_forcings - direct_forcings
#mean_error = np.mean(errors) #  3.72
#std_error  = np.std(errors)  # 29.04

## From astropy fit to error histogram.
## -----------------------------------
## noland72?
#mean_error = -3.321
#std_error  = 22.969
## noland50, 6 surface types
#mean_error = -2.571
#std_error  = 22.167

begin_date = '200504'
end_date   = '202009'
season     = 'sunlight'
minlat = 65
maxlat = 87
shawn_file = home_dir + '/Research/OMI/omi_shawn_daily_2005_2020_v2.hdf5'
#jz_file    = home_dir + '/Research/OMI/omi_VJZ211_daily_2005_2020.hdf5'


#NOTE: This is the format used to make the _v3 files
# Is this the "clean-sky" background version?
OMI_daily_VSJ4  = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -20.00, max_AI = 0.7, minlat = minlat, maxlat = maxlat)
# Is this the OMI AI monthly trend version?
#OMI_daily_VSJ4  = calcOMI_MonthAvg_FromDaily(shawn_file, min_AI = -0.10, max_AI = 20.0, minlat = minlat, maxlat = maxlat)

#OMI_daily_VJZ211 = calcOMI_MonthAvg_FromDaily(jz_file, min_AI = -0.10, minlat = minlat, maxlat = maxlat)
daily_VSJ4 = readOMI_daily_HDF5(shawn_file, minlat = minlat, maxlat = maxlat)

ai_thresh = 0.7
maxerr = 1.5


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Test the L2L3_err stuff for a whole study period
#
# NOTE: RUNNING WITH filter_bad_vals SET TO False
# 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Control daily values
# --------------------
#daily_filename = 'arctic_daily_est_forcing_v1.hdf5'        
if(sim_name == 'noland74'):
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

    # Daily values with ice modifiations
    # ----------------------------------
    #ice_filename = 'arctic_daily_est_forcing_iceerr_v1.hdf5'
    ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr.hdf5' # noland72
    ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v1.hdf5' # noland50
    ice_filename = 'arctic_daily_est_forcing_numsfcbins4_iceerr.hdf5' # noland50
    ice_filename = 'arctic_daily_est_forcing_numsfcbins4_iceerr_v1.hdf5' # noland72
    ice_filename = 'arctic_daily_est_forcing_numsfcbins6_iceerr_v2.hdf5' # noland74

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
else:
    print("INVALID SIM NAME")
    sys.exit()



#daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)
#ice_dict   = read_daily_month_force_L2L3_error_from_HDF5(ice_filename)
#cod_dict   = read_daily_month_force_L2L3_error_from_HDF5(cod_filename)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# This is how the daily, ice, and cod files were generated
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

## CONTROL RUN
## -----------
#cod_err_mean = None
#cod_err_std = None
#ice_err_mean = None
#ice_err_std = None
#all_month_vals_orig_alldaily = \
#    calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    filter_bad_vals = False, return_modis_nsidc = False, \
#    mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
#    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#    use_intercept = True, debug = False)
#
#write_daily_month_force_L2L3_error_to_HDF5(\
#    all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
#    minlat = minlat, maxlat = maxlat, \
#    maxerr = maxerr, ai_thresh = ai_thresh, \
#    L2L3_err_mean = None, L2L3_err_std = None, \
#    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#    dtype = None, \
#    overwrite_old_file = False, \
#    write_daily_values = True, \
#    OMI_daily_data = daily_VSJ4, \
#    save_path = './', name_add = '_numsfcbins6_refcomp')
#sys.exit()

## ICE RUN
## -----------
#cod_err_mean = None
#cod_err_std = None
#ice_err_mean = 0.
#ice_err_std = 15.
#all_month_vals_orig_alldaily = \
#    calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    filter_bad_vals = False, return_modis_nsidc = False, \
#    mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
#    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#    use_intercept = True, debug = False)
#
#write_daily_month_force_L2L3_error_to_HDF5(\
#    all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
#    minlat = minlat, maxlat = maxlat, \
#    maxerr = maxerr, ai_thresh = ai_thresh, \
#    L2L3_err_mean = None, L2L3_err_std = None, \
#    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#    dtype = None, \
#    overwrite_old_file = False, \
#    write_daily_values = True, \
#    OMI_daily_data = daily_VSJ4, \
#    save_path = './', name_add = '_numsfcbins6_iceerr')
#
## COD RUN
## -------
#cod_err_mean = 0
#cod_err_std = 5
#ice_err_mean = None
#ice_err_std = None
#all_month_vals_orig_alldaily = \
#    calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#    ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#    filter_bad_vals = False, return_modis_nsidc = False, \
#    mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
#    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#    use_intercept = True, debug = False)
#
#write_daily_month_force_L2L3_error_to_HDF5(\
#    all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
#    slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
#    minlat = minlat, maxlat = maxlat, \
#    maxerr = maxerr, ai_thresh = ai_thresh, \
#    L2L3_err_mean = None, L2L3_err_std = None, \
#    ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#    cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#    dtype = None, \
#    overwrite_old_file = False, \
#    write_daily_values = True, \
#    OMI_daily_data = daily_VSJ4, \
#    save_path = './', name_add = '_numsfcbins6_coderr')
#
#sys.exit()

# REF CLD 2005 RUN
# ----------------
#cod_err_mean = None
#cod_err_std = None
#ice_err_mean = None
#ice_err_std = None
#years = np.arange(2006, 2021)
#for year in years:
#    ref_cld = str(year)
#    all_month_vals_orig_alldaily = \
#        calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
#        slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#        ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#        reference_cld = ref_cld, reference_ice = None, \
#        filter_bad_vals = False, return_modis_nsidc = False, \
#        mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
#        ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#        cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#        use_intercept = True, debug = False)
#    
#    write_daily_month_force_L2L3_error_to_HDF5(\
#        all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
#        slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
#        minlat = minlat, maxlat = maxlat, \
#        maxerr = maxerr, ai_thresh = ai_thresh, \
#        reference_cld = ref_cld, reference_ice = None, \
#        L2L3_err_mean = None, L2L3_err_std = None, \
#        ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#        cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#        dtype = None, \
#        overwrite_old_file = False, \
#        write_daily_values = True, \
#        OMI_daily_data = daily_VSJ4, \
#        save_path = './', name_add = '_numsfcbins6_refcld' + ref_cld)
#sys.exit()

## REF ICE 2005 RUN
## ----------------
#cod_err_mean = None
#cod_err_std = None
#ice_err_mean = None
#ice_err_std = None
#years = np.arange(2006, 2021)
#for year in years:
#    #ref_ice = '2005'
#    ref_ice = str(year)
#    all_month_vals_orig_alldaily = \
#        calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
#        slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
#        ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
#        reference_cld = None, reference_ice = ref_ice, \
#        filter_bad_vals = False, return_modis_nsidc = False, \
#        mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
#        ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#        cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#        use_intercept = True, debug = False)
#    
#    write_daily_month_force_L2L3_error_to_HDF5(\
#        all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
#        slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
#        minlat = minlat, maxlat = maxlat, \
#        maxerr = maxerr, ai_thresh = ai_thresh, \
#        reference_cld = None, reference_ice = ref_ice, \
#        L2L3_err_mean = None, L2L3_err_std = None, \
#        ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
#        cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
#        dtype = None, \
#        overwrite_old_file = False, \
#        write_daily_values = True, \
#        OMI_daily_data = daily_VSJ4, \
#        save_path = './', name_add = '_numsfcbins6_refice' + ref_ice)
#sys.exit()

# REF BOTH RUNS
cod_err_mean = None
cod_err_std = None
ice_err_mean = None
ice_err_std = None
years = np.arange(2005, 2021, 2)
for year in years:
    #ref_ice = '2005'
    ref_both = str(year)
    all_month_vals_orig_alldaily = \
        calculate_type_forcing_v4_alldaily(daily_VSJ4, OMI_daily_VSJ4, \
        slope_dict_lin, bin_dict, 'all', minlat = minlat, maxlat = maxlat, \
        ai_thresh = ai_thresh, maxerr = maxerr, mod_slopes = None, \
        reference_cld = ref_both, reference_ice = ref_both, \
        filter_bad_vals = False, return_modis_nsidc = False, \
        mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
        ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
        cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
        use_intercept = True, debug = False)
    
    write_daily_month_force_L2L3_error_to_HDF5(\
        all_month_vals_orig_alldaily, OMI_daily_VSJ4, \
        slope_dict_lin, bin_dict, all_month_vals_orig_alldaily, \
        minlat = minlat, maxlat = maxlat, \
        maxerr = maxerr, ai_thresh = ai_thresh, \
        reference_cld = ref_both, reference_ice = ref_both, \
        L2L3_err_mean = None, L2L3_err_std = None, \
        ice_err_mean = ice_err_mean, ice_err_std = ice_err_std, \
        cod_err_mean = cod_err_mean, cod_err_std = cod_err_std, \
        dtype = None, \
        overwrite_old_file = False, \
        write_daily_values = True, \
        OMI_daily_data = daily_VSJ4, \
        save_path = './', name_add = '_numsfcbins6_refboth' + ref_both)
sys.exit()

#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Generate figure 10 for version 7 of the Arctic obs forcing paper. 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Functions for 4-panel plot with all analyzed errors
# - NN - obs
# - L2 - L3
# - Ice errors
# - COD errors
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
num_bins = 100
plot_error_components_combined(direct_forcings, calc_forcings, sim_name, \
    daily_filename, ice_filename, cod_filename, num_bins, \
    astrofit = True, log_scale = True, save = True)
sys.exit()


