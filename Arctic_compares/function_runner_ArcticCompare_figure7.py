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

sza_min = 40
sza_max = 90
ai_calc_bins = np.array([0,2,4,6,30])
slope_vals, intpt_vals, mean_force_vals, std_force_vals, ai_calc_bins = \
    plot_NN_scatter_combined_alltypes(test_dict, bin_dict, \
    ai_min_forslopes, sza_min, sza_max, trend_type = 'linregress', \
    show_specific_cod = None, min_ai_for_stats = 2.0, \
    plot_bounds = False, return_line_vals = True, ai_calc_bins = ai_calc_bins, \
    save = False)

# Print a table containing these output values
# --------------------------------------------
print('# # # # # # # # # # # #\n\n      Slopes\n\n# # # # # # # # # # # #')
title_str = '           '
for ii in range(slope_vals.shape[1]):
    title_str += '{0:>15s}'.format(str(bin_dict['cod_bin_edges'][ii]) + '-' + str(bin_dict['cod_bin_edges'][ii + 1]))
    
print(title_str)

for ii in range(slope_vals.shape[0]):
    out_str = '{0:11s}'.format(str(bin_dict['ice_bin_edges'][ii]) + '-' + str(bin_dict['ice_bin_edges'][ii + 1]))
    #out_str += ["{0:5.1f}".format(ff) for ff in slope_vals[ii,:]]
    for jj in range(slope_vals.shape[1]):
        out_str += "{0:15.1f}".format(slope_vals[ii,jj])
    print(out_str)

print('\n# # # # # # # # # # # #\n\n      Intercepts\n\n# # # # # # # # # # # #')
title_str = '           '
for ii in range(slope_vals.shape[1]):
    title_str += '{0:>15s}'.format(str(bin_dict['cod_bin_edges'][ii]) + '-' + str(bin_dict['cod_bin_edges'][ii + 1]))
    
print(title_str)

for ii in range(slope_vals.shape[0]):
    #out_str = ''
    out_str = '{0:11s}'.format(str(bin_dict['ice_bin_edges'][ii]) + '-' + str(bin_dict['ice_bin_edges'][ii + 1]))
    #out_str += ["{0:5.1f}".format(ff) for ff in slope_vals[ii,:]]
    for jj in range(slope_vals.shape[1]):
        out_str += "{0:15.1f}".format(intpt_vals[ii,jj])
    print(out_str)

# Print the averaged ADRF values in different AI ranges
#### mean_force_vals = np.full( (len(bin_dict['ice_bin_means']), \
####     len(bin_dict['cod_bin_means']), len(ai_calc_bins) - 1), np.nan)
#### std_force_vals = np.full( (len(bin_dict['ice_bin_means']), \
####     len(bin_dict['cod_bin_means']), len(ai_calc_bins) - 1), np.nan)

for ii in range(bin_dict['ice_bin_means'].shape[0]):
    print("\n# # # # # # # # # # # # # # # #\n    Ice bin = " + \
        str(bin_dict['ice_bin_edges'][ii]) + ' - ' + \
        str(bin_dict['ice_bin_edges'][ii + 1]) + \
        '\n# # # # # # # # # # # # # # # #')

    title_str = '     '
    for kk in range(ai_calc_bins.shape[0] - 1):
        title_str += '{0:>16s}'.format(str(ai_calc_bins[kk]) + '-' + \
            str(ai_calc_bins[kk + 1]))
    print(title_str)

    for jj in range(bin_dict['cod_bin_means'].shape[0]):
        #out_str = ''
        out_str = '{0:>11s}'.format(str(bin_dict['cod_bin_edges'][jj]) + \
            '-' + str(bin_dict['cod_bin_edges'][jj + 1]))
        for kk in range(ai_calc_bins.shape[0] - 1):
            out_str += "{0:6.1f} +/- {1:<5.1f}".format(mean_force_vals[ii,jj,kk], std_force_vals[ii,jj,kk])
        print(out_str)
        





