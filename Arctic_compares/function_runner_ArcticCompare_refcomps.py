#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *
import random
#from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

# = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = =

sim_name = 'noland105'
run_type = 'noback2'


if(sim_name == 'noland105'):
    if(run_type == 'noback2'):

        file_start = 'arctic_daily_est_forcing_numsfcbins6_noland105_noback2'
        
        # Remember: for the daily forcing values, noback1/2 are the same thing.
        # Noback2 only matters when calculating the sim values with the errors
        # because Noback2 uses the correct CERES - NN errors to perturb the
        # daily forcing values.
        daily_filename = 'arctic_daily_est_forcing_numsfcbins6_noland105_noback1.hdf5' 

else:
    file_start = 'arctic_daily_est_forcing_numsfcbins6'

    # Control daily values
    # --------------------
    #daily_filename = 'arctic_daily_est_forcing_v1.hdf5'        
    daily_filename = 'arctic_daily_est_forcing_numsfcbins6.hdf5' # noland72
    daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v1.hdf5' # noland50
    daily_filename = 'arctic_daily_est_forcing_numsfcbins4.hdf5' # noland50
    daily_filename = 'arctic_daily_est_forcing_numsfcbins4_v1.hdf5' # noland72
    daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v2.hdf5' # noland74
    #daily_filename = 'arctic_daily_est_forcing_numsfcbins6_v3.hdf5' # noland103
    #daily_filename = 'arctic_daily_est_forcing_numsfcbins6_refcomp.hdf5' # noland74
    daily_filename = 'arctic_daily_est_forcing_numsfcbins6_refcomp_v1.hdf5' # noland74, redone with max cod 50

    # Daily values with ref_cld
    # -------------------------
    refcld_filename = 'arctic_daily_est_forcing_numsfcbins6_refcld2005.hdf5' # noland74 , new error (doesn't matter)
    
    # Daily values with ref_ice
    # -------------------------
    refice_filename = 'arctic_daily_est_forcing_numsfcbins6_refice2005.hdf5' # noland74 , new error (doesn't matter)
    
    # Daily values with ref_ice
    # -------------------------
    refboth_filename = 'arctic_daily_est_forcing_numsfcbins6_refboth2005.hdf5' # noland74 , new error (doesn't matter)


daily_trends, cld_trends = plot_reficecld_comps_many_allregions(daily_filename, 'cld', \
    file_start = file_start, trend_type = 'linregress', \
    minlat = 65.5, maxlat = 90.5, save = False, return_trends = True)


daily_trends, ice_trends = plot_reficecld_comps_many_allregions(daily_filename, 'ice', \
    file_start = file_start, trend_type = 'linregress', \
    minlat = 65.5, maxlat = 90.5, save = False, return_trends = True)

#daily_trends, both_trends = plot_reficecld_comps_many_allregions(daily_filename, 'both', \
#    file_start = file_start, trend_type = 'linregress', \
#    minlat = 65.5, maxlat = 90.5, save = False, return_trends = True)


"""
daily_trends, cld_trends = plot_reficecld_comps_many_allregions(daily_filename, 'cld', \
    file_start = 'arctic_daily_est_forcing_numsfcbins6', trend_type = 'linregress', \
    minlat = 65.5, maxlat = 90.5, save = False, return_trends = True)

daily_trends, ice_trends = plot_reficecld_comps_many_allregions(daily_filename, 'ice', \
    file_start = 'arctic_daily_est_forcing_numsfcbins6', trend_type = 'linregress', \
    minlat = 65.5, maxlat = 90.5, save = False, return_trends = True)

daily_trends, both_trends = plot_reficecld_comps_many_allregions(daily_filename, 'both', \
    file_start = 'arctic_daily_est_forcing_numsfcbins6', trend_type = 'linregress', \
    minlat = 65.5, maxlat = 90.5, save = False, return_trends = True)
"""

raw_cld_errors = np.full( cld_trends.shape, np.nan)
raw_ice_errors = np.full( cld_trends.shape, np.nan)
pcnt_cld_errors = np.full( cld_trends.shape, np.nan)
pcnt_ice_errors = np.full( cld_trends.shape, np.nan)

# Dimensions: 6, 3, 16
# 6:  month idx
# 3:  region idx
# 16: refice/cld idx

# Month loop
for ii in range(daily_trends.shape[0]):
    # Region loop
    for jj in range(daily_trends.shape[1]):
        # Calculate the error statistics between the original value and ref val
        # ---------------------------------------------------------------------
        raw_cld_errors[ii,jj]  = cld_trends[ii,jj,:]  - daily_trends[ii,jj]
        pcnt_cld_errors[ii,jj] = abs(raw_cld_errors[ii,jj] / daily_trends[ii,jj]) * 100.

        raw_ice_errors[ii,jj]  = ice_trends[ii,jj,:]  - daily_trends[ii,jj]
        pcnt_ice_errors[ii,jj] = abs(raw_ice_errors[ii,jj] / daily_trends[ii,jj]) * 100.

        print(ii, jj, np.round(daily_trends[ii,jj], 4), \
            np.round(np.mean(raw_ice_errors[ii,jj]), 4), \
            np.round(np.mean(raw_cld_errors[ii,jj]), 4), \
            np.round(np.min(raw_ice_errors[ii,jj]), 4), \
            np.round(np.min(raw_cld_errors[ii,jj]), 4), \
            np.round(np.max(raw_ice_errors[ii,jj]), 4), \
            np.round(np.max(raw_cld_errors[ii,jj]), 4), \
            np.round(np.mean(pcnt_ice_errors[ii,jj]), 4), \
            np.round(np.mean(pcnt_cld_errors[ii,jj]), 4), \
            )
         
 

sys.exit()
# Combine refcld and refice values here
# -------------------------------------
daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)
refcld_dict = read_daily_month_force_L2L3_error_from_HDF5(refcld_filename)
refice_dict = read_daily_month_force_L2L3_error_from_HDF5(refice_filename)

# Find all the ref*** simulation files
# ------------------------------------
refcld_files = glob('arctic_daily_est_forcing_numsfcbins6_refcld*.hdf5')
refice_files = glob('arctic_daily_est_forcing_numsfcbins6_refice*.hdf5')

if( len(refcld_files) != len(refice_files) ):  
    print("ERROR: Not the same number of ref files between cloud and ice")
    sys.exit()

num_ref_sims = len(refcld_files)

# Calculate the values for the control run
# ----------------------------------------
daily_month_vals = calc_monthly_force_from_daily(daily_dict, minlat = 65.5, maxlat = 90.5, \
    return_std = False)

daily_vals = np.full( (3, daily_month_vals.shape[0]), np.nan)
refcld_vals = np.full( (num_ref_sims, 3, daily_month_vals.shape[0]), np.nan)
refice_vals = np.full( (num_ref_sims, 3, daily_month_vals.shape[0]), np.nan)

lat_mins = [65.5, 65.5, 75.5]
lat_maxs = [89.5, 75.5, 89.5]


for jj in range(len(lat_mins)):
    lat_idxs = np.where( (daily_dict['latitude'][:] >= lat_mins[jj]) & (daily_dict['latitude'][:] < lat_maxs[jj])) 
    lat_beg_idx = lat_idxs[0][0]
    lat_end_idx = lat_idxs[0][-1] + 1
    daily_vals[jj,:] = \
        np.nanmean(daily_month_vals[:,lat_beg_idx:lat_end_idx,:], axis = (1,2))

# Calculate the values for the refcld/ice sims
# --------------------------------------------
for ii in range(num_ref_sims):

    refcld_dict = read_daily_month_force_L2L3_error_from_HDF5(refcld_files[ii])
    refice_dict = read_daily_month_force_L2L3_error_from_HDF5(refice_files[ii])

    refcld_month_vals = calc_monthly_force_from_daily(refcld_dict, minlat = 65.5, maxlat = 90.5, \
        return_std = False)
    
    refice_month_vals = calc_monthly_force_from_daily(refice_dict, minlat = 65.5, maxlat = 90.5, \
        return_std = False)

    for jj in range(len(lat_mins)):
        lat_idxs = np.where( (daily_dict['latitude'][:] >= lat_mins[jj]) & (daily_dict['latitude'][:] < lat_maxs[jj])) 
        lat_beg_idx = lat_idxs[0][0]
        lat_end_idx = lat_idxs[0][-1] + 1
        #daily_regions[jj,:] = \
        #    np.nanmean(daily_month_vals[:,lat_beg_idx:lat_end_idx,:], axis = (1,2))
        refcld_vals[ii,jj,:] = \
            np.nanmean(refcld_month_vals[:,lat_beg_idx:lat_end_idx,:], axis = (1,2))
        refice_vals[ii,jj,:] = \
            np.nanmean(refice_month_vals[:,lat_beg_idx:lat_end_idx,:], axis = (1,2))

# Plot the results
# ----------------
fig = plt.figure(figsize = (8, 5))
axs = fig.subplots(2,3, sharex = True, sharey = True)
flat_axs = axs.flatten()

x_vals = np.arange(2005, 2021)   

def plot_ref_time_series(x_data, y_data, pax, trend_type, label, color = None):
    
    pax.plot(x_data, y_data, label = label, color = color)
    if((trend_type == 'standard') | (trend_type == 'linregress')): 
        result = stats.linregress(x_vals, y_data)
        forcing_trend = result.slope * len(x_vals)
        forcing_pval  = result.pvalue
        forcing_uncert = result.stderr * len(x_vals)
    else:
        res = stats.theilslopes(y_data, x_vals, 0.90)
        forcing_trend = res[0]*len(x_vals)


    return forcing_trend    

cld_trends = np.full( (6, num_ref_sims), np.nan)
ice_trends = np.full( (6, num_ref_sims), np.nan)

reg_idx = 1 
trend_type = 'linregress'
for ii in range(len(flat_axs)):

    for jj in range(refcld_vals.shape[0]):
        # Now, do ref cld trends
        # ----------------------
        cld_trends[ii,jj] = plot_ref_time_series(x_vals, \
            refcld_vals[jj,reg_idx,ii::6], flat_axs[ii], trend_type, 'Cld2005')

        # Now, do ref ice trends
        # ----------------------
        ice_trends[ii,jj] = plot_ref_time_series(x_vals, \
            refice_vals[jj,reg_idx,ii::6], flat_axs[ii], trend_type, 'Ice2005')

    # Do the daily trends first
    # -------------------------
    daily_trend = plot_ref_time_series(x_vals, \
        daily_vals[reg_idx,ii::6], flat_axs[ii], trend_type, 'Control', color = 'k')

    flat_axs[ii].axhline(0, color = 'k', linestyle = ':')

    print("DAILY:", np.round(daily_trend, 4), \
          "CLD:",   np.round(np.mean(cld_trends[ii,:]), 4), \
          "ICE:",   np.round(np.mean(cld_trends[ii,:]), 4))

flat_axs[5].legend()

fig.tight_layout()       
plt.show() 


sys.exit()

#    daily_month_vals = calc_monthly_force_from_daily(daily_dict, minlat = 65.5, maxlat = 90.5, \
#        return_std = False)
#    
#    refcld_month_vals = calc_monthly_force_from_daily(refcld_dict, minlat = 65.5, maxlat = 90.5, \
#        return_std = False)
#    
#    refice_month_vals = calc_monthly_force_from_daily(refice_dict, minlat = 65.5, maxlat = 90.5, \
#        return_std = False)
#    
#    daily_regions  = np.full( (3, daily_month_vals.shape[0]), np.nan)
#    refcld_regions = np.full( (3, daily_month_vals.shape[0]), np.nan)
#    refice_regions = np.full( (3, daily_month_vals.shape[0]), np.nan)
    

plot_reficecld_comps(daily_dict, refcld_dict, refice_dict, 1, \
    save = False)

sys.exit()

sys.exit()
# Plot time series of the many region-averaged monthly forcing values
# for a given region index (0 = entire Arctic, 1 = low Arctic, 
# 2 = high Arctic
# -------------------------------------------------------------------
plot_arctic_avg_region_trends(sim_values, 0)
sys.exit()    
#trends, pvals = calc_arctic_avg_region_trends(sim_values)

