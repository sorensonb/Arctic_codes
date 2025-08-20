#!/usr/bin/env python

"""
  NAME:
    auto_fuliou_runner_combined.py

  PURPOSE:
    Automate the running of FuLiou for each unique point in a passed
        CCPEXCV-HALO NAAPS/lidar file. This wrapper allows the user to
        use a combination of HALO or NAAPS water vapor or aerosol profiles 
        as input for FuLiou.

  NOTES:
    This script runs FuLiou for each individual HALO time index in the file,
        so it takes a while to run.  

    This script automatically generates an HDF5 output file. Unlike previous
        scripts, the user no longer needs to add a command line argument
        to create the HDF5 file.

    By default, the script uses the HALO water vapor and aerosol fields.
        Thus, when running the script with no optional command line 
        arguments, the output will contain heating rates calculated from
        the HALO water vapor and aerosol fields. 

    Running this script with NAAPS water vapor and aerosol should make 
        output that looks very similar to the auto_fuliou_runner_naaps.py
        output. The difference, though, is that this script performs the
        calculations at each HALO time index, so there will be thousands
        of calculations performed rather than the dozens performed in 
        auto_fuliou_runner_naaps.py.

    All NAAPS aerosol fields are currently applied when running the script
        with NAAPS aerosols (ABF, dust, smoke, and sea salt). For a 
        potentially better comparison with HALO, the code could be modified 
        so that only NAAPS dust and sea salt aerosols are used.

    When running with HALO water vapor, there are several things to note:
        - HALO provides water vapor mixing ratio, but FuLiou takes specific
            humidity and relative humidity as input. Specific humidity
            is calculated from only the HALO mixing ratio. However, to
            calculate relative humidity from mixing ratio, the NAAPS 
            temperature and air pressure fields must be used.
        - FuLiou crashes when there are input layers with no valid water
            vapor values. Thus, when running the script with HALO water
            vapor, if the script finds layers with missing HALO water 
            vapor values, it uses the NAAPS water vapor fields at that
            layer.
        - In (at least) the 20220909 CPEX file, there are sub-zero HALO mixing
            ratio values. These values are set to 0, and then are replaced
            with NAAPS water vapor fields.

  SYNTAX:
    python auto_fuliou_runner_combined.py CPEXCV-HALO_file -wv={halo,naaps} -aer={halo,naaps,none}
        wv : specifies which moisture variables to use in FuLiou")
             "-wv=halo"  - use HALO moisture vars
             "-wv=naaps" - use NAAPS moisture vars
        aer : specifies which aerosol variables to use in FuLiou
             "-aer=halo"  - use HALO aerosol vars
             "-aer=naaps" - use NAAPS aerosol vars
             "-aer=none"  - use no aerosol


  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2024/05/17:
      Written (based off of auto_fuliou_runner_lidar)

"""

import sys
from datetime import datetime, timedelta
import h5py
import numpy as np
import os
from scipy.stats import mode

# Use the given aerosol classification to calculate the total
# species-specific AODs
# -----------------------------------------------------------
def calc_species_aod(ext_val, aer_id, lidar_z, lidar_dz, abf_ids, dust_ids, \
        smoke_ids, salt_ids, ii):

    total_abf_aot = total_dust_aot = total_smoke_aot = total_salt_aot = \
        total_ice_aot = total_other_aot = 0.
    for jj in range(ext_val.shape[1]-1, -1, -1):
        if(ext_val[ii,jj] != -9.):
            if(aer_id[ii,jj] in abf_ids):
                total_abf_aot += ext_val[ii,jj] * lidar_dz[jj]  
            elif(aer_id[ii,jj] in dust_ids):
                total_dust_aot += ext_val[ii,jj] * lidar_dz[jj]
            elif(aer_id[ii,jj] in smoke_ids):
                total_smoke_aot += ext_val[ii,jj] * lidar_dz[jj]
            elif(aer_id[ii,jj] in salt_ids):
                total_salt_aot += ext_val[ii,jj] * lidar_dz[jj]
            elif(aer_id[ii,jj] == 1.):
                total_ice_aot += ext_val[ii,jj] * lidar_dz[jj]
            else:
                # If there is no aerosol classification for this bin, use the
                # following:
                # If the bin is below 500 m AND has extinction less than
                # or equal to 0.05 km-1, classify as sea salt. 
                # In any other case, classify as dust.
                # -----------------------------------------------------------
                if(lidar_z[jj] < 500):
                    if(ext_val[ii,jj] <= 5e-5):
                        total_salt_aot += ext_val[ii,jj] * lidar_dz[jj]
                    else:
                        total_dust_aot += ext_val[ii,jj] * lidar_dz[jj]
                else: 
                    total_dust_aot += ext_val[ii,jj] * lidar_dz[jj]

    spec_aods = np.array([total_abf_aot, total_dust_aot, total_smoke_aot, \
        total_salt_aot])
    other_aods = np.array([total_ice_aot, total_other_aot])
        
    return spec_aods, other_aods

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  INTERPOLATION FUNCTIONS
#
#  Used to interpolate the NAAPS fields to the lidar vertical resolution
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# identify_unique_idxs is my own version of the numpy unique function,
# which returns the indices of the first occurence of each unique value
# in a passed NAAPS field array. 
#
# For some reason, the numpy unique function doesn't work correctly
# in identifying the index of the first occurence of each NAAPS field.
# ---------------------------------------------------------------------
def identify_unique_idxs(field_var, xx = 0):

    h_idxs = np.zeros(25)
    h_idxs[:] = -9
    kk = 0
    check_val = -999.

    for jj in range(field_var.shape[1]):
        if(field_var[xx,jj] != check_val):
            h_idxs[kk] = int(jj)
            check_val = field_var[xx,jj]
            kk += 1
    
    kk += 1
    h_idxs[kk] = int(field_var.shape[1] - 1)

    h_idxs = h_idxs[h_idxs >= 0]
    h_idxs = np.array([int(hidx) for hidx in h_idxs])

    return h_idxs

# calc_regress_slope calculates the linear change in a NAAPS field value
# (i.e. temp, pressure, rh, etc.) as a function of the number of lidar
# layers within the NAAPS bin. 
# ---------------------------------------------------------------------
def calc_regress_slope(field_var, h_idxs, idx1, xx = 0):
    return ((field_var[xx, int(h_idxs[idx1 + 1])] - \
        field_var[xx, int(h_idxs[idx1])]) / \
        (h_idxs[idx1 + 1] - h_idxs[idx1]))

# interp_NAAPS interpolates the passed NAAPS fields through each NAAPS
# bin as a function of the number of lidar layers within the NAAPS bin.
# ---------------------------------------------------------------------
def interp_NAAPS(ii, naaps_press, naaps_temp, naaps_relhum, naaps_spechum, \
        naaps_dust_mass, naaps_abf_mass, naaps_salt_mass, naaps_smoke_mass, \
        naaps_dust_ext, naaps_abf_ext, naaps_salt_ext, naaps_smoke_ext, \
        ):

    # Identify the unique layer values
    # --------------------------------
    h_idxs = identify_unique_idxs(naaps_press, xx = ii)

    # Find the last index of the second lowest NAAPS pressure. Do not
    # interpolate through the highest NAAPS layer, because this results
    # in duplicate pressures
    # -----------------------------------------------------------------
    unique_vert_press = np.unique(naaps_press[ii,:])[::-1]
    keep_vert_idxs = np.array([\
        np.where(naaps_press[ii,:] == tup)[0][0] \
        for tup in unique_vert_press])

    new_min_press = naaps_press[ii,keep_vert_idxs[-2]]    
    highest_idx = np.where(naaps_press[ii,:] == new_min_press)[0][-1] + 1
    new_shape = naaps_press[ii,:highest_idx].shape[0]

    # Initialize slopes
    # -----------------
    press_slope   = 0.
    temp_slope    = 0.
    relhum_slope  = 0.
    spechum_slope = 0.

    dust_mass_slope  = 0.
    abf_mass_slope   = 0.
    salt_mass_slope  = 0.
    smoke_mass_slope = 0.
    
    dust_ext_slope  = 0.
    abf_ext_slope   = 0.
    salt_ext_slope  = 0.
    smoke_ext_slope = 0.
    
    # This is just used to keep track of which h_idxs value to use
    kk = 0
    
    interp_press   = np.zeros(new_shape)
    interp_temp    = np.zeros(new_shape)
    interp_relhum  = np.zeros(new_shape)
    interp_spechum = np.zeros(new_shape)
    
    interp_dust_mass  = np.zeros(new_shape)
    interp_abf_mass   = np.zeros(new_shape)
    interp_salt_mass  = np.zeros(new_shape)
    interp_smoke_mass = np.zeros(new_shape)
   
    interp_dust_ext   = np.zeros(new_shape)
    interp_abf_ext    = np.zeros(new_shape)
    interp_salt_ext   = np.zeros(new_shape)
    interp_smoke_ext  = np.zeros(new_shape)
   
    # This loops over the vertical indices. Can use this inside the
    # wrapper below.
    for jj in range(new_shape):
    
        # While looping, if the current loop index is one of the 
        # height bin indices (h_idxs), call calc_regress_slope to reset 
        # the regression equation
        # -------------------------------------------------------------
        if((jj in h_idxs) & (jj != (naaps_press.shape[1] - 1))):
            press_slope = calc_regress_slope(naaps_press, h_idxs, kk, \
                xx = ii)
            temp_slope = calc_regress_slope(naaps_temp, h_idxs, kk, \
                xx = ii)
            relhum_slope = calc_regress_slope(naaps_relhum, h_idxs, kk, \
                xx = ii)
            spechum_slope = calc_regress_slope(naaps_spechum, h_idxs, kk, \
                xx = ii)
    
            dust_mass_slope = calc_regress_slope(\
                naaps_dust_mass, h_idxs, kk, xx = ii)
            abf_mass_slope = calc_regress_slope(\
                naaps_abf_mass, h_idxs, kk, xx = ii)
            salt_mass_slope = calc_regress_slope(\
                naaps_salt_mass, h_idxs, kk, xx = ii)
            smoke_mass_slope = calc_regress_slope(\
                naaps_smoke_mass, h_idxs, kk, xx = ii)
    
            dust_ext_slope = calc_regress_slope(\
                naaps_dust_ext, h_idxs, kk, xx = ii)
            abf_ext_slope = calc_regress_slope(\
                naaps_abf_ext, h_idxs, kk, xx = ii)
            salt_ext_slope = calc_regress_slope(\
                naaps_salt_ext, h_idxs, kk, xx = ii)
            smoke_ext_slope = calc_regress_slope(\
                naaps_smoke_ext, h_idxs, kk, xx = ii)
    
            kk += 1
    
        # Otherwise, if the index is not one of the bins, just use
        # the pre-existing regress equation to find the new values
        # --------------------------------------------------------
        interp_press[jj]   = press_slope   * (jj - h_idxs[kk-1]) + \
            naaps_press[ii, h_idxs[kk-1]]
        interp_temp[jj]    = temp_slope    * (jj - h_idxs[kk-1]) + \
            naaps_temp[ii, h_idxs[kk-1]]
        interp_relhum[jj]  = relhum_slope  * (jj - h_idxs[kk-1]) + \
            naaps_relhum[ii, h_idxs[kk-1]]
        interp_spechum[jj] = spechum_slope * (jj - h_idxs[kk-1]) + \
            naaps_spechum[ii, h_idxs[kk-1]]
    
        interp_dust_mass[jj]   = dust_mass_slope  * (jj - h_idxs[kk-1]) + \
            naaps_dust_mass[ii, h_idxs[kk-1]]
        interp_abf_mass[jj]   = abf_mass_slope    * (jj - h_idxs[kk-1]) + \
            naaps_abf_mass[ii, h_idxs[kk-1]]
        interp_salt_mass[jj]   = salt_mass_slope  * (jj - h_idxs[kk-1]) + \
            naaps_salt_mass[ii, h_idxs[kk-1]]
        interp_smoke_mass[jj]   = smoke_mass_slope * (jj - h_idxs[kk-1]) + \
            naaps_smoke_mass[ii, h_idxs[kk-1]]
    
        interp_dust_ext[jj]   = dust_ext_slope  * (jj - h_idxs[kk-1]) + \
            naaps_dust_ext[ii, h_idxs[kk-1]]
        interp_abf_ext[jj]   = abf_ext_slope    * (jj - h_idxs[kk-1]) + \
            naaps_abf_ext[ii, h_idxs[kk-1]]
        interp_salt_ext[jj]   = salt_ext_slope  * (jj - h_idxs[kk-1]) + \
            naaps_salt_ext[ii, h_idxs[kk-1]]
        interp_smoke_ext[jj]   = smoke_ext_slope * (jj - h_idxs[kk-1]) + \
            naaps_smoke_ext[ii, h_idxs[kk-1]]
    
    return interp_press, interp_temp, interp_relhum, interp_spechum, \
        interp_dust_mass, interp_abf_mass, interp_salt_mass, interp_smoke_mass,\
        interp_dust_ext, interp_abf_ext, interp_salt_ext, interp_smoke_ext

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  END OF INTERPOLATION FUNCTIONS
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# This switch determines if the total AOD given in the first line of the
# FuLiou input file is either the straight lidar AOT file or is the
# hand-calculated value given by totalling the four hand-calculated
# species-specific AODs
# ----------------------------------------------------------------------
l_USE_lidar_AOD = False

# This switch determines if the lidar vertical bins (and interpolated
# NAAPS layers) are averaged to some pre-determined number of bins,
# currently set to 100. This is to address the fact that the FuLiou
# code doesn't run properly when the number of vertical layers
# is more than 100, so the calculation doesn't work when all 827 
# lidar layers are used. 
# ----------------------------------------------------------------------
l_BIN_LIDAR_VERTICAL = True

# This switch determines if an HDF5 file is generated from the
# FuLiou output files, in addition to the text file. The 
# HDF5 file is titled 'fuliou_out_YYYYMMDD.hdf5'. By default, this
# switch is set to False. The user may set this to "True" by
# adding "hdf" as a command line argument after the CPEX file.
# ----------------------------------------------------------------------
l_GEN_HDF5_OUTPUT = True

# This switch runs the FuLiou RTM code without any aerosol (or just dust?).
# Only water vapor is the key player then.
# -------------------------------------------------------------------------
l_NO_AEROSOL = False

# This switch uses the HALO mixing ratio to calculate specific humidity 
# and RH, and then these HALO water vapor values are printed to the 
# FuLiou input file
# ---------------------------------------------------------------------
l_USE_HALO_WV = True

# This switch determines if the NAAPS aerosol fields (at either the raw
# or thinned lidar resolution) are used instead of the HALO aerosol fields.
# Right now, when NAAPS aerosol fields are used, all species are included.
# -------------------------------------------------------------------------
l_USE_NAAPS_AER = False

# Check command line arguments
# ----------------------------
if(len(sys.argv) < 2):
    print("SYNTAX: python auto_fuliou_runner_combined.py CPEXCV-HALO_file "+\
            "-wv={halo,naaps} -aer={halo,naaps,none}")
    print("        wv : specifies which moisture variables to use in FuLiou")
    print("             \"-wv=halo\"   - use HALO moisture vars (DEFAULT)")
    print("             \"-wv=naaps\"  - use NAAPS moisture vars")
    print("        aer : specifies which aerosol variables to use in FuLiou")
    print("             \"-aer=halo\"  - use HALO aerosol vars (DEFAULT)")
    print("             \"-aer=naaps\" - use NAAPS aerosol vars")
    print("             \"-aer=none\"  - use no aerosol")
    sys.exit(1)

wv_adder  = '_wvHALO'
aer_adder = '_aerHALO'
for arg in sys.argv:
    if(arg[1:3] == 'wv'):
        arg_val = arg.strip().split('=')[1]
        if(arg_val == 'halo'):
            l_USE_HALO_WV = True
        elif(arg_val == 'naaps'):
            l_USE_HALO_WV = False
            wv_adder = '_wvNAAPS'
        else:
            print("WARNING: Invalid water vapor argument")
    elif(arg[1:4] == 'aer'):
        arg_val = arg.strip().split('=')[1]
        if(arg_val == 'halo'):
            l_USE_NAAPS_AER = False
        elif(arg_val == 'naaps'):
            l_USE_NAAPS_AER = True
            aer_adder = '_aerNAAPS'
        elif(arg_val == 'none'):
            l_USE_NAAPS_AER = False
            l_NO_AEROSOL = True
            aer_adder = '_aerNONE'
        else:
            print("WARNING: Invalid aerosol argument")

print("HALO WV argument:    ", l_USE_HALO_WV)
print("NAAPS AER argument:  ", l_USE_NAAPS_AER)
print("NO AEROSOL argument: ", l_NO_AEROSOL)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Set path variables
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# base_dir points to the 
# NOTE: Please modify this path accordingly for your system!!!
base_dir = "/Research/fuliou_lidarflux/naaps_lidar_fuliou_package_v2/"
mhome = base_dir + "test_run"
output_dir = mhome + '/data'

homedir = mhome

########################################################
#Parameter Group 2

########################################################
#locations for executables
naaps_fulioucode = base_dir + "sourcecodes"

tabledir  = naaps_fulioucode + "/table"             ##place to store .txt files needed for the program
srcdir  =  naaps_fulioucode 


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Write namelist file
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

namelist_file = 'namelist'
dtg = '2000010100'
with open(namelist_file, 'w') as fname:
    fname.write("    &nl\n" + \
                "    tape        = '" + dtg + "',\n" + \
                "    tabledir    = '" + tabledir + "',\n" + \
                "    sulfatefile = 'anthro_farop_v3.txt',\n" + \
                "    dustfile    = 'dust_farop_v3.txt',\n" + \
                "    smokefile   = 'smoke_farop_v3.txt',\n" + \
                "    saltfile    = 'salt_farop_v3.txt',\n" + \
                "    sceneidfile =  'global_scene_id.map',\n" + \
                "    &end\n")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Open the data file and extract the variables
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

infile = sys.argv[1]
hdf_data = h5py.File(infile)

naaps_press   = hdf_data['NAAPS_pressure'][:,:] / 100.

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Convert the mass concentrations from ug/m3 to kg/m3 by  multiplying
#  each value by 1e-9.
#
#  NOTE: All lidar bins are 14.989 m thick.  Assume this, or calculate
#        using a function?
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Extract metadata
# ----------------
split_date    = hdf_data['time'].attrs['units'].strip().split()
date_str      = ' '.join(split_date[2:])
base_date     = datetime.strptime(date_str, '%Y%m%d %H:%M:%S')
time          = hdf_data['time'][:]
lidar_z       = hdf_data['z'][:]
naaps_masl    = hdf_data['NAAPS_MASL'][:,:]  # Use NAAPS z and dz for interp
naaps_dz      = hdf_data['NAAPS_dz'][:,:]
lat           = hdf_data['lat'][:]
lon           = hdf_data['lon'][:]

lidar_dz      = np.ones(lidar_z.shape)
lidar_dz[:]   = 14.989

# Extract NAAPS metorological variables
# -------------------------------------
naaps_sfcpres = hdf_data['NAAPS_ps'][:]
naaps_temp    = hdf_data['NAAPS_temperature'][:,:]
naaps_spechum = hdf_data['NAAPS_q'][:,:] 
naaps_relhum  = hdf_data['NAAPS_rh'][:,:]

# Extract NAAPS aerosol variables
# -------------------------------
naaps_tot_aot   = hdf_data['NAAPS_total_aot'][:]

naaps_dust_ext  = hdf_data['NAAPS_dust_ext'][:,:]
naaps_dust_aot  = hdf_data['NAAPS_dust_aot'][:]
naaps_dust_mass = \
    hdf_data['NAAPS_dust_mass_concentration'][:,:] * 1e-9

naaps_abf_ext   = hdf_data['NAAPS_ABF_ext'][:,:]
naaps_abf_aot   = hdf_data['NAAPS_ABF_aot'][:]
naaps_abf_mass  = \
    hdf_data['NAAPS_ABF_mass_concentration'][:,:] * 1e-9

naaps_salt_ext  = hdf_data['NAAPS_SS_ext'][:,:]
naaps_salt_aot  = hdf_data['NAAPS_SS_aot'][:]
naaps_salt_mass = \
    hdf_data['NAAPS_SS_mass_concentration'][:,:] * 1e-9

naaps_smoke_ext  = hdf_data['NAAPS_smoke_ext'][:,:]
naaps_smoke_aot  = hdf_data['NAAPS_smoke_aot'][:]
naaps_smoke_mass = \
    hdf_data['NAAPS_smoke_mass_concentration'][:,:] * 1e-9

# Extract lidar aerosol variables
# 
# NOTE: For the lidar total AOT, am using the '532_AOT_lo' parameter from 
#       the files because it closely agrees with the AOT value that is 
#       calculated from the '532_ext' extinction profile data, using
#       
#           np.sum((data['532_ext'][1290,:] / 1e3) * 15)
#
#       From the CPEX subset file, the '532_AOT_lo' variable
#       is "532 nm Total Optical Thickness (extended ext. coeff integration)"
#
#       The other AOT option in the file is '532_AOT_hi', which is:
#       "532 nm Total Optical Thickness (using molecular signal)". However
#       from testing, it was found that these values do not agree as well
#       with the hand-calculated AOT values listed above.
#
# NOTE: Convert the lidar exts from in units of km-1 to units of m-1
#       Also, convert the HALO mixing ratio from g/kg to kg/kg by dividing by
#           1e3
# ---------------------------------------------------------------------------
ext_532       = hdf_data['532_ext'][:,:] / 1e3
if('h2o_mmr_v' in hdf_data.keys()):
    lidar_mixrat  = hdf_data['h2o_mmr_v'][:,:] / 1e3
elif('h2o_mass_mixing_ratio' in hdf_data.keys()):
    lidar_mixrat  = hdf_data['h2o_mass_mixing_ratio'][:,:] / 1e3
else:
    print("ERROR: No recognized HALO mixing ratio value in the file.")
    print("       Please check the variables in the file and modify ")
    print("       the code here.")
    hdf_data.close()
    sys.exit(2)

aot_532       = hdf_data['532_AOT_lo'][:]
aer_id        = hdf_data['Aerosol_ID'][:,:]

# Convert the HALO mixing ratio to specific humidity
# --------------------------------------------------
lidar_q = lidar_mixrat / (1 + lidar_mixrat)

# Replace masked extinctions or water vapors with 0
# -------------------------------------------------
ext_532 = np.ma.where(np.ma.masked_invalid(ext_532).mask == True, \
    0., ext_532)
lidar_q = np.ma.where(np.ma.masked_invalid(lidar_q).mask == True, \
    0., lidar_q)
lidar_mixrat = np.ma.where(np.ma.masked_invalid(lidar_mixrat).mask == True, \
    0., lidar_mixrat)

# NOTE: There may be negative HALO mixing ratios. Set them to 0?
# --------------------------------------------------------------
lidar_q = np.where(lidar_q < 0, 0, lidar_q)

# Close the file
# --------------
hdf_data.close()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Loop over the data and print the stuff
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Variable descriptions
#  0 - interp_p
#  1 - lidar_z
#  2 - lidar_dz
#  3 - interp_t
#  4 - interp_rh
#  5 - interp_sh
#  6 - sulext            
#  7 - dustext
#  8 - smokeext
#  9 - saltext
# 10 - sulmass           # THESE ARE NOT USED. KEEPING FOR NOW     
# 11 - dustmass          # THESE ARE NOT USED. KEEPING FOR NOW     
# 12 - smokemass         # THESE ARE NOT USED. KEEPING FOR NOW
# 13 - saltmass          # THESE ARE NOT USED. KEEPING FOR NOW

fmt_str = "{0:7.1f} {1:8.1f} {2:8.1f} {3:7.2f} " + \
    "{4:7.3f} {5:10.3e} " + \
    "{6:10.3e} {7:10.3e} {8:10.3e} {9:10.3e} " + \
    "{10:10.3e} {11:10.3e} {12:10.3e} {13:10.3e}\n"

# totaod = [aod1, aod2, aod3, aod4, aod5]
# totaod: [UNUSED, sulfate, dust, smoke, salt]

foutname = 'fuliou_input_naaps_file.txt'
out_data_name = 'fuliou_output_file_lidar.txt'

# If the outputfile exists in its current name, delete the old one
# ----------------------------------------------------------------
if(os.path.exists(out_data_name)):
    print("Removing old output data file")
    cmnd = 'rm ' + out_data_name
    os.system(cmnd)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# These are the aerosol IDs that will be binned into each NAAPS type
#
# lidar class values:
#  1 = Ice
#  2 = Dusty Mix
#  3 = Marine
#  4 = Urban/Pollution
#  5 = Smoke
#  6 = Fresh Smoke
#  7 = Pol. Marine
#  8 = Dust
#
# The lidar values are binned to the NAAPS species
# following:
#
#   abf:   4
#   dust:  2, 8
#   smoke: 5, 6
#   salt:  3, 7
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
abf_ids = [4.]
dust_ids = [2., 8.]
smoke_ids = [5., 6.]
salt_ids  = [3., 7.]


check_press = -999.
if(l_BIN_LIDAR_VERTICAL):
    # NOTE: num_layers gives the number of smoothed layers that will
    #       be in the FuLiou input file. 
    # --------------------------------------------------------------
    num_layers = 90
    split_dz = np.array(np.split(lidar_dz[:len(lidar_dz) - \
        (len(lidar_dz) % num_layers)], num_layers))
    avg_dz   = np.mean(np.split(lidar_dz[:len(lidar_dz) - \
        (len(lidar_dz) % num_layers)], num_layers), axis = 1)
    avg_z    = np.mean(np.array(np.split(lidar_z[:len(lidar_dz) - \
        (len(lidar_dz) % num_layers)], num_layers)), axis = 1)

    # Save a variable to hold the number of HALO vertical bins
    # in each of these split bins
    # -------------------------------------------------------
    num_layers_per_bin = split_dz.shape[1]
else:
    # Not sure yet if this works, but am inserting this line for now
    # --------------------------------------------------------------
    num_layers = ext_532.shape[1]
    num_layers_per_bin = 1

# Prep the output HDF5 file, if desired
# -------------------------------------
if(l_GEN_HDF5_OUTPUT):
    # NOTE: In the FuLiou code, the output printed to the screen begins at 
    # index 2 of the arrays in the code. This is mirrored here
    # --------------------------------------------------------------------
    num_pres_layers = num_layers - 2
    num_heat_layers = num_layers - 3

    output_time           = np.zeros((time.shape[0]))
    output_lat            = np.zeros((time.shape[0]))
    output_lon            = np.zeros((time.shape[0]))
    output_aod            = np.zeros((time.shape[0]))
    output_sfc_dwn_sw_flx = np.zeros((time.shape[0]))
    output_toa_up_sw_flx  = np.zeros((time.shape[0]))
    output_sfc_dwn_lw_flx = np.zeros((time.shape[0]))
    output_toa_up_lw_flx  = np.zeros((time.shape[0]))

    output_alt        = np.zeros((num_pres_layers))
    output_press      = np.zeros((time.shape[0], num_pres_layers))
    output_tmp        = np.zeros((time.shape[0], num_pres_layers))
    output_relhum     = np.zeros((time.shape[0], num_pres_layers))
    output_spchum     = np.zeros((time.shape[0], num_pres_layers))
    output_abf_ext    = np.zeros((time.shape[0], num_pres_layers))
    output_dust_ext   = np.zeros((time.shape[0], num_pres_layers))
    output_smoke_ext  = np.zeros((time.shape[0], num_pres_layers))
    output_salt_ext   = np.zeros((time.shape[0], num_pres_layers))

    output_heat       = np.zeros((time.shape[0], num_heat_layers))


    output_time[:]           = np.nan
    output_lat[:]            = np.nan
    output_lon[:]            = np.nan
    output_aod[:]            = np.nan
    output_sfc_dwn_sw_flx[:] = np.nan
    output_toa_up_sw_flx[:]  = np.nan
    output_sfc_dwn_lw_flx[:] = np.nan
    output_toa_up_lw_flx[:]  = np.nan
    output_press[:,:]        = np.nan
    output_heat[:,:]         = np.nan

    output_alt[:]          = np.nan
    output_tmp[:,:]        = np.nan
    output_relhum[:,:]     = np.nan
    output_spchum[:,:]     = np.nan
    output_abf_ext[:,:]    = np.nan
    output_dust_ext[:,:]   = np.nan
    output_smoke_ext[:,:]  = np.nan
    output_salt_ext[:,:]   = np.nan

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Now, run the FuLiou code for each point
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
for ii in range(0, time.shape[0]):

    # If the old single-run output file (fort.22) exists, remove it
    # -------------------------------------------------------------
    if(os.path.exists('fort.22')):
        print("Removing fort.22 file")
        cmnd = 'rm fort.22'
        os.system(cmnd)
    
    if(naaps_press[ii,0] != check_press):

        print("    interpolating")

        # Interpolate the values here
        # ----------------------------
        interp_press, interp_temp, interp_relhum, interp_spechum, \
            interp_dust_mass, interp_abf_mass, interp_salt_mass, \
            interp_smoke_mass, interp_dust_ext, interp_abf_ext, \
            interp_salt_ext, interp_smoke_ext = \
            interp_NAAPS(ii, naaps_press, naaps_temp, naaps_relhum, \
                naaps_spechum, \
                naaps_dust_mass, naaps_abf_mass, \
                naaps_salt_mass, naaps_smoke_mass, \
                naaps_dust_ext, naaps_abf_ext, \
                naaps_salt_ext, naaps_smoke_ext, \
                )

        if(l_BIN_LIDAR_VERTICAL):
            # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
            #
            # Thin the data, if desired
            #
            # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
            split_press       = np.mean(np.array(np.split(\
                interp_press[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_temp        = np.mean(np.array(np.split(\
                interp_temp[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_relhum      = np.mean(np.array(np.split(\
                interp_relhum[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_spechum     = np.mean(np.array(np.split(\
                interp_spechum[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)

            split_dust_mass   = np.mean(np.array(np.split(\
                interp_dust_mass[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_abf_mass    = np.mean(np.array(np.split(\
                interp_abf_mass[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_salt_mass   = np.mean(np.array(np.split(\
                interp_salt_mass[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_smoke_mass  = np.mean(np.array(np.split(\
                interp_smoke_mass[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)

            split_dust_ext   = np.mean(np.array(np.split(\
                interp_dust_ext[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_abf_ext    = np.mean(np.array(np.split(\
                interp_abf_ext[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_salt_ext   = np.mean(np.array(np.split(\
                interp_salt_ext[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_smoke_ext  = np.mean(np.array(np.split(\
                interp_smoke_ext[:len(interp_press) - \
                (len(interp_press) % num_layers)], num_layers)), axis = 1)

        check_press = naaps_press[ii,0]

    mode_vals = mode(split_press)
    if(mode_vals[1][0] != 1):
        print(ii, "BAD PRESS", mode_vals[0][0], mode_vals[1][0])

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Test the correction to the aerosol ID problem
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

    if(l_BIN_LIDAR_VERTICAL): 
        # Divide the lidar extinction, specific humidity, mixing ratio, and
        # aerosol IDs into even bins, depending on the number of layers
        # specified by the user
        # -----------------------------------------------------------------
        split_ext     = np.array(np.split(ext_532[ii,:len(interp_press) - \
            (len(interp_press) % num_layers)], num_layers))
        split_lidar_q = np.array(np.split(lidar_q[ii,:len(interp_press) - \
            (len(interp_press) % num_layers)], num_layers))
        split_lidar_mixrat = np.array(np.split(lidar_mixrat[ii,:len(interp_press) - \
            (len(interp_press) % num_layers)], num_layers))
        split_aerid   = np.array(np.split(aer_id[ii,:len(interp_press) - \
            (len(interp_press) % num_layers)], num_layers))
   
        # Calculate the bin averages for extinction, specific humidity, 
        # and mixing ratio
        # -------------------------------------------------------------
        avg_exts   = np.average(split_ext, axis = 1)
        avg_q      = np.average(split_lidar_q, axis = 1)
        avg_mixrat = np.average(split_lidar_mixrat, axis = 1)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Need to calculate RH based on the temp, press, and mixing ratio.
    # The HALO mixing ratio (and derived RH) will then be used in FuLiou
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    epsilon = 0.622   # molecular weight ratio of vapor to dry air
    if(l_BIN_LIDAR_VERTICAL):
        # First, calculate saturation vapor pressure from the NAAPS temp
        # Output is in hPa, and air temperature must be in degrees Celsius. 
        # NAAPS air temp is used here.
        # -----------------------------------------------------------------
        sat_vap_press = 6.112 * np.exp( (17.67 * (split_temp - 273.15)) / \
            ((split_temp - 273.15) + 243.5))        

        # Then, use the saturation vapor pressure and NAAPS air pressure to
        # calculate the saturation mixing ratio
        # ----------------------------------------------------------------- 
        sat_mix_rat = epsilon * (sat_vap_press / (split_press - sat_vap_press))

        # Finally, calculate the relative humidity given the HALO mixing
        # ratio and the derived saturation mixing ratio
        # -------------------------------------------------------------- 
        avg_lidar_rh = ((avg_mixrat / (epsilon + avg_mixrat)) * \
            ((epsilon + avg_mixrat) / sat_mix_rat)) * 100.
    else:
        # First, calculate saturation vapor pressure from the NAAPS temp
        # Output is in Pascals, and air temperature must be in degrees Celsius
        sat_vap_press = 6.112 * np.exp( (17.67 * (interp_temp - 273.15)) / \
            ((interp_temp - 273.15) + 243.5))        

        # Then, use the saturation vapor pressure and NAAPS air pressure to
        # calculate the saturation mixing ratio
        # ----------------------------------------------------------------- 
        sat_mix_rat = epsilon * (sat_vap_press / (interp_press - sat_vap_press))
   
        # Finally, calculate the relative humidity given the HALO mixing
        # ratio and the derived saturation mixing ratio
        # -------------------------------------------------------------- 
        lidar_rh = ((lidar_mixrat[ii,:interp_press.shape[0]] / \
            (epsilon + lidar_mixrat[ii,:interp_press.shape[0]])) * \
            ((epsilon + lidar_mixrat[ii,:interp_press.shape[0]]) / \
            sat_mix_rat)) * 100.

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Hand-calculate the species-specific AOTs here. These will replace the 
    #   'naaps_abf_aot'-like variables printed at the top of the output file
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    if(l_USE_NAAPS_AER):
        spec_aods = [0,0,0,0]
        other_aods = [0,0]
        spec_aods[0] = np.sum(split_abf_ext * (avg_dz * num_layers_per_bin))
        spec_aods[1] = np.sum(split_dust_ext * (avg_dz * num_layers_per_bin))
        spec_aods[2] = np.sum(split_smoke_ext * (avg_dz * num_layers_per_bin))
        spec_aods[3] = np.sum(split_salt_ext * (avg_dz * num_layers_per_bin))

    else:
        spec_aods, other_aods = \
            calc_species_aod(ext_532, aer_id, lidar_z, lidar_dz, abf_ids, \
                dust_ids, smoke_ids, salt_ids, ii)

    calc_total_aot = np.sum(spec_aods)


    with open(foutname, 'w') as fout:
    
        """
        
        Output format
        YYYYMMDDHHMM lat lon tau_aer
         sulfate_aod, dust_aod, smoke_aod, salt_aod
         press lidar_z lidar_dz temp rhum spechum abf_ext dust_ext smoke_ext salt_ext abf_mas dust_mas smoke_mas salt_mas ext_532
         ...
         ...
         ...
        
         
        """
   
        # Print the metadata information
        # ------------------------------
        local_dt_date = base_date + timedelta(seconds = time[ii]) 
        local_date = local_dt_date.strftime('%Y%m%d%H%M%S')
        local_fmt = "{0:7.4f} {1:7.4f} {2:7.4f} {3:7.4f} {4:7.4f}"
        if(l_NO_AEROSOL):
            print(ii, local_dt_date.strftime('%Y%m%d %H:%M:%S'),\
                local_fmt.format(0.0, 0.0, 0.0, 0.0, 0.0 ))
        else:
            print(ii, local_dt_date.strftime('%Y%m%d %H:%M:%S'),\
                local_fmt.format(aot_532[ii], \
                    calc_total_aot, other_aods[1] + calc_total_aot, \
                    other_aods[0], other_aods[1] ))

        if(l_NO_AEROSOL):
            fout.write("{0} {1:8.5f} {2:8.5f} {3:8.5f}\n".format(\
                local_date, lat[ii], lon[ii], 0.0))
        else:
            if(l_USE_lidar_AOD):
                fout.write("{0} {1:8.5f} {2:8.5f} {3:8.5f}\n".format(\
                    local_date, lat[ii], lon[ii], aot_532[ii]))
            else:
                fout.write("{0} {1:8.5f} {2:8.5f} {3:8.5f}\n".format(\
                    local_date, lat[ii], lon[ii], calc_total_aot))
        
        if(l_NO_AEROSOL):
            fout.write("{0:10.3e} {1:10.3e} {2:10.3e} {3:10.3e}\n".format(\
                0.0, 0.0, 0.0, 0.0))
        else: 
            fout.write("{0:10.3e} {1:10.3e} {2:10.3e} {3:10.3e}\n".format(\
                spec_aods[0], spec_aods[1], \
                spec_aods[2], spec_aods[3]))
        
        # Print the profile information in reverse order, with the lowest
        # data in the profile being at the last index.
        # ----------------------------------------------------------------------
        if(l_BIN_LIDAR_VERTICAL):
            vert_dim = avg_exts.shape[0] - 1
        else:
            vert_dim = ext_532.shape[1] - 1
        for jj in range(vert_dim, -1, -1):

            local_abf = local_dust = local_smoke = local_salt = 0.

            if(l_USE_NAAPS_AER):
                # If the "l_NO_AEROSOL" switch is True, then the 
                # local_**** aerosol values will remain at 0
                # ----------------------------------------------
                if(not l_NO_AEROSOL):
                    # ABF
                    local_abf = split_abf_ext[jj]
                    # Dust
                    local_dust = split_dust_ext[jj]
                    # Smoke
                    local_smoke = split_smoke_ext[jj]
                    # Salt 
                    local_salt = split_salt_ext[jj]
    
            else:

                # Use the lidar aerosol classification to put the lidar
                # extinction values into the correct type slot:
                # abf, dust, smoke, or salt
                #
                # lidar class values:
                #  1 = Ice
                #  2 = Dusty Mix
                #  3 = Marine
                #  4 = Urban/Pollution
                #  5 = Smoke
                #  6 = Fresh Smoke
                #  7 = Pol. Marine
                #  8 = Dust
                #
                # The lidar values are binned to the NAAPS species
                # following:
                #
                #   abf:   4
                #   dust:  2, 8
                #   smoke: 5, 6
                #   salt:  3, 7
                # -------------------------------------------------------

                if(l_BIN_LIDAR_VERTICAL):
                    # In this section, account for lidar bins that do not contain
                    # classification. Classify as "sea salt" if the bin is less
                    # than 500 m MSL and the extinction is less than 0.05 km-1.
                    # For any other case (greater than 500 m MSL or extinction
                    # greater than 0.05 km-1), classify as dust.
                    local_ids = np.ma.masked_invalid(split_aerid[jj])
                    local_ids = np.where(local_ids.mask, -9, local_ids)
                    if(avg_z[jj] > 500):
                        local_ids = np.where(local_ids == -9, 2, local_ids)
                    else:
                        local_ids = np.where((local_ids == -9) & (split_ext[jj] <= 5e-5), 3, local_ids)
                        local_ids = np.where((local_ids == -9) & (split_ext[jj] > 5e-5), 2, local_ids)

                    unique_ids = np.unique(local_ids)
                    work_id = unique_ids[np.argmax([\
                        len(np.where(local_ids == uq)[0]) for uq in unique_ids])]
                    work_ext = avg_exts[jj]
                else:
                    work_id = aer_id[ii,jj]
                    work_ext = ext_532[ii,jj]

                # If the "l_NO_AEROSOL" switch is True, then the 
                # local_**** aerosol values will remain at 0
                # ----------------------------------------------
                if(not l_NO_AEROSOL):
                    # ABF
                    if(work_id in abf_ids):
                        local_abf = work_ext
                    # Dust
                    elif(work_id in dust_ids):
                        local_dust = work_ext
                    # Smoke
                    elif(work_id in smoke_ids):
                        local_smoke = work_ext
                    # Salt 
                    elif(work_id in salt_ids):
                        local_salt = work_ext

            # 2024/05/14: Modifying the print statements to use the HALO
            #       moisture variables rather than NAAPS
            if(l_BIN_LIDAR_VERTICAL):
                if( (l_USE_HALO_WV) and (avg_q[jj] != 0.)):
                    local_relhum  = avg_lidar_rh[jj]
                    local_spechum = avg_q[jj]
                else:
                    local_relhum  = split_relhum[jj]
                    local_spechum = split_spechum[jj]

                if(local_spechum != 0):
                    fout.write(fmt_str.format(\
                        split_press[jj], float(avg_z[jj]), avg_dz[jj], \
                        split_temp[jj], local_relhum, \
                        local_spechum, \
                        local_abf, local_dust, local_smoke, local_salt, \
                        split_abf_mass[jj], split_dust_mass[jj], \
                        split_smoke_mass[jj], split_salt_mass[jj], \
                    ))
            else:
                if( (l_USE_HALO_WV) and (lidar_q[jj] != 0.)):
                    local_relhum  = lidar_rh[jj]
                    local_spechum = lidar_q[jj]
                else:
                    local_relhum  = interp_relhum[jj]
                    local_spechum = interp_spechum[jj]

                if(local_spechum != 0):
                    fout.write(fmt_str.format(\
                        interp_press[jj], float(lidar_z[jj]), lidar_dz[jj], \
                        interp_temp[jj], local_relhum, \
                        local_spechum,\
                        local_abf, local_dust, local_smoke, local_salt, \
                        interp_abf_mass[jj], interp_dust_mass[jj], \
                        interp_smoke_mass[jj], interp_salt_mass[jj], \
                    ))
   
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #  Call the FuLiou code
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    print(' Execute naaps_fuliou.exe')
    cmnd = srcdir + '/bin/naaps_fuliou.exe'
    print(cmnd)
    os.system(cmnd)
   
    # Concatenate the temp output file to the end of the actual
    # output file. Only do this if the code generated an output file, 
    # which only happens if the code successfully runs.
    # ---------------------------------------------------------
    if(os.path.exists('fort.22')):
        print("COPYING FILE OUTPUTS")
        cmnd = 'cat fort.22 >> ' + out_data_name
        print(cmnd)
        os.system(cmnd)

        if(l_GEN_HDF5_OUTPUT):
        
            output_time[ii]            = time[ii]

            # Read in the last output file
            # ----------------------------
            with open('fort.22','r') as fin:
                flines = fin.readlines()

            line1 = flines[0].strip().split()
            if(len(line1) == 7):
                output_lat[ii]            = float(line1[0])
                output_lon[ii]            = float(line1[1])
                output_aod[ii]            = float(line1[2])
                output_sfc_dwn_sw_flx[ii] = float(line1[3])
                output_toa_up_sw_flx[ii]  = float(line1[4])
                output_sfc_dwn_lw_flx[ii] = float(line1[5])
                output_toa_up_lw_flx[ii]  = float(line1[6])

                # Figure out how to properly divide the file
                # ------------------------------------------
                num_data_lines = (len(flines) - 1) / 2
               
                work_press = np.array([float(tvar) for tvar in \
                    ''.join(flines[1:num_data_lines+1]).strip().split()])
                work_heat = np.array([float(tvar) for tvar in \
                    ''.join(flines[num_data_lines+1:]).strip().split()])

                # Read in the input file data, containing met profiles and aerosol
                # ----------------------------------------------------------------
                indata = np.loadtxt('fuliou_input_naaps_file.txt', skiprows = 2)
                work_idx = num_layers - num_pres_layers

                # If the number of pressures here matches the expected
                # size, then just insert the work arrays into the output
                # arrays
                if(len(work_press) == num_pres_layers):
                    print("MATCH")
                    output_press[ii,:] = work_press
                    output_heat[ii,:] = work_heat

                    output_alt[:]           = indata[work_idx:,1]
                    output_tmp[ii,:]        = indata[work_idx:,3]
                    output_relhum[ii,:]     = indata[work_idx:,4]
                    output_spchum[ii,:]     = indata[work_idx:,5]
                    output_abf_ext[ii,:]    = indata[work_idx:,6]
                    output_dust_ext[ii,:]   = indata[work_idx:,7]
                    output_smoke_ext[ii,:]  = indata[work_idx:,8]
                    output_salt_ext[ii,:]   = indata[work_idx:,9]
                else:
                    print("MISMATCH")

            else:
                print("WARNING: BAD fort.22 FILE AT INDEX",ii)
    
# At the end, if desired, make the output HDF5 file
if(l_GEN_HDF5_OUTPUT):

    out_name = 'fuliou_output_lidar' + wv_adder + aer_adder + '_' + \
        infile.strip().split('/')[-1].split('_')[2] + '.h5'

    # Here, remove the indices that contain missing values?
    # -----------------------------------------------------
    keep_idxs = np.where(~np.isnan(output_lat[:]))

    dset = h5py.File(out_name, 'w')

    cdt = dset.create_dataset('time',   data = output_time[keep_idxs])
    cdt.attrs['long_name'] = 'UTC time of measurement instant'
    cdt.attrs['units'] = 'seconds after midnight on date provided in file name'
    cdt = dset.create_dataset('latitude',   data = output_lat[keep_idxs])
    cdt.attrs['long_name'] = 'latitude'
    cdt.attrs['units'] = 'degrees north (negative for degrees south)'
    cdt = dset.create_dataset('longitude',  data = output_lon[keep_idxs])
    cdt.attrs['long_name'] = 'longitude'
    cdt.attrs['units'] = 'degrees east (negative for degrees west)'
    cdt = dset.create_dataset('AOD_550',    data = output_aod[keep_idxs])
    cdt.attrs['long_name'] = 'lidar total AOD at 532 nm'
    cdt.attrs['units'] = 'none'
    cdt = dset.create_dataset('sfc_dwn_sw_flux',  data = output_sfc_dwn_sw_flx[keep_idxs])
    cdt.attrs['long_name'] = 'Surface downward SW flux'
    cdt.attrs['units'] = 'Wm^-2'
    cdt = dset.create_dataset('toa_up_sw_flux',   data = output_toa_up_sw_flx[keep_idxs])
    cdt.attrs['long_name'] = 'TOA upward SW flux'
    cdt.attrs['units'] = 'Wm^-2'
    cdt = dset.create_dataset('sfc_dwn_lw_flux',  data = output_sfc_dwn_lw_flx[keep_idxs])
    cdt.attrs['long_name'] = 'Surface downward LW flux'
    cdt.attrs['units'] = 'Wm^-2'
    cdt = dset.create_dataset('toa_up_lw_flux',   data = output_toa_up_lw_flx[keep_idxs])
    cdt.attrs['long_name'] = 'TOA upward LW flux'
    cdt.attrs['units'] = 'Wm^-2'
    cdt = dset.create_dataset('press_lvl',        data = output_press[keep_idxs])
    cdt.attrs['long_name'] = 'NAAPS Pressure'
    cdt.attrs['units'] = 'mb'
    cdt = dset.create_dataset('altitude',         data = output_alt[:])
    cdt.attrs['long_name'] = 'lidar bin altitude (above sea level)'
    cdt.attrs['units'] = 'm'
    cdt = dset.create_dataset('temperature',      data = output_tmp[keep_idxs])
    cdt.attrs['long_name'] = 'NAAPS air temperature'
    cdt.attrs['units'] = 'K'
    cdt = dset.create_dataset('rel_hum',          data = output_relhum[keep_idxs])
    if(l_USE_HALO_WV):
        cdt.attrs['long_name'] = 'HALO relative humidity (calculated using ' + \
            'HALO mixing ratio and NAAPS air pressure/temperature)'
    else:
        cdt.attrs['long_name'] = 'NAAPS relative humidity'
    cdt.attrs['units'] = '%'
    cdt = dset.create_dataset('spec_hum',         data = output_spchum[keep_idxs])
    if(l_USE_HALO_WV):
        cdt.attrs['long_name'] = 'HALO specific humidity (derived from ' + \
            'HALO mixing ratio)'
    else:
        cdt.attrs['long_name'] = 'NAAPS specific humidity'
    cdt.attrs['units'] = 'kg/kg'
    cdt = dset.create_dataset('abf_ext',          data = output_abf_ext[keep_idxs])
    cdt.attrs['long_name'] = 'lidar abf extinction'
    cdt.attrs['units'] = 'm^-1'
    cdt = dset.create_dataset('dust_ext',         data = output_dust_ext[keep_idxs])
    cdt.attrs['long_name'] = 'lidar dust extinction'
    cdt.attrs['units'] = 'm^-1'
    cdt = dset.create_dataset('smoke_ext',        data = output_smoke_ext[keep_idxs])
    cdt.attrs['long_name'] = 'lidar smoke extinction'
    cdt.attrs['units'] = 'm^-1'
    cdt = dset.create_dataset('salt_ext',         data = output_salt_ext[keep_idxs])
    cdt.attrs['long_name'] = 'lidar salt extinction'
    cdt.attrs['units'] = 'm^-1'
    cdt = dset.create_dataset('net_heating_rate', data = output_heat[keep_idxs])
    cdt.attrs['long_name'] = 'net daily heating rate'
    cdt.attrs['units'] = 'K/day'

    dset.close()
    print("Saved HDF5 output file: " + out_name)
