#!/usr/bin/env python

"""
  NAME:
    auto_fuliou_runner_lidar.py

  PURPOSE:
    Automate the running of FuLiou for each unique point in a passed
        CCPEXCV-HALO NAAPS/lidar file. This wrapper uses the
        lidar extinction variables rather than the NAAPS data.

  SYNTAX:
    python extract_naapslidar.py CCPEXCV-HALO_file

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2023/10/03:
      Written

"""

import sys
from datetime import datetime, timedelta
import h5py
import numpy as np
import os

# NOTES FROM JEFF
#
#       EXT up to 0.05 km-1 < 500 = sea salt, for now
#           anything above that (presume this means EXT > 0.05 km-1 and z < 500) = dust
#
#   Z < 500
#       EXT <= 0.05 km-1 = sea salt
#       EXT >  0.05 km-1 = dust
#   Z > 500
#       dust
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
                if(lidar_z[jj] < 500):
                    if(ext_val[ii,jj] <= 5e-5):
                        total_salt_aot += ext_val[ii,jj] * lidar_dz[jj]
                    else:
                        total_dust_aot += ext_val[ii,jj] * lidar_dz[jj]
                else: 
                    total_dust_aot += ext_val[ii,jj] * lidar_dz[jj]
                #total_other_aot += ext_val[ii,jj] * lidar_dz[jj]

    spec_aods = np.array([total_abf_aot, total_dust_aot, total_smoke_aot, \
        total_salt_aot])
    other_aods = np.array([total_ice_aot, total_other_aot])
        
    return spec_aods, other_aods


def test_aod_calc(ext_val, aer_id, aot_val, lidar_dz, abf_ids, dust_ids, \
        smoke_ids, salt_ids, ii):
    
    spec_aods, other_aods = \
        calc_species_aod(ext_val, aer_id, lidar_dz, abf_ids, dust_ids, \
        smoke_ids, salt_ids, ii)

    print(spec_aods, other_aods)

    calc_total_aod = np.sum(spec_aods)
    #calc_total_aod = total_abf_aot + total_dust_aot + total_smoke_aot + total_salt_aot

    calc_caliop_aod = np.sum((ext_val[ii,:][ext_val[ii,:] != 0.]) * 15)

    print('lidar AOD:              ',aot_val[ii])
    print('calc AOD from all bins   ',calc_caliop_aod)
    print('sum of 4-spec aods       ',calc_total_aod)

    return calc_caliop_aod, calc_total_aod 

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
        #naaps_dust_ext, naaps_abf_ext, naaps_salt_ext, naaps_smoke_ext, \
        naaps_dust_mass, naaps_abf_mass, naaps_salt_mass, naaps_smoke_mass, \
        ):

    # Identify the unique layer values
    # --------------------------------
    h_idxs = identify_unique_idxs(naaps_press)

    # Add a check for broken unique
    # -----------------------------
    #if(h_idxs[0] != 0): h_idxs[0] = 0

    res, ind = np.unique(naaps_press[:,0], return_index = True)
    v_idxs = np.append(ind[np.argsort(ind)], naaps_press[:,0].shape[0]-1)

    press_slope   = 0.
    temp_slope    = 0.
    relhum_slope  = 0.
    spechum_slope = 0.

    #dust_ext_slope   = 0.
    #abf_ext_slope    = 0.
    #salt_ext_slope   = 0.
    #smoke_ext_slope  = 0.

    dust_mass_slope  = 0.
    abf_mass_slope   = 0.
    salt_mass_slope  = 0.
    smoke_mass_slope = 0.
    
    # This is just used to keep track of which h_idxs value to use
    kk = 0
    
    interp_press   = np.zeros(naaps_press.shape[1])
    interp_temp    = np.zeros(naaps_press.shape[1])
    interp_relhum  = np.zeros(naaps_press.shape[1])
    interp_spechum = np.zeros(naaps_press.shape[1])

    #interp_dust_ext  = np.zeros(naaps_press.shape[1])
    #interp_abf_ext   = np.zeros(naaps_press.shape[1])
    #interp_salt_ext  = np.zeros(naaps_press.shape[1])
    #interp_smoke_ext = np.zeros(naaps_press.shape[1])
    
    interp_dust_mass  = np.zeros(naaps_press.shape[1])
    interp_abf_mass   = np.zeros(naaps_press.shape[1])
    interp_salt_mass  = np.zeros(naaps_press.shape[1])
    interp_smoke_mass = np.zeros(naaps_press.shape[1])
   
    # This loops over the vertical indices. Can use this inside the
    # wrapper below.
    for jj in range(naaps_press.shape[1]):
    
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

            #dust_ext_slope = calc_regress_slope(\
            #    naaps_dust_ext, h_idxs, kk, xx = ii)
            #abf_ext_slope = calc_regress_slope(\
            #    naaps_abf_ext, h_idxs, kk, xx = ii)
            #salt_ext_slope = calc_regress_slope(\
            #    naaps_salt_ext, h_idxs, kk, xx = ii)
            #smoke_ext_slope = calc_regress_slope(\
            #    naaps_smoke_ext, h_idxs, kk, xx = ii)
    
            dust_mass_slope = calc_regress_slope(\
                naaps_dust_mass, h_idxs, kk, xx = ii)
            abf_mass_slope = calc_regress_slope(\
                naaps_abf_mass, h_idxs, kk, xx = ii)
            salt_mass_slope = calc_regress_slope(\
                naaps_salt_mass, h_idxs, kk, xx = ii)
            smoke_mass_slope = calc_regress_slope(\
                naaps_smoke_mass, h_idxs, kk, xx = ii)
    
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

        #interp_dust_ext[jj]   = dust_ext_slope   * (jj - h_idxs[kk-1]) + \
        #    naaps_dust_ext[ii, h_idxs[kk-1]]
        #interp_abf_ext[jj]   = abf_ext_slope     * (jj - h_idxs[kk-1]) + \
        #    naaps_abf_ext[ii, h_idxs[kk-1]]
        #interp_salt_ext[jj]   = salt_ext_slope   * (jj - h_idxs[kk-1]) + \
        #    naaps_salt_ext[ii, h_idxs[kk-1]]
        #interp_smoke_ext[jj]   = smoke_ext_slope * (jj - h_idxs[kk-1]) + \
        #    naaps_smoke_ext[ii, h_idxs[kk-1]]
    
        interp_dust_mass[jj]   = dust_mass_slope  * (jj - h_idxs[kk-1]) + \
            naaps_dust_mass[ii, h_idxs[kk-1]]
        interp_abf_mass[jj]   = abf_mass_slope    * (jj - h_idxs[kk-1]) + \
            naaps_abf_mass[ii, h_idxs[kk-1]]
        interp_salt_mass[jj]   = salt_mass_slope  * (jj - h_idxs[kk-1]) + \
            naaps_salt_mass[ii, h_idxs[kk-1]]
        interp_smoke_mass[jj]   = smoke_mass_slope * (jj - h_idxs[kk-1]) + \
            naaps_smoke_mass[ii, h_idxs[kk-1]]
    
    return interp_press, interp_temp, interp_relhum, interp_spechum, \
        interp_dust_mass, interp_abf_mass, interp_salt_mass, interp_smoke_mass
        #interp_dust_ext, interp_abf_ext, interp_salt_ext, interp_smoke_ext, \

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

# Check command line arguments
# ----------------------------
if(len(sys.argv) != 2):
    print("SYNTAX: python extract_naapslidar.py CCPEXCV-HALO_file")
    sys.exit(1)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Set path variables
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#base_dir = "/Research/fuliou_lidarflux/fuliou_package/"
base_dir = "/Research/fuliou_lidarflux/test_package/"
data_dir = "/Research/NAAPS_NOGAPS_2023/"
mhome = base_dir + "test_run"
output_dir = mhome + '/data'

homedir = mhome

########################################################
#Parameter Group 2

########################################################
#locations for executables
#naaps_fuliou
#setenv NAAPS_FULIOUCODE "/home/jzhang/fuliou_package/sourcecodes"
naaps_fulioucode = base_dir + "sourcecodes"

tabledir  = naaps_fulioucode + "/table"             ##place to store .txt files needed for the program
srcdir  =  naaps_fulioucode 

#tools
mtools = naaps_fulioucode + "/tools"
#setenv MTOOLS      "$NAAPS_FULIOUCODE/tools"


########################################################
#locations for data

#NOGAPS data
#setenv NOGAPS_SOURCE "/home/jzhang/fuliou_package/test_data/test_nogaps"
#setenv NOGAPS_SOURCE "/Research/for_blake_from_zhang/fuliou_package/test_data/test_nogaps"
#setenv NOGAPS_SOURCE "/Research/NAAPS_NOGAPS_2023/NAAPS/"

#setenv NOGAPS_SOURCE "/Research/NAAPS_NOGAPS_2023/nog_cmorph"
nogaps_source = data_dir + "nog_cmorph"

#NAAPS DATA
#setenv naapsdir  "/home/jzhang/fuliou_package/test_data/test_naaps"
#setenv naapsdir  "/Research/for_blake_from_zhang/fuliou_package/test_run/test_naaps"

#setenv naapsdir  "/Research/NAAPS_NOGAPS_2023/NAAPS"
naapsdir = data_dir + "NAAPS"
concdir = os.getcwd()
nogapsdir = os.getcwd()
predir = os.getcwd()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Write namelist file
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

namelist_file = 'namelist'
dtg = '2000010100'
#set namelist = "${runhome}/namelist"
#/bin/rm -r $namelist >& /dev/null
with open(namelist_file, 'w') as fname:
    fname.write("    &nl\n" + \
                "    tape        = '" + dtg + "',\n" + \
                "    tabledir    = '" + tabledir + "',\n" + \
                "    sulfatefile = 'anthro_farop_v3.txt',\n" + \
                "    dustfile    = 'dust_farop_v3.txt',\n" + \
                "    smokefile   = 'smoke_farop_v3.txt',\n" + \
                "    saltfile    = 'salt_farop_v3.txt',\n" + \
                "    concdir     = '" + concdir + "',\n" + \
                "    nogapsdir   = '" + nogapsdir + "',\n" + \
                "    predir      = '" + predir + "',\n" + \
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

#naaps_dust_ext  = hdf_data['NAAPS_dust_ext'][:,:]
naaps_dust_aot  = hdf_data['NAAPS_dust_aot'][:]
naaps_dust_mass = \
    hdf_data['NAAPS_dust_mass_concentration'][:,:] * 1e-9

#naaps_abf_ext   = hdf_data['NAAPS_ABF_ext'][:,:]
naaps_abf_aot   = hdf_data['NAAPS_ABF_aot'][:]
naaps_abf_mass  = \
    hdf_data['NAAPS_ABF_mass_concentration'][:,:] * 1e-9

#naaps_salt_ext  = hdf_data['NAAPS_SS_ext'][:,:]
naaps_salt_aot  = hdf_data['NAAPS_SS_aot'][:]
naaps_salt_mass = \
    hdf_data['NAAPS_SS_mass_concentration'][:,:] - 1e-9

#naaps_smoke_ext  = hdf_data['NAAPS_smoke_ext'][:,:]
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
# ---------------------------------------------------------------------------
ext_532       = hdf_data['532_ext'][:,:] / 1e3
aot_532       = hdf_data['532_AOT_lo'][:]
aer_id        = hdf_data['Aerosol_ID'][:,:]

ext_532 = np.ma.where(np.ma.masked_invalid(ext_532).mask == True, \
    0., ext_532)


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

#fmt_str = "{0:7.1f} {1:8.1f} {2:8.1f} {3:7.2f} " + \
#    "{4:7.3f} {5:10.3e} " + \
#    "{6:10.3e} {7:10.3e} {8:10.3e} {9:10.3e} " + \
#    "{10:10.3e} {11:10.3e} {12:10.3e} {13:10.3e}\n"

# totaod = [aod1, aod2, aod3, aod4, aod5]
# totaod: [UNUSED, sulfate, dust, smoke, salt]

#foutname = base_dir + \
#            'test_data/test_naaps_new/test_naaps_file_lidar.txt' 
foutname = 'test_naaps_file.txt'

#out_data_name = output_dir + '/test_output_file_lidar.txt'
out_data_name = 'test_output_file_lidar.txt'

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
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

abf_ids = [4.]
dust_ids = [2., 8.]
smoke_ids = [5., 6.]
salt_ids  = [3., 7.]

#good_idxs = np.where(np.ma.masked_invalid(aot_532).mask == False)[0]

#calc_caliop_aod, calc_total_aod = \
#    test_aod_calc(ext_532, aer_id, aot_532, lidar_dz, abf_ids, dust_ids, \
#        smoke_ids, salt_ids, good_idxs[100])

##!#work_file = 'test_comp_work.txt'
##!#if(os.path.exists(work_file)):
##!#    print("Removing work file")
##!#    cmnd = 'rm ' + work_file
##!#    os.system(cmnd)

check_press = -999.
if(l_BIN_LIDAR_VERTICAL):
    num_layers = 100
    split_dz = np.array(np.split(lidar_dz[:len(lidar_dz) - (len(lidar_dz) % num_layers)], num_layers))
    avg_dz   = np.mean(np.split(lidar_dz[:len(lidar_dz) - (len(lidar_dz) % num_layers)], num_layers), axis = 1)
    avg_z    = np.mean(np.array(np.split(lidar_z[:len(lidar_dz) - (len(lidar_dz) % num_layers)], num_layers)), axis = 1)
#for ii in range(time.shape[0]):
#for ii in range(0,time.shape[0],100):
for ii in range(1223,1224):

    #print(ii, naaps_press[ii,0])
    if(naaps_press[ii,0] != check_press):

        print("    interpolating")
        # Interpolate the values here
        # ----------------------------
        interp_press, interp_temp, interp_relhum, interp_spechum, \
            interp_dust_mass, interp_abf_mass, interp_salt_mass, \
            interp_smoke_mass = \
            interp_NAAPS(ii, naaps_press, naaps_temp, naaps_relhum, \
                naaps_spechum, \
                #naaps_dust_ext, naaps_abf_ext, naaps_salt_ext, \
                #naaps_smoke_ext, 
                naaps_dust_mass, naaps_abf_mass, \
                naaps_salt_mass, naaps_smoke_mass)
            #interp_dust_ext, interp_abf_ext, interp_salt_ext, interp_smoke_ext,\

        if(l_BIN_LIDAR_VERTICAL):
            # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
            #
            # Play around with thinning the data to 100 vertical layers
            #
            # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
            split_press       = np.mean(np.array(np.split(interp_press[:len(interp_press) - (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_temp        = np.mean(np.array(np.split(interp_temp[:len(interp_press) - (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_relhum      = np.mean(np.array(np.split(interp_relhum[:len(interp_press) - (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_spechum     = np.mean(np.array(np.split(interp_spechum[:len(interp_press) - (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_dust_mass   = np.mean(np.array(np.split(interp_dust_mass[:len(interp_press) - (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_abf_mass    = np.mean(np.array(np.split(interp_abf_mass[:len(interp_press) - (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_salt_mass   = np.mean(np.array(np.split(interp_salt_mass[:len(interp_press) - (len(interp_press) % num_layers)], num_layers)), axis = 1)
            split_smoke_mass  = np.mean(np.array(np.split(interp_smoke_mass[:len(interp_press) - (len(interp_press) % num_layers)], num_layers)), axis = 1)
    

        check_press = naaps_press[ii,0]

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Test the correction to the aerosol ID problem
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

    if(l_BIN_LIDAR_VERTICAL): 
        split_ext     = np.array(np.split(ext_532[ii,:len(interp_press) - (len(interp_press) % num_layers)], num_layers))
        split_aerid   = np.array(np.split(aer_id[ii,:len(interp_press) - (len(interp_press) % num_layers)], num_layers))
   
        avg_exts = np.average(split_ext, axis = 1)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Hand-calculate the species-specific AOTs here. These will replace the 
    #   'naaps_abf_aot'-like variables printed at the top of the output file
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

    spec_aods, other_aods = \
        calc_species_aod(ext_532, aer_id, lidar_z, lidar_dz, abf_ids, \
            dust_ids, smoke_ids, salt_ids, ii)

    ##!#total_abf_aot = total_dust_aot = total_smoke_aot = total_salt_aot = \
    ##!#    total_other_aot = 0.
    ##!#for jj in range(ext_532.shape[1]-1, -1, -1):
    ##!#    if(ext_532[ii,jj] != -9.):
    ##!#    #if(ext_532[ii,jj] == 0):
    ##!#       if(aer_id[ii,jj] in abf_ids):
    ##!#           total_abf_aot += ext_532[ii,jj] * lidar_dz[jj]
    ##!#       # Dust
    ##!#       elif(aer_id[ii,jj] in dust_ids):
    ##!#           total_dust_aot += ext_532[ii,jj] * lidar_dz[jj]
    ##!#       # Smoke
    ##!#       elif(aer_id[ii,jj] in smoke_ids):
    ##!#           total_smoke_aot += ext_532[ii,jj] * lidar_dz[jj]
    ##!#       # Salt 
    ##!#       elif(aer_id[ii,jj] in salt_ids):
    ##!#           total_salt_aot += ext_532[ii,jj] * lidar_dz[jj]
    ##!#       # Other - either "ice" or unclassified
    ##!#       else:
    ##!#           total_other_aot += ext_532[ii,jj] * lidar_dz[jj]
          
    calc_total_aot = np.sum(spec_aods)
    #calc_total_aot += total_other_aot

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
        local_date = local_dt_date.strftime('%Y%m%d%H%M')
        local_fmt = "{0:7.4f} {1:7.4f} {2:7.4f} {3:7.4f} {4:7.4f}"
        print(ii, local_dt_date.strftime('%Y%m%d %H:%M:%S'),\
            local_fmt.format(aot_532[ii], \
                calc_total_aot, other_aods[1] + calc_total_aot, other_aods[0], other_aods[1] ))

        if(l_USE_lidar_AOD):
            fout.write("{0} {1:8.5f} {2:8.5f} {3:8.5f}\n".format(\
                local_date, lat[ii], lon[ii], aot_532[ii]))
                #local_date, lat[ii], lon[ii], naaps_tot_aot[ii]))
        else:
            fout.write("{0} {1:8.5f} {2:8.5f} {3:8.5f}\n".format(\
                local_date, lat[ii], lon[ii], calc_total_aot))
                #local_date, lat[ii], lon[ii], naaps_tot_aot[ii]))
        
        fout.write("{0:10.3e} {1:10.3e} {2:10.3e} {3:10.3e}\n".format(\
            spec_aods[0], spec_aods[1], \
            spec_aods[2], spec_aods[3]))
            #naaps_abf_aot[ii], naaps_dust_aot[ii], \
            #naaps_smoke_aot[ii], naaps_salt_aot[ii]))
        
        # Print the profile information in reverse order, with the lowest
        # data in the profile being at the last index.
        # ----------------------------------------------------------------------
        if(l_BIN_LIDAR_VERTICAL):
            vert_dim = avg_exts.shape[0] - 1
        else:
            vert_dim = ext_532.shape[1] - 1
        for jj in range(vert_dim, -1, -1):

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
            local_abf = local_dust = local_smoke = local_salt = 0.

            if(l_BIN_LIDAR_VERTICAL):
                local_ids = np.ma.masked_invalid(split_aerid[jj])
                local_ids = np.where(local_ids.mask, -9, local_ids)
                if(avg_z[jj] > 500):
                    local_ids = np.where(local_ids == -9, 2, local_ids)
                else:
                    local_ids = np.where((local_ids == -9) & (split_ext[jj] <= 5e-5), 3, local_ids)
                    local_ids = np.where((local_ids == -9) & (split_ext[jj] > 5e-5), 2, local_ids)
                ###if(lidar_z[jj] < 500):
                ###    if(ext_val[ii,jj] <= 5e-5):
                ###        total_salt_aot += ext_val[ii,jj] * lidar_dz[jj]
                ###    else:
                ###        total_dust_aot += ext_val[ii,jj] * lidar_dz[jj]
                ###else: 
                ###    total_dust_aot += ext_val[ii,jj] * lidar_dz[jj]

                unique_ids = np.unique(local_ids)
                work_id = unique_ids[np.argmax([\
                    len(np.where(local_ids == uq)[0]) for uq in unique_ids])]
                work_ext = avg_exts[jj]
            else:
                work_id = aer_id[ii,jj]
                work_ext = ext_532[ii,jj]

            # ABF
            #if(aer_id[ii,jj] == 4.):
            #if(aer_id[ii,jj] in abf_ids):
            #    local_abf = ext_532[ii,jj]
            if(work_id in abf_ids):
                local_abf = work_ext
                if(work_ext != 0.0):
                    print(ii, jj, local_abf, local_dust, local_smoke, local_salt, 'ABF')
            # Dust
            #elif( (aer_id[ii,jj] == 2.) | (aer_id[ii,jj] == 8.) ):
            #elif(aer_id[ii,jj] in dust_ids):
            #    local_dust = ext_532[ii,jj]
            elif(work_id in dust_ids):
                local_dust = work_ext
            # Smoke
            #elif(aer_id[ii,jj] in smoke_ids):
            #    local_smoke = ext_532[ii,jj]
            elif(work_id in smoke_ids):
                local_smoke = work_ext
                if(work_ext != 0.0):
                    print(ii, jj, local_abf, local_dust, local_smoke, local_salt, 'SMOKE')
            # Salt 
            #elif(aer_id[ii,jj] in salt_ids):
            #    local_salt = ext_532[ii,jj]
            elif(work_id in salt_ids):
                local_salt = work_ext
                if(work_ext != 0.0):
                    print(ii, jj, local_abf, local_dust, local_smoke, local_salt, 'SALT')


            if(l_BIN_LIDAR_VERTICAL):
                fout.write(fmt_str.format(\
                    split_press[jj], float(avg_z[jj]), avg_dz[jj], \
                    split_temp[jj], split_relhum[jj], \
                    split_spechum[jj], \
                    local_abf, local_dust, local_smoke, local_salt, \
                    #interp_abf_ext[jj], interp_dust_ext[jj], \
                    #interp_smoke_ext[jj], interp_salt_ext[jj], \
                    split_abf_mass[jj], split_dust_mass[jj], \
                    split_smoke_mass[jj], split_salt_mass[jj], \
                ))
            else:
                fout.write(fmt_str.format(\
                    interp_press[jj], float(lidar_z[jj]), lidar_dz[jj], \
                    interp_temp[jj], interp_relhum[jj], \
                    interp_spechum[jj], \
                    local_abf, local_dust, local_smoke, local_salt, \
                    #interp_abf_ext[jj], interp_dust_ext[jj], \
                    #interp_smoke_ext[jj], interp_salt_ext[jj], \
                    interp_abf_mass[jj], interp_dust_mass[jj], \
                    interp_smoke_mass[jj], interp_salt_mass[jj], \
                ))
   
    ###cmnd = 'cp ' + foutname + ' ' + base_dir + \
    ###        'test_naaps_file_lidar_' + str(int(ii)) + '.txt'
    ###print(cmnd)
    ###os.system(cmnd)
    #os.system('cp ' + foutname + ' ' + base_dir + \
    #        'test_data/test_naaps_new/test_naaps_file_lidar_' + str(int(ii)) + '.txt')

    # Copy the first couple of lines to the work file 
    # -----------------------------------------------
    cmnd = 'head -n 2 ' + foutname + ' >> ' + work_file
    print(cmnd)
    os.system(cmnd)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #  Call the FuLiou code
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #print("Calling FuLiou code") 
    print(' Execute naaps_fuliou.exe')
    cmnd = 'time ' + srcdir + '/bin/naaps_fuliou.exe'
    print(cmnd)
    os.system(cmnd)

    # Concatenate the temp output file to the end of the actual
    # output file
    # ---------------------------------------------------------
    cmnd = 'cat fort.22 >> ' + out_data_name
    print(cmnd)
    os.system(cmnd)

    cmnd = 'head -n 1 fort.22 >> ' + work_file
    print(cmnd)
    os.system(cmnd)
