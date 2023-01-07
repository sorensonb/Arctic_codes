#!/usr/bin/env python

"""


"""

import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import sys
import metpy.calc as mpcalc
from metpy.units import units
import scipy
from copy import deepcopy
from matplotlib.lines import Line2D

home_dir = os.environ['HOME']

sys.path.append(home_dir + '/')
from python_lib import *

h_const = 6.626e-34 #J*s
k_const = 1.3806e-23 #J/K
c_const = 3e8 # m/s
sb_path = home_dir + '/Research/SBDART/'

# wavel in microns
# radiance in (W/(m2*μm*Sr))
def inverse_planck(radiance, wavel):
    return ((h_const * c_const) / (wavel*1e-6 * k_const)) * \
        (np.log( ((2.*h_const * (c_const**2.)) / (wavel * (wavel * 1e-6)**4. \
        * radiance)) + 1.))**-1.

# radiance is a function that calculates the blackbody radiance for a given 
# wavelength (in meters) and temperature (in Kelvin)
# -------------------------------------------------------------------------
def radiance(wavelength,tmp):
    return ((2.*h_const*(c_const**2.)) / \
           ((wavelength**5.)*(np.exp((h_const*c_const)/\
            (wavelength*k_const*tmp)) - 1.))) * \
           1e-6

# NOTE: modify to use actual CrIS ozone profile?
def prep_cris_profile(CrIS_prof):
    
    # Read in the default data for oxygen data
    # ----------------------------------------
    atms_file = sb_path + 'src/atms_mdlat_sum.dat'
    temp_data = pd.read_csv(atms_file, names = \
        ['z','p','t','wv','o3'], delim_whitespace = True)

    # Mask any data that have nan temps (indicating below-sfc data)
    # -------------------------------------------------------------
    bottom_idx = np.where(np.isnan(CrIS_prof['temp']))[0][0]
    local_temp  = CrIS_prof['temp'][:bottom_idx]
    local_wv    = CrIS_prof['wv'][:bottom_idx]
    local_rh    = CrIS_prof['rh'][:bottom_idx]
    local_press = CrIS_prof['press'][:bottom_idx]

    # Convert the CrIS pressures to heights using standard atmosphere
    # ---------------------------------------------------------------
    calc_hgt = np.array(mpcalc.pressure_to_height_std(local_press * units('hPa'))) 

    # Convert the mixing ratio to water vapor density (g/m3)
    # ------------------------------------------------------
    vap_press = mpcalc.vapor_pressure(local_press * units('hPa'),\
        local_wv * units('g/kg'))
    abs_hum = (np.array(vap_press.to('Pa')) / \
        (461.5 * local_temp)) * 1e3

    # Interpolate the ozone values from the default profile
    # to match the new pressures here
    # -----------------------------------------------------
    def_o3 = pd.to_numeric(temp_data['o3']).values
    def_p  = pd.to_numeric(temp_data['p']).values

    o3_interp1d = scipy.interpolate.interp1d(def_p, def_o3,\
        kind = 'linear')
    new_o3 = o3_interp1d(local_press) 
 
    # Insert the wanted variables into a pandas dataframe
    # ---------------------------------------------------
    data = pd.DataFrame()
    data['z'] = calc_hgt
    data['p'] = local_press
    data['t'] = local_temp
    data['wv'] = abs_hum
    data['mix_rat'] = local_wv
    data['o2'] = new_o3
    data['rhum'] = local_rh

    return data

# Reads the SBDART model profile. If infile is '', the default
# profile will be used. If not, a model profile from 
# home_dir + /Research/SBDART/model/ and reformats into the 
# right format.
def read_atmos_profile(infile = ''):
    atms_file = sb_path + 'src/atms_mdlat_sum.dat'
    #atms_file = sb_path + 'src/atms_orig.dat'
    temp_data = pd.read_csv(atms_file, names = \
        ['z','p','t','wv','o2'], delim_whitespace = True)
    #temp_data = temp_data[1:]

    if((infile == '') | (infile.strip().split('/')[-1][:4] == 'atms')):
        infile = atms_file
        print('HERE: Reading atmosphere from ', atms_file)
        data = temp_data

        # Convert wv density from g/m3 to kg/m3
        # -------------------------------------
        wv_dense_kgm3 = data['wv'] * 1e-3       
 
        # Calculate vapor pressure from wv density
        # ----------------------------------------
        vap_press = wv_dense_kgm3 * 461.5 * data['t']

        # Calculate mixing ratio from vapor pressure and prssure, convert 
        # to g/kg
        # ---------------------------------------------------------------
        mix_rat = 0.622 * (vap_press / (data['p'] * 100. - vap_press)) * 1e3

        data['mix_rat'] = mix_rat

        # Calculate relative humidity and add to structure
        # ------------------------------------------------
        data['rhum'] = mpcalc.relative_humidity_from_mixing_ratio(\
            np.array(data['p']) * units('mbar'), np.array(data['t']) * units('K'), \
            np.array(mix_rat) / 1e3) * 100.

        ##!#infile = sb_path + 'src/atms_orig.dat'
        ##!#data = pd.read_csv(infile,skiprows = 0, names = \
        ##!#    ['z','p','t','wv','o2'], delim_whitespace = True)
    else:
        # Read the file lines with readlines
        # ----------------------------------
        print('Reading atmosphere from ', infile)
        with open(infile,'r') as fin:
            flines = fin.readlines()
        np_lines = np.array(flines)

        # Find where each of the data columns starts
        # ------------------------------------------
        starters = np.array([tline.startswith('Hmsl') for tline in flines])
        s_idx    = np.where(starters == True)[0] + 2

        # To account for the extra surface line at the beginning
        # of the profile, add one to this line
        s_idx[0] += 1

        # Find where each of the data columns ends
        # ----------------------------------------
        enders = np.array([tline.startswith(' ==') for tline in flines])
        e_idx  = np.where(enders == True)[0][2:] - 1
        e_idx = np.append(e_idx[:], len(flines))

        # Extract all the variables         
        # -------------------------
        hgt_z = np.array([float(tline.strip().split()[0])\
            for tline in np_lines[s_idx[0]:e_idx[0]]])      
        tmp_c = np.array([float(tline.strip().split()[1])\
            for tline in np_lines[s_idx[0]:e_idx[0]]])      
        rhum  = np.array([float(tline.strip().split()[1])\
            for tline in np_lines[s_idx[1]:e_idx[1]]])      
        wdir  = np.array([float(tline.strip().split()[1])\
            for tline in np_lines[s_idx[2]:e_idx[2]]])      
        wspd  = np.array([float(tline.strip().split()[1])\
            for tline in np_lines[s_idx[3]:e_idx[3]]])      

        # Convert the heights to pressure using the standard
        # atmosphere
        # --------------------------------------------------
        p_std = np.array(mpcalc.height_to_pressure_std(hgt_z * \
            units('m')))

        # Interpolate the ozone values from the default profile
        # to match the new pressures here
        # -----------------------------------------------------
        def_o2 = pd.to_numeric(temp_data['o2']).values
        def_p  = pd.to_numeric(temp_data['p']).values

        o2_interp1d = scipy.interpolate.interp1d(def_p, def_o2,\
            kind = 'linear')
        new_o2 = o2_interp1d(p_std) 

        # Calculate absolute humidity from the other variables
        # ----------------------------------------------------
        mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(\
            p_std * units('mbar'), tmp_c * units('degC'), \
            rhum / 100.)
        vap_press = mpcalc.vapor_pressure(p_std * units('mbar'),\
            mixing_ratio)
        abs_hum = (np.array(vap_press.to('Pa')) / \
            (461.5 * (tmp_c + 273.15))) * 1e3

        # Insert the wanted variables into a pandas dataframe
        # ---------------------------------------------------
        data = pd.DataFrame()
        data['z'] = hgt_z / 1e3
        data['p'] = p_std
        data['t'] = tmp_c + 273.15
        data['wv'] = abs_hum
        data['mix_rat'] = np.array(mixing_ratio.to('g/kg'))
        data['o2'] = new_o2
        data['rhum'] = rhum

    return data 

def prep_filter_function(satellite):
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Filter function prep
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    # Copy the filter function to working directory
    # ---------------------------------------------
    os.system('cp ' + sb_path + 'data/filter/' + satellite + '_filter.dat ./filter.dat')
    
    # Modify the filter file
    # ----------------------
    with open('filter.dat','r') as fin:
        indata = fin.readlines()
    with open('filter.dat','w') as fout:
        fout.writelines(indata[4:])
    
    # Read with pandas and convert to microns
    # ---------------------------------------
    ff_data = pd.read_csv('filter.dat', names = ['wavenum','filter'], \
        delim_whitespace = True)
    ff_data['wavenum'] = (1. / (ff_data['wavenum'] * 100.)) * 1e6
    ff_data = ff_data.reindex(index = ff_data.index[::-1])

    # If the filter function is greater than 1000 elements long, rescale
    # the pandas dataframe to have less than 1000 elements
    # ------------------------------------------------------------------
    if(len(ff_data) > 1000):
        ff_data = ff_data[ff_data.index % int(np.ceil(len(ff_data) / 1000)) == 0]
    
    ff_data.to_csv('filter.dat', \
            index = False, sep = ' ', header = None)

def prep_sbdart_input(isat = -1, iout = 5, nzen = None, uzen = [0, 60],\
        set_co2_mix = None, add_co2_mix = None, \
        set_ch4_mix = None, add_ch4_mix = None, \
        set_no2_mix = None, add_no2_mix = None, \
        set_co_mix = None,  add_co_mix = None, \
        set_nh3_mix = None, add_nh3_mix = None, \
        btemp = None):
 
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Input file prep
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   
    print('Prepping INPUT file')
 
    #iout = 20
    #iout = 0

    # Prep the CO2 mixing ratio
    # -------------------------
    if((set_co2_mix is not None) & (add_co2_mix is None)):
        xco2 = set_co2_mix
    elif((set_co2_mix is None) & (add_co2_mix is not None)):
        xco2 = 360 + add_co2_mix
    elif((set_co2_mix is None) & (add_co2_mix is None)):
        xco2 = -1
    else:
        print("WARNING: Invalid co2 options. Setting to default")
        xco2 = -1

    # Prep the CH4 mixing ratio
    # -------------------------
    if((set_ch4_mix is not None) & (add_ch4_mix is None)):
        xch4 = set_ch4_mix
    elif((set_ch4_mix is None) & (add_ch4_mix is not None)):
        xch4 = 1.74 + add_ch4_mix
    elif((set_ch4_mix is None) & (add_ch4_mix is None)):
        xch4 = -1
    else:
        print("WARNING: Invalid ch4 options. Setting to default")
        xch4 = -1

    # Prep the NO2 mixing ratio
    # -------------------------
    if((set_no2_mix is not None) & (add_no2_mix is None)):
        xno2 = set_no2_mix
    elif((set_no2_mix is None) & (add_no2_mix is not None)):
        xno2 = 2.3e-5 + add_no2_mix
    elif((set_no2_mix is None) & (add_no2_mix is None)):
        xno2 = -1
    else:
        print("WARNING: Invalid no2 options. Setting to default")
        xno2 = -1

    # Prep the CO mixing ratio
    # -------------------------
    if((set_co_mix is not None) & (add_co_mix is None)):
        xco = set_co_mix
    elif((set_co_mix is None) & (add_co_mix is not None)):
        xco = 0.15 + add_co_mix
    elif((set_co_mix is None) & (add_co_mix is None)):
        xco = -1
    else:
        print("WARNING: Invalid co options. Setting to default")
        xco = -1

    # Prep the NH3 mixing ratio
    # -------------------------
    if((set_nh3_mix is not None) & (add_nh3_mix is None)):
        xnh3 = set_nh3_mix
    elif((set_nh3_mix is None) & (add_nh3_mix is not None)):
        xnh3 = 5.0e-4 + add_nh3_mix
    elif((set_nh3_mix is None) & (add_nh3_mix is None)):
        xnh3 = -1
    else:
        print("WARNING: Invalid nh3 options. Setting to default")
        xnh3 = -1


 
    with open('./INPUT','w') as fout:
        #fout.write(" &INPUT\n" + \
        out_str = " &INPUT\n" + \
            "   idatm=0\n" + \
            "   isat={0}\n".format(isat) + \
            "   wlinc=0.01\n" + \
            "   iout={0}\n".format(iout) + \
            "   xco2={0}\n".format(xco2) + \
            "   xch4={0}\n".format(xch4) + \
            "   xno2={0}\n".format(xno2) + \
            "   xco={0}\n".format(xco) + \
            "   xnh3={0}\n".format(xnh3)
            #"   wlinf=8.00\n" + \
            #"   wlsup=12.00\n" + \
            #"   nzen={0}\n".format(nzen) + \
        #print(out_str)
        if(uzen is not None):
            if(isinstance(uzen, list)):
                if(nzen is not None):
                    out_str = out_str + \
                        "   nzen={0}\n".format(nzen)
                out_str = out_str + "   uzen=" + \
                    ','.join([str(uz) for uz in uzen]) + "\n"
#                out_str = out_str + \
#                   "   uzen={0},{1}\n".format(uzen[0], uzen[1])
            else:
                out_str = out_str + \
                   "   uzen={0}\n".format(uzen)
            #"   vzen=180\n" + \
        if(btemp is not None):
            out_str = out_str + \
                "   btemp={0}\n".format(int(btemp))
        out_str = out_str + \
            "   phi=0\n" + \
            "/\n"

        fout.write(out_str)

def check_max_rhum(data):
    # See if any regions of this profile exceed 80% relative humidity
    num_exceed = len(np.where(data['rhum'] > 80.)[0])
    if(num_exceed > 0):
        # If any are over 4 g/kg, set the maximum to 4 and continue
        data['rhum'][data['rhum'] > 80.] = 80.

        # Recalculate mixing ratios using the 80% max RH
        data['mix_rat'][data['rhum'] > 80.] = \
            mpcalc.mixing_ratio_from_relative_humidity(\
            np.array(data['p'][data['rhum'] > 80.]) * units('mbar'), \
            np.array(data['t'][data['rhum'] > 80.]) * units('K'), \
            np.array(data['rhum'][data['rhum'] > 80.]) * units('%')).to('g/kg')
    
    return data


def prep_atmos_file(data, \
        zmin = None, zmax = None, \
        set_tmp = None, add_tmp = None, add_tmp_sfc = None,\
        set_wv_mix = None, add_wv_mix = None, \
        set_rh = None, \
        out_path = './', saver = None):
    #print("processing for wv multiplier = ",np.round(mtp,2))

    # Make copies for editing
    # -----------------------
    new_data = data.copy(deep = True)

    if(set_tmp is not None):
        # Set the layer temps
        # -------------------
        new_data['t'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = set_tmp

        # Recalculate relative humidity in these layers
        # ---------------------------------------------
        new_rh = mpcalc.relative_humidity_from_mixing_ratio(\
            np.array(new_data['p'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('mbar'), \
            np.array(new_data['t'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('K'), \
            np.array(new_data['mix_rat'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) / 1e3) * 100.
        new_data['rhum'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = new_rh

    if(add_tmp_sfc is not None):
        new_data['t'][:2] += add_tmp_sfc

        # Recalculate relative humidity in these layers
        # ---------------------------------------------
        new_rh = mpcalc.relative_humidity_from_mixing_ratio(\
            np.array(new_data['p'][:2]) * units('mbar'), \
            np.array(new_data['t'][:2]) * units('K'), \
            np.array(new_data['mix_rat'][:2]) / 1e3) * 100.

        new_data['rhum'][:2] = new_rh

    elif(add_tmp is not None): 
        print("Adding temp of ", add_tmp)
        # Add the layer temps
        # -------------------
        new_data['t'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] += add_tmp

        # Recalculate relative humidity in these layers
        # ---------------------------------------------
        new_rh = mpcalc.relative_humidity_from_mixing_ratio(\
            np.array(new_data['p'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('mbar'), \
            np.array(new_data['t'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('K'), \
            np.array(new_data['mix_rat'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) / 1e3) * 100.
        new_data['rhum'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = new_rh
 
    if(set_wv_mix is not None):  
        new_data['mix_rat'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = set_wv_mix

        # Recalculate relative humidity in these layers
        # ---------------------------------------------
        new_rh = mpcalc.relative_humidity_from_mixing_ratio(\
            np.array(new_data['p'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('mbar'), \
            np.array(new_data['t'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('K'), \
            np.array(new_data['mix_rat'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) / 1e3) * 100.
        new_data['rhum'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = new_rh

        new_data = check_max_rhum(new_data)

        ##!## See if any regions of this profile exceed 80% relative humidity
        ##!#num_exceed = len(np.where(new_data['rhum'] > 80.)[0])
        ##!#if(num_exceed > 0):
        ##!#    # If any are over 4 g/kg, set the maximum to 4 and continue
        ##!#    new_data['rhum'][new_data['rhum'] > 80.] = 80.

        ##!#    # Recalculate mixing ratios using the 80% max RH
        ##!#    new_data['mix_rat'][new_data['rhum'] > 80.] = \
        ##!#        mpcalc.mixing_ratio_from_relative_humidity(\
        ##!#        np.array(new_data['p'][new_data['rhum'] > 80.]) * units('mbar'), \
        ##!#        np.array(new_data['t'][new_data['rhum'] > 80.]) * units('K'), \
        ##!#        np.array(new_data['rhum'][new_data['rhum'] > 80.]) * units('%')).to('g/kg')

        # Convert the added water vapor mixing ratio to added
        # absolute humidity
        # ----------------------------------------------------
        vap_press = mpcalc.vapor_pressure(\
            np.array(new_data['p']) * units('mbar'), \
            np.array(new_data['mix_rat']) * units('g/kg'))
        new_hum = (vap_press.to('Pa') / (461.5 * (new_data['t']))) * 1e3

        # Set the extra absolute humidity
        # -------------------------------
        new_data['wv'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = new_hum

    elif(add_wv_mix is not None):  
        new_data['mix_rat'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] += add_wv_mix

        # Recalculate relative humidity in these layers
        # ---------------------------------------------
        new_rh = mpcalc.relative_humidity_from_mixing_ratio(\
            np.array(new_data['p'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('mbar'), \
            np.array(new_data['t'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('K'), \
            np.array(new_data['mix_rat'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) / 1e3) * 100.
        new_data['rhum'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = new_rh

        new_data = check_max_rhum(new_data)

        # Convert the added water vapor mixing ratio to added
        # absolute humidity
        # Calculate the new water vapor mixing ratio to new 
        # absolute humidity.
        # ----------------------------------------------------
        vap_press = mpcalc.vapor_pressure(\
            np.array(new_data['p']) * units('mbar'), \
            np.array(new_data['mix_rat']) * units('g/kg'))
        new_hum = (vap_press.to('Pa') / (461.5 * (data['t']))) * 1e3

        # Add the extra absolute humidity to the original values
        # Put the new absolute humidity values in the profile 
        # in the correct locations.
        # ------------------------------------------------------
        new_data['wv'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = new_hum
            #(new_data['z'] <= zmax)] += new_hum

    if(set_rh is not None):  
        # Calculate mixing ratio from relative humidity
        # ---------------------------------------------
        mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(\
            np.array(new_data['p']) * units('mbar'), \
            np.array(new_data['t']) * units('K'), \
            set_rh * units('%')).to('g/kg') 

        # See if any regions of this profile exceed the 4 g/kg maximum
        exceed = np.where(mixing_ratio > 4.0)[0]
        if(len(exceed) > 0):
            # If any are over 4 g/kg, set the maximum to 4 and continue
            mixing_ratio[exceed] = 4.0

        new_data['mix_rat'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = np.array(mixing_ratio.to('g/kg'))

        # Convert the water vapor mixing ratio to added
        # absolute humidity
        # ----------------------------------------------------
        vap_press = mpcalc.vapor_pressure(\
            np.array(new_data['p']) * units('mbar'), mixing_ratio.to('dimensionless'))
        new_hum = (vap_press.to('Pa') / (461.5 * (new_data['t']))) * 1e3

        # Set the extra absolute humidity
        # -------------------------------
        new_data['wv'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = new_hum

        # Recalculate relative humidity in these layers
        # ---------------------------------------------
        new_rh = mpcalc.relative_humidity_from_mixing_ratio(\
            np.array(new_data['p'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('mbar'), \
            np.array(new_data['t'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) * units('K'), \
            np.array(new_data['mix_rat'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)]) / 1e3) * 100.
        new_data['rhum'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = new_rh
    #new_data['wv'][new_data['z'] <= 5] *= mtp
    ##!#for zz, xx, yy in zip(data['z'], data['wv'], new_data['wv']):
    ##!#    print(zz,xx,yy)
   
    header_line = [str(len(new_data['wv'])),' ',' ',' ',' ']

    # Save to outfile
    # ---------------
    outname = out_path + 'atms.dat'
    #outname = 'atms_' + str(int(mtp*100)) + '.dat'
    new_data.to_csv(outname, \
        index = False, sep = ' ', header = header_line, \
        columns = ['z','p','t','wv','o2'])
    print('Saved file', outname)
    
    return new_data

def run_sbdart(satellite, calc_radiance, run = True, vza = None, \
        nzen = None, atms_file = '', cris_profile = None, \
        z_mins = None, z_maxs = None, \
        set_tmp = None, add_tmp = None, add_tmp_sfc = None, \
        set_wv_mix = None, add_wv_mix = None, \
        set_co2_mix = None, add_co2_mix = None, \
        set_ch4_mix = None, add_ch4_mix = None, \
        set_no2_mix = None, add_no2_mix = None, \
        set_co_mix = None,  add_co_mix = None, \
        set_nh3_mix = None, add_nh3_mix = None, \
        set_rh = None, btemp = None):
    if(calc_radiance):
        file_adder = '_rad.dat'
    else:
        file_adder = '.dat'

    if(vza is None):
        # Determine viewing angle values
        # ------------------------------
        if(satellite.strip().split('_')[0] == 'modis'):
            nzen = None
            uzen = 3.15
        elif(satellite.strip().split('_')[0] == 'goes17'):
            nzen = None
            uzen = 49.61378
        else:
            print("INVALID SATELLITE OPTION")
            nzen = 9
            uzen = [0, 60]
    else:
        uzen = vza
    
    # Set up outdict stuff
    # -------------------- 
    names = ['wavelength','filter','down_sol_flux_toa','up_rad_flux_toa',\
             'dir_sol_flux_toa','down_rad_flux_sfc','up_rad_flux_sfc',\
             'dir_sol_flux_sfc']

    if(run): 
        # Prep the filter function for the given satellite
        # ------------------------------------------------
        prep_filter_function(satellite)

        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #
        # SBDART runs
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        
        # Read the atmos file
        # -------------------
        # NOTE: If reading in a model profile for the atmosphere, 
        #       read it in here.
        if(cris_profile is None):
            data = read_atmos_profile(infile = atms_file)
        else:
            print("Converting CrIS profile to SBDART format")
            data = prep_cris_profile(cris_profile)
            btemp = cris_profile['sfc_temp']
        ##!#data = pd.read_csv(infile,skiprows = 0, names = \
        ##!#    ['z','p','t','wv','o2'], delim_whitespace = True)
        
        # Added water vapor mixing ratio, must be converted to 
        # water vapor density.
        # ----------------------------------------------------
        ##!#adders  = np.array([0.0, 2.0, 4.0])
        ##!#z_mins = np.arange(0.,6.)
        ##!#z_maxs = np.arange(5.,11.)
        #z_mins = np.arange(0.,6.)
        #z_maxs = np.arange(2.,8.)
        ##!#set_wvs = np.full(z_mins.shape, 4.)
        #multipliers = np.array([1.0, 2.0, 3.0, 3.5])
        #multipliers = np.arange(0.1,2.1,0.1)
       
 
        outfile_list = []
      
        data_dict = {}
        data_dict['output'] = {}
        data_dict['info'] = {}
        data_dict['orig_atmos'] = data
        if(set_wv_mix is not None): 
            data_dict['info']['met_var_name'] = 'set_wv_mix'
            data_dict['info']['met_var_vals'] = set_wv_mix
        elif(add_wv_mix is not None):
            data_dict['info']['met_var_name'] = 'add_wv_mix'
            data_dict['info']['met_var_vals'] = add_wv_mix
        elif(set_tmp is not None):
            data_dict['info']['met_var_name'] = 'set_tmp'
            data_dict['info']['met_var_vals'] = set_tmp
        elif(add_tmp is not None):
            data_dict['info']['met_var_name'] = 'add_tmp'
            data_dict['info']['met_var_vals'] = add_tmp
        elif(add_tmp_sfc is not None):
            data_dict['info']['met_var_name'] = 'add_tmp_sfc'
            data_dict['info']['met_var_vals'] = add_tmp_sfc
        elif(set_co2_mix is not None):
            data_dict['info']['met_var_name'] = 'set_co2_mix'
            data_dict['info']['met_var_vals'] = set_co2_mix
        elif(add_co2_mix is not None):
            data_dict['info']['met_var_name'] = 'add_co2_mix'
            data_dict['info']['met_var_vals'] = add_co2_mix
        elif(set_ch4_mix is not None):
            data_dict['info']['met_var_name'] = 'set_ch4_mix'
            data_dict['info']['met_var_vals'] = set_ch4_mix
        elif(add_ch4_mix is not None):
            data_dict['info']['met_var_name'] = 'add_ch4_mix'
            data_dict['info']['met_var_vals'] = add_ch4_mix
        elif(set_no2_mix is not None):
            data_dict['info']['met_var_name'] = 'set_no2_mix'
            data_dict['info']['met_var_vals'] = set_no2_mix
        elif(add_no2_mix is not None):
            data_dict['info']['met_var_name'] = 'add_no2_mix'
            data_dict['info']['met_var_vals'] = add_no2_mix
        elif(set_co_mix is not None):
            data_dict['info']['met_var_name'] = 'set_co_mix'
            data_dict['info']['met_var_vals'] = set_co_mix
        elif(add_co_mix is not None):
            data_dict['info']['met_var_name'] = 'add_co_mix'
            data_dict['info']['met_var_vals'] = add_co_mix
        elif(set_nh3_mix is not None):
            data_dict['info']['met_var_name'] = 'set_nh3_mix'
            data_dict['info']['met_var_vals'] = set_nh3_mix
        elif(add_nh3_mix is not None):
            data_dict['info']['met_var_name'] = 'add_nh3_mix'
            data_dict['info']['met_var_vals'] = add_nh3_mix
        elif(set_rh  is not None):
            data_dict['info']['met_var_name'] = 'set_rh'
            data_dict['info']['met_var_vals'] = set_rh

        # POSSIBLE ADDED COMPLEXITY?
        # Add second dimension to z variables for ease of looping
        # -------------------------------------------------------
        if((len(z_mins.shape) == 1) & (len(z_maxs.shape) == 1)):
            # Expand both dimensions
            z_mins = np.expand_dims(z_mins, axis = 0)
            z_maxs = np.expand_dims(z_maxs, axis = 0)
        elif((len(z_maxs.shape) == 1) & (len(z_mins.shape) != 1)):
            # Handle case when the user wants runs with different lower
            # plume heights. 
            # Expand z_maxs
            z_maxs = np.expand_dims(z_maxs, axis = 0)
        elif((len(z_maxs.shape) != 1) & (len(z_mins.shape) == 1)):
            # Handle case when the user wants runs with different higher
            # plume heights. 
            # Expand z_maxs
            z_mins = np.expand_dims(z_mins, axis = 0)

        data_dict['info']['z_mins'] = z_mins
        data_dict['info']['z_maxs'] = z_maxs
        
        if(len(data_dict['info']['met_var_vals'].shape) == 1):
            # User only wants one set of runs. Expand the dimensions
            data_dict['info']['met_var_vals'] = \
                np.expand_dims(data_dict['info']['met_var_vals'], axis = 0)
        else:
            # User wants multiple runs      
            print(' multiple runs with met variables')
            new_zmins = np.full(data_dict['info']['met_var_vals'].shape, np.nan)
            new_zmaxs = np.full(data_dict['info']['met_var_vals'].shape, np.nan)
            for ii in range(data_dict['info']['met_var_vals'].shape[0]):
               new_zmins[ii,:] = z_mins[0, :]
               new_zmaxs[ii,:] = z_maxs[0, :]
            z_mins = new_zmins
            z_maxs = new_zmaxs

 
        print("SHAPES", z_mins.shape, z_maxs.shape, data_dict['info']['met_var_vals'].shape)
 
        ##!## Look at size of met variable and determine the shapes
        ##!## of the arrays
        ##!## -----------------------------------------------------
        ##!#if(z_mins.shape != data_dict['info']['met_var_vals'].shape):
        ##!#    # Shape mismatch. User wants two sets of runs with two
        ##!#    # met variables
        ##!#    new_vals = np.full((len(data_dict['info']['met_var_vals']), len(z_mins)), np.nan)
        ##!#    for ii in range(len(data_dict['info']['met_var_vals'])):
        ##!#        new_vals[ii,:] = data_dict['info']['met_var_vals'][ii]
        ##!#    data_dict['info']['met_var_vals'] = new_vals
        ##!#else:
        ##!#    # Shapes are the same. One set of runs
        ##!#    print("same shapes")

        for kk in range(z_mins.shape[0]):
            r1_key = str(int(kk+1))
            data_dict['output'][r1_key] = {}
            for ii in range(z_mins.shape[1]): 
            #for adr in adders:

                # Prep the INPUT namelist file
                # ----------------------------
                if(set_co2_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        set_co2_mix = set_co2_mix[kk,ii], btemp = btemp)
                elif(add_co2_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        add_co2_mix = add_co2_mix[kk,ii], btemp = btemp)
                elif(set_ch4_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        set_ch4_mix = set_ch4_mix[kk,ii], btemp = btemp)
                elif(add_ch4_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        add_ch4_mix = add_ch4_mix[kk,ii], btemp = btemp)
                elif(set_no2_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        set_no2_mix = set_no2_mix[kk,ii], btemp = btemp)
                elif(add_no2_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        add_no2_mix = add_n2o_mix[kk,ii], btemp = btemp)
                elif(set_co_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        set_co_mix = set_co_mix[kk,ii], btemp = btemp)
                elif(add_co_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        add_co_mix = add_co_mix[kk,ii], btemp = btemp)
                elif(set_nh3_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        set_nh3_mix = set_nh3_mix[kk,ii], btemp = btemp)
                elif(add_nh3_mix is not None):
                    prep_sbdart_input(nzen = nzen, uzen = uzen, \
                        add_nh3_mix = add_nh3_mix[kk,ii], btemp = btemp)
                else:
                    prep_sbdart_input(nzen = nzen, uzen = uzen, btemp = btemp)
                #prep_sbdart_input(uzen = uzen)
        
                # Prep the atmosphere file 
                # NOTE: for now, only allow one meteorological variable to change
                # ---------------------------------------------------------------
                new_data = None
                if(set_wv_mix is not None):
                    new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                        zmax = z_maxs[kk,ii], set_wv_mix = data_dict['info']['met_var_vals'][kk,ii])
                        #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'setwv'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(add_wv_mix is not None):
                    new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                        zmax = z_maxs[kk,ii], add_wv_mix = data_dict['info']['met_var_vals'][kk,ii])
                        #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'addwv'+str(data_dict['info']['met_var_vals'][kk,ii])

                elif(set_tmp is not None):
                    new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                        zmax = z_maxs[kk,ii], set_tmp = data_dict['info']['met_var_vals'][kk,ii])
                        #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'settmp'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(add_tmp is not None):
                    new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                        zmax = z_maxs[kk,ii], add_tmp = data_dict['info']['met_var_vals'][kk,ii])
                        #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'addtmp'+str(data_dict['info']['met_var_vals'][kk,ii])

                elif(set_co2_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], set_co2_mix = data_dict['info']['met_var_vals'][kk,ii])
                        #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'setco2mix'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(add_co2_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], add_co2_mix = data_dict['info']['met_var_vals'][kk,ii])
                    #    #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'addco2mix'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(set_ch4_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], set_ch4_mix = data_dict['info']['met_var_vals'][kk,ii])
                        #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'setch4mix'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(add_ch4_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], add_ch4_mix = data_dict['info']['met_var_vals'][kk,ii])
                    #    #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'addch4mix'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(set_no2_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], set_no2_mix = data_dict['info']['met_var_vals'][kk,ii])
                        #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'setno2mix'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(add_no2_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], add_no2_mix = data_dict['info']['met_var_vals'][kk,ii])
                    #    #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'addno2mix'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(set_co_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], set_co_mix = data_dict['info']['met_var_vals'][kk,ii])
                        #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'setcomix'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(add_co_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], add_co_mix = data_dict['info']['met_var_vals'][kk,ii])
                    #    #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'addcomix'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(set_nh3_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], set_nh3_mix = data_dict['info']['met_var_vals'][kk,ii])
                        #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'setnh3mix'+str(data_dict['info']['met_var_vals'][kk,ii])
                elif(add_nh3_mix is not None):
                    new_data = data
                    #new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                    #    zmax = z_maxs[kk,ii], add_nh3_mix = data_dict['info']['met_var_vals'][kk,ii])
                    #    #set_wv_mix = adr)
                    run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                'addnh3mix'+str(data_dict['info']['met_var_vals'][kk,ii])

                # Allows the user to modify the surface temperature in 
                # addition to other variables
                # ----------------------------------------------------
                if(add_tmp_sfc is not None):
                    if(new_data is None):
                        new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                            zmax = z_maxs[kk,ii], add_tmp_sfc = data_dict['info']['met_var_vals'][kk,ii])
                            #set_wv_mix = adr)
                        run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                    'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                    'addtmpssfc'+str(data_dict['info']['met_var_vals'][kk,ii])
                    else:
                        new_data = prep_atmos_file(new_data, zmin = z_mins[kk,ii], \
                            zmax = z_maxs[kk,ii], add_tmp_sfc = add_tmp_sfc)
                            #set_wv_mix = adr)
                        run_adder = run_adder + '_addtmpsfc' + str(int(add_tmp_sfc))

                elif(set_rh is not None):
                    if(np.isnan(data_dict['info']['met_var_vals'][kk,ii])):
                        local_rh = None
                        run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                    'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                    'setrhcontrol'
                    else:
                        run_adder = 'zmin'+str(int(z_mins[kk,ii])) + '_' + \
                                    'zmax'+str(int(z_maxs[kk,ii])) + '_' + \
                                    'setrh'+str(data_dict['info']['met_var_vals'][kk,ii])
                        local_rh = data_dict['info']['met_var_vals'][kk,ii]
                    new_data = prep_atmos_file(data, zmin = z_mins[kk,ii], \
                        zmax = z_maxs[kk,ii], set_rh = local_rh)
                        #set_wv_mix = adr)

                print(run_adder) 
                #prep_atmos_file(data, zmin = 0., zmax = 5., \
                #    add_wv_mix = adr)
 
                # Run SBDART with the new atmosphere profile
                # ------------------------------------------
                outfile = 'sbout_' + satellite + '_' + run_adder + '_'+\
                    file_adder
                cmnd = sb_path + 'bin/sbdart > ' + outfile
                print(cmnd)
                os.system(cmnd)
            
                # Modify the output file
                # ----------------------
                if(not calc_radiance):
                    with open(outfile,'r') as fin:
                        indata = fin.readlines()
                    with open(outfile,'w') as fout:
                        fout.writelines(indata[3:])
       
                # Read the output file contents into the structure
                # ------------------------------------------------
                #wv_pct = int(ofile.split('/')[-1].split('_')[1])
                d_key = run_adder
                r2_key = str(int(ii+1))
                data_dict['output'][r1_key][r2_key] = {}
                data_dict['output'][r1_key][r2_key]['data'] = {}
                data_dict['output'][r1_key][r2_key]['mod_str'] = d_key
                data_dict['output'][r1_key][r2_key]['atmos'] = new_data
                with open(outfile,'r') as fin: 
                    flines = fin.readlines()

                    # Extract the number of reports
                    # -----------------------------
                    num_rep = int(flines[2].strip().split()[0])

                    # Pull out the number of VZAs in the file
                    # ---------------------------------------
                    num_vzen = int(flines[4].strip().split()[1])

                    # Extract the VZAs
                    # ----------------
                    vza = np.array([float(tval) for tval in flines[6].strip().split()])
               
                    # Extract the wavelength and filter values
                    # ----------------------------------------
                    datalines = flines[3:][::4+num_vzen]
                    wavel  = np.array([float(dstr.strip().split()[0]) for dstr in datalines])
                    fvals  = np.array([float(dstr.strip().split()[1]) for dstr in datalines])
                    topdwn = np.array([float(dstr.strip().split()[2]) for dstr in datalines])
                    topup  = np.array([float(dstr.strip().split()[3]) for dstr in datalines])
                    topdir = np.array([float(dstr.strip().split()[4]) for dstr in datalines])
                    botdwn = np.array([float(dstr.strip().split()[5]) for dstr in datalines])
                    botup  = np.array([float(dstr.strip().split()[6]) for dstr in datalines])
                    botdir = np.array([float(dstr.strip().split()[7]) for dstr in datalines])

                    # Declare an array to hold all of the radiances
                    # ---------------------------------------------
                    total_radiances = np.full((num_rep, num_vzen), np.nan)

                    #vzen_lines = flines[6:][:
 
                    # Extract the radiances
                    # --------------------- 
                    for ii in range(num_vzen):
                        total_radiances[:,ii] = np.array([float(tstr.strip()) for tstr in flines[(7+ii):][::4+num_vzen]])
                    #radlines = flines[7:][::5]
                    #radiances = np.array([float(tstr.strip()) for tstr in radlines])
    
    
                    # Insert all the data into the dictionary
                    # ---------------------------------------
                    data_dict['output'][r1_key][r2_key]['vza'] = vza
                    data_dict['output'][r1_key][r2_key]['data']['wavelength'] = wavel
                    data_dict['output'][r1_key][r2_key]['data']['filter']     = fvals
                    data_dict['output'][r1_key][r2_key]['data']['radiance']   = total_radiances
                    data_dict['output'][r1_key][r2_key]['data']['down_sol_flux_toa'] = topdwn
                    data_dict['output'][r1_key][r2_key]['data']['up_rad_flux_toa']   = topup
                    data_dict['output'][r1_key][r2_key]['data']['dir_sol_flux_toa']  = topdwn
                    data_dict['output'][r1_key][r2_key]['data']['down_rad_flux_sfc'] = botdwn
                    data_dict['output'][r1_key][r2_key]['data']['up_rad_flux_sfc']   = botup 
                    data_dict['output'][r1_key][r2_key]['data']['dir_sol_flux_sfc']  = botdir
    
                # Calculated weighted average wavelength
                # --------------------------------------
                data_dict['output'][r1_key][r2_key]['avg_wavel'] = \
                    np.sum(data_dict['output'][r1_key][r2_key]['data']['wavelength'] * \
                    data_dict['output'][r1_key][r2_key]['data']['filter']) / \
                    np.sum(data_dict['output'][r1_key][r2_key]['data']['filter'])
    
                # Calculated weighted average radiance  
                # ------------------------------------
                data_dict['output'][r1_key][r2_key]['avg_rads'] = \
                    np.array([np.sum(data_dict['output'][r1_key][r2_key]['data']['radiance'][:,jj] * \
                    data_dict['output'][r1_key][r2_key]['data']['filter']) / \
                    np.sum(data_dict['output'][r1_key][r2_key]['data']['filter']) \
                    for jj in range(len(data_dict['output'][r1_key][r2_key]['vza']))])

                #data_dict[d_key]['avg_rad'] = \
                #    np.sum(data_dict[d_key]['data']['radiance'] * data_dict[d_key]['data']['filter']) / \
                #    np.sum(data_dict[d_key]['data']['filter'])

                # Calculate temperatures for all wavelengths
                # ------------------------------------------
                data_dict['output'][r1_key][r2_key]['data']['bght_tmp_all'] = \
                    np.array([inverse_planck(data_dict['output'][r1_key][r2_key]['data']['radiance'][:,jj], \
                    data_dict['output'][r1_key][r2_key]['data']['wavelength']) \
                    for jj in range(len(data_dict['output'][r1_key][r2_key]['vza']))]).transpose()

                #data_dict[r1_key][r2_key]['data']['bght_tmp_all'] = \
                #    inverse_planck(data_dict[r1_key][r2_key]['data']['radiance'], \
                #    data_dict[r1_key][r2_key]['data']['wavelength'])

                # Calculate temperature for filter averaged values
                # ------------------------------------------------
                data_dict['output'][r1_key][r2_key]['bght_tmp'] = \
                    inverse_planck(data_dict['output'][r1_key][r2_key]['avg_rads'], \
                    data_dict['output'][r1_key][r2_key]['avg_wavel'])
 
                outfile_list.append(outfile)

        # Insert the satellite name into the dictionary
        data_dict['info']['run_name'] = satellite
        data_dict['info']['satellite'] = satellite.split('_')[0]
        data_dict['info']['channel'] = satellite.split('_')[1][2:]

    else:
        search_str = sb_path + 'src/sbout_*' + satellite + \
            file_adder
        print(search_str)
        outfile_list = glob.glob(search_str)

        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #
        # Data reading
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        
        # Read in the SBDART output 
        # -------------------------
        #outfile_list = glob.glob(sb_path + 'src/sbout_*.dat')
        
        names = ['wavelength','filter','down_sol_flux_toa','up_rad_flux_toa',\
                 'dir_sol_flux_toa','down_rad_flux_sfc','up_rad_flux_sfc',\
                 'dir_sol_flux_sfc']
   
        print(outfile_list)
 
        data_dict = {}
        data_dict['adders'] = adders
        #data_dict['multipliers'] = multipliers
        for ofile in outfile_list:
            if(not calc_radiance):
                wv_pct = int(ofile.split('.')[0].split('_')[-1])
                data_dict[wv_pct] = {}
                data_dict[wv_pct]['data'] = pd.read_csv(ofile,skiprows = 0, names = names,\
                        delim_whitespace = True, header = None)
                data_dict[wv_pct]['data'] = data_dict[wv_pct]['data'].set_index('wavelength')
        
                # Calculated weighted average wavelength
                # --------------------------------------
                data_dict[wv_pct]['avg_wavel'] = \
                    np.sum(data_dict[wv_pct]['data'].index * data_dict[wv_pct]['data']['filter']) / \
                    np.sum(data_dict[wv_pct]['data']['filter'])
            else:
                wv_pct = int(ofile.split('/')[-1].split('_')[1])
                data_dict[wv_pct] = {}
                data_dict[wv_pct]['data'] = {}
                with open(ofile,'r') as fin: 
                    flines = fin.readlines()

                    # Extract the number of reports
                    # -----------------------------
                    num_rep = int(flines[2].strip().split()[0])

                    # Pull out the number of VZAs in the file
                    # ---------------------------------------
                    num_vzen = int(flines[4].strip().split()[1])

                    # Extract the VZAs
                    # ----------------
                    vza = np.array([float(tval) for tval in flines[6].strip().split()])
               
                    # Extract the wavelength and filter values
                    # ----------------------------------------
                    datalines = flines[3:][::4+num_vzen]
                    wavel  = np.array([float(dstr.strip().split()[0]) for dstr in datalines])
                    fvals  = np.array([float(dstr.strip().split()[1]) for dstr in datalines])
                    topdwn = np.array([float(dstr.strip().split()[2]) for dstr in datalines])
                    topup  = np.array([float(dstr.strip().split()[3]) for dstr in datalines])
                    topdir = np.array([float(dstr.strip().split()[4]) for dstr in datalines])
                    botdwn = np.array([float(dstr.strip().split()[5]) for dstr in datalines])
                    botup  = np.array([float(dstr.strip().split()[6]) for dstr in datalines])
                    botdir = np.array([float(dstr.strip().split()[7]) for dstr in datalines])

                    # Declare an array to hold all of the radiances
                    # ---------------------------------------------
                    total_radiances = np.full((num_rep, num_vzen), np.nan)

                    #vzen_lines = flines[6:][:
 
                    # Extract the radiances
                    # --------------------- 
                    for ii in range(num_vzen):
                        total_radiances[:,ii] = np.array([float(tstr.strip()) for tstr in flines[(7+ii):][::4+num_vzen]])
                    #radlines = flines[7:][::5]
                    #radiances = np.array([float(tstr.strip()) for tstr in radlines])
        
        
                    # Insert all the data into the dictionary
                    # ---------------------------------------
                    data_dict[wv_pct]['vza'] = vza
                    data_dict[wv_pct]['data']['wavelength'] = wavel
                    data_dict[wv_pct]['data']['filter']     = fvals
                    data_dict[wv_pct]['data']['radiance']   = total_radiances
                    data_dict[wv_pct]['data']['down_sol_flux_toa'] = topdwn
                    data_dict[wv_pct]['data']['up_rad_flux_toa']   = topup
                    data_dict[wv_pct]['data']['dir_sol_flux_toa']  = topdwn
                    data_dict[wv_pct]['data']['down_rad_flux_sfc'] = botdwn
                    data_dict[wv_pct]['data']['up_rad_flux_sfc']   = botup 
                    data_dict[wv_pct]['data']['dir_sol_flux_sfc']  = botdir
        
                # Calculated weighted average wavelength
                # --------------------------------------
                data_dict[wv_pct]['avg_wavel'] = \
                    np.sum(data_dict[wv_pct]['data']['wavelength'] * data_dict[wv_pct]['data']['filter']) / \
                    np.sum(data_dict[wv_pct]['data']['filter'])
        
                # Calculated weighted average radiance  
                # ------------------------------------
                data_dict[wv_pct]['avg_rads'] = \
                    np.array([np.sum(data_dict[wv_pct]['data']['radiance'][:,jj] * \
                    data_dict[wv_pct]['data']['filter']) / np.sum(data_dict[wv_pct]['data']['filter']) \
                    for jj in range(len(data_dict[wv_pct]['vza']))])

                #data_dict[wv_pct]['avg_rad'] = \
                #    np.sum(data_dict[wv_pct]['data']['radiance'] * data_dict[wv_pct]['data']['filter']) / \
                #    np.sum(data_dict[wv_pct]['data']['filter'])

                # Calculate temperatures for all wavelengths
                # ------------------------------------------
                data_dict[wv_pct]['data']['bght_tmp_all'] = \
                    np.array([inverse_planck(data_dict[wv_pct]['data']['radiance'][:,jj], data_dict[wv_pct]['data']['wavelength']) \
                    for jj in range(len(data_dict[wv_pct]['vza']))]).transpose()

                #data_dict[wv_pct]['data']['bght_tmp_all'] = \
                #    inverse_planck(data_dict[wv_pct]['data']['radiance'], \
                #    data_dict[wv_pct]['data']['wavelength'])

                # Calculate temperature for filter averaged values
                # ------------------------------------------------
                data_dict[wv_pct]['bght_tmp'] = \
                    inverse_planck(data_dict[wv_pct]['avg_rads'], \
                    data_dict[wv_pct]['avg_wavel'])

        # Insert the satellite name into the dictionary
        data_dict['run_name'] = satellite
        data_dict['satellite'] = satellite.split('_')[0]
        data_dict['channel'] = satellite.split('_')[1][2:]

    return data_dict

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Calculations
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def plot_bright(modis_data1, modis_data2, vza_idx = 0, save = False):

    # Pull out the brightness temperatures for each water vapor amount
    # ----------------------------------------------------------------
    m31_bghts = np.array([modis_data1[key]['bght_tmp'][vza_idx] for key in modis_data1.keys()])
    m32_bghts = np.array([modis_data2[key]['bght_tmp'][vza_idx] for key in modis_data2.keys()])
    diffs = m31_bghts - m32_bghts
    wvs = modis31.keys()
    
    plt.close('all')
    fig1 = plt.figure(figsize = (10,5))
    ax0 = fig1.add_subplot(1,2,1)
    ax1 = fig1.add_subplot(1,2,2)
    ax0.plot(wvs, m31_bghts, label = modis_data1['satellite'].upper() + ' Ch ' + modis_data1['channel'])
    ax0.plot(wvs, m32_bghts, label = modis_data2['satellite'].upper() + ' Ch ' + modis_data2['channel'])
    ax1.plot(wvs, diffs, label = 'MODIS Ch 32')
    ax0.set_xlabel('Percent of original WV in lowest 5 km')
    ax0.set_ylabel('Brightness temperature [K]')
    ax1.set_xlabel('Percent of original WV in lowest 5 km')
    ax1.set_ylabel('Brightness temperature Difference [K]')

    if(save):
        outname = 'modis_ch31vch32_btmp.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show()

add_label_dict = {
    'set_wv_mix': 'w = {0} g/kg', \
    'add_wv_mix': 'w = w$_{{o}}$ + {0} g/kg',\
    'set_tmp':    'T = {0} K',\
    'add_tmp':    'T = T$_{{o}}$ + {0}  K',
    'add_tmp_sfc':    'T$_{{sfc}}$ = T$_{{sfc_o}}$ + {0}  K',
    'set_co2_mix': 'CO$_2$ = {0} PPM', \
    'add_co2_mix': 'CO$_2$ = CO$_2$ + {0} PPM',\
    'set_ch4_mix': 'CH$_4$ = {0} PPM', \
    'add_ch4_mix': 'CH$_4$ = CH$_4$ + {0} PPM',\
    'set_no2_mix': 'NO$_2$ = {0} PPM', \
    'add_no2_mix': 'NO$_2$ = NO$_2$ + {0} PPM',\
    'set_co_mix': 'CO = {0} PPM', \
    'add_co_mix': 'CO = CO + {0} PPM',\
    'set_nh3_mix': 'NH$_3$ = {0} PPM',\
    'add_nh3_mix': 'NH$_3$ = NH$_3$ + {0} PPM',\
    'set_rh':     'RH = {0}%',\
}

def plot_bright_vza(sbout, pax = None, ptitle = '', plabelsize = 12, \
        legendsize = 9, titlesize = 12, save = False):
    in_pax = True
    if(pax is None):
        in_pax = False 
        plt.close('all')
        fig1 = plt.figure(figsize = (8,4))
        pax = fig1.add_subplot(1,1,1)

    ii = 0 
    for r1_key in sbout['output'].keys():
        for r2_key in sbout['output'][r1_key].keys():
            #for wv in sbout['adders']:
            #for wv in sbout['multipliers']:
            #int_wv = int(wv*100.)
            wavel = np.round(sbout['output'][r1_key][r2_key]['avg_wavel'],2)

            # Plot the correct wv amounts
            pax.plot(sbout['output'][r1_key][r2_key]['vza'],  \
                sbout['output'][r1_key][r2_key]['bght_tmp'], label = \
                    add_label_dict[sbout['info']['met_var_name']].format(\
                    sbout['info']['met_var_vals'][ii,0])) 
        ii += 1

    pax.set_xlabel('Viewing Zenith Angle [$^{o}$]', fontsize = plabelsize, \
        weight = 'bold')
    pax.set_ylabel('Brightness temperature [K]',    fontsize = plabelsize, \
        weight = 'bold')
    if(ptitle == ''):
        #ptitle = 'SBDART-simulated ' + sbout['info']['satellite'].upper() + \
        ptitle = sbout['info']['satellite'].upper() + \
        ' Channel ' + sbout['info']['channel'] + ' (' + str(wavel) + ' μm)'
    pax.set_title(ptitle, fontsize = titlesize)
    pax.legend(fontsize = legendsize, loc = 'lower left')
    if(not in_pax):
        if(save):
            outname = 'sbdart_'+sbout['run_name']+'_bright_vza.png'
            fig1.savefig(outname, dpi = 300)
            print("Saved image",outname)
        else:
            plt.show() 

def plot_bright_sfc_tmp(sbout, pax = None, relative = False, multi_sat = False, save = False):
    # Set up axis
    # -----------
    in_ax = True 
    if(pax is None):
        in_ax = False
        plt.close('all')   
        fig = plt.figure()
        pax = fig.add_subplot(1,1,1)

    if((sbout['info']['z_mins'].shape != sbout['info']['z_mins'].squeeze().shape) & \
       (sbout['info']['z_maxs'].shape != sbout['info']['z_maxs'].squeeze().shape)):
        local_z_mins = sbout['info']['z_mins'].squeeze()
        local_z_maxs = sbout['info']['z_maxs'].squeeze()
        # Check if user wants multiple meteorological values
        if(sbout['info']['met_var_vals'].shape != sbout['info']['met_var_vals'].squeeze().shape):
            # Only 1 critical dimension to each of the arrays. No dual runs here
            print("Single run")
        else:
            # Multiple critical dimension on met variable. Multiple runs
            print("Multiple met runs")
    else: 
        if(sbout['info']['z_mins'].shape != sbout['info']['z_mins'].squeeze().shape):
            print("Multiple lower plume heights")
        else:
            print("Multiple upper plume heights")

    # Extract the single brightness temperature for each run
    # ------------------------------------------------------
    bght_tmps = np.array([[sbout['output'][r1_key][r2_key]['bght_tmp'] for r2_key \
        in sbout['output'][r1_key].keys()] for r1_key in sbout['output'].keys()])

    if(relative):
        bght_tmps = bght_tmps - bght_tmps[0]
    
    print(sbout['info']['run_name'], bght_tmps.squeeze())

    # X axis is the low-level cooling. Y-axis is the brightness temp
    # --------------------------------------------------------------
    x_vals = sbout['info']['met_var_vals']
    if(np.mean(x_vals) < 0):
        label_str = 'Surface cooling [K]'
        x_vals = np.abs(x_vals)
    else:
        label_str = 'Surface warming [K]'

    if(multi_sat):
        pax.plot(x_vals, bght_tmps.squeeze(), label = title_dict[sbout['info']['run_name']])
        axis_fontsize = 10
    else:
        pax.plot(x_vals, bght_tmps.squeeze())
        pax.set_title(title_dict[sbout['info']['run_name']], fontsize = 10)
        if((pax.get_ylim()[1] - pax.get_ylim()[0]) < 0.5):
            mean_val = np.mean(bght_tmps.squeeze())
            pax.set_ylim([mean_val - 0.5, mean_val + 0.5])
        axis_fontsize = 10
    pax.set_xlabel(label_str, fontsize = axis_fontsize)
    if(relative):
        ylabel_str = 'Brightness temperature difference [K]'
    else:
        ylabel_str = 'Brightness temperature [K]'
    pax.set_ylabel(ylabel_str, fontsize = axis_fontsize)
    

    if(not in_ax):
        if(save):
            print("SAVE NAME NEEDED")
        else:
            plt.show()

def plot_bright_plume_height(sbout, pax = None, save = False):
    # Set up axis
    # -----------
    in_ax = True 
    if(pax is None):
        in_ax = False
        plt.close('all')   
        fig = plt.figure()
        pax = fig.add_subplot(1,1,1)

    if((sbout['info']['z_mins'].shape != sbout['info']['z_mins'].squeeze().shape) & \
       (sbout['info']['z_maxs'].shape != sbout['info']['z_maxs'].squeeze().shape)):
        local_z_mins = sbout['info']['z_mins'].squeeze()
        local_z_maxs = sbout['info']['z_maxs'].squeeze()
        # Check if user wants multiple meteorological values
        if(sbout['info']['met_var_vals'].shape != sbout['info']['met_var_vals'].squeeze().shape):
            # Only 1 critical dimension to each of the arrays. No dual runs here
            print("Single run")
        else:
            # Multiple critical dimension on met variable. Multiple runs
            print("Multiple met runs")
    else: 
        if(sbout['info']['z_mins'].shape != sbout['info']['z_mins'].squeeze().shape):
            print("Multiple lower plume heights")
        else:
            print("Multiple upper plume heights")

    # Extract the single brightness temperature for each run
    # ------------------------------------------------------
    bght_tmps = np.array([[sbout['output'][r1_key][r2_key]['bght_tmp'] for r2_key \
        in sbout['output'][r1_key].keys()] for r1_key in sbout['output'].keys()])

    # Look at the mod string in the sbout structure and figure out what's
    # edited. If the first two plume min heights are the same and the
    # first two plume max heights are different, then it's a plume 
    # thickness run. If both the first two plume min heights and max
    # heights are different, then it's a plume height run.
    # -------------------------------------------------------------------
    # If these two are the same, check the two max heights
    ptype = 'met'
    if(local_z_mins[0] == local_z_mins[1]):
        # If these two are the same, ...
        if(local_z_maxs[0] == local_z_maxs[1]):
            # None of the heights are changed, so must be other changes
            print("No height changes. Met param change")
        else:
            # The max heights vary, so it's a plume thickness run
            print("plume thickness run: top varying") 
            ptype = 'thicktop'
    else:
        # If these two are the same, ...
        if(local_z_maxs[0] == local_z_maxs[1]):
            # The max heights vary, so it's a plume thickness run
            print("plume thickness run: bottom varying") 
            ptype = 'thickbottom'
        else:
            # It's a plume height change
            print("Plume height run")
            ptype = 'height'
        
    print('local_z_mins = ', local_z_mins)
    print('local_z_maxs = ', local_z_maxs)

    if(ptype != 'met'):
        if(ptype == 'height'):
            # Calculate each plume height
            x_vals = (local_z_maxs + local_z_mins) / 2.
            pax.set_xlabel('Plume height [km]')
        else:
            # Calculate each plume thickness
            x_vals = local_z_maxs - local_z_mins
            pax.set_xlabel('Plume thickness [km]')
       
        for ii in range(bght_tmps.shape[0]): 
            if(np.isnan(sbout['info']['met_var_vals'][ii,0])):
                label_str = 'RH = original'
            else:
                label_str = \
                    add_label_dict[sbout['info']['met_var_name']].format(\
                    int(sbout['info']['met_var_vals'][ii,0])) 
            pax.plot(x_vals, bght_tmps[ii].squeeze(), label = label_str)
        pax.legend(fontsize = 8)
        pax.set_ylabel('Brightness temperature [K]')

    if(not in_ax):
        if(save):
            print("SAVE NAME NEEDED")
        else:
            plt.show()

title_dict = {
    'modis_ch31':  'Aqua MODIS Band 31 (11.0 μm)', \
    'goes17_ch08': 'GOES 17 Band 8 (6.2 μm)', \
    'goes17_ch09': 'GOES 17 Band 9 (6.94 μm)', \
    'goes17_ch10': 'GOES 17 Band 10 (7.34 μm)', \
    'goes17_ch13': 'GOES 17 Band 13 (10.35 μm)', \
}

def plot_SBDART_dual_height_thickness(data1, data2, ax1 = None, ax2 = None, save = False):

    in_ax = True
    if((ax1 is None) and (ax2 is None)):
        in_ax = False
        plt.close('all')
        fig = plt.figure(figsize = (8, 4))
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)
    
    plot_bright_plume_height(data1, pax = ax1)
    plot_bright_plume_height(data2, pax = ax2)
    if(in_ax):
        ax1.set_title(title_dict[data1['info']['run_name']])
        ax2.set_title(title_dict[data2['info']['run_name']])
    else: 
        plt.suptitle(title_dict[data1['info']['run_name']])
        fig.tight_layout()
        if(save):
            outname = 'sbdart_dual_height_thickness_' + \
                data1['info']['run_name'] + '_' + \
                data1['info']['met_var_name'] + '.png'
            fig.savefig(outname, dpi = 300)
            print("Saved image", outname)
        else:
            plt.show()

def plot_SBDART_atmos(sbout, r1_key = '1', ptype = 'met', close_plots = True, save = False):

    # Use the met_var_name to figure out how to make the plot
    if((sbout['info']['met_var_name'][4:] == 'wv_mix') | ((sbout['info']['met_var_name'][4:] == 'rh'))):
        ncols = 3
        plot_vars = ['wv','mix_rat','rhum']
        plot_names = ['Absolute Humidity [g/m$^{3}$]', \
            'Mixing Ratio [g/kg]', 'Relative Humidity [%]']
        tmp_fig = False
    else:
        ncols = 2
        plot_vars = ['t','mix_rat']
        plot_names = ['Temperature [K]', 'Mixing Ratio [g/kg]']
        tmp_fig = True 

    # Set up axis
    # -----------
    if(close_plots):
        plt.close('all')   
    figsize = (2 * ncols, 7)
    fig = plt.figure(figsize = figsize)
    r_num = len(sbout['output'].keys())
    axs = fig.subplots(nrows = r_num, ncols = ncols)

    xfontsize = 10
    for ii, r1_key in enumerate(sbout['output'].keys()):
        print(r1_key) 
        for kk in range(len(plot_names)):
            axs[ii,kk].plot(sbout['orig_atmos'][plot_vars[kk]], sbout['orig_atmos']['z'], \
                label = 'control') 
            for jj, r2_key in enumerate(sbout['output'][r1_key].keys()):
                axs[ii,kk].plot(sbout['output'][r1_key][r2_key]['atmos'][plot_vars[kk]], \
                    sbout['output'][r1_key][r2_key]['atmos']['z'], label = \
                    add_label_dict[sbout['info']['met_var_name']].format(\
                    sbout['info']['met_var_vals'][ii,0])) 
            axs[ii,kk].set_ylim(-1, 16)
            if(kk == 0):
                axs[ii,kk].set_ylabel('Height [km]')
            if(ii == r_num - 1):
                axs[ii,kk].set_xlabel(plot_names[kk])

    fig.tight_layout()

    if(save):
        outname = 'sbdart_atmos_profiles_' + sbout['info']['met_var_name'] + '_' + ptype + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def process_SBDART_multi_plume_height_thick(atms_file = '', plot_atmos = False, save = False):

    satellites = ['goes17_ch08','goes17_ch09','goes17_ch10','goes17_ch13', 'modis_ch31']

    # Set up the overall figures. Figure 1 is for the mixing ratio changes. 
    # Figure 2 is for the RH changes
    # ---------------------------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (6, 11))
    axs1 = fig1.subplots(nrows = len(satellites), ncols = 2)
    fig2 = plt.figure(figsize = (6, 11))
    axs2 = fig2.subplots(nrows = len(satellites), ncols = 2)
    
    for ii, satellite in enumerate(satellites):
    
        # ADD MIX RUNS
        
        # Run num 1
        z_maxs = np.arange(1.0, 6., 1.00)
        z_mins = np.full(z_maxs.shape, 0.)
        add_wv_mix = np.full((3, len(z_maxs)), np.nan)
        add_wv_mix[0,:] = 0
        add_wv_mix[1,:] = 2
        add_wv_mix[2,:] = 4
        
        data1 = run_sbdart(satellite, calc_radiance = True, atms_file = atms_file, \
            z_mins = z_mins, z_maxs = z_maxs, add_wv_mix = add_wv_mix)
        
        # Run num 2
        z_maxs = np.arange(2, 6.)
        z_mins = np.arange(0, 4.)
        add_wv_mix = np.full((3, len(z_maxs)), np.nan)
        add_wv_mix[0,:] = 0
        add_wv_mix[1,:] = 2
        add_wv_mix[2,:] = 4
        
        data2 = run_sbdart(satellite, calc_radiance = True, atms_file = atms_file, \
            z_mins = z_mins, z_maxs = z_maxs, add_wv_mix = add_wv_mix)
        
        
        plot_SBDART_dual_height_thickness(data1, data2, \
            ax1 = axs1[ii, 0], ax2 = axs1[ii, 1], save = True)
        if(plot_atmos):
            plot_SBDART_atmos(data1, ptype = 'thicktop', close_plots = False, save = True)
            plot_SBDART_atmos(data2, ptype = 'height', close_plots = False, save = True)
        
        # RHUM runs
        
        # Run num 1
        #z_maxs = np.arange(1.0, 11., 1.00)
        z_maxs = np.arange(1.0, 6., 1.00)
        z_mins = np.full(z_maxs.shape, 0.)
        set_rh = np.full((3, len(z_maxs)), np.nan)
        set_rh[0,:] = None
        set_rh[1,:] = 40
        set_rh[2,:] = 80
        
        data3 = run_sbdart(satellite, calc_radiance = True, atms_file = atms_file, \
            z_mins = z_mins, z_maxs = z_maxs, set_rh = set_rh)
        
        # Run num 2
        #z_maxs = np.arange(2, 11.)
        #z_mins = np.arange(0, 9.)
        z_maxs = np.arange(2, 6.)
        z_mins = np.arange(0, 4.)
        ##!#set_rh = np.full(len(z_maxs), 4.)
        set_rh = np.full((3, len(z_maxs)), np.nan)
        set_rh[0,:] = None
        set_rh[1,:] = 40
        set_rh[2,:] = 80
        
        data4 = run_sbdart(satellite, calc_radiance = True, atms_file = atms_file, \
            z_mins = z_mins, z_maxs = z_maxs, set_rh = set_rh)
        
        
        plot_SBDART_dual_height_thickness(data3, data4, \
            ax1 = axs2[ii, 0], ax2 = axs2[ii, 1], save = True)
        if(plot_atmos):
            plot_SBDART_atmos(data3, ptype = 'thicktop', close_plots = False, save = True)
            plot_SBDART_atmos(data4, ptype = 'height', close_plots = False, save = True)

    fig1.tight_layout()
    fig2.tight_layout()
    if(save):
        outname1 = 'sbdart_multi_height_thick_mixrat.png'
        outname2 = 'sbdart_multi_height_thick_RH.png'
        fig1.savefig(outname1, dpi = 300)
        fig2.savefig(outname2, dpi = 300)
        print("Saved image", outname1)
        print("Saved image", outname2)
    else:
        plt.show()

def process_SBDART_multi_lower_tmps(atms_file = '', save = False, multi_sat = True, relative = True):
    satellites = ['goes17_ch08','goes17_ch09','goes17_ch10','goes17_ch13', 'modis_ch31']
   
    plt.close('all')
    if(multi_sat):
        figsize = (5,5)
    else:
        figsize = (9,6)
    fig = plt.figure(figsize = figsize)
    if(multi_sat):
        pax = fig.add_subplot(1,1,1)
    #axs = fig.subplots(nrows = 2, ncols = 3)
 
    for jj, satellite in enumerate(satellites):
    
        # These max vals cool the surface in the lowest 1 km 
        z_maxs = np.arange(2.0, 3., 1.00)
        z_mins = np.full(z_maxs.shape, 0.)
        add_tmp = np.full((11, len(z_maxs)), np.nan)
        for ii in range(11):
            add_tmp[ii,:] = -ii
        
        data1 = run_sbdart(satellite, calc_radiance = True, atms_file = atms_file, \
            z_mins = z_mins, z_maxs = z_maxs, add_tmp = add_tmp)
  
        if(not multi_sat): 
            pax = fig.add_subplot(2,3,jj+1)    
        plot_bright_sfc_tmp(data1, pax = pax, save = False, \
            multi_sat = multi_sat, relative = relative)

    if(multi_sat):
        pax.legend(fontsize = 9)

    fig.tight_layout()
    if(save):
        rel_add = ''
        mult_add = ''
        if(relative):
            rel_add = '_relative'
        if(multi_sat):
            mult_add = '_multisat' 
        outname = 'sbdart_multi_sfc_cool'+rel_add+mult_add+'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This runs the SBDART simulations for several CrIS profiles
def process_SBDART_multi_sat_vza_CrIS(date_str, row_str = 'ml', \
        save = False):

    if('/home/bsorenson/Research/CrIS' not in sys.path):
        sys.path.append('/home/bsorenson/Research/CrIS')
    from CrISLib import cris_loc_dict, readCrIS_retrieval_profile

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M%S')

    z_maxs = np.array([7.])
    z_mins = np.array([0.])
    add_wv_mix = np.full((3, len(z_maxs)), np.nan)
    add_wv_mix[0,:] = 0. 
    add_wv_mix[1,:] = 2. 
    add_wv_mix[2,:] = 4. 

    # Get the lat/lon pairs for the CrIS data
    # ---------------------------------------    
    smoke_lat = cris_loc_dict[row_str]['smoke_lat']
    smoke_lon = cris_loc_dict[row_str]['smoke_lon']
    clear_lat1 = cris_loc_dict[row_str]['clear_lat1']
    clear_lon1 = cris_loc_dict[row_str]['clear_lon1']
    clear_lat2 = cris_loc_dict[row_str]['clear_lat2']
    clear_lon2 = cris_loc_dict[row_str]['clear_lon2']

    # Read in CrIS profiles for the clear and smoky points
    # ----------------------------------------------------
    CrIS_smoke  =  readCrIS_retrieval_profile(date_str, smoke_lat, smoke_lon)
    CrIS_clear1 =  readCrIS_retrieval_profile(date_str, clear_lat1, clear_lon1)
    CrIS_clear2 =  readCrIS_retrieval_profile(date_str, clear_lat2, clear_lon2)

    # Modify the clear 1 profile to contain the smoky moisture amounts
    # ----------------------------------------------------------------
    CrIS_clear1_wv = deepcopy(CrIS_clear1)
    bottom_idx = np.where(np.isnan(CrIS_smoke['temp']))[0][0]
    CrIS_clear1_wv['wv'][:bottom_idx] = CrIS_smoke['wv'][:bottom_idx] 

    # Modify the modified profile to now include the smoky temperatures
    # -----------------------------------------------------------------
    CrIS_clear1_both = deepcopy(CrIS_clear1_wv)
    CrIS_clear1_both['temp'][:bottom_idx] = CrIS_smoke['temp'][:bottom_idx] 
 

    channels = ['goes17_ch08','goes17_ch09','goes17_ch10','modis_ch31']

    plt.close('all')
    fig = plt.figure(figsize = (14, 7))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
   
    axs = [ax1, ax2, ax3, ax4]
    plabelsize = 10
    
    # Run SBDART for one of the clear points, but using the water vapor
    # profile from the smoky point. Does the altered water vapor profile
    # explain the changes observed in the GOES WV channels or MODIS/CrIS
    # skin temperatures?
    # -----------------------------------------------------------------
    for chl, tax in zip(channels, axs):
        data1_smoke = run_sbdart(chl, calc_radiance = True, \
            cris_profile = CrIS_smoke, z_mins = z_mins, z_maxs = z_maxs, \
            add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
        data1_clear1 = run_sbdart(chl, calc_radiance = True, \
            cris_profile = CrIS_clear1, z_mins = z_mins, z_maxs = z_maxs, \
            add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
        #data1_clear2 = run_sbdart(chl, calc_radiance = True, \
        #    cris_profile = CrIS_clear2, z_mins = z_mins, z_maxs = z_maxs, \
        #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
        data1_clear1_wv = run_sbdart(chl, calc_radiance = True, \
            cris_profile = CrIS_clear1_wv, z_mins = z_mins, z_maxs = z_maxs, \
            add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
        data1_clear1_both = run_sbdart(chl, calc_radiance = True, \
            cris_profile = CrIS_clear1_both, z_mins = z_mins, z_maxs = z_maxs, \
            add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])

        # Plot the data
        # -------------
        tax.plot(data1_smoke['output']['1']['1']['vza'], \
            data1_smoke['output']['1']['1']['bght_tmp'], \
            label = 'Smoky (Blue)')
        tax.plot(data1_clear1['output']['1']['1']['vza'], \
            data1_clear1['output']['1']['1']['bght_tmp'], \
            label = 'Clear 1 (Orange)')
        tax.plot(data1_clear1_wv['output']['1']['1']['vza'], \
            data1_clear1_wv['output']['1']['1']['bght_tmp'], \
            label = 'Clear 1 (Orange) w/Blue wv')
        tax.plot(data1_clear1_both['output']['1']['1']['vza'], \
            data1_clear1_both['output']['1']['1']['bght_tmp'], \
            label = 'Clear 1 (Orange) w/Blue wv and temps')
        tax.plot(data1_clear1['output']['2']['1']['vza'], \
            data1_clear1['output']['2']['1']['bght_tmp'], \
            label = 'Clear 1 (Orange) + 2 g/kg')
        tax.plot(data1_clear1['output']['3']['1']['vza'], \
            data1_clear1['output']['3']['1']['bght_tmp'], \
            label = 'Clear 1 (Orange) + 4 g/kg')
        #tax.plot(data1_clear1_wv['output']['1']['1']['vza'], \
        #    data1_clear1_wv['output']['1']['1']['bght_tmp'], \
        #    label = 'Clear 1 (Orange) w/Blue wv')
        #tax.plot(data1_clear2['output']['1']['1']['vza'], \
        #    data1_clear2['output']['1']['1']['bght_tmp'], \
        #    label = 'Clear 2 (Green)')

        tax.set_xlabel('Viewing Zenith Angle [$^{o}$]', fontsize = plabelsize, \
            weight = 'bold')
        tax.set_ylabel('Brightness temperature [K]',    fontsize = plabelsize, \
            weight = 'bold')
        wavel = np.round(data1_smoke['output']['1']['1']['avg_wavel'],2)
        ptitle = data1_smoke['info']['satellite'].upper() + \
        ' Channel ' + data1_smoke['info']['channel'] + ' (' + str(wavel) + ' μm)'
        tax.set_title(ptitle)


    #data2_smoke = run_sbdart('goes17_ch09', calc_radiance = True, \
    #    cris_profile = CrIS_smoke, z_mins = z_mins, z_maxs = z_maxs, \
    #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    #data3_smoke = run_sbdart('goes17_ch10', calc_radiance = True, \
    #    cris_profile = CrIS_smoke, z_mins = z_mins, z_maxs = z_maxs, \
    #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    #data4_smoke = run_sbdart('modis_ch31', calc_radiance = True, \
    #    cris_profile = CrIS_smoke, z_mins = z_mins, z_maxs = z_maxs, \
    #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])

    #data2_clear1 = run_sbdart('goes17_ch09', calc_radiance = True, \
    #    cris_profile = CrIS_clear1, z_mins = z_mins, z_maxs = z_maxs, \
    #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    #data3_clear1 = run_sbdart('goes17_ch10', calc_radiance = True, \
    #    cris_profile = CrIS_clear1, z_mins = z_mins, z_maxs = z_maxs, \
    #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    #data4_clear1 = run_sbdart('modis_ch31', calc_radiance = True, \
    #    cris_profile = CrIS_clear1, z_mins = z_mins, z_maxs = z_maxs, \
    #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])

    #data2_clear2 = run_sbdart('goes17_ch09', calc_radiance = True, \
    #    cris_profile = CrIS_clear2, z_mins = z_mins, z_maxs = z_maxs, \
    #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    #data3_clear2 = run_sbdart('goes17_ch10', calc_radiance = True, \
    #    cris_profile = CrIS_clear2, z_mins = z_mins, z_maxs = z_maxs, \
    #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    #data4_clear2 = run_sbdart('modis_ch31', calc_radiance = True, \
    #    cris_profile = CrIS_clear2, z_mins = z_mins, z_maxs = z_maxs, \
    #    add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])

    # Set up the figure
    
    #plot_bright_vza(data1_smoke, pax = ax1, plabelsize = 10)
    #plot_bright_vza(data2_smoke, pax = ax2, plabelsize = 10)
    #plot_bright_vza(data3_smoke, pax = ax3, plabelsize = 10)
    #plot_bright_vza(data4_smoke, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_upper_right', fontsize = 11)


    #ax1.legend()
    #ax2.legend()
    #ax3.legend()
    #ax4.legend()
    custom_lines = [Line2D([0], [0], color='tab:blue'),
                    Line2D([0], [0], color='tab:orange'),
                    Line2D([0], [0], color='tab:green'), 
                    Line2D([0], [0], color='tab:red'), 
                    Line2D([0], [0], color='tab:purple'), 
                    Line2D([0], [0], color='tab:brown')] 

    font_size = 10
    plt.legend(custom_lines, ['Smoky', 'Clear', \
            'Clear w/smoky mixing ratio', \
            'Clear w/smoky mixing ratio + temperature', \
            'Clear + 2 g/kg', \
            'Clear + 4 g/kg'],\
        fontsize = font_size, bbox_to_anchor = (2.01, 1.0), \
        ncol = 1, loc = 'upper right')

    fig.tight_layout()

    if(save):
        outname = 'modis_goes_sbdart_comps_' + row_str + '_' + date_str + '_cris.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This makes a re-creation of the SBDART GOES/MODIS figure from the paper
def process_SBDART_multi_sat_vza_tsfc(atms_file = '', save = False):
    #atms_file = 'home_dir + /Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
    # Run num 1
    z_maxs = np.array([5.])
    z_mins = np.array([0.])
    add_tmp_sfc = np.full((3, len(z_maxs)), np.nan)
    add_tmp_sfc[0,:] = 0. 
    add_tmp_sfc[1,:] = -10. 
    add_tmp_sfc[2,:] = -20. 
    
    data1 = run_sbdart('goes17_ch08', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    data2 = run_sbdart('goes17_ch09', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    data3 = run_sbdart('goes17_ch10', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    data4 = run_sbdart('modis_ch31', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_tmp_sfc = add_tmp_sfc, nzen = 9, vza = [0, 60])
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    plot_bright_vza(data4, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_right', fontsize = 11)

    fig.tight_layout()

    if(save):
        if(atms_file == ''):
            atms_add = 'mdlat_sum'
        else:
            atms_add = 'atms_file'
        outname = 'modis_goes_sbdart_comps_' + atms_add + '_tsfc.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()


# This makes a re-creation of the SBDART GOES/MODIS figure from the paper
def process_SBDART_multi_sat_vza_wv(atms_file = '', btemp = 315., save = False):
    #atms_file = 'home_dir + /Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
    # Run num 1
    z_maxs = np.array([5.])
    z_mins = np.array([0.])
    add_wv_mix = np.full((3, len(z_maxs)), np.nan)
    add_wv_mix[0,:] = 0. 
    add_wv_mix[1,:] = 2. 
    add_wv_mix[2,:] = 4. 
    
    data1 = run_sbdart('goes17_ch08', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60], \
        btemp = btemp)
    data2 = run_sbdart('goes17_ch09', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60], \
        btemp = btemp)
    data3 = run_sbdart('goes17_ch10', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60], \
        btemp = btemp)
    data4 = run_sbdart('modis_ch31', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60], \
        btemp = btemp)
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    plot_bright_vza(data4, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_right', fontsize = 11)

    fig.tight_layout()

    if(save):
        if(atms_file == ''):
            atms_add = 'mdlat_sum'
        else:
            atms_add = 'atms_file'
        outname = 'modis_goes_sbdart_comps_' + atms_add + '_wv.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This makes a re-creation of the SBDART GOES/MODIS figure from the paper
def process_SBDART_multi_sat_vza(atms_file = '', save = False):
    #atms_file = 'home_dir + /Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
    # Run num 1
    z_maxs = np.array([5.])
    z_mins = np.array([0.])
    #set_nh3_mix = np.full((3, len(z_maxs)), np.nan)
    #set_nh3_mix[0,:] = 5.0e-4
    #set_nh3_mix[1,:] = 1.0e-3
    #set_nh3_mix[2,:] = 1.5e-3
    #set_co_mix = np.full((3, len(z_maxs)), np.nan)
    #set_co_mix[0,:] = 0.15 
    #set_co_mix[1,:] = 0.30
    #set_co_mix[2,:] = 0.45
    #set_no2_mix = np.full((3, len(z_maxs)), np.nan)
    #set_no2_mix[0,:] = 2.3e-5
    #set_no2_mix[1,:] = 4.6e-5
    #set_no2_mix[2,:] = 6.9e-5
    #set_ch4_mix = np.full((3, len(z_maxs)), np.nan)
    #set_ch4_mix[0,:] = 1.74 
    #set_ch4_mix[1,:] = 2.74
    #set_ch4_mix[2,:] = 3.74
    #set_co2_mix = np.full((3, len(z_maxs)), np.nan)
    #set_co2_mix[0,:] = 360. 
    #set_co2_mix[1,:] = 500.
    #set_co2_mix[2,:] = 1000.
    add_wv_mix = np.full((3, len(z_maxs)), np.nan)
    add_wv_mix[0,:] = 0. 
    add_wv_mix[1,:] = 2. 
    add_wv_mix[2,:] = 4. 
    
    data1 = run_sbdart('goes17_ch08', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
    data2 = run_sbdart('goes17_ch09', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
    data3 = run_sbdart('goes17_ch10', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
    data4 = run_sbdart('modis_ch31', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    plot_bright_vza(data4, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_right', fontsize = 11)

    fig.tight_layout()

    if(save):
        if(atms_file == ''):
            atms_add = 'mdlat_sum'
        else:
            atms_add = 'atms_file'
        outname = 'modis_goes_sbdart_comps_' + atms_add + '_wv.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This makes a re-creation of the SBDART GOES/MODIS figure from the paper
def process_SBDART_multi_sat_vza_co2(atms_file = '', save = False):
    #atms_file = 'home_dir + /Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
    # Run num 1
    z_maxs = np.array([5.])
    z_mins = np.array([0.])
    set_co2_mix = np.full((3, len(z_maxs)), np.nan)
    set_co2_mix[0,:] = 360. 
    set_co2_mix[1,:] = 500.
    set_co2_mix[2,:] = 1000.
    #add_wv_mix = np.full((3, len(z_maxs)), np.nan)
    #add_wv_mix[0,:] = 0. 
    #add_wv_mix[1,:] = 2. 
    #add_wv_mix[2,:] = 4. 
    
    data1 = run_sbdart('goes17_ch08', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_co2_mix = set_co2_mix, nzen = 9, vza = [0, 60])
    data2 = run_sbdart('goes17_ch09', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_co2_mix = set_co2_mix, nzen = 9, vza = [0, 60])
    data3 = run_sbdart('goes17_ch10', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_co2_mix = set_co2_mix, nzen = 9, vza = [0, 60])
    data4 = run_sbdart('modis_ch31', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_co2_mix = set_co2_mix, nzen = 9, vza = [0, 60])
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    plot_bright_vza(data4, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_right', fontsize = 11)

    fig.tight_layout()

    if(save):
        if(atms_file == ''):
            atms_add = 'mdlat_sum'
        else:
            atms_add = 'atms_file'
        outname = 'modis_goes_sbdart_comps_' + atms_add + '_co2.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This makes a re-creation of the SBDART GOES/MODIS figure from the paper
def process_SBDART_multi_sat_vza_ch4(atms_file = '', save = False):
    #atms_file = 'home_dir + /Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
    # Run num 1
    z_maxs = np.array([5.])
    z_mins = np.array([0.])
    set_ch4_mix = np.full((3, len(z_maxs)), np.nan)
    set_ch4_mix[0,:] = 1.74 
    set_ch4_mix[1,:] = 2.74
    set_ch4_mix[2,:] = 3.74
    
    data1 = run_sbdart('goes17_ch08', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_ch4_mix = set_ch4_mix, nzen = 9, vza = [0, 60])
    data2 = run_sbdart('goes17_ch09', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_ch4_mix = set_ch4_mix, nzen = 9, vza = [0, 60])
    data3 = run_sbdart('goes17_ch10', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_ch4_mix = set_ch4_mix, nzen = 9, vza = [0, 60])
    data4 = run_sbdart('modis_ch31', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_ch4_mix = set_ch4_mix, nzen = 9, vza = [0, 60])
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    plot_bright_vza(data4, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_right', fontsize = 11)

    fig.tight_layout()

    if(save):
        if(atms_file == ''):
            atms_add = 'mdlat_sum'
        else:
            atms_add = 'atms_file'
        outname = 'modis_goes_sbdart_comps_' + atms_add + '_ch4.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This makes a re-creation of the SBDART GOES/MODIS figure from the paper
def process_SBDART_multi_sat_vza_no2(atms_file = '', save = False):
    #atms_file = 'home_dir + /Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
    # Run num 1
    z_maxs = np.array([5.])
    z_mins = np.array([0.])
    set_no2_mix = np.full((3, len(z_maxs)), np.nan)
    set_no2_mix[0,:] = 2.3e-5
    set_no2_mix[1,:] = 4.6e-5
    set_no2_mix[2,:] = 6.9e-5
    #add_wv_mix = np.full((3, len(z_maxs)), np.nan)
    #add_wv_mix[0,:] = 0. 
    #add_wv_mix[1,:] = 2. 
    #add_wv_mix[2,:] = 4. 
    
    data1 = run_sbdart('goes17_ch08', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_no2_mix = set_no2_mix, nzen = 9, vza = [0, 60])
    data2 = run_sbdart('goes17_ch09', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_no2_mix = set_no2_mix, nzen = 9, vza = [0, 60])
    data3 = run_sbdart('goes17_ch10', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_no2_mix = set_no2_mix, nzen = 9, vza = [0, 60])
    data4 = run_sbdart('modis_ch31', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_no2_mix = set_no2_mix, nzen = 9, vza = [0, 60])
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    plot_bright_vza(data4, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_right', fontsize = 11)

    fig.tight_layout()

    if(save):
        if(atms_file == ''):
            atms_add = 'mdlat_sum'
        else:
            atms_add = 'atms_file'
        outname = 'modis_goes_sbdart_comps_' + atms_add + '_no2.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This makes a re-creation of the SBDART GOES/MODIS figure from the paper
def process_SBDART_multi_sat_vza_co(atms_file = '', save = False):
    #atms_file = 'home_dir + /Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
    # Run num 1
    z_maxs = np.array([5.])
    z_mins = np.array([0.])
    set_co_mix = np.full((3, len(z_maxs)), np.nan)
    set_co_mix[0,:] = 0.15 
    set_co_mix[1,:] = 0.30
    set_co_mix[2,:] = 0.45
    #add_wv_mix = np.full((3, len(z_maxs)), np.nan)
    #add_wv_mix[0,:] = 0. 
    #add_wv_mix[1,:] = 2. 
    #add_wv_mix[2,:] = 4. 
    
    data1 = run_sbdart('goes17_ch08', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_co_mix = set_co_mix, nzen = 9, vza = [0, 60])
    data2 = run_sbdart('goes17_ch09', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_co_mix = set_co_mix, nzen = 9, vza = [0, 60])
    data3 = run_sbdart('goes17_ch10', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_co_mix = set_co_mix, nzen = 9, vza = [0, 60])
    data4 = run_sbdart('modis_ch31', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_co_mix = set_co_mix, nzen = 9, vza = [0, 60])
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    plot_bright_vza(data4, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_right', fontsize = 11)

    fig.tight_layout()

    if(save):
        if(atms_file == ''):
            atms_add = 'mdlat_sum'
        else:
            atms_add = 'atms_file'
        outname = 'modis_goes_sbdart_comps_' + atms_add + '_co.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This makes a re-creation of the SBDART GOES/MODIS figure from the paper
def process_SBDART_multi_sat_vza_nh3(atms_file = '', save = False):
    #atms_file = 'home_dir + /Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
    # Run num 1
    z_maxs = np.array([5.])
    z_mins = np.array([0.])
    set_nh3_mix = np.full((3, len(z_maxs)), np.nan)
    set_nh3_mix[0,:] = 5.0e-4
    set_nh3_mix[1,:] = 1.0e-3
    set_nh3_mix[2,:] = 1.5e-3
    #add_wv_mix = np.full((3, len(z_maxs)), np.nan)
    #add_wv_mix[0,:] = 0. 
    #add_wv_mix[1,:] = 2. 
    #add_wv_mix[2,:] = 4. 
    
    data1 = run_sbdart('goes17_ch08', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_nh3_mix = set_nh3_mix, nzen = 9, vza = [0, 60])
    data2 = run_sbdart('goes17_ch09', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_nh3_mix = set_nh3_mix, nzen = 9, vza = [0, 60])
    data3 = run_sbdart('goes17_ch10', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_nh3_mix = set_nh3_mix, nzen = 9, vza = [0, 60])
    data4 = run_sbdart('modis_ch31', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        set_nh3_mix = set_nh3_mix, nzen = 9, vza = [0, 60])
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    plot_bright_vza(data4, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_right', fontsize = 11)

    fig.tight_layout()

    if(save):
        if(atms_file == ''):
            atms_add = 'mdlat_sum'
        else:
            atms_add = 'atms_file'
        outname = 'modis_goes_sbdart_comps_' + atms_add + '_nh3.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

    
def process_SBDART_single_sat_wv_tmp_vza(satellite = 'modis_ch31',\
        atms_file = '', save = False):
    #atms_file = 'home_dir + /Research/SBDART/data/model/210722_220000_XXX_HRRR.txt'
    # Run num 1
    z_maxs = np.array([5.])
    z_mins = np.array([0.])
    add_wv_mix_vals = np.arange(0, 4.1, 0.5)
    add_wv_mix = np.full((len(add_wv_mix_vals), len(z_maxs)), np.nan)
    for jj in range(len(add_wv_mix_vals)):
        add_wv_mix[jj,:] = add_wv_mix_vals[jj]
    ##!#add_wv_mix[0,:] = 0. 
    ##!#add_wv_mix[1,:] = 2. 
    ##!#add_wv_mix[2,:] = 4. 

    satellites = ['goes17_ch08','goes17_ch09','goes17_ch10','goes17_ch13',\
        'modis_ch31']
    titles = ['GOES-17 6.18 μm','GOES-17 6.95 μm','GOES-17 7.34 μm',\
        'GOES-17 10.35 μm','MODIS 11 μm']

    plt.close('all')
    fig = plt.figure()
    #ax1 = fig.add_subplot(1,1,1)

    for kk, satellite in enumerate(satellites):

        print(kk, satellite)
        tax = fig.add_subplot(2,3,kk + 1)
        ##!#tmp_diffs = np.arange(1.)
        tmp_diffs = np.arange(0., 11., 5.)
        bght_tmps = np.full((len(add_wv_mix), len(tmp_diffs)), np.nan)

        for ii, tmp_diff in enumerate(tmp_diffs):
            data = run_sbdart(satellite, calc_radiance = True, \
                atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
                add_wv_mix = add_wv_mix, add_tmp_sfc = -tmp_diff)
                #add_wv_mix = add_wv_mix)

            for jj, wvmix in enumerate(add_wv_mix[:,0]):
                bght_tmps[jj,ii] = data['output'][str(jj+1)]['1']['bght_tmp']
            ##!#bght_tmps[0,ii] = data['output']['1']['1']['bght_tmp']
            ##!#bght_tmps[1,ii] = data['output']['2']['1']['bght_tmp']
            ##!#bght_tmps[2,ii] = data['output']['3']['1']['bght_tmp']
        
        ##!#ax1.plot(bght_tmps[2,:] - bght_tmps[0,0], label = titles[kk])
        #tax.plot(bght_tmps
            tax.plot(add_wv_mix, bght_tmps[:,ii], label = str(int(tmp_diff)))
        #tax.plot(bght_tmps[1,:], label = 'wo + 2')
        #tax.plot(bght_tmps[2,:], label = 'wo + 4')
        tax.set_xlabel('Water Vapor Mixing Ratio Increase')
        #tax.set_xlabel('Surface temperature reduction [K]')
        tax.set_ylabel('BT')
        #tax.set_ylabel('TOA ΔBT')
        tax.set_title(satellite)
        tax.legend()
    ##!#ax1.set_xlabel('Surface temperature reduction [K]')
    ##!#ax1.set_ylabel('TOA ΔBT')
    ##!#ax1.set_title('SBDART-simulated TOA BT Cooling\n' + \
    ##!#    'w = w$_{o}$ + 4 g kg$^{-1}$ in lowest 5 km')
    ##!#ax1.legend()
    fig.tight_layout()
    plt.show()
    return 
    ##!#ax1 = fig.add_subplot(1,3,1)
    ##!#ax2 = fig.add_subplot(1,3,2)
    ##!#ax3 = fig.add_subplot(1,3,3)
    ##!#data1 = run_sbdart(satellite, calc_radiance = True, \
    ##!#    atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
    ##!#    add_wv_mix = add_wv_mix,  \
    ##!#    nzen = 9, vza = [0, 60])
    ##!#data2 = run_sbdart(satellite, calc_radiance = True, \
    ##!#    atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
    ##!#    add_wv_mix = add_wv_mix, add_tmp_sfc = add_tmp[1],\
    ##!#    nzen = 9, vza = [0, 60])
    ##!#data3 = run_sbdart(satellite, calc_radiance = True, \
    ##!#    atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
    ##!#    add_wv_mix = add_wv_mix, add_tmp_sfc = add_tmp[2],\
    ##!#    nzen = 9, vza = [0, 60])
    ##!#plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    ##!#plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    ##!#plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    ##!#fig.tight_layout()
    ##!#plt.show()
    ##!#return

    kk = 0
    for ii in range(len(add_wv_mix)):
        print(ii, add_wv_mix[ii])
        for jj in range(len(add_tmp)):
            print('  ', jj, add_tmp[jj])
            tax = fig.add_subplot(3,3,kk)

            data1 = run_sbdart(satellite, calc_radiance = True, \
                atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
                add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])

            kk += 1

    return
 
    data1 = run_sbdart('goes17_ch08', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
    data2 = run_sbdart('goes17_ch09', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
    data3 = run_sbdart('goes17_ch10', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
    data4 = run_sbdart('modis_ch31', calc_radiance = True, \
        atms_file = atms_file, z_mins = z_mins, z_maxs = z_maxs, \
        add_wv_mix = add_wv_mix, nzen = 9, vza = [0, 60])
    
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plot_bright_vza(data1, pax = ax1, plabelsize = 10)
    plot_bright_vza(data2, pax = ax2, plabelsize = 10)
    plot_bright_vza(data3, pax = ax3, plabelsize = 10)
    plot_bright_vza(data4, pax = ax4, plabelsize = 10)

    # Add dashed lines at the correct VZAs
    # ------------------------------------
    ax1.axvline(49.61378, linestyle = '--', color = 'black')
    ax2.axvline(49.61378, linestyle = '--', color = 'black')
    ax3.axvline(49.61378, linestyle = '--', color = 'black')
    ax4.axvline(3.15, linestyle = '--', color = 'black')

    plot_subplot_label(ax1, '(a)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'upper_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'upper_right', fontsize = 11)

    fig.tight_layout()

    if(save):
        if(atms_file == ''):
            atms_add = 'mdlat_sum'
        else:
            atms_add = 'atms_file'
        outname = 'modis_goes_sbdart_comps_' + atms_add + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()
