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

h_const = 6.626e-34 #J*s
k_const = 1.3806e-23 #J/K
c_const = 3e8 # m/s
sb_path = '/home/bsorenson/Research/SBDART/'

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

# Reads the SBDART model profile. If infile is '', the default
# profile will be used. If not, a model profile from 
# /home/bsorenson/Research/SBDART/model/ and reformats into the 
# right format.
def read_atmos_profile(infile = ''):
    atms_file = sb_path + 'src/atms_mdlat_sum.dat'
    #atms_file = sb_path + 'src/atms_orig.dat'
    temp_data = pd.read_csv(atms_file,skiprows = 0, names = \
        ['z','p','t','wv','o2'], delim_whitespace = True)
    temp_data = temp_data[1:]

    if((infile == '') | (infile.strip().split('/')[-1][:4] == 'atms')):
        atms_file = infile
        print('Reading atmosphere from ', atms_file)
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

def prep_sbdart_input(isat = -1, iout = 5, nzen = None, uzen = [0, 60]):
 
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Input file prep
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   
    print('Prepping INPUT file')
 
    #iout = 20
    #iout = 0
    
    with open('./INPUT','w') as fout:
        #fout.write(" &INPUT\n" + \
        out_str = " &INPUT\n" + \
            "   idatm=0\n" + \
            "   isat={0}\n".format(isat) + \
            "   wlinc=0.01\n" + \
            "   iout={0}\n".format(iout)
            #"   wlinf=8.00\n" + \
            #"   wlsup=12.00\n" + \
            #"   nzen={0}\n".format(nzen) + \
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
        out_str = out_str + \
            "   phi=0\n" + \
            "/\n"

        fout.write(out_str)

def prep_atmos_file(data, \
        zmin = None, zmax = None, \
        set_tmp = None, add_tmp = None, \
        set_wv_mix = None, add_wv_mix = None, \
        out_path = './'):
    #print("processing for wv multiplier = ",np.round(mtp,2))

    # Make copies for editing
    # -----------------------
    new_data = data.copy(deep = True)

    if(set_tmp is not None):
        # Set the layer temps
        # -------------------
        new_data['t'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = set_tmp

    elif(add_tmp is not None): 
        print("Adding temp of ", add_tmp)
        # Add the layer temps
        # -------------------
        new_data['t'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] += add_tmp
 
    if(set_wv_mix is not None):  
        # Convert the added water vapor mixing ratio to added
        # absolute humidity
        # ----------------------------------------------------
        vap_press = mpcalc.vapor_pressure(\
            np.array(data['p']) * units('mbar'), set_wv_mix * units('g/kg'))
        new_hum = (vap_press.to('Pa') / (461.5 * (data['t']))) * 1e3

        # Set the extra absolute humidity
        # -------------------------------
        new_data['wv'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] = new_hum

    elif(add_wv_mix is not None):  
        # Convert the added water vapor mixing ratio to added
        # absolute humidity
        # ----------------------------------------------------
        vap_press = mpcalc.vapor_pressure(\
            np.array(data['p']) * units('mbar'), add_wv_mix * \
            units('g/kg'))
        new_hum = (vap_press.to('Pa') / (461.5 * (data['t']))) * 1e3

        # Add the extra absolute humidity to the original values
        # ------------------------------------------------------
        new_data['wv'][(new_data['z'] >= zmin) & \
            (new_data['z'] <= zmax)] += new_hum
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

def run_sbdart(satellite, calc_radiance, run = True, vza = None, atms_file = ''):
    if(calc_radiance):
        file_adder = '_rad.dat'
    else:
        file_adder = '.dat'

    if(vza is None):
        # Determine viewing angle values
        # ------------------------------
        if(satellite.strip().split('_')[0] == 'modis'):
            uzen = 3.15
        elif(satellite.strip().split('_')[0] == 'goes17'):
            uzen = 49.61378
        else:
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

        # Prep the INPUT namelist file
        # ----------------------------
        #prep_sbdart_input(nzen = 9, uzen = [0,60])
        prep_sbdart_input(uzen = uzen)
        
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #
        # SBDART runs
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        
        # Read the atmos file
        # -------------------
        # NOTE: If reading in a model profile for the atmosphere, 
        #       read it in here.
        data = read_atmos_profile(infile = atms_file)
        ##!#data = pd.read_csv(infile,skiprows = 0, names = \
        ##!#    ['z','p','t','wv','o2'], delim_whitespace = True)
        
        # Added water vapor mixing ratio, must be converted to 
        # water vapor density.
        # ----------------------------------------------------
        adders  = np.array([0.0, 2.0, 4.0])
        z_mins = np.arange(0.,6.)
        z_maxs = np.arange(5.,11.)
        #z_mins = np.arange(0.,6.)
        #z_maxs = np.arange(2.,8.)
        set_wvs = np.full(z_mins.shape, 4.)
        #multipliers = np.array([1.0, 2.0, 3.0, 3.5])
        #multipliers = np.arange(0.1,2.1,0.1)
        
        outfile_list = []
      
        data_dict = {}
        data_dict['output'] = {}
        data_dict['info'] = {}
 
        for ii in range(len(z_mins)): 
        #for adr in adders:
            # Prep the atmosphere file 
            # ------------------------
            prep_atmos_file(data, zmin = z_mins[ii], \
                zmax = z_maxs[ii], add_wv_mix = set_wvs[ii])
                #set_wv_mix = adr)
            run_adder = 'zmin'+str(int(z_mins[ii])) + '_' + \
                        'zmax'+str(int(z_maxs[ii])) + '_' + \
                        'setwv'+str(int(set_wvs[ii]))
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
            data_dict['output'][d_key] = {}
            data_dict['output'][d_key]['data'] = {}
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
                data_dict['output'][d_key]['vza'] = vza
                data_dict['output'][d_key]['data']['wavelength'] = wavel
                data_dict['output'][d_key]['data']['filter']     = fvals
                data_dict['output'][d_key]['data']['radiance']   = total_radiances
                data_dict['output'][d_key]['data']['down_sol_flux_toa'] = topdwn
                data_dict['output'][d_key]['data']['up_rad_flux_toa']   = topup
                data_dict['output'][d_key]['data']['dir_sol_flux_toa']  = topdwn
                data_dict['output'][d_key]['data']['down_rad_flux_sfc'] = botdwn
                data_dict['output'][d_key]['data']['up_rad_flux_sfc']   = botup 
                data_dict['output'][d_key]['data']['dir_sol_flux_sfc']  = botdir
    
            # Calculated weighted average wavelength
            # --------------------------------------
            data_dict['output'][d_key]['avg_wavel'] = \
                np.sum(data_dict['output'][d_key]['data']['wavelength'] * \
                data_dict['output'][d_key]['data']['filter']) / \
                np.sum(data_dict['output'][d_key]['data']['filter'])
    
            # Calculated weighted average radiance  
            # ------------------------------------
            data_dict['output'][d_key]['avg_rads'] = \
                np.array([np.sum(data_dict['output'][d_key]['data']['radiance'][:,jj] * \
                data_dict['output'][d_key]['data']['filter']) / \
                np.sum(data_dict['output'][d_key]['data']['filter']) \
                for jj in range(len(data_dict['output'][d_key]['vza']))])

            #data_dict[d_key]['avg_rad'] = \
            #    np.sum(data_dict[d_key]['data']['radiance'] * data_dict[d_key]['data']['filter']) / \
            #    np.sum(data_dict[d_key]['data']['filter'])

            # Calculate temperatures for all wavelengths
            # ------------------------------------------
            data_dict['output'][d_key]['data']['bght_tmp_all'] = \
                np.array([inverse_planck(data_dict['output'][d_key]['data']['radiance'][:,jj], \
                data_dict['output'][d_key]['data']['wavelength']) \
                for jj in range(len(data_dict['output'][d_key]['vza']))]).transpose()

            #data_dict[d_key]['data']['bght_tmp_all'] = \
            #    inverse_planck(data_dict[d_key]['data']['radiance'], \
            #    data_dict[d_key]['data']['wavelength'])

            # Calculate temperature for filter averaged values
            # ------------------------------------------------
            data_dict['output'][d_key]['bght_tmp'] = \
                inverse_planck(data_dict['output'][d_key]['avg_rads'], \
                data_dict['output'][d_key]['avg_wavel'])
 
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

def plot_bright_vza(modis_data1, pax = None, save = False):
    in_pax = True
    if(pax is None):
        in_pax = False 
        plt.close('all')
        fig1 = plt.figure(figsize = (8,4))
        pax = fig1.add_subplot(1,1,1)
  
    #wvs = [20, 100, 150, 200, 250]

    # Determine the plot type by the keys in the structure
    # ----------------------------------------------------
    key1 = list(modis_data1['output'].keys())[0]
    key2 = list(modis_data1['output'].keys())[1]

    
 
    for dkey in modis_data1['output'].keys():
    #for wv in modis_data1['adders']:
    #for wv in modis_data1['multipliers']:
        int_wv = int(wv)
        #int_wv = int(wv*100.)
        wavel = np.round(modis_data1[int_wv]['avg_wavel'],2)
        if(wv == 0.):
        #if(wv == 1.):
            label = 'w$_{0}$'
        else:
            label = 'w$_{0}$ + ' + str(int(int_wv))+' g/kg'
            #label = str(int(int_wv - 100.))+'% higher WV'
        # Plot the correct wv amounts
        pax.plot(modis_data1[int_wv]['vza'],  \
            modis_data1[int_wv]['bght_tmp'], label = label)

    pax.set_xlabel('Viewing Zenith Angle [$^{o}$]', fontsize = 12, weight = 'bold')
    pax.set_ylabel('Brightness temperature [K]',    fontsize = 12, weight = 'bold')
    pax.set_title('SBDART-simulated ' + modis_data1['satellite'].upper() + \
        ' Channel ' + modis_data1['channel'] + ' (' + str(wavel) + ' μm)')
    pax.legend(fontsize = 9, loc = 'lower left')
    if(not in_pax):
        if(save):
            outname = 'sbdart_'+modis_data1['run_name']+'_bright_vza.png'
            fig1.savefig(outname, dpi = 300)
            print("Saved image",outname)
        else:
            plt.show() 

#plot_bright(modis31, modis32)
