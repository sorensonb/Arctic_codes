#!/usr/bin/env python

"""



"""

import numpy as np
import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# 1  SOYB    "soybean"
# 2  CO30    39.0    0.53    25.0     8.0     6.0    0.50   15.05   50.95    0.10
# 3  CO36    39.0    0.53    25.0     8.0     6.0    0.50   15.05   50.95    0.10
# 4  CO38    39.0    0.53    25.0     8.0     6.0    0.50   15.05   50.95    0.10
# 5  CO40    39.0    0.53    25.0     8.0     6.0    0.50   15.05   50.95    0.10
# 6  SO30    37.2    0.45    30.0     8.0     3.5    0.64   20.05   60.95    0.50
# 7  SO36    37.2    0.45    30.0     8.0     3.5    0.64   20.05   60.95    0.50
# 8  SO38    37.2    0.45    30.0     8.0     3.5    0.64   20.05   60.95    0.50
# 9  SO40    37.2    0.45    30.0     8.0     3.5    0.64   20.05   60.95    0.50
#10  SHCN    31.0    0.50    30.0     8.0     6.0    0.63   15.05   50.95    1.00
#11  WHET    "wheat"
#12  TXWT    30.0    0.40    15.0     0.0     5.0    0.50    5.05   45.95    1.00
#13  OATS    "oats"
#14  SUNF    "sunflower"
#15  SUGC    "sugar cane"
#16  ALFA    "alfalfa?"
#17  TIMO    20.0    0.25    20.0     4.0     3.0    0.90   15.01   50.95    0.50
#18  JHGR    35.0    0.15    30.0    11.0     5.0    0.50   15.05   57.95    1.00
#19  GTFX    37.0    0.10    25.0    10.0     5.0    0.60   15.05   54.95    2.00
#20  GRFX    37.0    0.15    25.0    10.0     3.5    0.60   15.05   54.95    1.00
#21  COCB    33.0    0.40    25.0    10.0     2.5    0.60   15.05   50.95    0.50
#22  VELV    27.0    0.01    25.0     4.0     3.8    0.90   15.05   50.95    1.00
#23  CEAT    18.0    0.01    15.0     0.0     1.0    0.60    5.05   50.95    1.00
#24  RYEA    30.0    0.42    15.0     0.0     3.0    0.80   30.01   50.95    1.00
#25  SWCH    49.0    0.01    25.0    12.0    12.0    0.70   10.20   20.95    1.00
#26  PAST    30.0    0.90    25.0     8.0     5.0    0.85    5.05   49.95    2.00
#27  SIDE    18.0    0.01    25.0    12.0     2.5    0.35    5.05   30.70    0.10
#28  BUFF    20.0    0.01    25.0    12.0     2.0    0.50   10.20   20.95    0.20
#29  BLUG    18.0    0.01    25.0    12.0     2.0    0.35    5.05   30.70    0.10
#30  LBL2    34.0    0.01    25.0    12.0     3.0    0.67    6.16   64.99    0.01
#31  BBLS    34.0    0.01    25.0    12.0     7.6    0.35    5.10   25.70    0.50
#32  EGAM    50.0    0.01    25.0    12.0     4.8    0.40    5.18   25.90    0.20
#33  INDN    34.0    0.01    25.0    12.0     5.0    0.35    5.10   25.70    0.50
#34  BLAG    18.0    0.01    25.0    12.0     1.5    0.35    5.05   30.70    0.10
#35  COAS    16.0    0.01    25.0    12.0     5.9    0.70   30.10   95.71    1.00
#36  BAHI    15.0    0.01    25.0    12.0     2.0    0.70   22.12   54.62    1.00
#37  ORCH    20.0    0.25    20.0     4.0     2.5    0.90    5.05   50.95    0.50
#38  OBLU    47.0    0.01    25.0    12.0     7.6    0.69   13.14   22.58    0.10
#39  CWGR    35.0    0.10    25.0     6.0     5.0    0.85   35.02   62.95    2.00
#40  POPL    "poplar"
#41  PINE    "pine"
#42   OAK    "oak"
#43  MSQT    "mesquite"
#44  ASPN    "aspen"
#45  WSPR    15.0    0.76    30.0     0.0     5.0    0.99   15.70   25.99    1.00
#46  BSPR    15.0    0.76    30.0     0.0     5.0    0.99   15.70   25.99    1.00
#47  LDGP    15.0    0.76    30.0     0.0     5.0    0.99   15.70   25.99    1.00
#48  SALT    30.0    0.05    30.0    10.0     5.0    0.50    5.05   40.95    0.05
#49  SAGE    11.5    0.76    30.0    10.0     2.0    0.90   15.05   50.80    1.00
#50  BSQU    18.0    0.01    15.0     0.0     0.5    0.60   31.07   57.95    1.00
#51  TXWG    18.0    0.01    15.0     0.0     1.5    0.60    5.05   50.95    1.00
#52  COTS    17.5    0.12    27.5    12.0     1.1    0.85   15.01   50.95    2.00
#53  SNLL    49.0    0.01    25.0    12.0     5.5    0.70   10.20   20.95    1.00
#54  SNUL    49.0    0.01    25.0    12.0     4.0    0.70   10.20   20.95    1.00
#55  SSLL    49.0    0.01    25.0    12.0     6.0    0.70   10.20   20.95    1.00
#56  SSUL    49.0    0.01    25.0    12.0     5.0    0.70   10.20   20.95    1.00
#57  TFES    14.0    0.01    16.0     5.0     3.0    0.50   10.20   20.95    0.01
#58  CEDR    16.0    0.01    30.0    10.0    12.0    0.99   20.20   99.99    0.01
#59  SQUT    26.1    0.01    15.0     0.0     1.3    0.60   31.07   57.95    1.00
#60  BBWG    16.4    0.01    15.0     0.0     1.3    0.60   31.07   57.95    1.00
#61  INRG    14.4    0.01    15.0     0.0     0.6    0.60   31.07   57.95    1.00
#62  BHMU    18.0    0.01    25.0    12.0     1.0    0.35    5.05   30.70    0.10
#63  GLLG    15.0    0.01    25.0    12.0     2.3    0.35    5.05   30.70    0.10
#64  SCTN    40.0    0.01    25.0    12.0    12.0    0.70   10.20   20.95    1.00
#65  VNMQ    16.0    0.01    25.0    12.0     5.9    0.70   30.10   95.71    1.00
#66  ESOG    49.0    0.90    27.5    10.0     5.0    0.55   15.01   50.95    0.01
#67  CANP    34.0    0.30    21.0     5.0     5.5    0.78   51.10   59.56    0.20
#68  CANA    34.0    0.23    21.0     5.0     5.5    0.78   51.10   59.56    0.20

#############################################################################
#
#                           Model running functions
#
#############################################################################
    
# Set up dictionary holding the conversion factors for metric
# tons per hectare to bushels per acre.
# ---------------------------------------------------------------
conv_dict = {
    'WHET': 14.87, \
    'SOYB': 14.87, \
    'CORN': 15.93, \
}

# Defaults:
# - crop_type     = 11 (wheat)
# - solar_units   = 1800 
# - plant_month   = 4
# - plant_day     = 10
# - harvest_month = 9
# - harvest_day   = 30
# - output_unit   = 'TPH' (BPA = bushels per acre, \
#                          TPH = metric tons per hectare)
def run_MSPAM_yearly_yields(year = None, county = None, \
        crop_type = 11, solar_units = 1800, \
        plant_month = 4, plant_day = 10, \
        harvest_month = 9, harvest_day = 30, \
        output_unit = 'TPH'):

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #                Set parameters and create temporary directories
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    os.system('mkdir ./tmp_outzip_dir')
    
    # Loop over the yearly yields file, extracting the mukeys and weather
    # files for the counties
    obs_dict = pd.read_csv('../run_batch/yearly_yield_2001_2015.csv')
    
    if(year is not None):
        years = [year]
    else:
        years = np.arange(2001, 2016)

    if(county is not None):
        counties = [county]
    else:
        counties = obs_dict['County'].values
    
    
    output_dict = {}
    for cty in counties:
        output_dict[cty] = {}        
    
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #                        Begin processing loop
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    crop_name = "XXXX"
    for year in years:
    
        #output_dict[year] = {}    
    
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #
        #                        Prep the almanac.input file
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        infile = 'almanac_original.input'
        outfile = 'almanac.input'
        
        with open(outfile, 'w') as fout:
            with open(infile, 'r') as fin:
                # Read a line 
                for ii, line in enumerate(fin):
                    #print(ii, line.strip())
        
                    if(ii == 0):
                        # Line idx 0 : date line
                        line = '    ' + str(year) + '   1   1    0   0   1\n'
                    elif(ii == 63):
                        # Line idx 63 : planting / solar line
                        '   4  10  17  11   0    1800       0      50'
                        line = '{0:4d}{1:4d}  17{2:4d}   0    {3:4d}       0      50\n'.format(\
                            plant_month, plant_day, crop_type,solar_units)
                    elif(ii == 64):
                        # Line idx 64 : harvesting line
                        '   9  30  46  11   0       0       0       0'
                        line = '{0:4d}{1:4d}  46{2:4d}   0       0       0       0\n'.format(\
                            harvest_month, harvest_day, crop_type)
        
        
                    fout.write(line)
                              
           
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #
        #                        Prep the cropALNC.dat file ?
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #
        #                     Select the county and weather files
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        for jj, cty in enumerate(counties):
            output_dict[cty][str(year)] = {}
    
            mukey   = obs_dict[obs_dict['County'] == cty]['MUKEY'].values[0]
            wx_file = obs_dict[obs_dict['County'] == cty]['WTH'].values[0]
    
            # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
            #
            #                           Run the model
            #
            # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
            cmnd = '../bin/model_almanac ' + str(mukey) + ' ' + \
                str(wx_file) + ' 1'
            print(cmnd)
            os.system(cmnd)
            
            # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
            #
            #             Parse the output file and extract yield
            #
            # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
            
            # Move the output zip file into a new, burner directory
            # -----------------------------------------------------
            outzip_name = 'almanac_Mu' + str(mukey) + '_Wx' + \
                str(wx_file) + '.tar.gz'
            outzip_path = '../output/' + outzip_name
            
            cmnd = 'mv ' + outzip_path + ' ./tmp_outzip_dir/'
            print(cmnd)
            os.system(cmnd)
            
            # Open the zip file and extract contents
            # --------------------------------------
            cmnd = 'tar -xvzf tmp_outzip_dir/' + outzip_name + \
                ' -C ./tmp_outzip_dir/'
            print(cmnd)
            os.system(cmnd)
            
            # Open the output file
            # --------------------
            infile = './tmp_outzip_dir/ALNC_Run.out'
            
            with open(infile, 'r') as fin:
                flines = fin.readlines()
    
            # Figure out where the output is located in the file
            # --------------------------------------------------
            good_output = False
            for ii, line in enumerate(flines):
                tline = line.strip()
                if(tline[:11] == 'CROP    YLD'):
                    start_idx = ii
                    good_output = True
    
    
            if(good_output):
                out_split = flines[start_idx + 2].strip().split()
                print(out_split)
                
                # Parse the desired output from the output line. Also, convert
                # the yields from metric tons per hectare to bushels per acre.
                #
                # https://www.extension.iastate.edu/agdm/wholefarm/pdf/c6-80.pdf
                #
                # 1 acre = 2.471 hectares
                #
                # BPA = bushels per acre
                # MTH = metric tons per hectare
                #
                # Wheat:    60 lbs (27.22 kg) per bushel
                # Soybeans: 60 lbs (27.22 kg) per bushel
                # Corn:     56 lbs (25.40 kg) per bushel
                #
                # Example: corn yield of 200 BPA
                #     
                #   200 bu * 56 lb/bu = 11,200 lbs
                #   11,200 lbs * 0.4536 kg/lb = 5080 kg
                #   5080 kg / acre * 2.471 acres / ha = 12553 kg / ha
                #   = 12.55 metric tons per hectare 
                #
                #
                # Reversed:
                #
                #   12.55 metric tons per hectare
                #   12.55 MTH * (1000 kg / ton) = 12553 kg / ha
                #   12553 kg/ha * (0.4047 ha / ac) = 5080 kg / ac
                #   5080 kg/ac * (2.205 lb / kg) = 11,200 lb / ac
                #   11,200 lb/ac / (56 lb/bu) = 200 bu
                #
                # --------------------------------------------------------------- 
                crop_name = out_split[0]
                cp_yield = float(out_split[1])
                if(crop_name in conv_dict):
                    if(output_unit == 'BPA'):
                        cp_yield *= conv_dict[crop_name]
                else:
                    print("ERROR: Unknown crop type. Printing output in ")
                    print("       tons per hectare.")
                    output_unit = 'TPH'

                output_dict[cty][str(year)]['YLD'] = cp_yield
                output_dict[cty][str(year)]['BIOM'] = float(out_split[2])
                output_dict[cty][str(year)]['RAD'] = float(out_split[3])
                output_dict[cty][str(year)]['HU'] = float(out_split[4])
                output_dict[cty][str(year)]['RD'] = float(out_split[5])
            else:
                print("ERROR: No output for ", cty)
                output_dict[cty][str(year)]['YLD'] =  -9999.
                output_dict[cty][str(year)]['BIOM'] = -9999.
                output_dict[cty][str(year)]['RAD'] = -9999.
                output_dict[cty][str(year)]['HU'] = -9999.
                output_dict[cty][str(year)]['RD'] = -9999.
            
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #                       Write yield data to output
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    outfile_name = 'test_out_' + crop_name + '_HU' + str(solar_units) + '_' \
        + output_unit + '.txt'
    
    fmt_str = '{0:15s}{1:<10d}{2:19s}{3:10.1f}{4:10.1f}{5:10.1f}{6:10.1f}' + \
        '{7:10.1f}{8:10.1f}{9:10.1f}{10:10.1f}{11:10.1f}{12:10.1f}' + \
        '{13:10.1f}{14:10.1f}{15:10.1f}{16:10.1f}{17:10.1f}\n'
    
    title_fmt = '{0:15s}{1:10s}{2:20s}{3:10s}{4:10s}{5:10s}{6:10s}' + \
        '{7:10s}{8:10s}{9:10s}{10:10s}{11:10s}{12:10s}{13:10s}{14:10s}' + \
        '{15:10s}{16:10s}{17:10s}\n'
    
    mukeys = obs_dict['MUKEY'].values
    wx_files = obs_dict['WTH'].values
    
    with open(outfile_name, 'w') as fout:
    
        fout.write(title_fmt.format('County','MUKEY','WTH','Yield2001',\
            'Yield2002','Yield2003','Yield2004','Yield2005','Yield2006',\
            'Yield2007','Yield2008','Yield2009','Yield2010','Yield2011',\
            'Yield2012','Yield2013','Yield2014','Yield2015'))
    
        for ii, cty in enumerate(counties):
            # Determine if spaces are in this county name. If there are
            # spaces, replace with underscores. This will be reversed
            # in the reader.
            # ---------------------------------------------------------
            if(' ' in cty):
                out_cty = re.sub(' ','_',cty)
            else:
                out_cty = cty

            fout.write(fmt_str.format(out_cty, mukeys[ii], wx_files[ii], \
                output_dict[cty]['2001']['YLD'], 
                output_dict[cty]['2002']['YLD'], 
                output_dict[cty]['2003']['YLD'], 
                output_dict[cty]['2004']['YLD'], 
                output_dict[cty]['2005']['YLD'], 
                output_dict[cty]['2006']['YLD'], 
                output_dict[cty]['2007']['YLD'], 
                output_dict[cty]['2008']['YLD'], 
                output_dict[cty]['2009']['YLD'], 
                output_dict[cty]['2010']['YLD'], 
                output_dict[cty]['2011']['YLD'], 
                output_dict[cty]['2012']['YLD'], 
                output_dict[cty]['2013']['YLD'], 
                output_dict[cty]['2014']['YLD'], 
                output_dict[cty]['2015']['YLD']))
    
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #                               Clean up
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    os.system('rm -r ./tmp_outzip_dir')

#############################################################################
#
#                            Reading functions
#
#############################################################################
def read_yield_observations(infile): 

    # Read the observations file
    # --------------------------
    obs_dict = pd.read_csv(infile)

    return obs_dict

def read_MSPAM_output(infile):

    data = pd.read_csv(infile, delim_whitespace = True)
    
    # Fix the county names and switch underscores for spaces
    # ------------------------------------------------------
    counties = data['County']
    for ii, cty in enumerate(data['County']):
        if('_' in cty):
            counties[ii] = re.sub('_',' ', cty)

    data['County'] = counties

    return data

#############################################################################
#
#                            Calculating functions
#
#############################################################################

# NOTE: Only set up to work with wheat yield output files now. Can modify
# to work with corn or soybeans later
def calc_model_errors():

    heat_units = np.arange(1000, 2100, 100)
    years = np.arange(2001, 2016)
    
    infile = "../run_batch/yearly_yield_2001_2015.csv"
    obs = read_yield_observations(infile) 
    counties = obs['County'].values
    
    pcnt_errors = np.full( (len(heat_units), len(counties), len(years)), np.nan)
    obs_yields = np.full( (len(heat_units), len(counties), len(years)), np.nan)
    model_yields = np.full( (len(heat_units), len(counties), len(years)), np.nan)
   
    # Loop over the heat units
    # ------------------------ 
    for ii, hu in enumerate(heat_units):

        # Loop over the counties
        # ---------------------- 
        for jj, cty in enumerate(counties):
    
            print(hu, cty)
   
            # Read the MSPAM modeled wheat yields 
            # -----------------------------------
            infile = 'test_out_WHET_HU' + str(int(hu)) + '_BPA.txt'
            data = read_MSPAM_output(infile)
           
            # Extract the observed and modeled yields from 2001 - 2015 for 
            # the current county
            # ------------------------------------------------------------ 
            yield_keys = data.keys()[3:]
            yearly_model = np.array([data[data['County'] == cty][tkey].values for tkey in yield_keys]).squeeze()
            yearly_obs   = np.array([obs[obs['County'] == cty][tkey].values for tkey in yield_keys]).squeeze()
           
            # Calculate the percent errors between the modeled and observed
            # yields
            # ------------------------------------------------------------- 
            errors = yearly_obs - yearly_model
            pcnt_err = (np.abs(errors) / yearly_obs) * 100.
   
            # Mask percent errors for years, heat units, or counties 
            # with missing obs or modeled yields
            # ------------------------------------------------------
            pcnt_err = np.where( (yearly_obs < 0) | (yearly_model < 0), np.nan, pcnt_err)
   
            # Insert the data into the arrays
            # ------------------------------- 
            pcnt_errors[ii,jj,:] = pcnt_err
            obs_yields[ii,jj,:]  = yearly_obs    
            model_yields[ii,jj,:]  = yearly_model  
    
    # Mask any invalid or missing values in the arrays
    # ------------------------------------------------    
    pcnt_errors = np.ma.masked_invalid(pcnt_errors)
    obs_yields = np.ma.masked_where( (obs_yields < 0) | (np.isnan(obs_yields)), obs_yields)
    model_yields = np.ma.masked_where( (model_yields < 0) | (np.isnan(model_yields)), model_yields)

    # Combine the observed yields, modeled yields, and percent errors
    # into a dictionary for output
    # ---------------------------------------------------------------
    out_dict = {}
    out_dict['dims'] = ['heat_units','counties','years']
    out_dict['heat_units'] = heat_units
    out_dict['counties'] = counties
    out_dict['years'] = years
    out_dict['obs_yields']   = obs_yields
    out_dict['model_yields'] = model_yields
    out_dict['pcnt_errors']  = pcnt_errors   

    return out_dict 

#############################################################################
#
#                            Plotting functions
#
#############################################################################

# Plot yearly-averaged errors for a single county as a function of
# heat units
# -----------------------------------------------------------------
def plot_MSPAM_errors_yearavg_single(mspam_dict, county, ax = None, \
        ptitle = None, save = False, show_xlabel = True, show_ylabel = True):

    # Extract the desired errors for the given county
    # -----------------------------------------------
    cty_idx = np.where(mspam_dict['counties'] == county)[0][0]
    cty_errs = mspam_dict['pcnt_errors'][:,cty_idx,:]

    # Average the errors across the years
    # -----------------------------------
    avg_errs = np.nanmean(cty_errs, axis = 1)

    # Plot the errors
    # ---------------
    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    ax.plot(mspam_dict['heat_units'], avg_errs, color = 'k')
    if(show_xlabel):
        ax.set_xlabel('Heat units')
    if(show_ylabel):
        ax.set_ylabel('Yield Percent Error')
    if(ptitle is None):
        ptitle = 'M-SPAM Wheat Yield Errors - Averaged 2001-2015\n' + county
    ax.set_title(ptitle)
    
    if(not in_ax):
        fig.tight_layout()

        if(save):
            outname = re.sub(' ','_', 'mspam_errors_' + county + '.png')
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show()

def plot_MSPAM_errors_yearavg_fourpanel(mspam_dict, cty1, cty2, cty3, cty4, \
        save = False):

    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    axs = fig.subplots(2,2, sharex = True, sharey = True)
    ax1 = axs[0,0]
    ax2 = axs[0,1]
    ax3 = axs[1,0]
    ax4 = axs[1,1]
     
    plot_MSPAM_errors_yearavg_single(mspam_dict, cty1, ax = ax1, \
        save = False, ptitle = cty1, show_xlabel = False)
    plot_MSPAM_errors_yearavg_single(mspam_dict, cty2, ax = ax2, \
        save = False, ptitle = cty2, show_xlabel = False, \
        show_ylabel = False)
    plot_MSPAM_errors_yearavg_single(mspam_dict, cty3, ax = ax3, \
        save = False, ptitle = cty3)
    plot_MSPAM_errors_yearavg_single(mspam_dict, cty4, ax = ax4, \
        save = False, ptitle = cty4, show_ylabel = False)

    plt.suptitle('M-SPAM Wheat Yield Errors\nAveraged 2001-2015')

    if(save):
        outname = re.sub(' ','_', 'mspam_errors_fourpanel_' + cty1 + \
            '_' + cty2 + '_' + cty3 + '_' + cty4 + '.png')
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_MSPAM_errors_yearavg_2D(mspam_dict, save = False):
    
    # Calculate the yearly averages of the percent errors
    # ---------------------------------------------------
    yearly_avg_errs = np.nanmean(mspam_dict['pcnt_errors'], axis = 2)

    # Plot the yearly averaged errors
    # -------------------------------
    plt.close('all')
    fig = plt.figure(figsize = (6, 8))
    ax = fig.add_subplot(1,1,1)
    mesh = ax.pcolormesh(mspam_dict['heat_units'], \
        np.arange(len(mspam_dict['counties'])), \
        yearly_avg_errs.T, \
        cmap = 'viridis', vmin = 15, vmax = 75)
    cbar = fig.colorbar(mesh, ax = ax, label = 'Pcnt Error')
    ax.set_xlabel('Heat units')
    ax.set_ylabel('County')
    ax.set_title('M-SPAM Wheat Yield Errors\nAveraged 2001-2015')
    fig.tight_layout()
    if(save):
        outname = 'mspam_errors_yearavg_2d.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()
