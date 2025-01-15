#!/usr/bin/env python

"""


"""

import numpy as np
import sys
import os
import pandas as pd

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#                Set parameters and create temporary directories
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

os.system('mkdir ./tmp_outzip_dir')

crop_type = 11      # wheat
solar_units = 1800  
plant_month = 4     
plant_day = 10
harvest_month = 9
harvest_day = 30

# Loop over the yearly yields file, extracting the mukeys and weather
# files for the counties
obs_dict = pd.read_csv('../run_batch/yearly_yield_2001_2015.csv')
counties = obs_dict['County'].values

#years = [2004]
years = np.arange(2001, 2016)
#counties = ['GRAND FORKS']


output_dict = {}
for cty in counties:
    output_dict[cty] = {}        


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#                        Begin processing loop
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

for year in years:

    #output_dict[year] = {}    

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #                        Prep the almanac.input file
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
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
                          
       
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #                        Prep the cropALNC.dat file ?
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    #                     Select the county and weather files
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    for jj, cty in enumerate(counties):
        output_dict[cty][str(year)] = {}

        mukey   = obs_dict[obs_dict['County'] == cty]['MUKEY'].values[0]
        wx_file = obs_dict[obs_dict['County'] == cty]['WTH'].values[0]

        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #
        #                               Run the model
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        cmnd = '../bin/model_almanac ' + str(mukey) + ' ' + str(wx_file) + ' 1'
        print(cmnd)
        os.system(cmnd)
        
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        #
        #                Parse the output file and extract yield
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        
        # Move the output zip file into a new, burner directory
        # -----------------------------------------------------
        outzip_name = 'almanac_Mu' + str(mukey) + '_Wx' + str(wx_file) + '.tar.gz'
        outzip_path = '../output/' + outzip_name
        
        cmnd = 'mv ' + outzip_path + ' ./tmp_outzip_dir/'
        print(cmnd)
        os.system(cmnd)
        
        # Open the zip file and extract contents
        # --------------------------------------
        cmnd = 'tar -xvzf tmp_outzip_dir/' + outzip_name + ' -C ./tmp_outzip_dir/'
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
        #if(len(flines) > 4115):
            # Grab only the necessary line
            #output_line = flines[4115].strip()
            #out_split = output_line.split()
            print(out_split)

            output_dict[cty][str(year)]['YLD'] = float(out_split[1])
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
        
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#                         Write yield data to output
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

outfile_name = 'test_out.txt'

fmt_str = '{0:15s}{1:<10d}{2:19s}{3:10.1f}{4:10.1f}{5:10.1f}{6:10.1f}{7:10.1f}{8:10.1f}' + \
    '{9:10.1f}{10:10.1f}{11:10.1f}{12:10.1f}{13:10.1f}{14:10.1f}{15:10.1f}' + \
    '{16:10.1f}{17:10.1f}\n'

title_fmt = '{0:15s}{1:10s}{2:19s}{3:10s}{4:10s}{5:10s}{6:10s}{7:10s}{8:10s}' + \
    '{9:10s}{10:10s}{11:10s}{12:10s}{13:10s}{14:10s}{15:10s}' + \
    '{16:10s}{17:10s}\n'

mukeys = obs_dict['MUKEY'].values
wx_files = obs_dict['WTH'].values

with open(outfile_name, 'w') as fout:

    fout.write(title_fmt.format('County','MUKEY','WTH','Yield2001','Yield2002','Yield2003','Yield2004',\
        'Yield2005','Yield2006','Yield2007','Yield2008','Yield2009',\
        'Yield2010','Yield2011','Yield2012','Yield2013','Yield2014','Yield2015'))

    for ii, cty in enumerate(counties):
        fout.write(fmt_str.format(cty, mukeys[ii], wx_files[ii], \
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


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#                               Clean up
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
os.system('rm -r ./tmp_outzip_dir')
