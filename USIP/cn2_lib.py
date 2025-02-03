"""


"""

# Import modules
import numpy as np 
import numpy.ma as ma
import sys
import matplotlib.pyplot as plt
import os
import warnings
from datetime import datetime
import pandas as pd
import subprocess
import random
from scipy import interpolate
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from glob import glob
from adpaa import ADPAA
from readsounding import *
from compare_cn2 import cn2_calc, cn2_calc_thermo, smooth, plot_cn2

sys.path.append('/home/bsorenson')
from sounding_lib import *
from python_lib import *

colors=['tab:cyan','tab:red','black','green','olive','blue','purple','yellow']

file_dict = {
    '2018050505': {
        'model_file':       '2018050505/2018050505/HRRR_2018050505_ANALYSIS_CKN_ARL',
        'thermo_file':      '2018050505/2018050505/original_data/18_05_05_05_11_03.GRAW.tempdiff.1Hz',
        'radio_file':       '2018050505/2018050505/180505_051104_CKN_GRAW.txt',
        'radio_file_orig':  '2018050505/2018050505/original_data/180505_051104_CKN_GRAW.txt',
        'bis_files': ['2018050505/2018050505/180505_000000_BIS_UWYO.txt',\
            '2018050505/2018050505/180505_120000_BIS_UWYO.txt']
    },
    '2019050403': {
        'model_file':       '2019050403/2019050403/HRRR_2019050403_ANALYSIS_MVL_ARL',
        'thermo_file':      '2019050403/2019050403/19_05_04_02_54_49.GRAW.tempdiff.1Hz',
        'radio_file':       '2019050403/2019050403/190504_030000_MVL_GRAW.txt',
        'radio_file_orig':  '2019050403/2019050403/190504_030000_MVL_GRAW.txt',
        'bis_files': ['2019050403/2019050403/190504_000000_BIS_UWYO.txt',\
            '2019050403/2019050403/190504_120000_BIS_UWYO.txt']
    }
}

# infile must be the 1Hz file
def read_temp_diffs(radio_file, thermo_file, save = False):
    if(radio_file.strip().split('/')[-1][:2] == '19'):
        pre2019 = False
    else:
        pre2019 = True
    
    radio = readsounding(radio_file,allGraw=True,keepMissing=True)
    thermo = ADPAA()
    thermo.ReadFile(thermo_file)
    # Mask any missing data in the thermosonde data
    thermo.mask_MVC()
    #k = kts2ms(a['UWIND'])
    #v = kts2ms(a['VWIND'])
    
    # Delete any constant altitude data or descent data from the radiosonde
    # data
    
    ###############################################################################
    #
    # TRANSPLANTED FROM comparison_checker.py
    #
    # Stretch the thermosonde times to match the radiosonde's times. 
    #
    ###############################################################################
    
    # Delete radiosonde data where the altitudes are missing (only at the end of the file)
    radio['Time'] = np.delete(radio['Time'],np.where(radio['ALT']==999999.9999))
    radio['UTCTime'] = np.delete(radio['UTCTime'],np.where(radio['ALT']==999999.9999))
    radio['TEMP'] = np.delete(radio['TEMP'],np.where(radio['ALT']==999999.9999))
    radio['PRESS'] = np.delete(radio['PRESS'],np.where(radio['ALT']==999999.9999))
    radio['UWIND'] = np.delete(radio['UWIND'],np.where(radio['ALT']==999999.9999))
    radio['VWIND'] = np.delete(radio['VWIND'],np.where(radio['ALT']==999999.9999))
    radio['SPD'] = np.delete(radio['SPD'],np.where(radio['ALT']==999999.9999))
    radio['DIR'] = np.delete(radio['DIR'],np.where(radio['ALT']==999999.9999))
    radio['ALT'] = np.delete(radio['ALT'],np.where(radio['ALT']==999999.9999))
    
    if(pre2019 == True):
        # Thermo data
        start_time = 23331.0
        stop_time = 27380.0 
        # Radio data
        radio_ascent_start = 20154.5
        radio_ascent_stop  = 24628.5
        radio_last_contact = 26398.5
         
        radio_last_index = 5964 # burst
    
        thermo_start_index = 4626
        thermo_last_index = 8675
        thermo_final_last_index = 10232  # last index in the file
    else:
        # Thermo data
        start_time = 13152.0
        stop_time = 14220.0 
        # After 2019/05/03
        radio_ascent_start = 13334.7
        radio_ascent_stop  = 14410.7
        radio_last_contact = 14410.7
    
        # After 2019/05/03
        radio_last_index = 1076 # burst
    
        thermo_start_index= 2622
        thermo_last_index = 3690
        thermo_final_last_index = 3918  # last index in the file
    
    # After 2019/05/03
    #start_time = 13152.0
    #stop_time = 14220.0 
    # Before 2019/05/03
    #start_time = 23331.0
    #stop_time = 27380.0 
    # Radio data
    # After 2019/05/03
    #radio_ascent_start = 13334.7
    #radio_ascent_stop  = 14410.7
    #radio_last_contact = 14410.7
    # Before 2019/05/03
    #radio_ascent_start = 20154.5
    #radio_ascent_stop  = 24628.5
    #radio_last_contact = 26398.5
    #radio_last_index = 7734 # last contact
    # After 2019/05/03
    #radio_last_index = 1076 # burst
    # Before 2019/05/03
    #radio_last_index = 5964 # burst
    
    #thermo_last_index = 10226 # last contact 
    # After 2019/05/03
    #thermo_start_index= 2622
    #thermo_last_index = 3690
    #thermo_final_last_index = 3918  # last index in the file
    # Before 2019/05/03
    #thermo_last_index = 8675
    #thermo_final_last_index = 10232  # last index in the file
    
    
    
    #radio_last_contact = 26383.0
    radio_start_time  = radio_ascent_start 
    radio_last_time  = radio_last_contact 
    #radio_diff = (radio_last_time-radio_start_time)/(7054.0-1490.0)
    radio_diff = (radio_last_time-radio_start_time)/(np.where(radio['UTCTime']==radio_last_time)[0][0]-np.where(radio['UTCTime']==radio_start_time)[0][0])
    
    thermo_diff = 1.0/radio_diff
    
    diff = start_time-radio_start_time # first launch 2018 05 05 
    #plotTime = radio['UTCTime']+diff
    
    # Match the launch times between the thermosonde and radiosonde
    plotTime = thermo.data['Time']-diff
    thermo.data['Time'] = thermo.data['Time']-diff
    
    if(pre2019 == True):
        # Before 2019/05/03
        ascent_rtime = radio['UTCTime'][1490:radio_last_index]
        ascent_ralt = radio['ALT'][1490:radio_last_index]
        ascent_ttime = plotTime[thermo_start_index:thermo_last_index]
        ascent_talt = thermo.data['Alt'][thermo_start_index:thermo_last_index]
    else:
        # After 2019/05/03
        ascent_rtime = radio['UTCTime'][0:radio_last_index]
        ascent_ralt = radio['ALT'][0:radio_last_index]
        ascent_ttime = plotTime[thermo_start_index:thermo_last_index]
        ascent_talt = thermo.data['Alt'][thermo_start_index:thermo_last_index]
    
    descent_rtime = radio['UTCTime'][radio_last_index:len(radio['UTCTime'])]
    descent_ralt = radio['ALT'][radio_last_index:len(radio['UTCTime'])]
    descent_ttime = plotTime[thermo_last_index:thermo_final_last_index]
    
    # Find the ratio of the number of radiosonde ascent times to the number of
    # thermosonde ascent times
    ascent_diff_ratio = float(len(ascent_rtime))/float(len(ascent_ttime))
    # Use that ratio to create new thermosonde ascent times that match with the
    # radiosonde's ascent times. Now, the radiosonde and thermosonde ascent times
    # have launches at the same time and bursts at the same time, although the 
    # burst altitudes differ by a few dozen meters.
    ascent_ttime = np.arange(ascent_ttime[0],ascent_rtime[-1]+1,ascent_diff_ratio)
    
    # Unneeded for the 2019/05/03 launch
    
    rplot = np.arange(0,len(ascent_rtime))
    tplot = np.arange(0,len(ascent_ttime))
    
    if(pre2019 == True):
        # Repeat the process for the descent
        descent_rtime = radio['UTCTime'][radio_last_index:len(radio['UTCTime'])]
        descent_ralt = radio['ALT'][radio_last_index:len(radio['UTCTime'])]
        #descent_ttime = plotTime[thermo_last_index:len(thermo.data['Time'])]
        descent_talt = thermo.data['Alt'][thermo_last_index:thermo_final_last_index]
        
        descent_diff_ratio = float(len(descent_rtime))/float(len(descent_ttime)-1.0)
        descent_ttime = np.arange(ascent_ttime[-1]+1,descent_rtime[-1]+1,descent_diff_ratio)
        # Combine the new thermosonde times into a single array
        ad_ttime = np.concatenate([ascent_ttime,descent_ttime])
        ad_talt  = np.concatenate([ascent_talt,descent_talt])
    
        # Before 2019/05/03
        thermo.data['Time'][thermo_start_index:thermo_final_last_index] = ad_ttime
        thermo.data['Alt'][thermo_start_index:thermo_final_last_index] = ad_talt
    else:
        # After 2019/05/03
        thermo.data['Time'][thermo_start_index:thermo_last_index] = ascent_ttime
        thermo.data['Alt'][thermo_start_index:thermo_last_index] = ascent_talt
    
    #thermo.data['Time'][2622:thermo_final_last_index] = ascent_ttime
    #thermo.data['Alt'][2622:thermo_final_last_index] = ascent_talt
    
    thermo.mask_MVC()
    
    ###############################################################################
    #
    # END OF comparison_checker.py TRANSPLANT
    #
    ###############################################################################
    
    # Account for the fact that the 
    ###time_offset = thermo.data['Time'][-1]-radio['UTCTime'][-1]
    #thermo_start_time = 86332.0
    #thermo_stop_time = 86606.0   # these three are for the tethered test
    #radio_start_time  = 83093.0
    #thermo_start_time = 23331.0
    #thermo_stop_time = 27380.0   # these three are for the tethered test
    #radio_start_time  = 20154.0
    
    # These three are for the second full launch
    thermo_start_time = 13152.0
    thermo_stop_time = 14220.0
    radio_start_time = 13334.7
    
    if(pre2019 == True):
        # Before 2019/05/03
        combined_start_time = 20154.5
        combined_stop_time = 24629.5
    else:
        # After 2019/05/03
        combined_start_time = 13334.7
        combined_stop_time = 14410.7
    
    # This part is made obselete by the comparison_checker.py transplant
    ####time_offset = thermo_start_time-radio_start_time-23.0 # tethered test
    ###time_offset = thermo_start_time-radio_start_time
    ####radio['UTCTime']+=time_offset
    ###thermo.data['Time']-=time_offset
    
    # Grab the "matching" times to plot a time series of the data later
    closetime_thermo = np.array([])
    matchalt_thermo = np.array([])
    closetime_radio = np.array([])
    
    # Get rid of data that is outside the ascent time
    #closeindex_radio = np.where((radio['UTCTime']>=86322.0) & (radio['UTCTime']<=86606.0))[0] # tethered test
    closeindex_radio = np.where((radio['UTCTime']>=combined_start_time) & (radio['UTCTime']<=combined_stop_time))[0]
    for key in radio.keys():
        if((key != 'UNITS') & (type(radio[key]) is not str)):
            radio[key] = radio[key][closeindex_radio]
    
    tempdiff = np.array([])
    for time in radio['UTCTime']:
        closetime=thermo.data['Time'][:].flat[np.abs(thermo.data['Time'][:]-time).argmin()]
        closetime_thermo = np.append(closetime_thermo[:],closetime)
        close_index = np.where(thermo.data['Time']==closetime)[0]
        matchalt_thermo = np.append(matchalt_thermo[:],thermo.data['Alt'][close_index])
    #    print radio['ALT'][np.where(radio['UTCTime']==time)],thermo.data['Alt'][close_index]
    #    print thermo.data['TempDiff'][np.where(thermo.data['Time']==closetime)]
        tempdiff = np.append(tempdiff[:],thermo.data['TempDiff'][np.where(thermo.data['Time']==closetime)])
    
    
    thermo_cn2 = dict(radio)
    thermo_cn2 = cn2_calc_thermo(tempdiff,thermo_cn2)
    # Adding the method='thermo' flag causes the radiosonde Cn2 to be several
    # of magnitude higher than without
    sradio = dict(radio)
    sradio = smooth(sradio)
    sradio = cn2_calc(sradio,method='thermo')
    radio = cn2_calc(radio)
    
    # Ignore the missing thermosonde_data
    masked_tempdiff = ma.masked_values(tempdiff,1e6)
    masked_indices = np.where(masked_tempdiff!=ma.masked)
    for key in thermo_cn2:
        if type(thermo_cn2[key]) is np.ndarray:
            thermo_cn2[key] = thermo_cn2[key][masked_indices]
    
    
    # SMOOTHER
    # Put an 11-point smoother on the thermosonde data
    # Declare arrays to hold smoothed data
    avg_t = np.array([])
    avg_u = np.array([])
    avg_v = np.array([])
    mid_alt = np.array([])
    mid_press = np.array([])
    avg_cn2t = np.array([])
    temp_cn2t = np.zeros(11)
    j=0
    # Set the totals equal to the first elements of the t, u, and v. When
    # the for loop started at 0, an extra 0 was showing up at the beginning
    # of the averaged data arrays. Starting the loop at 1 removes the 0, 
    # but requires the first elements of the data arrays to be added
    # to the total before the start.
    
    # Convert thermo_cn2['CN2T'] to logarithmic for smoothing
    thermo_scn2 = dict(thermo_cn2)
    thermo_scn2['LABEL'] = 'Smoothed Thermosonde'
    thermo_scn2['CN2T'] = np.log10(thermo_scn2['CN2T'])
    total_t=thermo_scn2['TEMP'][0]
    total_cn2t=thermo_scn2['CN2T'][0]
    total_u=thermo_scn2['UWIND'][0]
    total_v=thermo_scn2['VWIND'][0]
    # Loop through the t, u, and v data
    for i in range(1, len(thermo_cn2['CN2T'])):
        # If 11 elements have been summed, average the current total and 
        # append the averages to arrays. 
        if(i%11==0):
            avg_t = np.append(avg_t[:],total_t/11)
            avg_u = np.append(avg_u[:],total_u/11)
            avg_v = np.append(avg_v[:],total_v/11)
            mid_alt = np.append(mid_alt[:],thermo_cn2['ALT'][i-6])
            mid_press = np.append(mid_press[:],thermo_cn2['PRESS'][i-6])
            #avg_cn2t = np.append(avg_cn2t[:],total_cn2t/11)
            avg_cn2t = np.append(avg_cn2t[:],np.average(temp_cn2t))
            j=0
            total_t=0
            total_u=0
            total_v=0
        # Add the current data to the totals
        total_t+=thermo_scn2['TEMP'][i]
        #total_cn2t+=thermo_cn2['CN2T'][i]
        temp_cn2t[j] = thermo_scn2['CN2T'][i]
        j+=1
        total_u+=thermo_scn2['UWIND'][i]
        total_v+=thermo_scn2['VWIND'][i]
    
    # REMOVE to prevent resetting the data with the smoothed data
    thermo_scn2['CN2T']=avg_cn2t
    thermo_scn2['ALT']=mid_alt
    
    # Convert thermo_cn2['CN2T'] back to linear
    thermo_scn2['CN2T'] = 10.**(thermo_scn2['CN2T'])
   
    thermo_scn2['tempdiff'] = tempdiff
    thermo_scn2['masked_indices'] = masked_indices
    thermo_scn2['closetime_thermo'] = closetime_thermo
    thermo_scn2['matchalt_thermo'] = matchalt_thermo

    return thermo_scn2

# A module version of sounding_splitter.py
def split_soundings(archive_file):
    start_indices = []
    end_indices = []

    infile = archive_file   
 
    # Get the year from the archive file's name
    year = infile.split('/')[-1].split('_')[1]
    # Read in all the lines from the archive file with readlines
    lines = open(infile,'r').readlines()
    
    total_lines = 0
    with open(infile) as g:
        num_lines = sum(1 for line in g)
        total_lines = num_lines
    
    # Find the indices of the beginning line of each sounding in the large file
    with open(infile) as f:
        for i, line in enumerate(f):
            if(line[0]=='7'):
                start_indices.append(i)
    
    # Find the indices of the end of each sounding
    for i in range(1,len(start_indices)):
        end_indices.append(start_indices[i]-1)
    
    # Add the last index in the file to the end list 
    end_indices.append(total_lines-2)
    
    # Convert the lists into numpy arrays
    start_indices = np.array(start_indices)
    end_indices = np.array(end_indices)
    
    file_dump = '/home/bsorenson/Research/thermosonde/monthly_analysis/sounding_analysis/'+year+'/'
    #file_dump = '/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/sounding_archive/'+year+'/'
    # Loop over all the files and make extra files for each sounding

    for i in range(0, len(start_indices)):
        # Grab the first line of each sounding that contains date and time info
        first_line = lines[start_indices[i]:end_indices[i]][0].split(' ')
        # Get the sounding 
        location = first_line[1]
        str_year  = first_line[-1].strip()[2:]
        int_month = datetime.strptime(first_line[-2],'%b').month
        if(int_month<10):
            str_month = '0'+str(int_month)
        else:
            str_month = str(int_month)
        str_day   = first_line[-3]
        str_time  = first_line[-4][:2]
        filename = file_dump+location+'/soundings/'+str_year+str_month+str_day+'_'+str_time+'0000_'+location+'_UWYO.txt'
        
        outfile = open(filename,'w')
        outfile.writelines(lines[start_indices[i]:end_indices[i]])
        outfile.close()    
        print("Saved file: "+filename)

# date_str: YYYYMM
#    time:  00 or 12
def calc_monthly_climatology(date_str, time, location = 'BIS', METER = False):

    #work_path='/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/'
    work_path='/home/bsorenson/Research/thermosonde/monthly_analysis/'

    dt_date_str = datetime.strptime(date_str + ' ' + time, '%Y%m %H')

    # Build the climatology file filename based on the input information
    # ------------------------------------------------------------------
    #archive_file = work_path + 'sounding_analysis/
    
    # Extract the individual soundings from the archive file
    #infile = sys.argv[2]
    # Grab the month information from the monthly archive file name
    #date = infile.split('/')[-1].split('_')
    #month = date[0]
    #month_no = strptime(month[:3],'%b').tm_mon
    #str_monthno = str(month_no)
    #if(dt_date_str.month < 10):
    #    str_monthno = '0'+str_monthno
    #tmonth = month_name[month_no]
    #year = date[1]
    #location = date[2]
    year = dt_date_str.strftime('%Y')
    file_path = '/home/bsorenson/Research/thermosonde/monthly_analysis/sounding_analysis/'+year+'/'+location+'/soundings/'
    archive_path = '/home/bsorenson/Research/thermosonde/monthly_analysis/sounding_analysis/'+year+'/'+location+'/'
    # If the files are already extracted and located in 
    #   /nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/2017/***/soundings/,
    # ignore this step
    files_there = glob(dt_date_str.strftime(file_path + '%y%m*.txt'))
    if(len(files_there) == 0):
        print("Extracting soundings from archive file")
        month_name = dt_date_str.strftime('%B').lower()
        archive_file = dt_date_str.strftime(work_path + 'sounding_analysis/%Y/' + \
            location + '/soundings/' + month_name + '_%Y_' + location + '_archive.txt')
        split_soundings(archive_file)
    else:
        print("Soundings already extracted from archive file. Proceeding")
    
    # Get the file names from the recently split archive file
    files = glob(dt_date_str.strftime(file_path + '%y%m*_%H0000_' + location + '_UWYO.txt'))
    list.sort(files)

    # Calculate Cn2 for each sounding and calculate the average and standard deviation
    # for each sounding
    avgs = []
    stdevs = []
    cn2_sfcto925 = []
    cn2_925to850 = []
    cn2_850to700 = []
    cn2_700to500 = []
    cn2_500to400 = []
    cn2_400to300 = []
    cn2_300to200 = []
    cn2_200to150 = []
    cn2_150to100 = []
    cn2_100to50 = []
    cn2_50totop = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for f in files:
            data = readsounding(f)
            data = cn2_calc(data) 
        #    log_cn2 = np.log10(data['CN2'])
            cn2_sfcto925.append(np.average(data['CN2'][np.where(data['PRESS']>925)]))
            cn2_925to850.append(np.average(data['CN2'][np.where((data['PRESS']<=925) & (data['PRESS']>850))]))
            cn2_850to700.append(np.average(data['CN2'][np.where((data['PRESS']<=850) & (data['PRESS']>700))]))
            cn2_700to500.append(np.average(data['CN2'][np.where((data['PRESS']<=700) & (data['PRESS']>500))]))
            cn2_500to400.append(np.average(data['CN2'][np.where((data['PRESS']<=500) & (data['PRESS']>400))]))
            cn2_400to300.append(np.average(data['CN2'][np.where((data['PRESS']<=400) & (data['PRESS']>300))]))
            cn2_300to200.append(np.average(data['CN2'][np.where((data['PRESS']<=300) & (data['PRESS']>200))]))
            cn2_200to150.append(np.average(data['CN2'][np.where((data['PRESS']<=200) & (data['PRESS']>150))]))
            cn2_150to100.append(np.average(data['CN2'][np.where((data['PRESS']<=150) & (data['PRESS']>100))]))
            cn2_100to50.append(np.average(data['CN2'][np.where((data['PRESS']<=100) & (data['PRESS']>50))]))
            cn2_50totop.append(np.average(data['CN2'][np.where(data['PRESS']<=50)]))
    
    log_data_list = [np.log10(cn2_sfcto925),np.log10(cn2_925to850),\
        np.log10(cn2_850to700),np.log10(cn2_700to500),\
        np.log10(cn2_500to400),np.log10(cn2_400to300),\
        np.log10(cn2_300to200),np.log10(cn2_200to150),\
        np.log10(cn2_150to100),np.log10(cn2_100to50),\
        np.log10(cn2_50totop)]
    
    #avg_data_list = [np.average(x) for x in log_data_list]
    #stdev_data_list = [np.std(y) for y in log_data_list]
    
    # Get rid of missing values
    for i,s in enumerate(log_data_list):
        log_data_list[i] = log_data_list[i][np.logical_not(np.isnan(log_data_list[i]))]

    return log_data_list

# date_str = 'YYYYmm HH'
def read_synthetic(date_str, location):
    work_path='/home/bsorenson/Research/thermosonde/'
    # Add check for compare_cn2.py, plot_cn2.py, and readsounding.py
    #if (os.path.isfile(work_path+'sounding_splitter.py') is False):
    #    print "ERROR: "+work_path+"sounding_splitter.py not found"
    #    print "       Please ensure the necessary program is in the correct location"
    #    sys.exit(1)
   
    dt_date_str = datetime.strptime(date_str, '%Y%m %H')
     
    # Extract the individual soundings from the archive file
    infile = dt_date_str.strftime(work_path + 'monthly_analysis/sounding_analysis/%Y/' + \
        location + '/') + dt_date_str.strftime('%b').lower() + \
        dt_date_str.strftime('_%Y_' + location + '_archive.txt')

    #infile = sys.argv[2]
    # Grab the month information from the monthly archive file name
    archive_path = dt_date_str.strftime(\
        '/home/bsorenson/Research/thermosonde/monthly_analysis/sounding_analysis/%Y/'+location+'/')
    file_path    = archive_path + 'soundings/'
    # If the files are already extracted and located in 
    #   /nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/2017/***/soundings/,
    # ignore this step
    files_there = subprocess.check_output(dt_date_str.strftime(\
        'ls '+ file_path + '%y%m*.txt'), \
        shell = True).decode('utf-8').strip().split('\n')

    if(len(files_there)==1):
        print("Extracting soundings from archive file")
        os.system(work_path+'sounding_splitter.py '+ archive_path + infile.strip().split('/')[-1])
    else:
        print("Soundings already extracted from archive file. Proceeding")

    # Get the file names from the recently split archive file
    output = sorted(subprocess.check_output(dt_date_str.strftime(\
        'ls ' + file_path + '%y%m*_%H*_' +location+'_UWYO.txt'), \
        shell = True).decode('utf-8').strip().split('\n'))
   
    # Calculate Cn2 for each sounding and calculate the average and standard deviation
    # for each sounding
    avgs = []
    stdevs = []
    cn2_sfcto925 = []
    cn2_925to850 = []
    cn2_850to700 = []
    cn2_700to500 = []
    cn2_500to400 = []
    cn2_400to300 = []
    cn2_300to200 = []
    cn2_200to150 = []
    cn2_150to100 = []
    cn2_100to50 = []
    cn2_50totop = []
    
    temp_sfcto925 = []
    temp_925to850 = []
    temp_850to700 = []
    temp_700to500 = []
    temp_500to400 = []
    temp_400to300 = []
    temp_300to200 = []
    temp_200to150 = []
    temp_150to100 = []
    temp_100to50 = []
    temp_50totop = []
    
    height_sfc = []
    height_925 = []
    height_850 = []
    height_700 = []
    height_500 = []
    height_400 = []
    height_300 = []
    height_200 = []
    height_150 = []
    height_100 = []
    height_50 = []
    height_top = []
    
    sfc_press = []
    top_press = []
    j=0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for f in output:
            data = readsounding(f)
            data = cn2_calc(data) 
        #    log_cn2 = np.log10(data['CN2'])
            cn2_sfcto925.append(np.average(data['CN2'][np.where(data['PRESS']>925)]))
            cn2_925to850.append(np.average(data['CN2'][np.where((data['PRESS']<=925) & (data['PRESS']>850))]))
            cn2_850to700.append(np.average(data['CN2'][np.where((data['PRESS']<=850) & (data['PRESS']>700))]))
            cn2_700to500.append(np.average(data['CN2'][np.where((data['PRESS']<=700) & (data['PRESS']>500))]))
            cn2_500to400.append(np.average(data['CN2'][np.where((data['PRESS']<=500) & (data['PRESS']>400))]))
            cn2_400to300.append(np.average(data['CN2'][np.where((data['PRESS']<=400) & (data['PRESS']>300))]))
            cn2_300to200.append(np.average(data['CN2'][np.where((data['PRESS']<=300) & (data['PRESS']>200))]))
            cn2_200to150.append(np.average(data['CN2'][np.where((data['PRESS']<=200) & (data['PRESS']>150))]))
            cn2_150to100.append(np.average(data['CN2'][np.where((data['PRESS']<=150) & (data['PRESS']>100))]))
            cn2_100to50.append(np.average(data['CN2'][np.where((data['PRESS']<=100) & (data['PRESS']>50))]))
            cn2_50totop.append(np.average(data['CN2'][np.where(data['PRESS']<=50)]))
    
            temp_sfcto925.append(np.average(data['TEMP'][np.where(data['PRESS']>925)]))
            temp_925to850.append(np.average(data['TEMP'][np.where((data['PRESS']<=925) & (data['PRESS']>850))]))
            temp_850to700.append(np.average(data['TEMP'][np.where((data['PRESS']<=850) & (data['PRESS']>700))]))
            temp_700to500.append(np.average(data['TEMP'][np.where((data['PRESS']<=700) & (data['PRESS']>500))]))
            temp_500to400.append(np.average(data['TEMP'][np.where((data['PRESS']<=500) & (data['PRESS']>400))]))
            temp_400to300.append(np.average(data['TEMP'][np.where((data['PRESS']<=400) & (data['PRESS']>300))]))
            temp_300to200.append(np.average(data['TEMP'][np.where((data['PRESS']<=300) & (data['PRESS']>200))]))
            temp_200to150.append(np.average(data['TEMP'][np.where((data['PRESS']<=200) & (data['PRESS']>150))]))
            temp_150to100.append(np.average(data['TEMP'][np.where((data['PRESS']<=150) & (data['PRESS']>100))]))
            temp_100to50.append(np.average(data['TEMP'][np.where((data['PRESS']<=100) & (data['PRESS']>50))]))
            temp_50totop.append(np.average(data['TEMP'][np.where(data['PRESS']<=50)]))
    
            height_sfc.append(data['ALT'][0]) 
            if(len(np.where(data['PRESS']==925)[0]) != 0):
                height_925.append(data['ALT'][np.where(data['PRESS']==925)])
            if(len(np.where(data['PRESS']==850)[0]) != 0):
                height_850.append(data['ALT'][np.where(data['PRESS']==850)])
            if(len(np.where(data['PRESS']==700)[0]) != 0):
                height_700.append(data['ALT'][np.where(data['PRESS']==700)])
            if(len(np.where(data['PRESS']==500)[0]) != 0):
                height_500.append(data['ALT'][np.where(data['PRESS']==500)])
            if(len(np.where(data['PRESS']==400)[0]) != 0):
                height_400.append(data['ALT'][np.where(data['PRESS']==400)])
            if(len(np.where(data['PRESS']==300)[0]) != 0):
                height_300.append(data['ALT'][np.where(data['PRESS']==300)])
            if(len(np.where(data['PRESS']==200)[0]) != 0):
                height_200.append(data['ALT'][np.where(data['PRESS']==200)])
            if(len(np.where(data['PRESS']==150)[0]) != 0):
                height_150.append(data['ALT'][np.where(data['PRESS']==150)])
            if(len(np.where(data['PRESS']==100)[0]) != 0):
                height_100.append(data['ALT'][np.where(data['PRESS']==100)])
            if(len(np.where(data['PRESS']==50)[0]) != 0):
                height_50.append(data['ALT'][np.where(data['PRESS']==50)])
            height_top.append(data['ALT'][-1])
    
            sfc_press.append(data['PRESS'][0])
            top_press.append(data['PRESS'][-1])
    
    log_data_list = [np.log10(cn2_sfcto925),np.log10(cn2_925to850),\
        np.log10(cn2_850to700),np.log10(cn2_700to500),\
        np.log10(cn2_500to400),np.log10(cn2_400to300),\
        np.log10(cn2_300to200),np.log10(cn2_200to150),\
        np.log10(cn2_150to100),np.log10(cn2_100to50),\
        np.log10(cn2_50totop)]
    
    temp_list = [np.array(temp_sfcto925),np.array(temp_925to850),\
        np.array(temp_850to700),np.array(temp_700to500),\
        np.array(temp_500to400),np.array(temp_400to300),\
        np.array(temp_300to200),np.array(temp_200to150),\
        np.array(temp_150to100),np.array(temp_100to50),\
        np.array(temp_50totop)]
    
    avg_heights = [np.average(height_sfc),np.average(height_925),\
        np.average(height_850),np.average(height_700),np.average(height_500),\
        np.average(height_400),np.average(height_300),np.average(height_200),\
        np.average(height_150),np.average(height_100),np.average(height_50),\
        np.average(height_top)]
    #avg_heights = [np.average(height_sfc),np.average(height_top)]
    
    avg_sfc_press = np.average(sfc_press)
    avg_top_press = np.average(top_press)
    
    # Get rid of missing values
    for i,s in enumerate(log_data_list):
        log_data_list[i] = log_data_list[i][np.logical_not(np.isnan(log_data_list[i]))]
    for i,s in enumerate(temp_list):
        temp_list[i] = temp_list[i][np.logical_not(np.isnan(temp_list[i]))]
    
    avg_cn2 = [np.average(x) for x in log_data_list]
    std_cn2 = [np.std(y) for y in log_data_list]
    
    avg_temp = [np.average(z) for z in temp_list]
    std_temp = [np.std(a) for a in temp_list]
    
    press_list = [avg_sfc_press,925.,850.,700.,500.,400.,300.,200.,150.,100.,50.,avg_top_press]
    
    # Create the synthetic pressure and altitude profiles
    interp_alt = [0]*12
    interp_press = [0]*11
    avg_heights = np.array(avg_heights)
    for i in range(1,len(avg_heights)):
        talt = np.arange(avg_heights[i-1],avg_heights[i],5.0)
        tpress = interpolate.interp1d([talt[0],talt[-1]],[press_list[i-1],press_list[i]])
        ipress = tpress(talt)
        interp_press[i-1]=ipress
        interp_alt[i-1]=talt
    
    # Concatenate the pressure and altitude data into single arrays
    interp_alt = interp_alt[:-1]
    for i in range(1,len(interp_alt)):
        interp_alt[i] = np.delete(interp_alt[i],0)
        interp_press[i] = np.delete(interp_press[i],0)
    all_press = np.concatenate(interp_press)
    all_alt = np.concatenate(interp_alt)
    
    # Create synthetic Cn2 and temp datasets
    # s_temp and s_cn2 are lists of lists that will be concatenated later
    s_temp = [0]*11
    s_cn2 = [0]*11
    # Loop over each average and standard deviation and create synthetic cn2 and 
    # temperature data for each layer. After the loop, s_temp and s_cn2 will
    # contain eleven lists containing synthetic data for each layer
            # # # # # # # # # # # # # # 
            # RANDOM WALK STUFF
            # # # # # # # # # # # # # # 
    ###st1 = avg_temp[0]
    ###test = avg_temp[0]
    ###sc1 = avg_cn2[0]
            # # # # # # # # # # # # # # 
            # END RANDOM WALK STUFF
            # # # # # # # # # # # # # # 
    for i in range(0,len(avg_temp)):
        ###s_temp[i] = [0]*len(interp_press[i])
        ###s_cn2[i]  = [0]*len(interp_press[i])
        ###s_temp[i][0] = st1
        ###s_cn2[i][0]  = sc1
    
        # Interpolate the cn2 average values above and below
        cn2_bottom = avg_cn2[i]
        if(i==(len(avg_temp)-1)):
            cn2_top = avg_cn2[i]
        else:
            cn2_top = avg_cn2[i+1]
        interp_cn2 = interpolate.interp1d([interp_press[i][0],interp_press[i][-1]],[cn2_bottom,cn2_top])
        cn2_interp = interp_cn2(interp_press[i])
    
        # Interpolate the cn2 average values above and below
        temp_bottom = avg_temp[i]
        if(i==(len(avg_temp)-1)):
            temp_top = avg_temp[i]
        else:
            temp_top = avg_temp[i+1]
        interp_temp = interpolate.interp1d([interp_press[i][0],interp_press[i][-1]],[temp_bottom,temp_top])
        temp_interp = interp_temp(interp_press[i])
    
        #for j in range(0,len(interp_press[i])):
            # Use the distance from the mean to determine how
            # much to change
            # # # # # # # # # # # # # # 
            # RANDOM WALK STUFF
            # # # # # # # # # # # # # # 
            ###t_distance_factor = -((s_temp[i][j]-avg_temp[i])/std_temp[i])*(1.0/len(interp_press[i]))
            ###c_distance_factor = -((s_cn2[i][j]-avg_cn2[i])/std_cn2[i])*(1.0/len(interp_press[i]))
            ###if(s_temp[i][j]==avg_temp[i]):
            ###    t_distance_factor = 1.
            ###    c_distance_factor = 1.
    
            #### Create the new random value
            ###r_temp = (2*std_temp[i]*random.random())-std_temp[i]
            ###r_cn2 = (2*std_cn2[i]*random.random())-std_cn2[i]
    
            ###if(test==avg_temp[i]): test += r_temp
            ###else: test += r_temp*((test-avg_temp[i])/std_temp[i])
            ###print avg_temp[i], std_temp[i], ((s_temp[i][j]-avg_temp[i])/std_temp[i]), test
    
            ###s_temp[i][j] = s_temp[i][j]+r_temp*t_distance_factor
            ###s_cn2[i][j] = s_cn2[i][j]+r_cn2*c_distance_factor
    
            #print s_temp[i][j],avg_temp[i],std_temp[i],t_distance_factor
    
            # # # # # # # # # # # # # # 
            # END RANDOM WALK STUFF
            # # # # # # # # # # # # # # 
    
        s_temp[i] = [temp_interp[j]+(2*std_temp[i]*random.random())-std_temp[i] for j in range(0,len(interp_press[i]))]
        s_cn2[i] = [cn2_interp[j]+(2*std_cn2[i]*random.random())-std_cn2[i] for j in range(0,len(interp_press[i]))]
            #s_temp[i] = [avg_temp[i]+(2*std_temp[i]*random.random())-std_temp[i] for j in range(0,len(interp_press[i]))]
            #s_cn2[i] = [avg_cn2[i]+(2*std_cn2[i]*random.random())-std_cn2[i] for j in range(0,len(interp_press[i]))]
            ###st1 = s_temp[i][j]
            ###sc1 = s_cn2[i][j]
    
    # Concatenate the lists in s_temp and s_cn2 into single arrays
    synthetic_temp = np.concatenate(s_temp)
    synthetic_cn2 = np.concatenate(s_cn2)
    
    # Convert the synthetic cn2 from logarithmic to linear units
    linear_cn2 = 10.**synthetic_cn2
    
    # Convert pressure to Pascals and temperature to Kelvin
    t_all_press=all_press*100.
    t_synthetic_temp=synthetic_temp+273.15
    # Use the pressure, temp, and cn2 data to find a synthetic CT2 dataset
    synthetic_ct2 = (linear_cn2/(((79.*10.**-8.)*(t_all_press/(t_synthetic_temp**2.)))**2.))**0.5
    
    recalculated_cn2 = cn2_calc_thermo(synthetic_ct2,all_press,synthetic_temp)
    rcn2 = np.log10(recalculated_cn2['CN2T'])
    
    
    # Smooth the scn2, rcn2, and tempdiff data
    smooth_scn2 = np.array([])
    smooth_rcn2 = np.array([])
    smooth_tempdiff = np.array([])
    total_s = synthetic_cn2[0]
    total_r = rcn2[0]
    total_t = synthetic_ct2[0]
    mid_alt = np.array([])
    mid_press = np.array([])
    for i in range(1,len(synthetic_cn2)):
        if(i%11==0.0):
            smooth_scn2 = np.append(smooth_scn2[:],total_s/11.0)
            smooth_rcn2 = np.append(smooth_rcn2[:],total_r/11.0)
            smooth_tempdiff = np.append(smooth_tempdiff[:],total_t/11.0)
            mid_alt = np.append(mid_alt[:],all_alt[i-6])
            mid_press = np.append(mid_press[:],all_press[i-6])
            total_s=0
            total_r=0
            total_t=0

        # Add the current data to the totals
        total_s+=synthetic_cn2[i]
        total_r+=rcn2[i]
        total_t+=synthetic_ct2[i]

    synth_data = {}
    synth_data['date_str'] = date_str
    synth_data['dt_date_str'] = dt_date_str
    synth_data['location'] = location
    synth_data['synthetic_cn2'] = synthetic_cn2
    synth_data['all_alt'] = all_alt
    synth_data['smooth_cn2'] = smooth_scn2
    synth_data['mid_alt'] = mid_alt    
    synth_data['synthetic_tempdiff'] = synthetic_ct2
    synth_data['smooth_tempdiff'] = smooth_tempdiff
   
    return synth_data 

def plot_synthetic(synth_data, pvar, location = 'BIS', plot_lines = 'both', \
        ptitle = None, ax = None, save = True, pcolor = 'red'):
    
    if(isinstance(synth_data, str)):
        synth_data = read_synthetic(synth_data, location)
   
    dt_date_str = synth_data['dt_date_str']
 
    in_ax = True
    if(ax is None):
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    if(pvar == 'cn2'):
        label_add = '$C_n^2$'
        axis_label = '$log_{10}$[$C_n^2$ [$m^{2/3}$]]'
    elif(pvar == 'tempdiff'):
        label_add = 'ΔT'
        axis_label = 'Temperature Difference [K]'

    if((plot_lines == 'raw') | (plot_lines == 'both')):
        ax.plot(synth_data['synthetic_' + pvar],synth_data['all_alt']/1000.,\
            color='black',label = 'Synthetic ' + label_add)
    if((plot_lines == 'smooth') | (plot_lines == 'both')):
        ax.plot(synth_data['smooth_' + pvar],synth_data['mid_alt']/1000.,\
            color=pcolor,label='Smooth synth. ' + label_add)

    if(not in_ax):
        ax.legend()
        ax.set_xlabel(axis_label,fontsize=11)
    #plt.xlim(-20,-13)
        ax.set_ylabel('Altitude [km]',fontsize=11)
    #plt.ylim(0,35)
    if(ptitle is None):
        ptitle = dt_date_str.strftime('Synthetic ' + label_add + \
        ' Profile for ' + synth_data['location'] + ' ' + '%HZ %b %Y')
    ax.set_title(ptitle)

    if(not in_ax):
        if(save):
            newname = dt_date_str.strftime('synthetic_' + pvar + '_%HZ_%b_%Y_interp.png')
            fig.savefig(newname,dpi=300)
            print("Saved image: "+newname)
        else:
            plt.show()


def plot_synthetic_figure(date_str, save = False):
   
    title_dict = {
        '201805 12': 'Synthetic Data from the 12Z May 2018 Climatology',
        '201905 12': 'Synthetic Data from the 12Z May 2019 Climatology'
    }
 
    plt.close('all')
    fig = plt.figure(figsize = (10,5))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    plt.suptitle(title_dict[date_str])
    plot_synthetic(date_str, 'cn2', ax = ax1, ptitle = '')
    plot_synthetic(date_str, 'tempdiff', ax = ax2, ptitle = '')

    plot_subplot_label(ax1, '(a)', fontsize = 11, location = 'upper_upper_left')
    plot_subplot_label(ax2, '(b)', fontsize = 11, location = 'upper_upper_left')
    
    if(save):
        outname = 'synthetic_combined_' + date_str.split()[0] + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_tempdiffs(thermo_scn2, ax = None, save = False): 
    #print(thermo_scn2['NAME'])
    if(thermo_scn2['NAME'][:2] == '19'):
        pre2019 = False
    else:
        pre2019 = True
    ## Plot the thermosonde Cn2 and radiosonde Cn2 on a graph
    #plot_cn2(thermo_scn2,sradio,'save')
    ##plot_cn2()
    ##plt.plot(np.log10(radio['CN2']),radio['ALT']/1000.)
    ##plt.plot(np.log10(sradio['CN2']),sradio['ALT']/1000.)
    ##plt.plot(np.log10(thermo_cn2['CN2T']),thermo_cn2['ALT']/1000.)
    #plt.show()
    
    # Make a plot of the temperature differences versus height
    in_ax = True
    if(ax is None):
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    ax.plot(thermo_scn2['tempdiff'][thermo_scn2['masked_indices']],\
        thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']]/1000.,color='black', \
        label = 'Thermosonde')
    if(pre2019 == True):
        ax.set_xlim(0,0.18)
    else:
        ax.set_xlim(0,0.70)
    ##!#if(pre2019):
    ##!#    ax.set_title('Free-Flight Launch 2018/05/05 05:00 UTC',fontsize=9)
    ##!#else:
    ##!#    ax.set_title('Free-Flight Launch 2019/05/04 03:00 UTC',fontsize=9)
    ax.set_title('Thermosonde Temperature Differences', fontsize = 8)
    ax.set_xlabel('$\Delta$T [K]', fontsize = 9)
    ax.set_ylabel('Altitude [km]', fontsize = 9)
    ax.tick_params(axis = 'both', labelsize=9)
    ylims = ax.get_ylim()
    ax.set_ylim(0,ylims[1])
    ##!#if(pre2019 == True):
    ##!#    ax.set_ylim(0,thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']]/1000. + 2)
    ##!#    #ax.set_ylim(0,30)
    ##!#else:
    ##!#    ax.set_ylim(0,8)
    if(not in_ax):
        if(save):
            filename = radio['NAME']+'_tempdiff_vs_alt_newRange.png'
            plt.savefig(filename,dpi=300)
            print("Saved image",filename)
        else:
            plt.show()
    
   
def plot_vert_tempdiffs(date_str, ax = None, save = False, calc_abs = False,\
        fontsize = None):
 
    in_data = [file_dict[date_str]['radio_file'], \
        file_dict[date_str]['model_file'], file_dict[date_str]['thermo_file']]

    # Use readsounding to read the current sounding file into a dictionary
    radio = readsounding(in_data[0])

    # Read the thermosonde data
    thermo_scn2 = read_temp_diffs(file_dict[date_str]['radio_file_orig'], in_data[2])

    # Plot the radiosonde temperature differences versus the thermosonde
    # temperature differences.
    in_ax = True
    if(ax is None):
        in_ax = False
        fig = plt.figure(figsize = (8,6))
        ax = fig.add_subplot(1,1,1)

    good_indices = np.where(np.diff(radio['TEMP'])!=0.)
    radio_diffs = np.diff(radio['TEMP'])[good_indices]/np.diff(radio['ALT'])[good_indices]
    if(calc_abs == True):
        radio_diffs = abs(radio_diffs)
    #radio_diffs = abs(np.diff(radio['TEMP'])/np.diff(radio['ALT']))
    
    ax.plot(thermo_scn2['closetime_thermo'][thermo_scn2['masked_indices']],\
        thermo_scn2['tempdiff'][thermo_scn2['masked_indices']],label='Thermosonde [horizontal]')
    ax.plot(radio['UTCTime'][:-1][good_indices],radio_diffs,label='Radiosonde [vertical]')
    #ax.set_xticks(fontsize=11)
    #ax.set_yticks(fontsize=11)
    ax.set_xlabel('Seconds [UTC]',fontsize=fontsize)
    if(calc_abs):
        ax.set_ylabel('Absolute Value of\nTemperature Difference [K]',fontsize=fontsize)
    else:
        ax.set_ylabel('Temperature Difference [K]',fontsize=fontsize)
    #plt.ylabel('Temperature Difference [K]',fontsize=12)
    #if(pre2019 == True):
    #    plt.title('Free-Flight Launch 2018/05/05 05:00 UTC',fontsize=12)
    #else:
    #    plt.title('Free-Flight Launch 2019/05/04 03:00 UTC',fontsize=12)
    if(not in_ax):
        ax.legend()

        if(save):
            filename = radio['NAME']+'_tempdiff_radio_vs_thermo.png'
            if(calc_abs == True):
                filename = radio['NAME']+'_tempdiff_radio_vs_thermo_abs.png'
            plt.savefig(filename,dpi=300)
            print("Saved image",filename)
        else:
            plt.show()

def plot_dual_vert_tempdiffs(date_str, save = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    fig = plt.figure(figsize = (9, 4))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    
    plot_vert_tempdiffs(date_str, ax = ax1, save = False, calc_abs = False)
    plot_vert_tempdiffs(date_str, ax = ax2, save = False, calc_abs = True)
    
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d Thermosonde Flight'))
    
    plot_subplot_label(ax1, 'a)')
    plot_subplot_label(ax2, 'b)')
    
    ax1.legend()
    ax2.legend()
    
    fig.tight_layout()
    
    
    if(save):    
        filename = 'tempdiff_radio_vs_thermo_abs_' + date_str + '.png'
        fig.savefig(filename, dpi = 300)
        print("Saved image", filename)
    else:
        plt.show()


#def plot_cn2_figure(in_data, ax = None, raw = False, old = False, \
def plot_cn2_figure(date_str, ax = None, raw = False, old = False, \
        no_3km = True, ymax = None, show_xlabel = True, \
        show_ylabel = True, ptitle = None, save = False):

    in_data = [file_dict[date_str]['radio_file'], \
        file_dict[date_str]['model_file'], file_dict[date_str]['thermo_file']]

    if(isinstance(in_data, str)):
        in_data = [in_data]
 
    method="thermo"
    if(old):
        method = 'model'

    in_ax = True
    if(ax == None):
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    ax.set_ylim(0, ymax)
    leg_loc='upper right'
    # To view cn2 data, axis must be logarithmic
    #plt.xscale('log')
    #ax.set_xlim(-20, -13)
    #ax.set_title('$C_{n}^{2}$ Profile Comparison')
    if(ptitle is None):
        ptitle = 'Estimated $C_{n}^{2}$'
    ax.set_title(ptitle, fontsize = 10)
    if(show_xlabel):
        ax.set_xlabel("log$_{10}$ [$C_{n}^{2}$ [m$^{-2/3}$]]")
    if(show_ylabel):
        ax.set_ylabel('Altitude [km]', fontsize = 9)
    #ax.tick_params(axis = 'both', labelsize=9)

    max_alt = -99.

    for ii, ifile in enumerate(in_data):

        if(ifile.strip().split('.')[-1] == '1Hz'):
            # Read the thermosonde data
            thermo_scn2 = read_temp_diffs(file_dict[date_str]['radio_file_orig'], ifile)

            if(no_3km):
                thermo_scn2['CN2T'] = np.ma.masked_where(thermo_scn2['ALT'] / 1000. < 3., thermo_scn2['CN2T'])
                thermo_scn2['ALT']  = np.ma.masked_where(thermo_scn2['ALT'] / 1000. < 3., thermo_scn2['ALT'])

            ax.plot(np.log10(thermo_scn2['CN2T']),thermo_scn2['ALT']/1000., \
                color = colors[ii], label='Thermosonde') 

        else:
            # Use readsounding to read the current sounding file into a dictionary
            data = readsounding(ifile)


            # Convert wind speeds from knots to meters/second
            #a['UWIND'] = kts2ms(a['UWIND'])
            #a['VWIND'] = kts2ms(a['VWIND'])
            # Apply an 11-point smoother to the T, u, and v data in the sounding
            # if it is over 900 elements long 
            if((len(data['TEMP'])>900) & (raw==False)):
                data = smooth(data) 
                data = cn2_calc(data,method=method)
            else:
                data = cn2_calc(data,method=method)

            if(no_3km):
                data['CN2']  = np.ma.masked_where(data['ALT']/1000. < 3., data['CN2'])
                data['ALT']  = np.ma.masked_where(data['ALT']/1000. < 3., data['ALT'])

            # Plot cn2 with height
            #olabel=a['TITLE'].split(" ")[0]+" "+a['TITLE'].split(" ")[1]+" "+\
            #       a['TITLE'].split(" ")[-1]
            olabel=data['LABEL'].split()[0]
            if(olabel == 'GRAW'):
                olabel = 'RAOB'
            savename = data['NAME']+'_CN2.png'
            ax.plot(np.log10(data['CN2']),(data['ALT']/1000.),color=colors[ii],label=olabel)

            if(np.max(data['ALT']) / 1000. > max_alt):
                max_alt = np.max(data['ALT']) / 1000.

    ax.axhline(3.0, color = 'black', linestyle = ':')

    #ax.legend(loc='upper right', fontsize=10)
    ax.legend(loc='upper right', prop={'size': 8}, framealpha = 1)
    #ax.set_ylim(0,max_alt + 2)
    if(not in_ax):
        if save is True:
            plt.savefig(savename,dpi=300)
            print("Saved image: "+savename)
        else:
            plt.show()
    
##!#    # Plot a scatter plot of thermosonde temp diffs to radiosonde temp diffs
##!#    #fig6 = plt.figure()
##!#    #plt.scatter(tempdiff[masked_indices][:-3],radio_diffs[masked_indices[0][:-3]])
##!#    #plt.xlabel('Thermosonde Temperature Differences [K]')
##!#    #plt.xlabel('Radiosonde Temperature Differences [K]')
##!#    #plt.show()
##!#    
##!#    if(len(sys.argv)==4):
##!#        model = readsounding(sys.argv[3])
##!#        model = cn2_calc(model,method='thermo')
##!#        plot_cn2(sradio,thermo_scn2,model)
##!#    
##!#    sys.exit(1)
##!#    
##!#    plot_cn2(thermo_cn2,radio)
##!#    
##!#    fig2 = plt.figure()
##!#    plt.plot(radio['TEMP'],radio['ALT'])
##!#    plt.show()
##!#    
##!#    
##!#    press, temp, alt, u, v = smooth(a['PRESS'],a['TEMP'],a['ALT'],u,v) 
##!#    cn2 = cn2_calc(press,temp,alt,u,v,"thermo")
##!#    #cn2 = cn2_calc(a['PRESS'],a['TEMP'],a['ALT'],u,v,"thermo")
##!#    #alt = np.delete(a['ALT'],-1)
##!#    alt = np.delete(alt,-1)
##!#    alt = np.delete(alt,0)
##!#    
##!#    # Calculate test thermosonde cn2
##!#    at = np.array([0.004,0.05,0.001,0.003,0.001,0.005,0.001,0.005,0.015,0.003,0.002,0.0015,0.005,0.0078,0.003,0.01,0.08,0.0375,0.01567,0.01009791])
##!#    t = np.tile(at,int(len(alt)/len(at)))
##!#    diff = len(cn2)-len(t)
##!#    i=0
##!#    #while(len(t)!=len(a['PRESS'])):
##!#    while(len(t)!=len(press)):
##!#        t = np.append(t[:],at[i%len(at)])
##!#        i+=1
##!#    #cn2t = cn2_calc_thermo(a['PRESS'],a['TEMP'],t)
##!#    cn2t = cn2_calc_thermo(press,temp,t)
##!#    cn2t = np.delete(cn2t,-1)
##!#    cn2t = np.delete(cn2t,0)
##!#    cn2 = np.delete(cn2, np.argwhere(alt<3000))
##!#    cn2t = np.delete(cn2t,np.argwhere(alt<3000))
##!#    alt = np.delete(alt, np.argwhere(alt<3000))
##!#    generate_plot()
##!#    plt.plot(cn2,alt,color='black')
##!#    plt.plot(cn2t,alt,color='red')
##!#    plt.show()
##!#    #plt.savefig("17_11_17_23_02_35_CKN_corrected_compared.png",dpi=300)

def plot_combined_figure(date_str, show_synth_data = False, save = False):
   
    if(date_str == '2018050505'):
        synth_date = '201805 00'
        ylim = [0, 30]
    else:
        synth_date = '201905 00'
        ylim = [0, 6.5]

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H')
 
    in_data = [file_dict[date_str]['radio_file'], \
        file_dict[date_str]['model_file']]
    
    fig = plt.figure(figsize = (11,4))
    #fig = plt.figure(figsize = (14,4))
    gs = fig.add_gridspec(1,3)
    #gs = fig.add_gridspec(1,3, hspace = 0.3)
    #ax0  = fig.add_subplot(gs[0,0])   # true color    
    ax1  = fig.add_subplot(gs[0,1])   # true color    
    ax2  = fig.add_subplot(gs[0,2])   # true color 
    plot_sounding_figure(in_data, fig = fig, skew = gs[0,0], \
        save = False, color = None, \
        ptitle = dt_date_str.strftime('%Y-%m-%d %H:00 UTC'))
    #plot_sounding_figure(file_dict[date_str]['bis_files'], fig = fig, skew = gs[0,1], \
    #    save = False)
    
    thermo_scn2 = read_temp_diffs(file_dict[date_str]['radio_file_orig'], \
        file_dict[date_str]['thermo_file'])
    if(show_synth_data):
        plot_synthetic(synth_date, 'tempdiff', ax = ax1, \
            plot_lines = 'smooth', pcolor = 'tab:purple') 
    plot_tempdiffs(thermo_scn2, ax = ax1, save = False)
    ax1.set_ylim(ylim)
 
    #ax1_lims = ax1.get_ylim()
    #plot_cn2_figure(in_data, ax = ax2, raw = False, old = False, \
    if(show_synth_data):
        plot_synthetic(synth_date, 'cn2', ax = ax2, plot_lines = 'smooth', \
            pcolor = 'tab:purple') 
    plot_cn2_figure(date_str, ax = ax2, raw = False, old = False, \
            save = False, no_3km = False, ymax = ylim[1])
            #save = False, no_3km = False, ymax = ax1_lims[1])
    #print(thermo_scn2['CN2T'].shape, thermo_scn2['ALT'].shape)
    #ax2.plot(np.log10(thermo_scn2['CN2T']),thermo_scn2['ALT']/1000. ,label='Thermosonde') 
    
    #ax2.plot(np.log10(thermo_scn2['CN2']),(data['ALT']/1000.),color=colors[ii],label=olabel)
    #gs[0,1].set_title('  b)', loc = 'left', weight = 'bold')
    ax1.legend()
    ax1.set_title('Temperature Differences', fontsize = 10)
    ax1.set_title('  b)', loc = 'left', weight = 'bold')
    ax2.set_title('  c)', loc = 'left', weight = 'bold')

    #plt.suptitle(date_str)

    fig.tight_layout()
    if(save):
        if(show_synth_data):
            synth_add = '_synth'
        else:
            synth_add = ''
        outname = 'combined_cn2_'+date_str+'_singlerow' + synth_add + '.png'
        #outname = 'combined_cn2_'+date_str+'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def func(x, a, c):
    return a * np.log(x) + c

def plot_combined_figure_v2(date_str, show_synth_data = False, \
        save = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H')

    thermo_scn2 = read_temp_diffs(file_dict[date_str]['radio_file_orig'], \
        file_dict[date_str]['thermo_file'])
    
    # Make a plot of the temperature differences versus height
    ax = None
    if(int(date_str[:4]) < 2019):
        pre2019 = True
    else:
        pre2019 = False
    in_ax = True
    if(ax is None):
        in_ax = False
        fig = plt.figure(figsize = (11, 4))
        axs = fig.subplots(1, 4, sharey = True)
        ax1 = axs[0]
        ax2 = axs[1]
        ax3 = axs[2]
        ax4 = axs[3]
        ax5 = ax2.twiny()
    
    # Plot wind data
    # --------------
    ax1.plot(thermo_scn2['TEMP'], \
        thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']] / 1000., \
        color = 'tab:blue')
    ax1.set_xlabel('Temperature [$^{o}$C]')
    ax1.set_ylabel('Altitude [km]')
    ax2.plot(thermo_scn2['SPD'], \
        thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']] / 1000., \
        color = 'tab:blue')
    ax2.set_xlabel('Wind Speed [kts]', color = 'tab:blue')
    ax5.scatter(thermo_scn2['DIR'], \
        thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']] / 1000., \
        color = 'tab:orange', s = 3)
    ax5.set_xlabel('Wind Direction [degrees]', color = 'tab:orange')
    
    if(date_str == '2018050505'):
        synth_date = '201805 00'
        ylim = [0, 30]
    else:
        synth_date = '201905 00'
        ylim = [0, 6.5]
    
    # Plot thermosonde temperature differences
    # ----------------------------------------
    ax3.plot(thermo_scn2['tempdiff'][thermo_scn2['masked_indices']],\
        thermo_scn2['matchalt_thermo'][thermo_scn2['masked_indices']]/1000.,color='tab:blue', \
        label = 'Thermosonde')
    if(pre2019 == True):
        ax3.set_xlim(0,0.18)
    else:
        ax3.set_xlim(0,0.18)
    ##!#if(pre2019):
    ##!#    ax.set_title('Free-Flight Launch 2018/05/05 05:00 UTC',fontsize=9)
    ##!#else:
    ##!#    ax.set_title('Free-Flight Launch 2019/05/04 03:00 UTC',fontsize=9)
    #ax3.set_title('Thermosonde Temperature Differences', fontsize = 8)
    ax3.set_xlabel('Thermosonde $\Delta$T [K]', fontsize = 10)
    #ax3.set_ylabel('Altitude [km]', fontsize = 9)
    ax3.tick_params(axis = 'both', labelsize=10)
    
    
    plot_cn2_figure(date_str, ax = ax4, raw = False, old = False, \
            save = False, no_3km = False, ymax = ylim[1], ptitle = '', show_ylabel = False)
    
    ax1.grid(alpha = 0.50)
    ax2.grid(alpha = 0.50)
    ax3.grid(alpha = 0.50)
    ax4.grid(alpha = 0.50)
    
    
    ax1.set_ylim(ylim)
    ax2.set_ylim(ylim)
    ax3.set_ylim(ylim)
    ax4.set_ylim(ylim)

    plot_subplot_label(ax1, '(a)', location = 'middle_right', fontsize = 11)
    plot_subplot_label(ax2, '(b)', location = 'middle_right', fontsize = 11)
    plot_subplot_label(ax3, '(c)', location = 'middle_right', fontsize = 11)
    plot_subplot_label(ax4, '(d)', location = 'middle_right', fontsize = 11)
    
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    
    fig.tight_layout()
    
    if(save):
        if(show_synth_data):
            synth_add = '_synth'
        else:
            synth_add = ''
        outname = 'combined_cn2_v2_'+date_str+'_singlerow' + synth_add + '.png'
        #outname = 'combined_cn2_'+date_str+'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

#def plot_calibration_curves(data, split_val, outname):
def plot_calibration_curves(date_str, save = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')   
 
    # Read the calibration curve data
    # -------------------------------
    data = pd.read_csv(dt_date_str.strftime('calibration_values_%m_%d_%Y.csv'))
    #data2 = pd.read_csv('calibration_values_05_04_2019.csv')

    if(date_str == '20171117'):
        split_val = data['Vrms2'][5]
    else:
        split_val = data['Vrms2'][4]

    xdata1 = data['Vrms2']
    ydata1 = data['Vcorr2']

    xdata2 = data['Vcorr2']
    ydata2 = data['DeltaT']

    xdata1_1 = xdata1[xdata1 <= split_val]
    ydata1_1 = ydata1[xdata1 <= split_val]
    xdata1_2 = xdata1[xdata1 >= split_val]
    ydata1_2 = ydata1[xdata1 >= split_val]

    print('xdata1_1 = ', xdata1_1)

    popt, pcov = curve_fit(func, xdata1_1, ydata1_1)

    print('popt = ', popt)

    z1_2 = np.polyfit(xdata1_2, ydata1_2, 1)
    p1_2 = np.poly1d(z1_2)

    z2 = np.polyfit(xdata2, ydata2, 1)
    p2 = np.poly1d(z2)

    # Calculate R2
    # ------------
    r2_11 = r2_score(ydata1_1, func(xdata1_1, *popt))
    print("R2_11 = ", r2_11)
    r2_12 = r2_score(ydata1_2, p1_2(xdata1_2))
    print("R2_12 = ", r2_12)
    r2_2 = r2_score(ydata2, p2(xdata2))
    print("R2_2 = ", r2_2)
    

    # Plot everything
    # ---------------
    fig = plt.figure(figsize = (9,4))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    ax1.plot(xdata1_1, func(xdata1_1, *popt), linestyle = '-', color = 'tab:olive')
    ax1.text(xdata1_1[2] + 0.5, ydata1_1[2], \
        "V$_{{COR}}$= {0:.3f}*ln(V$_{{RMS}}$) + {1:.3f}\nR$^2$ = {2:}\nV$_{{RMS}}$ < {3:.3f}".format(*popt, \
        np.round(r2_11,3), split_val), \
        color = 'tab:olive', backgroundcolor = 'white', weight = 'bold')
    ##!#out = plot_trend_line(ax1, xdata1_1, ydata1_1, color = 'tab:red', \
    ##!#    linestyle = '-', slope = 'lin-regress', linewidth = 0.5)
    ##!#print(out)
    #out = plot_trend_line(ax1, xdata1_2, ydata1_2, color = 'tab:orange', \
    #    linestyle = '-', slope = 'lin-regress', linewidth = 0.5)
    ax1.plot(xdata1_2,p1_2(xdata1_2), color = 'tab:orange', linestyle = '-')
    if(z1_2[1] < 0):
        ax1.text(xdata1[2], ydata1[10], \
            "V$_{{COR}}$ = {0:.3f}*V$_{{RMS}}$ - {1:.3f}\nR$^2$ = {2}\nV$_{{RMS}}$ >= {3:.3f}".format(\
            z1_2[0], abs(z1_2[1]), np.round(r2_12), split_val), \
            color = 'tab:orange', backgroundcolor = 'white', weight = 'bold')
    else:       
        ax1.text(xdata1[2], ydata1[10], \
            "V$_{{COR}}$ = {0:.3f}*V$_{{RMS}}$ + {1:.3f}\nR$^2$ = {2}\nV$_{{RMS}}$ >= {3:.3f}".format(\
            z1_2[0], abs(z1_2[1]), np.round(r2_12), split_val), \
            color = 'tab:orange', backgroundcolor = 'white', weight = 'bold')
    print(z1_2)
    ax1.plot(xdata1,ydata1, 'o', markersize = 4, color = 'k')
    ax1.set_xlabel('V$_{RMS}$ (V)', fontsize = 12)
    ax1.set_ylabel('V$_{Corrected}$ (V)', fontsize = 12)
    ax1.set_title('Voltage Calibration')
    ax1.grid()

    ax2.plot(xdata2,p2(xdata2), color = 'tab:orange', linestyle = '-')
    #out = plot_trend_line(ax2, xdata2, ydata2, color = 'tab:orange', linestyle = '-', \
    #    slope = 'lin-regress', linewidth = 0.5)
    print(z2)
    if(z2[1] < 0):
        ax2.text(xdata2[2], ydata2[10], \
            "ΔT = {0:.3f}*V$_{{COR}}$ - {1:.3f}\nR$^2$ = {2:}".format(z2[0], abs(z2[1]), np.round(r2_2)), \
            color = 'k', backgroundcolor = 'white', weight = 'bold')
    else:       
        ax2.text(xdata2[2], ydata2[10], \
            "ΔT = {0:.3f}*V$_{{COR}}$ + {1:.3f}\nR$^2$ = {2:}".format(z2, np.round(r2_2)), \
            color = 'k', backgroundcolor = 'white', weight = 'bold')
 
    ax2.plot(xdata2,ydata2, 'o', markersize = 4, color = 'k')
    #ax2.plot(xdata2,p2(xdata2),'r--')
    ax2.set_xlabel('V$_{Corrected}$ (V)', fontsize = 12)
    ax2.set_ylabel('ΔT (K)', fontsize = 12)
    ax2.set_title('Temperature Conversion')
    ax2.grid()

    plot_subplot_label(ax1, '(a)', fontsize = 11)
    plot_subplot_label(ax2, '(b)', fontsize = 11)

    plt.suptitle(dt_date_str.strftime('%Y-%m-%d'))

    fig.tight_layout()

    if(save):
        outname = 'calib_curves_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_calibration_curves_both(save = False):

    lab_floor = 0.23
    site_floor = 0.25
    lab_floor2 = 0.255

    plot_calibration_curves('20171117', save = save)
    plot_calibration_curves('20190504', save = save)
    #plot_calibration_curves(data1, 0.667)
    #plot_calibration_curves(data2, 0.4)

    ##!#xdata1 = data['Vrms2']
    ##!#ydata1 = data['Vcorr2']

    ##!#xdata2 = data['Vcorr2']
    ##!#ydata2 = data['DeltaT']

    ##!## Divide the voltage calibration data
    ##!## -----------------------------------
    ##!#split_val = 0.6
    ##!#xdata1_1 = xdata1[xdata1 < split_val]
    ##!#ydata1_1 = ydata1[xdata1 < split_val]
    ##!#xdata1_2 = xdata1[xdata1 >= split_val]
    ##!#ydata1_2 = ydata1[xdata1 >= split_val]

    ##!#z2 = np.polyfit(xdata2, ydata2, 1)
    ##!#p2 = np.poly1d(z2)

    ##!## Plot everything
    ##!## ---------------
    ##!#fig = plt.figure(figsize = (9,3))
    ##!#ax1 = fig.add_subplot(1,2,1)
    ##!#ax2 = fig.add_subplot(1,2,2)

    ##!#out = plot_trend_line(ax1, xdata1_2, ydata1_2, color = 'tab:red', \
    ##!#    linestyle = '-', slope = 'lin-regress', linewidth = 0.5)
    ##!#print(out)
    ##!#ax1.plot(xdata1,ydata1, 'o', markersize = 4)
    ##!#ax1.set_xlabel('V$_{RMS}$ (V)')
    ##!#ax1.set_ylabel('V$_{Corrected}$ (V)')
    ##!#ax1.set_title('Voltage Calibration')
    ##!#ax1.grid()

    ##!#out = plot_trend_line(ax2, xdata2, ydata2, color = 'tab:red', linestyle = '-', \
    ##!#    slope = 'lin-regress', linewidth = 0.5)
    ##!#print(out)
    ##!#ax2.plot(xdata2,ydata2, 'o', markersize = 4)
    ##!##ax2.plot(xdata2,p2(xdata2),'r--')
    ##!#ax2.set_xlabel('V$_{Corrected}$ (V)')
    ##!#ax2.set_ylabel('ΔT (K)')
    ##!#ax2.set_title('Temperature Conversion')
    ##!#ax2.grid()

    ##!#plot_subplot_label(ax1, '(a)', fontsize = 11)
    ##!#plot_subplot_label(ax2, '(b)', fontsize = 11)

    ##!#fig.tight_layout()

    ##!#plt.show()

def plot_monthly_cn2_boxplots(date_str, time, location = 'BIS', ax = None, \
        save = False):

    dt_date_str = datetime.strptime(date_str + ' ' + time, '%Y%m %H')

    log_data_list = calc_monthly_climatology(date_str, time, \
        location = location, METER = False)

    plot_title=dt_date_str.strftime('Daily %HZ Averaged C$_n^2$ ' + \
        'Profiles: '+location+' %b %Y')
   
    in_ax = True
    if(ax is None): 
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    bx = ax.boxplot(log_data_list,vert=0)
    
    ax.set_yticklabels(['Sfc-925','925-850','850-700','700-500','500-400',\
        '400-300','300-200','200-150','150-100','100-50','50-top'])
    ax.tick_params(axis="y", labelsize = 8, rotation = 45)
    #ax.yticks(rotation=45,fontsize=6)
    ax.set_title(plot_title)
    ax.set_xlabel('$log_{10}$ $C_n^2$')
    ax.set_ylabel('Pressure levels [mb]',fontsize=10)
    ax.set_xlim(-20,-13)
  
    if(not in_ax): 
        if(save): 
            fig.tight_layout()
            #if(METER is True):
            #    image_name = month+'_'+year+'_'+location+'_'+stime+'Z_temp_meter_boxplots.png'
            #else:
            outname = dt_date_str.strftime('monthly_cn2_boxplots_' + location + '%Y%m_%HZ.png')
            fig.savefig(outname,dpi=200)
            print("Saved image: " + outname)
        else:
            plt.show()

def plot_monthly_cn2_boxplots_multi(date_str1, date_str2, date_str3, \
        date_str4, save = False):
    
    dt_date_str1 = datetime.strptime(date_str1, '%Y%m %H')
    dt_date_str2 = datetime.strptime(date_str2, '%Y%m %H')
    dt_date_str3 = datetime.strptime(date_str3, '%Y%m %H')
    dt_date_str4 = datetime.strptime(date_str4, '%Y%m %H')

    fig = plt.figure(figsize = (8, 8))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    plot_monthly_cn2_boxplots(dt_date_str1.strftime('%Y%m'), \
        dt_date_str1.strftime('%H'), location = 'BIS', ax = ax1)
    plot_monthly_cn2_boxplots(dt_date_str2.strftime('%Y%m'), \
        dt_date_str2.strftime('%H'), location = 'BIS', ax = ax2)
    plot_monthly_cn2_boxplots(dt_date_str3.strftime('%Y%m'), \
        dt_date_str3.strftime('%H'), location = 'BIS', ax = ax3)
    plot_monthly_cn2_boxplots(dt_date_str4.strftime('%Y%m'), \
        dt_date_str4.strftime('%H'), location = 'BIS', ax = ax4)
    
    plot_subplot_label(ax1, '(a)', location = 'upper_right')
    plot_subplot_label(ax2, '(b)', location = 'upper_right')
    plot_subplot_label(ax3, '(c)', location = 'upper_right')
    plot_subplot_label(ax4, '(d)', location = 'upper_right')
    
    plt.suptitle('Monthly C$_{n}^{2}$ Estimated Climatology\nBismarck, ND')
    ax1.set_title(dt_date_str1.strftime('%B %Y\nGenerated with %HZ Soundings'), fontsize = 10)
    ax2.set_title(dt_date_str2.strftime('%B %Y\nGenerated with %HZ Soundings'), fontsize = 10)
    ax3.set_title(dt_date_str3.strftime('%B %Y\nGenerated with %HZ Soundings'), fontsize = 10)
    ax4.set_title(dt_date_str4.strftime('%B %Y\nGenerated with %HZ Soundings'), fontsize = 10)
    
    fig.tight_layout()
    if(save):
        outname = 'monthly_cn2_boxplots_multi_201705_201711.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()
