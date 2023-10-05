"""
  NAME:

  PURPOSE:


"""
import os
home_dir = os.environ['HOME']
import sys
import json
from scipy.interpolate import interp1d
sys.path.append(home_dir)
import python_lib
import importlib
from python_lib import circle, plot_trend_line, nearest_gridpoint, \
    aerosol_event_dict, init_proj, plot_lat_circles, plot_figure_text, \
    plot_subplot_label
sys.path.append(home_dir + '/Research/OMI')
from OMILib import *
sys.path.append(home_dir + '/Research/CERES')
from gridCERESLib import *
sys.path.append(home_dir + '/Research/TROPOMI')
from TROPOMI_Lib import *
sys.path.append(home_dir + '/Research/MODIS/obs_smoke_forcing')
from MODISLib import *
sys.path.append(home_dir + '/Research/NSIDC')
from NSIDCLib import *
sys.path.append(home_dir + '/Research/FuLiou')
from FuLiouLib import *

data_dir = home_dir + '/Research/Arctic_compares/comp_data/'
datacrs = ccrs.PlateCarree()
mapcrs = ccrs.NorthPolarStereo()

json_time_database = home_dir + '/Research/Arctic_compares/json_comp_times.txt'
json_file_database = home_dir + '/Research/Arctic_compares/json_comp_files.txt'

case_dates = ['200607240029', # GOOD
              #'200607240208', # GOOD / CERES mismatch
              '200607240347', # GOOD
              '200607240526', # GOOD
              '200607240844', # GOOD
              '200607242155', # GOOD
              '200607242334', # GOOD
              '200607250112', # GOOD
              '200607250251', # GOOD
              '200607250748', # GOOD?
              '200607252238', # GOOD
              '200607260017', # GOOD
              '200607260156', # GOOD
              '200607260335', # GOOD
              '200607260513', # GOOD?
              '200607260831', # GOOD
              '200607262142', # GOOD
              '200607270100', # GOOD
              '200607270239', # GOOD?
              '200607270418', # GOOD?
              '200607270557', # GOOD?
              '200607270736', # GOOD?
              '200607272226', # GOOD
              '200804221841',  # GOOD
              '200804222020',  # GOOD
              '200804222159',  # GOOD
              '201408110046',
              '201408110404',
              '201408112032',
              '201408112211',
              '201408112350',
              #'201408120129',
              #'201408120308',
              #'201408122115',
              #'201408122254',
              '201506271538',
              '201506271717',
              '201506271856',
              #'201507061353',  # shadowed?
              #'201507061532',  # shadowed?
              '201507061711',        
              '201507061850',
              '201507062028',
              '201507062207',
              '201507062346',
              '201507070125',
              #'201507071436',  # shadowed?
              '201507071615',
              '201507071754',
              '201507071933',
              '201507072112',
              '201507080347',
              '201507080526',
              '201507082016',
              #'201507082115', # doesn't work, don't know why
              '201507090113',
              '201507090252',
              '201507090430',
              '201507090609',
              '201507090748',
              '201507091424',
              '201507100017',
              '201507100156',
              '201507100335',
              '201507100514',
              '201507100653',
              '201507101010',
              '201708161504',  # GOOD
              '201708161643',  # GOOD
              '201708161821',  # GOOD
              '201708171408',  # GOOD
              '201708171547',  # GOOD
              '201708171726',  # GOOD
              '201708171905',  # GOOD
              '201708172043',  # GOOD
              '201708181312',  # GOOD
              '201708181451',  # GOOD
              '201708181630',  # GOOD
              '201708181809',  # GOOD
              '201708181948',  # GOOD
              '201708191355',  # GOOD
              '201708191534',  # GOOD
              '201708191713',  # GOOD
              '201807051856',  # GOOD
              '201807052034',  # GOOD
              '201807052213',  # GOOD
              '201908102115',  # GOOD
              '201908102254',  # GOOD
              '201908110033',  # GOOD
              '201908110351',  # GOOD
             ]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Automation
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# This function automates everything for a single OMI date string:
# the data downloading/matching, the first-look comparison-making, the 
# hdf5 file-making, and the data copying to raindrop.
# Flags:
#   - include_tropomi: adds the colocated TROPOMI data to the colocation
#         NOTE: still need to implement this
#   - file package
def automate_all_preprocess(date_str, download = True, images = True, \
        process = True, omi_dtype = 'ltc3', include_tropomi = True, \
        copy_to_raindrop = False, remove_empty_scans = True, \
        minlat = 65., remove_ch2_file = False, remove_large_files = True):
    if(isinstance(date_str, str)):
        date_str = [date_str]

    if(download):
        # If desired, download all data and/or match up the files in the JSON
        # -------------------------------------------------------------------
        auto_all_download(date_str, download = download, rewrite_json = True, \
            include_tropomi = include_tropomi)

        # Reload the json file here
        with open(json_file_database, 'r') as fin:
            aerosol_event_dict = json.load(fin)

    if(omi_dtype == 'shawn'):
        omi_dtype2 = 'shawn'
        shawn_path = home_dir + '/data/OMI/shawn_files/'
    elif(omi_dtype == 'ltc3'):
        omi_dtype = 'shawn'
        omi_dtype2 = 'ltc3'
        shawn_path = home_dir+ '/data/OMI/shawn_files/ltc3/'

    # Load in the JSON file string relations
    # --------------------------------------
    with open(json_time_database, 'r') as fin:
        file_date_dict = json.load(fin)

    # date_str should consist of an array of OMI swath times
    for dstr in date_str:

        dt_date_str = datetime.strptime(dstr, '%Y%m%d%H%M')

        print(dstr)

        # Select the matching MODIS swath time from the time JSON
        # -------------------------------------------------------

        if(images):
            # Plot the first-look 6-panel comparison plot
            # -------------------------------------------
            #modis_date_str = file_date_dict[dstr]['MODIS'][0]
            if(dt_date_str.year > 2020):
                omi_dtype2 = 'control'
            plot_compare_OMI_CERES_MODIS_NSIDC(dstr, 7, \
                omi_dtype = omi_dtype2, minlat = minlat, zoom = True, \
                save = True)

        if(process):

            # Make a data storage directory, if one does not already exist
            # ------------------------------------------------------------
            save_dir = dt_date_str.strftime(data_dir + '%Y%m%d/')
            print(save_dir)

            if(not os.path.exists(save_dir)):
                print('Making ', save_dir)
                os.system('mkdir ' +  save_dir)

            save_dir2 = dt_date_str.strftime(save_dir + '%Y%m%d%H%M/')
            if(not os.path.exists(save_dir2)):
                print('Making ', save_dir2)
                os.system('mkdir ' +  save_dir2)
   
            short_save_dir2 = dt_date_str.strftime('%Y%m%d/%Y%m%d%H%M/')

            if(remove_ch2_file):
                ch2_file = glob(save_dir2 + '*ch2*')
                if(len(ch2_file) != 0):
                    print("Removing old CH2 file")
                    cmnd = 'rm ' + ch2_file[0]
                    print(cmnd)
                    os.system(cmnd)
                else:
                    print("No CH2 file to remove")
 
            # Run the data subsetter codes
            # ----------------------------
            if(dt_date_str.year <= 2020):
                omi_type = omi_dtype
            else:
                omi_type = 'control'
            write_swath_to_HDF5(dstr, omi_type, save_path = save_dir2, \
                minlat = minlat, shawn_path = shawn_path, \
                remove_empty_scans = remove_empty_scans)

            write_CERES_hrly_grid_to_HDF5(file_date_dict[dstr]['CERES'], \
                save_path = save_dir2, minlat = minlat, \
                remove_empty_scans = remove_empty_scans)

            MODIS_date = file_date_dict[dstr]['MODIS'][0]
            write_MODIS_to_HDF5(MODIS_date, channel = 1, swath = True, \
                save_path = save_dir2, minlat = minlat, \
                remove_empty_scans = remove_empty_scans)
            write_MODIS_to_HDF5(MODIS_date, channel = 7, swath = True, \
                save_path = save_dir2, minlat = minlat, \
                remove_empty_scans = remove_empty_scans)

            NSIDC_date = file_date_dict[dstr]['NSIDC'][:8]
            print(NSIDC_date)
            writeNSIDC_to_HDF5(NSIDC_date, save_path = save_dir2, \
                minlat = minlat, remove_empty_scans = remove_empty_scans)

            if(include_tropomi & (dt_date_str.year > 2017)):                    
                generate_TROPOMI_prep_data(dstr, copy_to_raindrop = \
                    copy_to_raindrop, minlat = minlat, \
                    trop_time = file_date_dict[dstr]['TROPOMI'], \
                    remove_large_files = remove_large_files)

            # Finally, gzip the data
            # ---------------------
            os.chdir(data_dir)
            cmnd = dt_date_str.strftime('tar -cvzf ' + \
                'combined_subsets_%Y%m%d%H%M.tar.gz ' + short_save_dir2)
            print(cmnd)
            os.system(cmnd)

            if(copy_to_raindrop):
                # Secure copy the gzipped file to Raindrop
                # ----------------------------------------
                cmnd = dt_date_str.strftime('scp ' + \
                    'combined_subsets_%Y%m%d%H%M.tar.gz ' + \
                    'bsorenson@134.129.222.68:' + \
                    '/home/bsorenson/OMI/arctic_comp/comp_data/')
                #cmnd = dt_date_str.strftime('scp ' + \
                #    'combined_subsets_%Y%m%d%H%M.tar.gz ' + \
                #    'bsorenson@raindrop.atmos.und.edu:' + \
                #    '/home/bsorenson/OMI/arctic_comp/comp_data/')
                print(cmnd)
                os.system(cmnd)

            if(remove_large_files):
                print("Removing MODIS MYD06 and TROPOMI files")
                

def single_wrap_function(date_str, minlat, min_AI, out_time_dict, download, \
        images, process, include_tropomi, new_only, copy_to_raindrop, \
        remove_empty_scans, remove_ch2_file, skiprows = [52]):

    print(date_str)

    good_list = download_identify_OMI_swaths(date_str, \
        minlat = minlat, min_AI = min_AI, remove_bad = True, \
        skiprows = skiprows, screen_SZA = False)

    #final_good_list = final_good_list + good_list

    for ttime in good_list:
        if(new_only and (ttime not in out_time_dict.keys())):
        #if(ttime not in out_time_dict.keys()):
            print("NEW TIME: " + ttime + " not in json database. Run preprocessor")
            automate_all_preprocess(ttime, download = download, \
                images = images, process = process, \
                omi_dtype = 'ltc3', include_tropomi = include_tropomi,\
                copy_to_raindrop = copy_to_raindrop, \
                minlat = minlat, \
                remove_empty_scans = remove_empty_scans, \
                remove_ch2_file = remove_ch2_file)

        else:
            print(ttime + " in json database. Reprocessing")
            automate_all_preprocess(ttime, download = download, \
                images = images, process = process, \
                omi_dtype = 'ltc3', include_tropomi = include_tropomi, \
                copy_to_raindrop = copy_to_raindrop, \
                minlat = minlat, \
                remove_empty_scans = remove_empty_scans, \
                remove_ch2_file = remove_ch2_file)

def entire_wrapper(min_AI = 1.0, minlat = 70., new_only = True, \
        download = True, images = False, process = False, run_list = None, \
        include_tropomi = True, copy_to_raindrop = True, \
        remove_empty_scans = True, remove_ch2_file = False, \
        skiprows = [52]):

    if(home_dir + '/Research/OMI/' not in sys.path):
        sys.path.append(home_dir + '/Research/OMI/')
    from OMILib import plotOMI_single_swath_figure

    # See which of the good list files are in the json database
    # ---------------------------------------------------------
    with open(json_time_database,'r') as fin:
        out_time_dict = json.load(fin)

    # Prep the run list variables
    # ---------------------------
    no_run_list = True
    if(run_list is not None):
        no_run_list = False

        for date_str in run_list:

            single_wrap_function(date_str, minlat, min_AI, out_time_dict, \
                download, images, process, include_tropomi, new_only, \
                copy_to_raindrop, remove_empty_scans, remove_ch2_file, \
                skiprows = skiprows)

    else:
    
        # Open the omi event list
        # -----------------------
        omi_event_file = home_dir + '/Research/OMI/omi_area_event_dates_minlat70.txt'
        fin = open(omi_event_file, 'r')
        flines = fin.readlines()
        fin.close()

        final_good_list = []
        for fline in flines[:]:
            date_str = fline.strip().split()[0]
            ai_val   = fline.strip().split()[1]

            if(no_run_list):
                run_list = [date_str]

            if((no_run_list) | \
               ((not no_run_list) & (date_str in run_list))):

                print(date_str, ai_val)

                single_wrap_function(date_str, minlat, min_AI, out_time_dict, \
                    download, images, process, include_tropomi, new_only,\
                    copy_to_raindrop, remove_empty_scans, remove_ch2_file, \
                    skiprows = skiprows)
                
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Downloading functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

##!## For a given date str (YYYYMMDDHH), downloads the corresponding OMI
##!## H5 and shawn files and determines which Shawn swaths have significant
##!## aerosol (AI_pert > 1.0)
##!#def download_OMI_files(date_str, omi_dtype = 'ltc3', minlat = 70., \
##!#        min_AI = 1.0, remove_bad = False):
##!#
##!#    dt_date_str = datetime.strptime(date_str, '%Y%m%d')
##!#
##!#    h5_path = home_dir + '/data/OMI/H5_files/'
##!#    if(omi_dtype == 'shawn'):
##!#        omi_dtype2 = 'shawn'
##!#        shawn_path = home_dir + '/data/OMI/shawn_files/'
##!#    elif(omi_dtype == 'ltc3'):
##!#        omi_dtype = 'shawn'
##!#        omi_dtype2 = 'ltc3'
##!#        shawn_path = home_dir+ '/data/OMI/shawn_files/ltc3/'
##!#
##!#    # For a date string, determine if H5 and shawn OMI files are downloaded
##!#    # ---------------------------------------------------------------------
##!#    h5_files    = glob(dt_date_str.strftime(h5_path + 'OMI*_%Ym%m%dt*.he5')) 
##!#    shawn_files = glob(dt_date_str.strftime(shawn_path + '%Y%m%d*')) 
##!#
##!#    # = = = = = = = = = = = = = = = = = = 
##!#    # 
##!#    # Download files for each day
##!#    #
##!#    # = = = = = = = = = = = = = = = = = = 
##!#
##!#    # if not, download here
##!#    if(len(h5_files) == 0):
##!#        # download H5 files
##!#        cmnd = dt_date_str.strftime(\
##!#            'scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/'+\
##!#            'H5_files/OMI*_%Ym%m%dt*.he5 ') + h5_path
##!#        print(cmnd)
##!#        os.system(cmnd)
##!#
##!#        h5_files    = glob(dt_date_str.strftime(h5_path + \
##!#            'OMI*_%Ym%m%dt*.he5')) 
##!#    
##!#    if(len(shawn_files) == 0):
##!#        # download shawn files
##!#        cmnd = dt_date_str.strftime(\
##!#            'scp -r bsorenson@raindrop.atmos.und.edu:/Research/OMI/'+\
##!#            'out_files-ltc3/%Y%m%d* ') + shawn_path
##!#        print(cmnd)
##!#        os.system(cmnd)
##!#
##!#        shawn_files = glob(dt_date_str.strftime(shawn_path + '%Y%m%d*')) 
##!#
##!#    # Loop over each of the files, reading in each file
##!#    file_times = [tfile.strip().split('/')[-1] for tfile in shawn_files]
##!#    
##!#    # = = = = = = = = = = = = = = = = = = = = = = = 
##!#    # 
##!#    # Determine which ones have significant aerosol
##!#    #
##!#    # = = = = = = = = = = = = = = = = = = = = = = = 
##!#    good_list = []
##!#    bad_list  = []
##!#    for ttime in file_times:
##!#        print(ttime)
##!#
##!#        OMI_base = readOMI_swath_shawn(ttime, latmin = minlat, \
##!#            shawn_path = shawn_path)
##!#
##!#        ##!## (Plot the file here?)
##!#        ##!#plt.close('all')
##!#        ##!#fig = plt.figure()
##!#        ##!#ax = fig.add_subplot(1,1,1, projection = mapcrs)
##!#        ##!#plotOMI_single_swath(ax, OMI_base, title = None, \
##!#        ##!#    circle_bound = False, gridlines = False, vmax = 3.0)
##!#        ##!##filename = 'omi_single_swath_' + ttime + '.png'
##!#        ##!#plt.show()
##!#
##!#
##!#
##!#        # Determine if the max AI pert. value north of 65 N
##!#        # is greater than 1 (or some threshold)
##!#        #if(np.nanmax(OMI_base['UVAI_pert']) >= min_AI):
##!#        if(num_above_threshold >= 100):
##!#            # If the value is over the threshold, add the YYYYMMDDHHMM
##!#            # date string to a keep list
##!#            good_list.append(ttime)
##!#        else:
##!#            # If not, add the date string to a throw away list
##!#            bad_list.append(ttime)
##!#
##!#    if(remove_bad):
##!#        # At the end, loop over the throw away list and delete those files.
##!#        for tstr in bad_list:
##!#            local_dstr = datetime.strptime(tstr, '%Y%m%d%H%M')
##!#            cmnd1 = local_dstr.strftime('rm ' + shawn_path + '%Y%m%d%H%M')
##!#            cmnd2 = local_dstr.strftime('rm ' + h5_path + \
##!#                'OMI*_%Ym%m%dt%H%M*.he5')
##!#        
##!#            print(cmnd1)
##!#            print(cmnd2)
##!#            os.system(cmnd1)
##!#            os.system(cmnd2)
##!#
##!#    # Return or print the keep list.
##!#    return good_list
##!#    
##!#    # Then here, loop over the keep list and run the preprocessor?
    
def auto_all_download(date_str, download = True, rewrite_json = False, \
        omi_dtype = 'ltc3', include_tropomi = True):
    if(isinstance(date_str, str)):
        date_str = [date_str]

    out_time_dict = {}
    out_file_dict = {}

    time_file_exists = False
    file_file_exists = False

    # if the json file exists, read it in
    # -----------------------------------
    if(os.path.exists(json_time_database)):
        print("JSON time already exists. Reading in")
        if(not rewrite_json):
            time_file_exists = True
        with open(json_time_database,'r') as fin:
            out_time_dict = json.load(fin)

    if(os.path.exists(json_file_database)):
        print("JSON file already exists. Reading in")
        if(not rewrite_json):
            file_file_exists = True
        with open(json_file_database,'r') as fin:
            out_file_dict = json.load(fin)

    if(omi_dtype == 'shawn'):
        omi_dtype2 = 'shawn'
        shawn_path = home_dir + '/data/OMI/shawn_files/'
    elif(omi_dtype == 'ltc3'):
        omi_dtype = 'shawn'
        omi_dtype2 = 'ltc3'
        shawn_path = home_dir+ '/data/OMI/shawn_files/ltc3/'

    for dstr in date_str:
        print(dstr)
        if(dstr not in out_time_dict.keys() or rewrite_json):

            dt_date_str = datetime.strptime(dstr, '%Y%m%d%H%M')

            # NOTE: 04/17/2023 - modifying to read the "control" data here,
            #       rather than the shawn data. The actual UVAI data doesn't
            #       matter here, only the time.
            #OMI_base = readOMI_swath_shawn(dstr, latmin = 65., \
            #    shawn_path = shawn_path)
            OMI_base = readOMI_swath_hdf(dstr, "control", latmin = 65.)

            omi_file = subprocess.check_output('ls ' + \
                home_dir + '/data/OMI/H5_files/' + \
                'OMI-*_' + OMI_base['date'][:4] + 'm' + OMI_base['date'][4:8] + 't' + \
                    OMI_base['date'][8:] + '*.he5', shell = True).decode(\
                'utf-8').strip().split('\n')[0]


            omi_shawn = shawn_path + OMI_base['date']
            #omi_shawn = home_dir + '/data/OMI/shawn_files/' + OMI_base['date']

            if(dstr == '200607270100'):
                min_time = np.min(OMI_base['TIME'][~OMI_base['UVAI'].mask]) - timedelta(minutes = 5)
                #min_time = np.min(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]) - timedelta(minutes = 5)
            else:
                min_time = np.min(OMI_base['TIME'][~OMI_base['UVAI'].mask])
                #min_time = np.min(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask])
            CERES_date_str = min_time.strftime('%Y%m%d%H')
    
            # Look for the CERES file
            # -----------------------
            try:
                ceres_file = subprocess.check_output('ls ' + \
                    home_dir + '/data/CERES/SSF_Level2/Aqua/' + \
                    'CERES_SSF_Aqua-XTRK_Edition4A_Subset_' + \
                    CERES_date_str[:8] + '*', shell = True).decode(\
                    'utf-8').strip().split('\n')[0]
                print(ceres_file)
            except subprocess.CalledProcessError:
                print("CERES file not found. Must download")
                print(" Cycling the loop.")
                continue

            try:
                #modis_date_list, modis_file_list = download_MODIS_swath(CERES_date_str, \
                        #dest_dir = home_dir + '/data/MODIS/Aqua/', download = download)
                download_dict = download_MODIS_swath(CERES_date_str, \
                        dest_dir = modis_dir, download = download)
            except TypeError: 
                print("ERROR: problems determining MODIS swath times " + \
                    "using CERES data")
                print(" Cycling the loop.")
                continue
 
            if(download): 
                download_NSIDC_daily(CERES_date_str[:8])

            nsidc_file = \
                home_dir + '/data/NSIDC/NSIDC0051_SEAICE_PS_N25km_' + \
                CERES_date_str[:8] + '_v2.0.nc'
   
            out_time_dict[dstr] = {}
            out_time_dict[dstr]['CERES'] = CERES_date_str
            out_time_dict[dstr]['MODIS'] = download_dict['modis_date_list']
            out_time_dict[dstr]['NSIDC'] = CERES_date_str[:8]

            # Download the TROPOMI file here (and/or get its name)
            # ----------------------------------------------------
            if(include_tropomi & (dt_date_str.year > 2017)):                    
                trop_name = download_TROPOMI_file(dstr)
                trop_time = datetime.strptime(\
                    trop_name.strip().split('/')[-1].split('_')[2][:16],\
                    '%Ym%m%dt%H%M%S').strftime('%Y%m%d%H%M')
                out_time_dict[dstr]['TROPOMI'] = trop_time

            dt_date_str = datetime.strptime(\
                download_dict['modis_date_list'][0], '%Y%m%d%H%M')
            local_date = dt_date_str.strftime('%Y-%m-%d')
            if(local_date not in out_file_dict.keys()):
                out_file_dict[local_date] = {}
            if('Lat' not in out_file_dict[local_date].keys()):
                out_file_dict[local_date]['Lat'] = [60., 90.]
            if('Lon' not in out_file_dict[local_date].keys()):
                out_file_dict[local_date]['Lon'] = [-180., 180.]
       
            #for ttime, tfile in zip(modis_date_list, file_list):
            for ii in range(len(download_dict['modis_date_list'])):

                ttime = download_dict['modis_date_list'][ii]
                mfile = download_dict['modis_file_list'][ii]
                cfile = download_dict['cldmk_file_list'][ii]               
                my6fl = download_dict['myd06_file_list'][ii]               
 
                dt_date_str = datetime.strptime(ttime, '%Y%m%d%H%M')

                local_time = dt_date_str.strftime('%H%M')

                if(local_time not in out_file_dict[local_date].keys()):
                    out_file_dict[local_date][local_time] = {}

                # Make sure the MODIS files have the path attached
                # ------------------------------------------------
                if(len(mfile.split('/')) == 1):
                    mfile = modis_dir + 'MYD/' + mfile
                if(len(cfile.split('/')) == 1):
                    cfile = modis_dir + 'CLDMSK/' + cfile
                if(len(my6fl.split('/')) == 1):
                    my6fl = modis_dir + 'MYD06/' + my6fl

                out_file_dict[local_date][local_time]['omi']   = omi_file
                if(include_tropomi & (dt_date_str.year > 2017)):                    
                    out_file_dict[local_date][local_time]['tropomi'] = \
                        home_dir + '/data/TROPOMI/' + trop_name
                out_file_dict[local_date][local_time]['ceres'] = ceres_file
                out_file_dict[local_date][local_time]['modis'] = mfile
                out_file_dict[local_date][local_time]['modis_cloud'] = cfile
                out_file_dict[local_date][local_time]['modis_myd06'] = my6fl
                out_file_dict[local_date][local_time]['ceres_time'] = CERES_date_str[8:10]
                out_file_dict[local_date][local_time]['swath']      = \
                    download_dict['modis_date_list']
                if('Lat' not in out_file_dict[local_date][local_time].keys()):
                    out_file_dict[local_date][local_time]['Lat'] = [60., 90.]
                elif( (out_file_dict[local_date]['Lat'] != [60., 90.]) & 
                      (out_file_dict[local_date][local_time]['Lat'] != out_file_dict[local_date]['Lat'])):
                    out_file_dict[local_date][local_time]['Lat'] = \
                      out_file_dict[local_date]['Lat']
                    print("Reset Lat value")
                if('Lon' not in out_file_dict[local_date][local_time].keys()):
                    out_file_dict[local_date][local_time]['Lon'] = [-180., 180.]
                elif( (out_file_dict[local_date]['Lon'] != [-180., 180.]) & 
                      (out_file_dict[local_date][local_time]['Lon'] != out_file_dict[local_date]['Lon'])):
                    out_file_dict[local_date][local_time]['Lon'] = \
                      out_file_dict[local_date]['Lon']
                if('modis_Lat' not in out_file_dict[local_date][local_time].keys()):
                    out_file_dict[local_date][local_time]['modis_Lat'] = [60., 90.]
                elif( (out_file_dict[local_date]['Lat'] != [60., 90.]) & 
                      (out_file_dict[local_date][local_time]['modis_Lat'] != out_file_dict[local_date]['Lat'])):
                    out_file_dict[local_date][local_time]['modis_Lat'] = \
                      out_file_dict[local_date]['Lat']
                if('modis_Lon' not in out_file_dict[local_date][local_time].keys()):
                    out_file_dict[local_date][local_time]['modis_Lon'] = [-180., 180.]
                elif( (out_file_dict[local_date]['Lon'] != [-180., 180.]) & 
                      (out_file_dict[local_date][local_time]['modis_Lon'] != out_file_dict[local_date]['Lon'])):
                    out_file_dict[local_date][local_time]['modis_Lon'] = \
                      out_file_dict[local_date]['Lon']
    

    with open(json_time_database,'w') as fout:
        json.dump(out_time_dict, fout, indent = 4)

    with open(json_file_database,'w') as fout:
        json.dump(out_file_dict, fout, indent = 4, sort_keys = True)

    return out_time_dict, out_file_dict

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Reading functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def read_colocated(date_str, minlat = 70., zoom = True, \
        compare_tropomi = False):
   
    filename =  data_dir + 'colocated_subset_' + date_str + '.hdf5'
    #filename =  data_dir + date_str[:8] + '/colocated_subset_' + date_str + '.hdf5'
    try:
        print(filename)
        data = h5py.File(filename,'r')
    except FileNotFoundError:
        print("WARNING: Unable to find colocation file with OMI dtg.")
        print("  Using the time database to find the MODIS dtg.")
        print("  Files now use the OMI dtg. Using the JSON comp times")
        print('  file to correct')

        with open(json_time_database, 'r') as fin:
            file_date_dict = json.load(fin)

        if(date_str in file_date_dict.keys()):
            date_str2 = file_date_dict[date_str]['MODIS'][0]
            filename =  data_dir + 'colocated_subset_' + date_str2 + '.hdf5'
            #filename =  data_dir + date_str2[:8] + '/colocated_subset_' + date_str2 + '.hdf5'
            print(filename)
            date_str = date_str2
            data = h5py.File(filename,'r')
        
    if('modis_ch2' in data.keys()):
        print("ERROR: Using old version of colocated file that uses MODIS CH2")
        print("       instead of MODIS CH1. Must reprocess file for ",date_str)
        data.close()
        return

    coloc_data = {}
    coloc_data['date_str'] = date_str
    coloc_data['LAT'] = data['omi_lat'][:,:]
    coloc_data['LON'] = data['omi_lon'][:,:]
    coloc_data['OMI_PERT'] = np.ma.masked_invalid(data['omi_uvai_pert'][:,:])
    coloc_data['OMI_PERT'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['OMI_PERT'])
    coloc_data['OMI_RAW']  = np.ma.masked_invalid(data['omi_uvai_raw'][:,:])
    #coloc_data['OMI_RAW'][:,23:53] = -9e9
    coloc_data['OMI_RAW'] = np.ma.masked_where(coloc_data['OMI_RAW'] < -15, \
        coloc_data['OMI_RAW'])
    coloc_data['OMI_RAW'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['OMI_RAW'])
    if('trop_uvai' in data.keys()):
        coloc_data['TROP_AI'] = np.ma.masked_invalid(data['trop_uvai'][:,:])
        coloc_data['TROP_AI'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
            coloc_data['TROP_AI'])
        if(compare_tropomi):
            coloc_data['TROP_AI'] = np.ma.masked_where(\
                coloc_data['OMI_RAW'].mask, coloc_data['TROP_AI'])
        coloc_data['TROP_AI'] = np.ma.masked_where(\
            coloc_data['TROP_AI'] < -900., coloc_data['TROP_AI'])
    if('trop_ssa0' in data.keys()):
        coloc_data['TROP_SSA0'] = np.ma.masked_invalid(data['trop_ssa0'][:,:])
        coloc_data['TROP_SSA0'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
            coloc_data['TROP_SSA0'])
        coloc_data['TROP_SSA0'] = np.ma.masked_where(\
            coloc_data['TROP_SSA0'] == -999., coloc_data['TROP_SSA0'])
    if('trop_ssa1' in data.keys()):
        coloc_data['TROP_SSA1'] = np.ma.masked_invalid(data['trop_ssa1'][:,:])
        coloc_data['TROP_SSA1'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
            coloc_data['TROP_SSA1'])
        coloc_data['TROP_SSA1'] = np.ma.masked_where(\
            coloc_data['TROP_SSA1'] == -999., coloc_data['TROP_SSA1'])
    if('trop_ssa2' in data.keys()):
        coloc_data['TROP_SSA2'] = np.ma.masked_invalid(data['trop_ssa2'][:,:])
        coloc_data['TROP_SSA2'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
            coloc_data['TROP_SSA2'])
        coloc_data['TROP_SSA2'] = np.ma.masked_where(\
            coloc_data['TROP_SSA2'] == -999., coloc_data['TROP_SSA2'])
    coloc_data['MODIS_CH1'] = np.ma.masked_where((\
        data['modis_ch1'][:,:] == -999.) | (data['modis_ch1'][:,:] > 1.0), \
        data['modis_ch1'][:,:])
    coloc_data['MODIS_CH1'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['MODIS_CH1'])
    coloc_data['MODIS_CH7'] = np.ma.masked_where((\
        data['modis_ch7'][:,:] == -999.) | (data['modis_ch7'][:,:] > 1.0), \
        data['modis_ch7'])
    coloc_data['MODIS_CH7'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['MODIS_CH7'])
    if('modis_cod' in data.keys()):
        coloc_data['MODIS_COD'] = np.ma.masked_where((\
            data['modis_cod'][:,:] == -999.) | (data['modis_cod'][:,:] > 1.0), \
            data['modis_cod'])
        coloc_data['MODIS_COD'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
            coloc_data['MODIS_COD'])
    coloc_data['MODIS_CLD'] = np.ma.masked_where(\
        (data['modis_cld'][:,:] < 0.), data['modis_cld'])
    coloc_data['MODIS_CLD'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['MODIS_CLD'])
    coloc_data['NSIDC_ICE'] = np.ma.masked_where((\
        #data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.), \
        data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.) | \
        (data['nsidc_ice'][:,:] < 80.), \
        data['nsidc_ice'])
    coloc_data['NSIDC_ICE'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['NSIDC_ICE'])
    coloc_data['NSIDC_IO_MIX'] = np.ma.masked_where((\
        #data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.), \
        data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] >= 80.) | \
        (data['nsidc_ice'][:,:] <= 0.), \
        data['nsidc_ice'])
    coloc_data['NSIDC_IOMIX'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['NSIDC_IO_MIX'])
    coloc_data['NSIDC_OCEAN'] = np.ma.masked_where((\
        data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 0.), \
        data['nsidc_ice'])
    coloc_data['NSIDC_OCEAN'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['NSIDC_OCEAN'])
    coloc_data['NSIDC_LAND'] = np.ma.masked_where((\
        data['nsidc_ice'][:,:] == -999.) |  (data['nsidc_ice'][:,:] != 254.), \
        data['nsidc_ice'])
    coloc_data['NSIDC_LAND'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['NSIDC_LAND'])
    coloc_data['NSIDC_OTHER'] = np.ma.masked_where((\
        data['nsidc_ice'][:,:] == -999.) |  (data['nsidc_ice'][:,:] < 251.) | \
        (data['nsidc_ice'][:,:] > 253.), \
        data['nsidc_ice'])
    coloc_data['NSIDC_OTHER'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['NSIDC_OTHER'])
    coloc_data['CERES_SWF'] = np.ma.masked_where((\
        data['ceres_swf'][:,:] == -999.) | (data['ceres_swf'][:,:] > 5000.), \
        data['ceres_swf'])
    coloc_data['CERES_SWF'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['CERES_SWF'])
    coloc_data['CERES_LWF'] = np.ma.masked_where((\
        data['ceres_lwf'][:,:] == -999.) | (data['ceres_lwf'][:,:] > 5000.), \
        data['ceres_lwf'])
    coloc_data['CERES_LWF'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['CERES_LWF'])
    coloc_data['SZA'] = data['omi_sza'][:,:]
    #coloc_data['VZA'] = data['omi_vza'][:,:]
    #coloc_data['AZM'] = data['omi_azm'][:,:]
     
    data.close()
 
    return coloc_data

def read_colocated_combined(date_str, zoom = True):

    if(len(date_str) == 6):
        fmt = '%Y%m'
    elif(len(date_str) == 8):
        fmt = '%Y%m%d'
    elif(date_str == 'all'):
        fmt = 'NONE' 
   
    dt_date_str = datetime.strptime(date_str, fmt)

    # Locate all matching files
    # -------------------------
    files = subprocess.check_output(dt_date_str.strftime('ls ' + data_dir + 'colocated*' + fmt + \
        '*.hdf5'), shell = True).decode('utf-8').strip().split('\n')

    for ii, tfile in enumerate(files):
        ttime = tfile.strip().split('/')[-1][-17:-5] 
        print(ttime)
        coloc_data = read_colocated(ttime)

        if(ii == 0):
            temp_data = coloc_data
        else:
            for key in coloc_data.keys():
                if(key != 'date_str'):
                    temp_data[key] = np.concatenate((temp_data[key], coloc_data[key]))
        #test_data = np.concatenate((test_data, temp_data['omi_uvai_pert'][:,:]))
        #print(tfile, temp_data['omi_uvai_pert'].shape) 
        #temp_data.close()

    temp_data['date_str'] = date_str

    total_data = {}
    total_data['date_str'] = date_str

    total_data['OMI'] = np.ma.masked_invalid(temp_data['OMI'][:,:])
    total_data['MODIS_CH1'] = np.ma.masked_where((\
        temp_data['MODIS_CH1'][:,:] == -999.) | (temp_data['MODIS_CH1'][:,:] > 1.0), \
        temp_data['MODIS_CH1'][:,:])
    total_data['MODIS_CH7'] = np.ma.masked_where((\
        temp_data['MODIS_CH7'][:,:] == -999.) | (temp_data['MODIS_CH7'][:,:] > 1.0), \
        temp_data['MODIS_CH7'])
    total_data['NSIDC_ICE'] = np.ma.masked_where((\
        #data['NSIDC_ICE'][:,:] == -999.) | (data['NSIDC_ICE'][:,:] > 100.), \
        temp_data['NSIDC_ICE'][:,:] == -999.) | (temp_data['NSIDC_ICE'][:,:] > 100.) | \
        (temp_data['NSIDC_ICE'][:,:] < 80.), \
        temp_data['NSIDC_ICE'])
    total_data['NSIDC_OCEAN'] = np.ma.masked_where((\
        temp_data['NSIDC_ICE'][:,:] == -999.) | (temp_data['NSIDC_ICE'][:,:] > 0.), \
        temp_data['NSIDC_ICE'])
    total_data['NSIDC_LAND'] = np.ma.masked_where((\
        temp_data['NSIDC_ICE'][:,:] == -999.) |  (temp_data['NSIDC_ICE'][:,:] != 254.), \
        temp_data['NSIDC_ICE'])
    total_data['CERES_SWF'] = np.ma.masked_where((\
        temp_data['CERES_SWF'][:,:] == -999.) | (temp_data['CERES_SWF'][:,:] > 5000.), \
        temp_data['CERES_SWF'])
    total_data['CERES_LWF'] = np.ma.masked_where((\
        temp_data['CERES_LWF'][:,:] == -999.) | (temp_data['CERES_LWF'][:,:] > 5000.), \
        temp_data['CERES_LWF'])

    return total_data
 
    #print(test_data.shape)    

def read_comp_grid_climo(filename):

    data = h5py.File(filename,'r')

    comp_grid_data = {}
    comp_grid_data['swf_climo'] = np.ma.masked_where(\
        (data['ceres_swf_count'][:] == -9) | \
        (data['ceres_swf_climo'][:] > 2e3), data['ceres_swf_climo'][:])
    comp_grid_data['swf_count'] = np.ma.masked_where(\
        (data['ceres_swf_count'][:] == -9) | \
        (data['ceres_swf_climo'][:] > 2e3), data['ceres_swf_count'][:])
   

 
    comp_grid_data['ai_bins'] = data['omi_ai_bins'][:]
    comp_grid_data['ai_edges'] = data['omi_ai_edges'][:]
    comp_grid_data['sza_bins'] = data['omi_sza_bins'][:]
    comp_grid_data['sza_edges'] = data['omi_sza_edges'][:]
    comp_grid_data['ice_bins'] = data['nsidc_ice_bins'][:]
    comp_grid_data['ice_edges'] = data['nsidc_ice_edges'][:]
    if('modis_cod_bins' in data.keys()):
        comp_grid_data['cod_bins'] = data['modis_cod_bins'][:]
        comp_grid_data['cod_edges'] = data['modis_cod_edges'][:]
    else:
        comp_grid_data['ch7_bins'] = data['modis_ch7_bins'][:]
        comp_grid_data['ch7_edges'] = data['modis_ch7_edges'][:]
    #comp_grid_data['cld_bins'] = data['ceres_cld_bins'][:]
    #comp_grid_data['cld_edges'] = data['ceres_cld_edges'][:]

    data.close()

    return comp_grid_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Processing functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def select_category(coloc_data, var, cat):

    ch7_lim = 0.05

    if(cat == 'ICE'):
        return np.ma.masked_where((coloc_data['NSIDC_ICE'].mask ), coloc_data[var])
    elif(cat == 'OCEAN'):
        return np.ma.masked_where((coloc_data['NSIDC_OCEAN'].mask ), coloc_data[var])
    elif(cat == 'LAND'):
        return np.ma.masked_where((coloc_data['NSIDC_LAND'].mask ), coloc_data[var])
    elif(cat == 'ICE_CLEAR'):
        return np.ma.masked_where((coloc_data['NSIDC_ICE'].mask) | \
            (coloc_data['MODIS_CH7'] > ch7_lim), coloc_data[var])
    elif(cat == 'ICE_CLOUD'):
        return np.ma.masked_where((coloc_data['NSIDC_ICE'].mask) | \
            (coloc_data['MODIS_CH7'] < ch7_lim), coloc_data[var])
    elif(cat == 'OCEAN_CLEAR'):
        return np.ma.masked_where((coloc_data['NSIDC_OCEAN'].mask) | \
            (coloc_data['MODIS_CH7'] > ch7_lim), coloc_data[var])
    elif(cat == 'OCEAN_CLOUD'):
        return np.ma.masked_where((coloc_data['NSIDC_OCEAN'].mask) | \
            (coloc_data['MODIS_CH7'] < ch7_lim), coloc_data[var])
    elif(cat == 'LAND_CLEAR'):
        return np.ma.masked_where((coloc_data['NSIDC_LAND'].mask) | \
            (coloc_data['MODIS_CH7'] > ch7_lim), coloc_data[var])
    elif(cat == 'LAND_CLOUD'):
        return np.ma.masked_where((coloc_data['NSIDC_LAND'].mask) | \
            (coloc_data['MODIS_CH7'] < ch7_lim), coloc_data[var])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Compare the OMI, CERES, and MODIS data over the Arctic
# ------------------------------------------------------
def plot_compare_OMI_CERES_MODIS_NSIDC(date_str, ch1, \
        omi_dtype = 'shawn', minlat = 65., zoom = False, save = False):

    if(home_dir + '/Research/OMI/' not in sys.path):
        sys.path.append(home_dir + '/Research/OMI/')
    if(home_dir + '/Research/CERES/' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES/')
    if(home_dir + '/Research/MODIS/obs_smoke_forcing/' not in sys.path):
        sys.path.append(home_dir + '/Research/MODIS/obs_smoke_forcing/')
    if(home_dir + '/Research/NSIDC/' not in sys.path):
        sys.path.append(home_dir + '/Research/NSIDC/')

    from OMILib import plotOMI_single_swath_figure
    from gridCERESLib import plotCERES_hrly_figure
    from MODISLib import plot_MODIS_channel
    from NSIDCLib import plotNSIDC_daily_figure

    if(omi_dtype == 'shawn'):
        shawn_path = home_dir + '/data/OMI/shawn_files/'
    elif(omi_dtype == 'ltc3'):
        omi_dtype = 'shawn'
        shawn_path = home_dir + '/data/OMI/shawn_files/ltc3/'
        print('In plotter, looking for ltc3')
    else:
        shawn_path = home_dir + '/data/OMI/shawn_files/ltc3/'

    # Load in the JSON file string relations
    # --------------------------------------
    with open(json_time_database, 'r') as fin:
        file_date_dict = json.load(fin)

    if((date_str[:8] == '20080422') | \
       (date_str[:8] == '20180705')\
        ):
        size = (14, 5)
    elif(date_str[:8] == '20190811'\
        ):
        size = (12, 6.5)
    elif(date_str[:7] == '2006072'\
        ):
        size = (12, 8)
    else:
        size = (10, 9)
    ##!#file_date_dict = {
    ##!#    '200804221935': {
    ##!#        'size': (14,5), 
    ##!#        'OMI': '200804221841', 
    ##!#        'CERES': '2008042219', 
    ##!#    },
    ##!#    '200804222110': {
    ##!#        'size': (14,5), 
    ##!#        'OMI': '200804222020', 
    ##!#        'CERES': '2008042221', 
    ##!#    },
    ##!#    '200804222250': {
    ##!#        'size': (14,5), 
    ##!#        'OMI': '200804222159', 
    ##!#        'CERES': '2008042222', 
    ##!#    },
    ##!#    '201807051950': {
    ##!#        'size': (14,5), 
    ##!#        'OMI': '201807051856', 
    ##!#        'CERES': '2018070519', 
    ##!#    },
    ##!#    '201807052125': {
    ##!#        'size': (14,5), 
    ##!#        'OMI': '201807052034', 
    ##!#        'CERES': '2018070521', 
    ##!#    },
    ##!#    '201807052305': {
    ##!#        'size': (14,5), 
    ##!#        'OMI': '201807052213',
    ##!#        'CERES': '2018070523', 
    ##!#    },
    ##!#    '201808241435': {
    ##!#        'OMI': '201808241343', 
    ##!#        'CERES': '2018082414', 
    ##!#    },
    ##!#    '201908110125': {
    ##!#        'size': (12,6.5), 
    ##!#        'OMI': '201908110033', 
    ##!#        'CERES': '2019081101', 
    ##!#    },
    ##!#    '201908110440': {
    ##!#        'size': (12,6.5), 
    ##!#        'OMI': '201908110351', 
    ##!#        'CERES': '2019081104', 
    ##!#    },
    ##!#}

    # Plot the MODIS true-color and channel data
    # ------------------------------------------
    modis_date = file_date_dict[date_str]['MODIS'][0]
    ceres_date = file_date_dict[date_str]['CERES']
    nsidc_date = file_date_dict[date_str]['NSIDC'][:8]
    omi_date   = date_str

    # Read in data for value printing
    # -------------------------------
    MODIS_ch7 = read_MODIS_channel(modis_date, 7, swath = True)
    CERES_data = readgridCERES_hrly_grid(ceres_date, 'SWF', minlat = 65. )

    # Make the overall figure
    # -----------------------
    plt.close('all')
    #fig1 = plt.figure(figsize = file_date_dict[modis_date_str]['size'])
    fig1 = plt.figure(figsize = size)
    ax1 = fig1.add_subplot(2,3,1, projection = mapcrs)
    ax2 = fig1.add_subplot(2,3,2, projection = mapcrs)
    ax3 = fig1.add_subplot(2,3,3, projection = mapcrs)
    ax4 = fig1.add_subplot(2,3,4, projection = mapcrs)
    ax5 = fig1.add_subplot(2,3,5, projection = mapcrs)
    ax6 = fig1.add_subplot(2,3,6, projection = mapcrs)

    def mouse_event(event):
        ix, iy = event.xdata, event.ydata
        event_list = [ix, iy]
        # convert from display coordinates to data coordinates
        p_a_cart = datacrs.transform_point(event_list[0], event_list[1], src_crs=mapcrs)

        match_idx = nearest_gridpoint(p_a_cart[1], p_a_cart[0], \
            MODIS_ch7['lat'], MODIS_ch7['lon'])
        c_match_idx = nearest_gridpoint(p_a_cart[1], p_a_cart[0], \
            CERES_data['lat'], CERES_data['lon'])

        print('x: {} and y: {}'.format(event.xdata, event.ydata))
        print('lon: {} and lat: {}'.format(p_a_cart[0], p_a_cart[1]))
        print('MODIS ch7 = ', MODIS_ch7['data'][match_idx])
        print('CERES alb = ', CERES_data['alb'][c_match_idx])
        print('CERES SWF = ', CERES_data['swf'][c_match_idx])
        print('CERES sza = ', CERES_data['sza'][c_match_idx])
        print('CERES vza = ', CERES_data['vza'][c_match_idx])

    cid = fig1.canvas.mpl_connect('button_press_event', mouse_event)
   
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = zoom, ax = ax1)
    plot_MODIS_channel(modis_date, ch1, swath = True, \
        zoom = zoom, ax = ax2, vmax = 0.4)
    #plot_MODIS_channel(modis_date_str, ch1, swath = True, \
    #    zoom = zoom, ax = ax3)

    # Plot the NSIDC data
    # -------------------
    plotNSIDC_daily_figure(nsidc_date, minlat = minlat, \
        zoom = zoom, ax = ax3, gridlines = False, save = False)

    # Plot the OMI data
    # -----------------
    plotOMI_single_swath_figure(omi_date, \
            dtype = omi_dtype, only_sea_ice = False, minlat = minlat, \
            ax = ax4, skiprows = [52], lat_circles = None, save = False, \
            zoom = zoom, shawn_path = shawn_path)
    
    # Plot the CERES data
    # -------------------
    plotCERES_hrly_figure(ceres_date, 'SWF',  \
        minlat = minlat, lat_circles = None, ax = ax5, title = 'SWF',\
        grid_data = True, zoom = zoom, vmax = 450, vmin = None, save = False)
    plotCERES_hrly_figure(ceres_date, 'ALB',  \
        minlat = minlat, lat_circles = None, ax = ax6, title = 'ALB',\
        grid_data = True, zoom = zoom, vmax = None, vmin = None, save = False)

    fig1.tight_layout()

    if(save):
        outname = 'omi_ceres_modis_nsidc_compare_' + omi_date + '.png'
        fig1.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_spatial(ax, lon, lat, data, date_str, cmap = 'Greys_r', \
        plabel = '', colorbar = True, zoom = False, vmin = None, \
        vmax = None, minlat = 65):
  
    if(vmin is not None):
        data = np.ma.masked_where(data < xmin, data) 
    if(vmax is not None):
        data = np.ma.masked_where(data > vmax, data) 
     
    mesh = ax.pcolormesh(lon, lat, data, \
                   transform = datacrs, shading = 'auto', cmap = cmap)
    ax.coastlines()

    if(colorbar):
        cbar = plt.colorbar(mesh,\
            ax = ax, orientation='vertical',shrink = 0.8, extend = 'both')
        #cbar10.set_label('UV Aerosol Index',weight='bold',fontsize=colorbar_label_size)
        #cbar.ax.tick_params(labelsize=14)
        cbar.set_label(plabel, weight = 'bold')
    
    if(zoom):
        #ax.set_extent([-180, 180, 60, 90], datacrs)
        ax.set_extent([aerosol_event_dict[date_str]['Lon'][0], \
                       aerosol_event_dict[date_str]['Lon'][1], \
                       aerosol_event_dict[date_str]['Lat'][0], \
                       aerosol_event_dict[date_str]['Lat'][1]],\
                       datacrs)
    else:
        ax.set_extent([-180, 180, minlat, 90], datacrs)

def plot_compare_colocate_spatial(coloc_data, minlat = 65., zoom = False, \
        save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data)
    else:
        dt_date_str = datetime.strptime(coloc_data['date_str'], '%Y%m%d%H%M')

    ## Read the colocated data
    ## -----------------------
    #data = h5py.File('colocated_subset_' + date_str + '.hdf5','r')

    # Make the overall figure
    # -----------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (10,8))
    ax1 = fig1.add_subplot(2,3,1, projection = mapcrs)
    ax2 = fig1.add_subplot(2,3,2, projection = mapcrs)
    ax3 = fig1.add_subplot(2,3,3, projection = mapcrs)
    ax4 = fig1.add_subplot(2,3,4, projection = mapcrs)
    ax5 = fig1.add_subplot(2,3,5, projection = mapcrs)
    ax6 = fig1.add_subplot(2,3,6, projection = mapcrs)
  
    pdate = dt_date_str.strftime('%Y-%m-%d') 
    # Plot the MODIS true-color and channel data
    # ------------------------------------------
    plot_spatial(ax1, coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH1'], \
        pdate, cmap = 'Greys_r', zoom = zoom)
    ##!#ax1.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH1'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'Greys_r')
    ##!#ax1.coastlines()
    ##!#ax1.set_extent([-180, 180, 60, 90], datacrs)
    plot_spatial(ax2, coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH7'], \
        pdate, cmap = 'Greys_r', zoom = zoom)
    ##!#ax2.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH7'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'Greys_r')
    ##!#ax2.coastlines()
    ##!#ax2.set_extent([-180, 180, 60, 90], datacrs)

    # Plot the NSIDC coloc_data
    # -------------------
    plot_spatial(ax3, coloc_data['LON'], coloc_data['LAT'], coloc_data['NSIDC_ICE'], \
        pdate, cmap = 'ocean', zoom = zoom)
    ##!#ax3.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['NSIDC'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'ocean')
    ##!#ax3.coastlines()
    ##!#ax3.set_extent([-180, 180, 60, 90], datacrs)

    # Plot the OMI coloc_data
    # -----------------
    plot_spatial(ax4, coloc_data['LON'], coloc_data['LAT'], coloc_data['OMI_RAW'], \
        pdate, cmap = 'jet', zoom = zoom)
    ##!#ax4.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['OMI'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'jet')
    ##!#ax4.coastlines()
    ##!#ax4.set_extent([-180, 180, 60, 90], datacrs)
    
    # Plot the CERES coloc_data
    # -------------------
    plot_spatial(ax5, coloc_data['LON'], coloc_data['LAT'], coloc_data['CERES_SWF'], \
        pdate, cmap = 'plasma', zoom = zoom)
    ##!#ax5.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['CERES_SWF'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'plasma')
    ##!#ax5.coastlines()
    ##!#ax5.set_extent([-180, 180, 60, 90], datacrs)
    plot_spatial(ax6, coloc_data['LON'], coloc_data['LAT'], coloc_data['CERES_LWF'], \
        pdate, cmap = 'plasma', zoom = zoom)
    ##!#ax6.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['CERES_LWF'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'plasma')
    ##!#ax6.coastlines()
    ##!#ax6.set_extent([-180, 180, 60, 90], datacrs)


    if(save):
        outname = 'arctic_compare_spatial_' + \
            dt_date_str.strftime('%Y%m%d%H%M') + '.png'
        fig1.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_compare_scatter(coloc_data, var1, var2, var3 = None, minlat = 65., \
        xmin = None, xmax = None, ymin = None, ymax = None, \
        colorbar = True, trend = False, zoom = False, save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data)
    else:
        dt_date_str = datetime.strptime(coloc_data['date_str'], '%Y%m%d%H%M')

    if(xmin is None):
        xmin = np.nanmin(coloc_data[var1])
    if(xmax is None):
        xmax = np.nanmax(coloc_data[var1])
    if(ymin is None):
        ymin = np.nanmin(coloc_data[var2])
    if(ymax is None):
        ymax = np.nanmax(coloc_data[var2])
    
    in_data = np.where((coloc_data[var1] > xmin) & \
        (coloc_data[var1] < xmax) & (coloc_data[var2] > ymin) &
        (coloc_data[var2] < ymax))

    xdata = coloc_data[var1][in_data]
    ydata = coloc_data[var2][in_data]

    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
   
    if(var3 is None):
        ax.scatter(xdata, ydata, s = 6)
    else: 
        scat = ax.scatter(xdata, ydata, c = coloc_data[var3][in_data], \
            s = 6)
        if(colorbar):
            cbar = plt.colorbar(scat, ax = ax, pad = 0.03, fraction = 0.052, \
                extend = 'both')
            cbar.set_label(var3)

    if(trend):
        plot_trend_line(ax, xdata, ydata, color='tab:red', linestyle = '-', \
            slope = 'thiel-sen')

    ax.set_title(coloc_data['date_str'])
 
    #ax.set_xlim([xmin, xmax])
    #ax.set_ylim([ymin, ymax])
    ax.set_xlabel(var1)
    ax.set_ylabel(var2)

    if(save):
        outname = 'arctic_compare_scatter_' + coloc_data['date_str'] + '_' + \
            var1 + '_' + var2 + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname) 
    else:
        plt.show()

def plot_compare_colocate_cloud(coloc_data, save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data)
    else:
        dt_date_str = datetime.strptime(coloc_data['date_str'], '%Y%m%d%H%M')

    pdate = dt_date_str.strftime('%Y-%m-%d') 

    plt.close('all')
    mapcrs = ccrs.NorthPolarStereo()
    fig = plt.figure(figsize = (9, 9))
    ax1 = fig.add_subplot(2,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(2,2,2, projection = mapcrs)
    ax3 = fig.add_subplot(2,2,3, projection = mapcrs)
    ax4 = fig.add_subplot(2,2,4, projection = mapcrs)
    
    zoom = False
    plot_spatial(ax1, coloc_data['LON'], coloc_data['LAT'], \
        coloc_data['MODIS_CH1'], \
        pdate, cmap = 'Greys_r', zoom = zoom)
    plot_spatial(ax2, coloc_data['LON'], coloc_data['LAT'], \
        coloc_data['MODIS_CH7'], \
        pdate, cmap = 'Greys_r', vmax = 0.4, zoom = zoom)
    plot_spatial(ax3, coloc_data['LON'], coloc_data['LAT'], \
        coloc_data['MODIS_CLD'], \
        pdate, cmap = 'jet', zoom = zoom)
    plot_spatial(ax4, coloc_data['LON'], coloc_data['LAT'], \
        coloc_data['OMI_RAW'], \
        pdate, cmap = 'jet', zoom = zoom)

    ax1.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax2.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax3.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax4.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax1.coastlines()
    ax2.coastlines()
    ax3.coastlines()
    ax4.coastlines()

    fig.tight_layout()

    if(save):
        outname = 'arctic_comp_cloud_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image",outname)
    else:
        plt.show()

def plot_compare_colocate_spatial_category(coloc_data, cat = 'ALL', minlat = 65., \
        zoom = False, save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data)
    else:
        dt_date_str = datetime.strptime(coloc_data['date_str'], '%Y%m%d%H%M')

    ## Read the colocated data
    ## -----------------------
    #data = h5py.File('colocated_subset_' + date_str + '.hdf5','r')

    if(cat == 'ALL'):
        plot_var = 'NSIDC_ICE'
    else:
        plot_var = 'NSIDC_' + cat 

    # Make the overall figure
    # -----------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (10,8))
    ax1 = fig1.add_subplot(2,3,1, projection = mapcrs)
    ax2 = fig1.add_subplot(2,3,2, projection = mapcrs)
    ax3 = fig1.add_subplot(2,3,3, projection = mapcrs)
    ax4 = fig1.add_subplot(2,3,4, projection = mapcrs)
    ax5 = fig1.add_subplot(2,3,5, projection = mapcrs)
    ax6 = fig1.add_subplot(2,3,6, projection = mapcrs)
  
    pdate = dt_date_str.strftime('%Y-%m-%d') 
    # Plot the MODIS true-color and channel data
    # ------------------------------------------
    plot_spatial(ax1, coloc_data['LON'], coloc_data['LAT'], \
        select_category(coloc_data, 'MODIS_CH1', cat), \
        #np.ma.masked_where(coloc_data[plot_var].mask, coloc_data['MODIS_CH1']), \
        pdate, cmap = 'Greys_r', zoom = zoom)
    ##!#ax1.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH1'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'Greys_r')
    ##!#ax1.coastlines()
    ##!#ax1.set_extent([-180, 180, 60, 90], datacrs)
    plot_spatial(ax2, coloc_data['LON'], coloc_data['LAT'], \
        select_category(coloc_data, 'MODIS_CH7', cat), \
        #np.ma.masked_where(coloc_data[plot_var].mask, coloc_data['MODIS_CH7']), \
        pdate, cmap = 'Greys_r', zoom = zoom)
    ##!#ax2.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH7'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'Greys_r')
    ##!#ax2.coastlines()
    ##!#ax2.set_extent([-180, 180, 60, 90], datacrs)

    # Plot the NSIDC coloc_data
    # -------------------
    plot_spatial(ax3, coloc_data['LON'], coloc_data['LAT'], \
        select_category(coloc_data, plot_var, cat), \
        #np.ma.masked_where(coloc_data[plot_var].mask, coloc_data[plot_var]), \
        pdate, cmap = 'ocean', zoom = zoom)
    ##!#ax3.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['NSIDC'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'ocean')
    ##!#ax3.coastlines()
    ##!#ax3.set_extent([-180, 180, 60, 90], datacrs)

    # Plot the OMI coloc_data
    # -----------------
    plot_spatial(ax4, coloc_data['LON'], coloc_data['LAT'], \
        select_category(coloc_data, 'OMI', cat), \
        #np.ma.masked_where(coloc_data[plot_var].mask, coloc_data['OMI']), \
        pdate, cmap = 'jet', zoom = zoom)
    ##!#ax4.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['OMI'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'jet')
    ##!#ax4.coastlines()
    ##!#ax4.set_extent([-180, 180, 60, 90], datacrs)
    
    # Plot the CERES coloc_data
    # -------------------
    plot_spatial(ax5, coloc_data['LON'], coloc_data['LAT'], \
        select_category(coloc_data, 'CERES_SWF', cat), \
        #np.ma.masked_where(coloc_data[plot_var].mask, coloc_data['CERES_SWF']), \
        pdate, cmap = 'plasma', zoom = zoom)
    ##!#ax5.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['CERES_SWF'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'plasma')
    ##!#ax5.coastlines()
    ##!#ax5.set_extent([-180, 180, 60, 90], datacrs)
    plot_spatial(ax6, coloc_data['LON'], coloc_data['LAT'], \
        select_category(coloc_data, 'CERES_LWF', cat), \
        #np.ma.masked_where(coloc_data[plot_var].mask, coloc_data['CERES_LWF']), \
        pdate, cmap = 'plasma', zoom = zoom)
    ##!#ax6.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['CERES_LWF'], \
    ##!#    transform = datacrs, shading = 'auto', cmap = 'plasma')
    ##!#ax6.coastlines()
    ##!#ax6.set_extent([-180, 180, 60, 90], datacrs)


    if(save):
        outname = 'arctic_compare_spatial_' + \
            dt_date_str.strftime('%Y%m%d%H%M') + '.png'
        fig1.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def event_category_slopes_all(coloc_data, var1, var2, var3 = None, \
    cat = "ALL", minlat = 65., xmin = None, xmax = None, ymin = None, \
    ymax = None, trend = False, num_points = 10, \
    restrict_sza = False, color = None, save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data)
    else:
        if(len(coloc_data['date_str']) == 6):
            fmt = '%Y%m'
        elif(len(coloc_data['date_str']) == 8):
            fmt = '%Y%m%d'
        elif(len(coloc_data['date_str']) == 12):
            fmt = '%Y%m%d%H%M'
        dt_date_str = datetime.strptime(coloc_data['date_str'], fmt)

    cat = 'ICE_CLOUD'
    cat = 'ICE_CLEAR'
    cat = 'OCEAN_CLOUD'
    cat = 'OCEAN_CLEAR'
    cat = 'LAND_CLOUD'
    cat = 'LAND_CLEAR'
    begin_cats = ['ICE','OCEAN','LAND']
    end_cats   = ['CLOUD','CLEAR']

    return_dict = {}
    for bc in begin_cats:
        for ec in end_cats:
            total_cat = bc + '_' + ec
            return_dict[total_cat] = {}
            thiel_slope, lin_slope, numx, numy, lin_pval = \
                calc_category_slopes(coloc_data, var1, \
                var2, cat = total_cat, xmin = xmin, xmax = xmax, \
                ymin = ymin, ymax = ymax, num_points = num_points)

            return_dict[total_cat]['Thiel'] = thiel_slope
            return_dict[total_cat]['Linear'] = lin_slope
            return_dict[total_cat]['num-x'] = numx
            return_dict[total_cat]['num-y'] = numy
            return_dict[total_cat]['lin_pval'] = lin_pval

    return return_dict 

def calc_category_slopes(coloc_data, var1, var2, var3 = None, \
        cat = "ALL", minlat = 65., xmin = None, xmax = None, ymin = None, \
        ymax = None, num_points = 10, restrict_sza = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data)
    else:
        if(len(coloc_data['date_str']) == 6):
            fmt = '%Y%m'
        elif(len(coloc_data['date_str']) == 8):
            fmt = '%Y%m%d'
        elif(len(coloc_data['date_str']) == 12):
            fmt = '%Y%m%d%H%M'
        dt_date_str = datetime.strptime(coloc_data['date_str'], fmt)

    # Pull out the category data for each variable
    # --------------------------------------------
    xdata = select_category(coloc_data, var1, cat)
    ydata = select_category(coloc_data, var2, cat)

    if(xmin is None):
        xmin = np.nanmin(xdata)
    if(xmax is None):
        xmax = np.nanmax(xdata)
    if(ymin is None):
        ymin = np.nanmin(ydata)
    if(ymax is None):
        ymax = np.nanmax(ydata)

    mask_xdata = np.ma.masked_where((xdata < xmin) | \
        (xdata > xmax) | (ydata < ymin) |
        (ydata > ymax), xdata).compressed()
    mask_ydata = np.ma.masked_where((xdata < xmin) | \
        (xdata > xmax) | (ydata < ymin) |
        (ydata > ymax), ydata).compressed()

    #print(cat, len(mask_xdata), len(mask_ydata))

    #plot_trend_line(ax, mask_xdata.compressed(), mask_ydata.compressed(), color=trend_color, linestyle = '-', \
    #    slope = 'thiel-sen')

    if((len(mask_xdata) < num_points) | (len(mask_ydata) < num_points)):
        return np.nan, np.nan, len(mask_xdata), len(mask_ydata), np.nan
    
    else:
        res = stats.theilslopes(mask_ydata, mask_xdata, 0.95)
        #print("\tTheil-Sen: {0}x + {1}".format(res[0], res[1]))

        #zdata = np.polyfit(mask_xdata, mask_ydata, 1)
        slope,intercept,r_value,p_value,std_err = \
            stats.linregress(mask_xdata,mask_ydata)
        #print("\tLin Regress: {0}x + {1}".format(*zdata))

        # Thiel slope, Lin regress slope
        return res[0], slope, len(mask_xdata), len(mask_ydata), p_value
        #return res[0], zdata[0], len(mask_xdata), len(mask_ydata), p_value


# cat is either "LAND", "ICE", "OCEAN", or "ALL"
def plot_compare_scatter_category(coloc_data, var1, var2, var3 = None, \
        cat = "ALL", minlat = 65., xmin = None, xmax = None, ymin = None, \
        ymax = None, ax = None, colorbar = True, trend = None, zoom = False, \
        restrict_sza = False, color = None, save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data)
    else:
        if(len(coloc_data['date_str']) == 6):
            fmt = '%Y%m'
        elif(len(coloc_data['date_str']) == 8):
            fmt = '%Y%m%d'
        elif(len(coloc_data['date_str']) == 12):
            fmt = '%Y%m%d%H%M'
        dt_date_str = datetime.strptime(coloc_data['date_str'], fmt)

    if(cat == 'ALL'):
        plot_var = 'NSIDC_ICE'
    else:
        plot_var = 'NSIDC_' + cat 

    # Pull out the category data for each variable
    # --------------------------------------------
    xdata = select_category(coloc_data, var1, cat)
    ydata = select_category(coloc_data, var2, cat)

    if(xmin is None):
        xmin = np.nanmin(xdata)
    if(xmax is None):
        xmax = np.nanmax(xdata)
    if(ymin is None):
        ymin = np.nanmin(ydata)
    if(ymax is None):
        ymax = np.nanmax(ydata)
  
    # If desired, remove data on other size of terminator
    # ---------------------------------------------------
    #if(restrict_sza):
    #    xdata = np.ma.masked_where(coloc_data[' 
 
    #in_data = np.where((coloc_data[var1] > xmin) & \
    #    (coloc_data[var1] < xmax) & (coloc_data[var2] > ymin) &
    #    (coloc_data[var2] < ymax))

    mask_xdata = np.ma.masked_where((xdata < xmin) | \
        (xdata > xmax) | (ydata < ymin) |
        (ydata > ymax), xdata)
    mask_ydata = np.ma.masked_where((xdata < xmin) | \
        (xdata > xmax) | (ydata < ymin) |
        (ydata > ymax), ydata)

    print("Here:", cat, np.nanmax(mask_xdata), np.nanmax(mask_ydata))

    ## Now mask using the surface type
    ## -------------------------------
    #xdata = np.ma.masked_where(coloc_data[plot_var].mask, xdata)
    #ydata = np.ma.masked_where(coloc_data[plot_var].mask, ydata)


    #xdata = coloc_data[var1][in_data]
    #ydata = coloc_data[var2][in_data]

    # Now, mask the data 
    if(color is None):
        pcolor = 'tab:blue'
        trend_color = 'tab:red'
    else:
        pcolor = color
        trend_color = color

    # Set up the figure
    # -----------------
    in_ax = True 
    if(ax is None): 
        plt.close('all')
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    
    print(plot_var, len(mask_xdata), len(mask_ydata), np.nanmax(mask_xdata), np.nanmax(mask_ydata))

    if(var3 is None):
        ax.scatter(mask_xdata, mask_ydata, s = 6, color = pcolor, label = cat)
    else: 
        zdata = select_category(coloc_data, var3, cat)
        mask_zdata = np.ma.masked_where((xdata < xmin) | \
            (xdata > xmax) | (ydata < ymin) |
            (ydata > ymax), zdata)

        scat = ax.scatter(mask_xdata, mask_ydata, c = mask_zdata, \
            s = 6)
        #scat = ax.scatter(xdata, ydata, c = coloc_data[var3][in_data], \

        if(colorbar):
            cbar = plt.colorbar(scat, ax = ax, pad = 0.03, fraction = 0.052, \
                extend = 'both')
            cbar.set_label(var3)

    if((trend is not None) & \
       (len(mask_xdata.compressed()) != 0) & \
       (len(mask_ydata.compressed()) != 0)):
        print(cat, len(mask_xdata.compressed()), len(mask_ydata.compressed()))
        plot_trend_line(ax, mask_xdata.compressed(), mask_ydata.compressed(), color=trend_color, linestyle = '-', \
            slope = trend)

    ax.set_title(coloc_data['date_str'] + '\n' + plot_var)
 
    #ax.set_xlim([xmin, xmax])
    #ax.set_ylim([ymin, ymax])
    ax.set_xlabel(var1)
    ax.set_ylabel(var2)

    if(not in_ax):
        if(save):
            outname = 'arctic_compare_scatter_' + cat + '_' + coloc_data['date_str'] + '_' + \
                var1 + '_' + var2 + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname) 
        else:
            plt.show()

def plot_compare_AI_combined_category(coloc_data, var2 = 'CERES_SWF', \
        cat = "ALL", minlat = 65., xmin = 1.0, xmax = None, ymin = None, \
        ymax = None, ax = None, \
        colorbar = True, trend = False, zoom = False, color = None, \
        compare_tropomi = False, mask_spatial = False, save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data, minlat = minlat, \
            compare_tropomi = compare_tropomi)
    else:
        dt_date_str = datetime.strptime(coloc_data['date_str'], '%Y%m%d%H%M')

    raw_min = xmin
    pert_min = xmin
    trop_min = xmin

    # Make the overall figure
    # -----------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (12,12))
    ax1  = fig1.add_subplot(4,3,1, projection = mapcrs)
    ax2  = fig1.add_subplot(4,3,2, projection = mapcrs)
    ax3  = fig1.add_subplot(4,3,3, projection = mapcrs)
    ax4  = fig1.add_subplot(4,3,4)
    ax5  = fig1.add_subplot(4,3,5)
    ax6  = fig1.add_subplot(4,3,6)
    ax7  = fig1.add_subplot(4,3,7)
    ax8  = fig1.add_subplot(4,3,8)
    ax9  = fig1.add_subplot(4,3,9)
    ax10 = fig1.add_subplot(4,3,10)
    ax11 = fig1.add_subplot(4,3,11)
    ax12 = fig1.add_subplot(4,3,12)
  
    pdate = dt_date_str.strftime('%Y-%m-%d') 
    # Plot the OMI coloc_data
    # -----------------
    raw_plot_min = None
    pert_plot_min = None
    trop_plot_min = None
    if(mask_spatial):
        raw_plot_min = raw_min 
        pert_plot_min = raw_min 
        trop_plot_min = raw_min 
    plot_spatial(ax1, coloc_data['LON'], coloc_data['LAT'], coloc_data['OMI_RAW'], \
        pdate, cmap = 'jet', zoom = zoom, xmin = raw_plot_min)
    ax1.set_title('OMI_RAW')
    plot_spatial(ax2, coloc_data['LON'], coloc_data['LAT'], coloc_data['OMI_PERT'], \
        pdate, cmap = 'jet', zoom = zoom, xmin = pert_plot_min)
    ax2.set_title('OMI_PERT')
    plot_spatial(ax3, coloc_data['LON'], coloc_data['LAT'], coloc_data['TROP_AI'], \
        pdate, cmap = 'jet', zoom = zoom, xmin = trop_plot_min)
    ax3.set_title('TROP_AI')
    
    # Plot the scatter data
    date_str = dt_date_str.strftime('%Y%m%d%H%M')

    # Plot the OMI RAW comparison stuff
    plot_compare_scatter_category(coloc_data, 'OMI_RAW', var2, var3 = None, \
        cat = 'ICE_CLOUD', minlat = 65., xmin = raw_min, xmax = None, ymin = None, ymax = None, \
        ax = ax4, colorbar = True, trend = trend, zoom = zoom, save = False,\
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, 'OMI_RAW', var2, var3 = None, \
        cat = 'ICE_CLEAR', minlat = 65., xmin = raw_min, xmax = None, ymin = None, ymax = None, \
        ax = ax4, colorbar = True, trend = trend, zoom = zoom, save = False,\
        color = 'tab:orange')
    plot_compare_scatter_category(coloc_data, 'OMI_RAW', var2, var3 = None, \
        cat = 'OCEAN_CLOUD', minlat = 65., xmin = raw_min, xmax = None, ymin = None, ymax = None, \
        ax = ax5, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, 'OMI_RAW', var2, var3 = None, \
        cat = 'OCEAN_CLEAR', minlat = 65., xmin = raw_min, xmax = None, ymin = None, ymax = None, \
        ax = ax5, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:orange')
    plot_compare_scatter_category(coloc_data, 'OMI_RAW', var2, var3 = None, \
        cat = 'LAND_CLOUD', minlat = 65., xmin = raw_min, xmax = None, ymin = None, ymax = None, \
        ax = ax6, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, 'OMI_RAW', var2, var3 = None, \
        cat = 'LAND_CLEAR', minlat = 65., xmin = raw_min, xmax = None, ymin = None, ymax = None, \
        ax = ax6, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:orange')

    # Plot the OMI PERT comparison stuff
    plot_compare_scatter_category(coloc_data, 'OMI_PERT', var2, var3 = None, \
        cat = 'ICE_CLOUD', minlat = 65., xmin = pert_min, xmax = None, ymin = None, ymax = None, \
        ax = ax7, colorbar = True, trend = trend, zoom = zoom, save = False,\
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, 'OMI_PERT', var2, var3 = None, \
        cat = 'ICE_CLEAR', minlat = 65., xmin = pert_min, xmax = None, ymin = None, ymax = None, \
        ax = ax7, colorbar = True, trend = trend, zoom = zoom, save = False,\
        color = 'tab:orange')
    plot_compare_scatter_category(coloc_data, 'OMI_PERT', var2, var3 = None, \
        cat = 'OCEAN_CLOUD', minlat = 65., xmin = pert_min, xmax = None, ymin = None, ymax = None, \
        ax = ax8, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, 'OMI_PERT', var2, var3 = None, \
        cat = 'OCEAN_CLEAR', minlat = 65., xmin = pert_min, xmax = None, ymin = None, ymax = None, \
        ax = ax8, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:orange')
    plot_compare_scatter_category(coloc_data, 'OMI_PERT', var2, var3 = None, \
        cat = 'LAND_CLOUD', minlat = 65., xmin = pert_min, xmax = None, ymin = None, ymax = None, \
        ax = ax9, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, 'OMI_PERT', var2, var3 = None, \
        cat = 'LAND_CLEAR', minlat = 65., xmin = pert_min, xmax = None, ymin = None, ymax = None, \
        ax = ax9, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:orange')

    # Plot the OMI RAW comparison stuff
    plot_compare_scatter_category(coloc_data, 'TROP_AI', var2, var3 = None, \
        cat = 'ICE_CLOUD', minlat = 65., xmin = trop_min, xmax = None, ymin = None, ymax = None, \
        ax = ax10, colorbar = True, trend = trend, zoom = zoom, save = False,\
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, 'TROP_AI', var2, var3 = None, \
        cat = 'ICE_CLEAR', minlat = 65., xmin = trop_min, xmax = None, ymin = None, ymax = None, \
        ax = ax10, colorbar = True, trend = trend, zoom = zoom, save = False,\
        color = 'tab:orange')
    plot_compare_scatter_category(coloc_data, 'TROP_AI', var2, var3 = None, \
        cat = 'OCEAN_CLOUD', minlat = 65., xmin = trop_min, xmax = None, ymin = None, ymax = None, \
        ax = ax11, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, 'TROP_AI', var2, var3 = None, \
        cat = 'OCEAN_CLEAR', minlat = 65., xmin = trop_min, xmax = None, ymin = None, ymax = None, \
        ax = ax11, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:orange')
    plot_compare_scatter_category(coloc_data, 'TROP_AI', var2, var3 = None, \
        cat = 'LAND_CLOUD', minlat = 65., xmin = trop_min, xmax = None, ymin = None, ymax = None, \
        ax = ax12, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, 'TROP_AI', var2, var3 = None, \
        cat = 'LAND_CLEAR', minlat = 65., xmin = trop_min, xmax = None, ymin = None, ymax = None, \
        ax = ax12, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:orange')


    # Add legend to last box
    # ----------------------
    custom_lines = [Line2D([0], [0], color='tab:blue', linestyle = ':'),
                    Line2D([0], [0], color='tab:orange', linestyle = ':')]

    ax12.legend(custom_lines, ['Cloud', 'Clear'],\
        fontsize = 10, loc = 1)

    fig1.tight_layout()

    if(save):
        date_str = dt_date_str.strftime('%Y%m%d%H%M')
        outname = 'arctic_compare_AI_combined_category_nomin_' + date_str + '.png'
        fig1.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


def plot_compare_combined_category(coloc_data, var1 = 'OMI', \
        var2 = 'CERES_SWF', var3 = None, cat = "ALL", minlat = 65., \
        xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
        colorbar = True, trend = False, zoom = False, color = None, \
        compare_tropomi = False, \
        save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data, minlat = minlat, \
            compare_tropomi = compare_tropomi)
    else:
        dt_date_str = datetime.strptime(coloc_data['date_str'], '%Y%m%d%H%M')

    # Make the overall figure
    # -----------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (12,9))
    ax1 = fig1.add_subplot(3,3,1, projection = mapcrs)
    ax2 = fig1.add_subplot(3,3,2, projection = mapcrs)
    ax3 = fig1.add_subplot(3,3,3, projection = mapcrs)
    ax4 = fig1.add_subplot(3,3,4, projection = mapcrs)
    ax5 = fig1.add_subplot(3,3,5, projection = mapcrs)
    ax6 = fig1.add_subplot(3,3,6, projection = mapcrs)
    ax7 = fig1.add_subplot(3,3,7)
    ax8 = fig1.add_subplot(3,3,8)
    ax9 = fig1.add_subplot(3,3,9)
  
    pdate = dt_date_str.strftime('%Y-%m-%d') 
    # Plot the MODIS true-color and channel data
    # ------------------------------------------
    plot_spatial(ax1, coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH1'], \
        pdate, cmap = 'Greys_r', zoom = zoom)
    plot_spatial(ax2, coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH7'], \
        pdate, cmap = 'Greys_r', zoom = zoom)

    ax1.set_title('MODIS 0.64 μm')
    ax2.set_title('MODIS 2.13 μm')

    # Plot the NSIDC coloc_data
    # -------------------
    plot_spatial(ax3, coloc_data['LON'], coloc_data['LAT'], coloc_data['NSIDC_ICE'], \
        pdate, cmap = 'ocean', zoom = zoom)

    ax3.set_title('NSIDC Ice Conc')

    # Plot the OMI coloc_data
    # -----------------
    plot_spatial(ax4, coloc_data['LON'], coloc_data['LAT'], coloc_data[var1], \
        pdate, cmap = 'jet', zoom = zoom)
    ax4.set_title(var1)
    
    # Plot the CERES coloc_data
    # -------------------
    plot_spatial(ax5, coloc_data['LON'], coloc_data['LAT'], coloc_data['CERES_SWF'], \
        pdate, cmap = 'plasma', zoom = zoom)
    plot_spatial(ax6, coloc_data['LON'], coloc_data['LAT'], coloc_data['CERES_LWF'], \
        pdate, cmap = 'plasma', zoom = zoom)

    ax5.set_title('CERES SWF')
    ax6.set_title('CERES LWF')

    # Plot the box plots
    # ------------------
    #plot_FuLiou_boxes(axs = [ax7, ax8, ax9])

    # Plot the scatter data
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
        cat = 'ICE_CLOUD', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax7, colorbar = True, trend = trend, zoom = zoom, save = False,\
        color = 'tab:blue')
    plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
        cat = 'ICE_CLEAR', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax7, colorbar = True, trend = trend, zoom = zoom, save = False,\
        color = 'tab:orange')
    plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
        cat = 'OCEAN_CLOUD', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax8, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
        cat = 'OCEAN_CLEAR', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax8, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:orange')
    plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
        cat = 'LAND_CLOUD', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax9, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
        cat = 'LAND_CLEAR', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax9, colorbar = True, trend = trend, zoom = zoom, save = False, \
        color = 'tab:orange')

    # Add legend to last box
    # ----------------------
    custom_lines = [Line2D([0], [0], color='tab:blue', linestyle = ':'),
                    Line2D([0], [0], color='tab:orange', linestyle = ':')]

    ax9.legend(custom_lines, ['Cloud', 'Clear'],\
        fontsize = 10, loc = 1)

    fig1.tight_layout()

    if(save):
        date_str = dt_date_str.strftime('%Y%m%d%H%M')
        outname = 'arctic_compare_combined_category_nomin_' + date_str + '.png'
        fig1.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_compare_combined_category_event(coloc_data, var1 = 'OMI', \
        var2 = 'CERES_SWF', cat = "ALL", minlat = 65., \
        xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
        colorbar = True, trend = True, zoom = False, color = None, \
        save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d')
        coloc_data = read_colocated_combined(coloc_data)
    else:
        dt_date_str = datetime.strptime(coloc_data['date_str'], '%Y%m%d')

    fig = plt.figure(figsize = (12,4))
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    ##ax1 = fig.add_subplot(2,3,1)
    ##ax2 = fig.add_subplot(2,3,4)
    ##ax3 = fig.add_subplot(2,3,2)
    ##ax4 = fig.add_subplot(2,3,5)
    ##ax5 = fig.add_subplot(2,3,3)
    ##ax6 = fig.add_subplot(2,3,6)
    
    plot_compare_scatter_category(coloc_data, var1, var2, var3 = None, \
        cat = 'ICE_CLOUD', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax1, colorbar = True, trend = trend, zoom = False, save = False,\
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, var1, var2, var3 = None, \
        cat = 'ICE_CLEAR', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax1, colorbar = True, trend = trend, zoom = False, save = False,\
        color = 'tab:orange')
    plot_compare_scatter_category(coloc_data, var1, var2, var3 = None, \
        cat = 'OCEAN_CLOUD', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax2, colorbar = True, trend = trend, zoom = False, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, var1, var2, var3 = None, \
        cat = 'OCEAN_CLEAR', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax2, colorbar = True, trend = trend, zoom = False, save = False, \
        color = 'tab:orange')
    plot_compare_scatter_category(coloc_data, var1, var2, var3 = None, \
        cat = 'LAND_CLOUD', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax3, colorbar = True, trend = trend, zoom = False, save = False, \
        color = 'tab:blue')
    plot_compare_scatter_category(coloc_data, var1, var2, var3 = None, \
        cat = 'LAND_CLEAR', minlat = 65., xmin = xmin, xmax = None, ymin = None, ymax = None, \
        ax = ax3, colorbar = True, trend = trend, zoom = False, save = False, \
        color = 'tab:orange')
    ax1.legend()
    ax2.legend()
    ax3.legend()
    #plt.suptitle(coloc_data['date_str'])
    #plt.suptitle(coloc_data)
    fig.tight_layout()

    if(save):    
        date_str = dt_date_str.strftime('%Y%m%d')
        outname = 'arctic_daily_scatter_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:    
        plt.show()
    
def plot_compare_all_slopes(date_strs = None, save = False, \
        return_slope_dict = False, return_slope_means = False, \
        plot_comp_figs = True):

    if(date_strs is None):
        date_strs = case_dates
    ##             ##!#'201605151925',  # MEDIOCRE
    ##             ##!#'201605152104',  # MEDIOCRE
    ##             ##!#'201605152243',  # MEDIOCRE
    ##             ##!#'201605162148',  # MEDIOCRE
    ##             ##!#'200607260017',  # GOOD
    ##             ##!#'200607252238',  # GOOD
    ##             ##!#'200607260156',  # GOOD
    ##             ##!#'200607260335',  # GOOD
    ##             ##!#'200607260513',  # GOOD
    ##             '201808241343',
    ##            ]
    
    slope_dict = {}
    
    #coloc_data = '201708161504'
    num_points = 2
    for date_str in date_strs:
        print(date_str)
        slope_dict[date_str] = event_category_slopes_all(date_str, 'OMI', 'CERES_SWF', var3 = None, \
            cat = "ALL", minlat = 65., xmin = 1.0, xmax = None, ymin = None, \
            ymax = None, trend = False, num_points = num_points, \
            restrict_sza = False, color = None, save = False)
    
    
    # Extract just the lin regress slopes
    # -----------------------------------
    extract_dict = {}
    dkeys = slope_dict['201708181451'].keys()
    
    for dkey in dkeys:
        extract_dict[dkey] = {}
    
        lin_slopes    = np.ma.masked_invalid(np.array([slope_dict[tkey][dkey]['Linear'] for tkey in slope_dict.keys()]))
        lin_pvals     = np.ma.masked_invalid(np.array([slope_dict[tkey][dkey]['lin_pval'] for tkey in slope_dict.keys()]))
        thiel_slopes  = np.ma.masked_invalid(np.array([slope_dict[tkey][dkey]['Thiel'] for tkey in slope_dict.keys()]))
    
        extract_dict[dkey]['Linear']   = np.ma.masked_where((lin_slopes > 200) | (lin_slopes < -200) | (lin_pvals > 0.1), lin_slopes)
        #extract_dict[dkey]['Linear']   = np.ma.masked_where((lin_slopes > 500), lin_slopes)
        extract_dict[dkey]['Thiel']    = np.ma.masked_where((thiel_slopes > 200) | (thiel_slopes < -200), thiel_slopes)
        extract_dict[dkey]['lin_pval'] = lin_pvals

    # Calculate statistics
    # --------------------   
    mean_lin_ICD = np.nanmean(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Linear']))
    mean_lin_ICR = np.nanmean(np.ma.masked_invalid(extract_dict['ICE_CLEAR']['Linear']))
    mean_lin_OCD = np.nanmean(np.ma.masked_invalid(extract_dict['OCEAN_CLOUD']['Linear']))
    mean_lin_OCR = np.nanmean(np.ma.masked_invalid(extract_dict['OCEAN_CLEAR']['Linear']))
    mean_lin_LCD = np.nanmean(np.ma.masked_invalid(extract_dict['LAND_CLOUD']['Linear']))
    mean_lin_LCR = np.nanmean(np.ma.masked_invalid(extract_dict['LAND_CLEAR']['Linear']))

    mean_thl_ICD = np.nanmean(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Thiel']))
    mean_thl_ICR = np.nanmean(np.ma.masked_invalid(extract_dict['ICE_CLEAR']['Thiel']))
    mean_thl_OCD = np.nanmean(np.ma.masked_invalid(extract_dict['OCEAN_CLOUD']['Thiel']))
    mean_thl_OCR = np.nanmean(np.ma.masked_invalid(extract_dict['OCEAN_CLEAR']['Thiel']))
    mean_thl_LCD = np.nanmean(np.ma.masked_invalid(extract_dict['LAND_CLOUD']['Thiel']))
    mean_thl_LCR = np.nanmean(np.ma.masked_invalid(extract_dict['LAND_CLEAR']['Thiel']))

    std_lin_ICD = np.nanstd(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Linear']))
    std_lin_ICR = np.nanstd(np.ma.masked_invalid(extract_dict['ICE_CLEAR']['Linear']))
    std_lin_OCD = np.nanstd(np.ma.masked_invalid(extract_dict['OCEAN_CLOUD']['Linear']))
    std_lin_OCR = np.nanstd(np.ma.masked_invalid(extract_dict['OCEAN_CLEAR']['Linear']))
    std_lin_LCD = np.nanstd(np.ma.masked_invalid(extract_dict['LAND_CLOUD']['Linear']))
    std_lin_LCR = np.nanstd(np.ma.masked_invalid(extract_dict['LAND_CLEAR']['Linear']))

    std_thl_ICD = np.nanstd(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Thiel']))
    std_thl_ICR = np.nanstd(np.ma.masked_invalid(extract_dict['ICE_CLEAR']['Thiel']))
    std_thl_OCD = np.nanstd(np.ma.masked_invalid(extract_dict['OCEAN_CLOUD']['Thiel']))
    std_thl_OCR = np.nanstd(np.ma.masked_invalid(extract_dict['OCEAN_CLEAR']['Thiel']))
    std_thl_LCD = np.nanstd(np.ma.masked_invalid(extract_dict['LAND_CLOUD']['Thiel']))
    std_thl_LCR = np.nanstd(np.ma.masked_invalid(extract_dict['LAND_CLEAR']['Thiel']))
    
    fig = plt.figure(figsize = (9, 11))
    ax1 = fig.add_subplot(3,2,1)
    ax2 = fig.add_subplot(3,2,2)
    ax3 = fig.add_subplot(3,2,3)
    ax4 = fig.add_subplot(3,2,4)
    ax5 = fig.add_subplot(3,2,5)
    ax6 = fig.add_subplot(3,2,6)
    
    num_bins = 20
    ax1.hist(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Linear']), bins = num_bins, alpha = 0.5, label = 'Linear')
    ax1.hist(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Thiel']), bins = num_bins, alpha = 0.5, label = 'Thiel')
    print('ICE_CLOUD')
    print('\tLinear:', mean_lin_ICD, std_lin_ICD)
    print('\tThiel: ', mean_thl_ICD, std_thl_ICD)
    ax2.hist(extract_dict['ICE_CLEAR']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
    ax2.hist(extract_dict['ICE_CLEAR']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
    print('ICE_CLEAR')
    print('\tLinear:', mean_lin_ICR, std_lin_ICR)
    print('\tThiel: ', mean_thl_ICR, std_thl_ICR)
    ax3.hist(extract_dict['OCEAN_CLOUD']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
    ax3.hist(extract_dict['OCEAN_CLOUD']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
    print('OCEAN_CLOUD')
    print('\tLinear:', mean_lin_OCD, std_lin_OCD)
    print('\tThiel: ', mean_thl_OCD, std_thl_OCD)
    ax4.hist(extract_dict['OCEAN_CLEAR']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
    ax4.hist(extract_dict['OCEAN_CLEAR']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
    print('OCEAN_CLEAR')
    print('\tLinear:', mean_lin_OCR, std_lin_OCR)
    print('\tThiel: ', mean_thl_OCR, std_thl_OCR)
    ax5.hist(extract_dict['LAND_CLOUD']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
    ax5.hist(extract_dict['LAND_CLOUD']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
    print('LAND_CLOUD')
    print('\tLinear:', mean_lin_LCD, std_lin_LCD)
    print('\tThiel: ', mean_thl_LCD, std_thl_LCD)
    ax6.hist(extract_dict['LAND_CLEAR']['Linear'], bins = num_bins, alpha = 0.5, label = 'Linear')
    ax6.hist(extract_dict['LAND_CLEAR']['Thiel'], bins = num_bins, alpha = 0.5, label = 'Thiel')
    print('LAND_CLEAR')
    print('\tLinear:', mean_lin_LCR, std_lin_LCR)
    print('\tThiel: ', mean_thl_LCR, std_thl_LCR)
    #ax2.scatter(np.ma.masked_invalid(extract_dict[var]['Linear']), np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Thiel']))
   
    # Calculate max bin center values
    # -------------------------------
    counts, hist = np.histogram(np.ma.masked_invalid(extract_dict['ICE_CLOUD']['Linear']).compressed(), bins = num_bins)
    max_bin_lin_ICD = 0.5*(hist[np.argmax(counts)] + hist[np.argmax(counts)+1])
    counts, hist = np.histogram(np.ma.masked_invalid(extract_dict['ICE_CLEAR']['Linear']).compressed(), bins = num_bins)
    max_bin_lin_ICR = 0.5*(hist[np.argmax(counts)] + hist[np.argmax(counts)+1])
    counts, hist = np.histogram(np.ma.masked_invalid(extract_dict['OCEAN_CLOUD']['Linear']).compressed(), bins = num_bins)
    max_bin_lin_OCD = 0.5*(hist[np.argmax(counts)] + hist[np.argmax(counts)+1])
    counts, hist = np.histogram(np.ma.masked_invalid(extract_dict['OCEAN_CLEAR']['Linear']).compressed(), bins = num_bins)
    max_bin_lin_OCR = 0.5*(hist[np.argmax(counts)] + hist[np.argmax(counts)+1])
    counts, hist = np.histogram(np.ma.masked_invalid(extract_dict['LAND_CLOUD']['Linear']).compressed(), bins = num_bins)
    max_bin_lin_LCD = 0.5*(hist[np.argmax(counts)] + hist[np.argmax(counts)+1])
    counts, hist = np.histogram(np.ma.masked_invalid(extract_dict['LAND_CLEAR']['Linear']).compressed(), bins = num_bins)
    max_bin_lin_LCR = 0.5*(hist[np.argmax(counts)] + hist[np.argmax(counts)+1])
 
    ax1.set_title('ICE_CLOUD\nLinear Max Bin Value: {0}\nLinear: SWF = {1}*AI +/- {2}\nThiel:  SWF = {3}*AI +/- {4}'.format(\
        np.round(max_bin_lin_ICD,2), \
        np.round(mean_lin_ICD,2), np.round(std_lin_ICD,2), \
        np.round(mean_thl_ICD,2), np.round(std_thl_ICD,2)))
    ax2.set_title('ICE_CLEAR\nLinear Max Bin Value: {0}\nLinear: SWF = {1}*AI +/- {2}\nThiel:  SWF = {3}*AI +/- {4}'.format(\
        np.round(max_bin_lin_ICR,2), \
        np.round(mean_lin_ICR,2), np.round(std_lin_ICR,2), \
        np.round(mean_thl_ICR,2), np.round(std_thl_ICR,2)))
    ax3.set_title('OCEAN_CLOUD\nLinear Max Bin Value: {0}\nLinear: SWF = {1}*AI +/- {2}\nThiel:  SWF = {3}*AI +/- {4}'.format(\
        np.round(max_bin_lin_OCD,2), \
        np.round(mean_lin_OCD,2), np.round(std_lin_OCD,2), \
        np.round(mean_thl_OCD,2), np.round(std_thl_OCD,2)))
    ax4.set_title('OCEAN_CLEAR\nLinear Max Bin Value: {0}\nLinear: SWF = {1}*AI +/- {2}\nThiel:  SWF = {3}*AI +/- {4}'.format(\
        np.round(max_bin_lin_OCR,2), \
        np.round(mean_lin_OCR,2), np.round(std_lin_OCR,2), \
        np.round(mean_thl_OCR,2), np.round(std_thl_OCR,2)))
    ax5.set_title('LAND_CLOUD\nLinear Max Bin Value: {0}\nLinear: SWF = {1}*AI +/- {2}\nThiel:  SWF = {3}*AI +/- {4}'.format(\
        np.round(max_bin_lin_LCD,2), \
        np.round(mean_lin_LCD,2), np.round(std_lin_LCD,2), \
        np.round(mean_thl_LCD,2), np.round(std_thl_LCD,2)))
    ax6.set_title('LAND_CLEAR\nLinear Max Bin Value: {0}\nLinear: SWF = {1}*AI +/- {2}\nThiel:  SWF = {3}*AI +/- {4}'.format(\
        np.round(max_bin_lin_LCR,2), \
        np.round(mean_lin_LCR,2), np.round(std_lin_LCR,2), \
        np.round(mean_thl_LCR,2), np.round(std_thl_LCR,2)))
   
    ax1.set_xlabel('SWF-AI slope (Wm-2/AI)')
    ax2.set_xlabel('SWF-AI slope (Wm-2/AI)')
    ax3.set_xlabel('SWF-AI slope (Wm-2/AI)')
    ax4.set_xlabel('SWF-AI slope (Wm-2/AI)')
    ax5.set_xlabel('SWF-AI slope (Wm-2/AI)')
    ax6.set_xlabel('SWF-AI slope (Wm-2/AI)')

    ax1.set_ylabel('Counts')
    ax2.set_ylabel('Counts')
    ax3.set_ylabel('Counts')
    ax4.set_ylabel('Counts')
    ax5.set_ylabel('Counts')
    ax6.set_ylabel('Counts')
 
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax5.legend()
    ax6.legend()

    fig.tight_layout()

    if(save):
        outname = 'omi_ceres_slope_histogram.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

    if(return_slope_dict):   
        return slope_dict, extract_dict
    elif(return_slope_means):
        return_dict = {}
        return_dict['ICE_CLOUD']    = np.nanmean(extract_dict['ICE_CLOUD']['Linear'])
        return_dict['ICE_CLEAR']    = np.nanmean(extract_dict['ICE_CLEAR']['Linear'])
        return_dict['OCEAN_CLOUD']  = np.nanmean(extract_dict['OCEAN_CLOUD']['Linear'])
        return_dict['OCEAN_CLEAR']  = np.nanmean(extract_dict['OCEAN_CLEAR']['Linear'])
        return_dict['LAND_CLOUD']   = np.nanmean(extract_dict['LAND_CLOUD']['Linear'])
        return_dict['LAND_CLEAR']   = np.nanmean(extract_dict['LAND_CLEAR']['Linear'])
        
        return return_dict

# dtype: 'clear', or 'cloud'
def calculate_interp_forcings(coloc_dict, month_idx, minlat, maxlat, \
        dtype, cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw'):

    ocean_slopes = calc_slope_clear_clean_sfctype(coloc_dict, 0, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
    ice_slopes   = calc_slope_clear_clean_sfctype(coloc_dict, 1, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
    land_slopes  = calc_slope_clear_clean_sfctype(coloc_dict, 2, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
    if(coloc_dict['raw_slopes'].shape[0] == 4):
        include_mix = True
        mix_slopes  = calc_slope_clear_clean_sfctype(coloc_dict, 3, cld_idx, \
            maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
    else:
        include_mix = False

    # Calculate the average declination angles for each month.
    # The resulting array is of shape (num_months, num_latitudes)
    # -----------------------------------------------------------
    days = np.arange(1, 366)
    del_angles = -23.45 * np.cos( np.radians((360 / 365) * (days + 10)))
    base_date = datetime(2005, 1, 1)
    months = np.array(\
        [(base_date + timedelta(days = int(ii - 1))).month \
        for ii in days])
    latitudes = np.arange(minlat + 0.5, maxlat + 0.5, 1.0)
    lat_szas = np.array([[\
        tlat - np.mean(del_angles[np.where(months == ii)]) \
        for tlat in latitudes] \
        for ii in range(1,13)])[4:10]

    # Calculate the SZA bins from the combined arctic comp output   
    # dictionary. This dictionary contains three sub-dictionaries, 
    # which are the 
    # --------------------------------------------------------------
    sza_bins = np.nanmean(np.array([coloc_dict['sza_mins'], \
                                    coloc_dict['sza_maxs']]), axis = 0)
    
    # NOTE: FOR NOW, SETTING ANY LATITUDE-CALCULATED SZA VALUES THAT
    #       LIE OUTSIDE OF THE SZA BIN VALUES TO EITHER THE MINIMUM
    #       OR MAXIMUM SZA VALUE. ANOTHER OPTION IS TO MASK LATITUDE
    #       SZA VALUES THAT ARE OUTSIDE OF THIS RANGE.
    # --------------------------------------------------------------
    lat_szas = np.where( lat_szas < np.min(sza_bins), np.min(sza_bins), \
        lat_szas)
    lat_szas = np.where( lat_szas > np.max(sza_bins), np.max(sza_bins), \
        lat_szas)

    # Interpolate the bin SZAs to the calculated solar zenith angles
    # Extract the interpolated bin forcing efficiency values
    # ------------------------------------------------------
    land_forcing_values  = interp1d(sza_bins, \
        land_slopes['sza_' + dtype + '_means'])(lat_szas[month_idx])
    ocean_forcing_values = interp1d(sza_bins, \
        ocean_slopes['sza_' + dtype + '_means'])(lat_szas[month_idx])
    ice_forcing_values   = interp1d(sza_bins, \
        ice_slopes['sza_' + dtype + '_means'])(lat_szas[month_idx])
    if(include_mix):
        mix_forcing_values  = interp1d(sza_bins, \
            mix_slopes['sza_clear_means'])(lat_szas[month_idx])
      
    interp_dict = {}
    interp_dict['land_forcing']  = land_forcing_values
    interp_dict['ocean_forcing'] = ocean_forcing_values 
    interp_dict['ice_forcing']  = ice_forcing_values 
    if(include_mix):
        interp_dict['mix_forcing'] = mix_forcing_values

    return interp_dict

# NSIDC_dict: also output from NSIDC library
# dtype: 'cloud', or 'clear'
# data_type: 'raw' or 'grid'
def calculate_type_forcing(OMI_data, NSIDC_data, month_idx, dtype, minlat = 70.,\
        maxlat = 87., debug = False, use_szas = False, coloc_dict = None, \
        cld_idx = 0, min_cloud = 0.95, maxerr = 2, data_type = 'raw'): 

    # Calculate the type values for this month
    # ----------------------------------------
    type_vals = np.array([[\
        check_type_change(NSIDC_data, month_idx, ii, jj, \
        min_ice_for_ice = 80., \
        max_ice_for_ocean = 20., 
        bin_size = 6) \
        for jj in range(NSIDC_data['grid_lon'].shape[1])] \
        for ii in range(NSIDC_data['grid_lon'].shape[0])])

    ai_trends, ai_pvals, ai_uncert = \
        calcOMI_grid_trend(OMI_data, month_idx, 'standard', \
        minlat)

    include_mix = True
    if(coloc_dict is not None):

        interp_dict = calculate_interp_forcings(coloc_dict, month_idx, \
            minlat, maxlat, dtype, cld_idx = cld_idx, maxerr = maxerr, \
            min_cloud = min_cloud, data_type = data_type)

        land_forcing_values  = interp_dict['land_forcing']
        ocean_forcing_values = interp_dict['ocean_forcing']
        ice_forcing_values   = interp_dict['ice_forcing']
        if('mix_forcing' in interp_dict.keys()):
            include_mix = True
            mix_forcing_values = interp_dict['mix_forcing']
        else:
            mix_forcing_values = interp_dict['ice_forcing']
            include_mix = False

    else:
        # NOTE: THESE VALUES ARE ROUGH ESTIMATES!
        # 
        # NOTE: ALSO, USING THE ICE FORCING VALUE FOR BOTH 'Ice' AND 'Mix'
        # ----------------------------------------------------------------
        if(dtype == 'clear'):
            land_forcing = 5
            ocean_forcing = 5
            mix_forcing = -10
            mix_forcing = -4
            ice_forcing = -10
        elif(dtype == 'cloud'): 
            # cloud values
            land_forcing = -5
            ocean_forcing = -5
            mix_forcing = -12
            mix_forcing = -10
            ice_forcing = -12
        else:
            print("WARNING: INCORRECT DATA TYPE")
            return 0

        land_forcing_values  = np.full(OMI_data['LAT'].shape[0], land_forcing)
        ocean_forcing_values = np.full(OMI_data['LAT'].shape[0], ocean_forcing)
        ice_forcing_values   = np.full(OMI_data['LAT'].shape[0], ice_forcing)
        mix_forcing_values   = np.full(OMI_data['LAT'].shape[0], mix_forcing)
 
    estimate_forcings = np.full(ai_trends.shape, np.nan)
    
    for ii in range(OMI_data['LAT'].shape[0]):
        local_land_forcing  = land_forcing_values[ii]
        local_ocean_forcing = ocean_forcing_values[ii]
        local_ice_forcing   = ice_forcing_values[ii]
        if(include_mix):
            local_mix_forcing   = mix_forcing_values[ii]
        for jj in range(OMI_data['LAT'].shape[1]):
            # Grab the curent NSIDC type flag
            # -------------------------------
            local_type = type_vals[ii,jj]
    
            # Land or coastline
            # -----------------
            if((local_type == 15) | (local_type == 14)):
                if(debug): print(ii,jj,"Land",local_land_forcing, \
                    np.round(ai_trends[ii,jj] * local_land_forcing, 3))
                estimate_forcings[ii,jj] = \
                    np.round(ai_trends[ii,jj] * local_land_forcing, 3)
     
            # Pure Ocean
            # ----------
            elif(((local_type == 0) | (local_type == 3))):
                if(debug): print(ii,jj,"Ocean",local_ocean_forcing, \
                    np.round(ai_trends[ii,jj] * local_ocean_forcing, 3))
                estimate_forcings[ii,jj] = \
                    np.round(ai_trends[ii,jj] * local_ocean_forcing, 3)
    
            # Pure Ice
            # --------
            elif((local_type == 2) | (local_type == 9)):
                if(debug): print(ii,jj,"Ice",local_ice_forcing, \
                    np.round(ai_trends[ii,jj] * local_ice_forcing, 3))
                estimate_forcings[ii,jj] = \
                    np.round(ai_trends[ii,jj] * local_ice_forcing, 3)

            # Pure Mix   
            # --------
            elif((local_type == 1) | (local_type == 6)):
                if(debug): print(ii,jj,"Mix",local_mix_forcing, \
                    np.round(ai_trends[ii,jj] * local_mix_forcing, 3))
                estimate_forcings[ii,jj] = \
                    np.round(ai_trends[ii,jj] * local_mix_forcing, 3)

    #'Unchanged ocean':       0,   # 0
    #'Unchanged mix':         1,   # 1
    #'Unchanged ice':         2,   # 2 
    #'Still primarily ocean': 3,   # 0
    #'Mix to ocean':          4,   # 3
    #'Ice to ocean':          5,   # 3
    #'Still primarily mix':   6,   # 1
    #'Ocean to mix':          7,   # 4
    #'Ice to mix':            8,   # 4
    #'Still primarily ice':   9,   # 2
    #'Ocean to ice':         10,   # 5
    #'Mix to ice':           11,   #
    #'Pole Hole':            12,
    #'Unused':               13,
    #'Coastline':            14,
    #'Land':                 15,
    #'Other':                -1,

            # Change: ice to ocean (or ocean to ice)
            # --------------------------------------
            elif((local_type == 5) | (local_type == 10)):
                num_ice   = len(np.where(NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj] >= 80)[0])
                num_ocean = len(np.where(NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj] < 20)[0])
 
                weight_forcing = local_ice_forcing * (num_ice / (num_ocean + num_ice)) + \
                                 local_ocean_forcing * (num_ocean / (num_ocean + num_ice))

                if(debug): print(ii,jj,"Change Ice/Ocean", weight_forcing, num_ice,\
                    num_ocean, np.round(ai_trends[ii,jj] * weight_forcing, 3))
                estimate_forcings[ii,jj] = np.round(ai_trends[ii,jj] * weight_forcing, 3)

            # Change: mix to ocean (or ocean to mix)
            # --------------------------------------
            elif((local_type == 4) | (local_type == 7)):
                num_mix   = len(np.where( (NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj] >= 20) & \
                                          (NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj] < 80))[0])
                num_ocean = len(np.where(NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj] < 20)[0])
 
                weight_forcing = local_mix_forcing * (num_mix / (num_ocean + num_mix)) + \
                                 local_ocean_forcing * (num_ocean / (num_ocean + num_mix))

                if(debug): print(ii,jj,"Change Mix/Ocean", weight_forcing, num_mix,\
                    num_ocean, np.round(ai_trends[ii,jj] * weight_forcing, 3))
                estimate_forcings[ii,jj] = np.round(ai_trends[ii,jj] * weight_forcing, 3)

            # Change: ice to mix (or mix to ice)
            # --------------------------------------
            elif((local_type == 8) | (local_type == 11)):
                num_mix   = len(np.where( (NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj] >= 20) & \
                                          (NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj] < 80))[0])
                num_ice = len(np.where(NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj] >= 80)[0])
 
                weight_forcing = local_mix_forcing * (num_mix / (num_ice + num_mix)) + \
                                 local_ice_forcing * (num_ice / (num_ice + num_mix))

                if(debug): print(ii,jj,"Change Ice/Mix", weight_forcing, num_ice,\
                    num_mix, np.round(ai_trends[ii,jj] * weight_forcing, 3))
                estimate_forcings[ii,jj] = np.round(ai_trends[ii,jj] * weight_forcing, 3)
    
            else:
                if(debug): print(ii,jj, "NOT HANDLED", local_type)

    return estimate_forcings

# 
def calculate_type_forcing_v2(OMI_data, NSIDC_data, MYD08_data, coloc_dict, \
        month_idx, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw',\
        filter_bad_vals = True):
    
    # Necessary pieces:
    # 2. Monthly AI values
    # 4. NSIDC surface type classification for each individual month
    # 5. MODIS monthly cloud fraction 
    #    (/home/bsorenson/data/MODIS/Aqua/MYD08/modis_MYD08_subset_YYYYMM.nc)
    # 6. Arctic comp colocation output dictionary
    # 1. Baseline OMI AI clear-sky climatology (for each month (April, May, etc.))
    # 3. Threshold monthly AI value below which is clear-sky (no forcing) and
    #       above which is aerosol-sky 

    clear_sky_AI = np.array([np.nanmean(\
        np.ma.masked_where(OMI_data['AI'][tidx::6,:,:] > ai_thresh,\
        OMI_data['AI'][tidx::6,:,:]), axis = 0) for tidx in range(6)])
    #clear_sky_AI = np.array([np.nanmean(\
    #    np.ma.masked_where(OMI_data['AI'][tidx::6,:,:] > ai_thresh,\
    #    OMI_data['AI'][tidx::6,:,:]), axis = 0) for tidx in range(6)])
    
    clear_dict = calculate_interp_forcings(coloc_dict, month_idx, \
        minlat, maxlat, 'clear', cld_idx = cld_idx, maxerr = maxerr, \
        min_cloud = min_cloud, data_type = data_type)
    cloud_dict = calculate_interp_forcings(coloc_dict, month_idx, \
        minlat, maxlat, 'cloud', cld_idx = cld_idx, maxerr = maxerr, \
        min_cloud = min_cloud, data_type = data_type)

    land_mask   = NSIDC_data['grid_land'][month_idx::6,:,:].mask
    coast_mask  = NSIDC_data['grid_coastline'][month_idx::6,:,:].mask
    pole_mask   = NSIDC_data['grid_pole_hole'][month_idx::6,:,:].mask
    unused_mask = NSIDC_data['grid_unused'][month_idx::6,:,:].mask
    ice_mask    = NSIDC_data['grid_ice_conc'][month_idx::6,:,:].mask

    estimate_forcings = np.full(OMI_data['AI'][month_idx::6,:,:].shape, np.nan)

    # 2.  Loop over each individual month
    # -----------------------------------
    for nn in range(OMI_data['AI'][month_idx::6,:,:].shape[0]):
        for ii in range(OMI_data['AI'].shape[1]):
            for jj in range(OMI_data['AI'].shape[2]):
        
                # 3.  Determine if a single gridbox has AI above the threshold
                if(OMI_data['AI'][month_idx::6,ii,jj][nn] < ai_thresh): 
                    # 3a. If not above the threshold, set forcing to 0 and continue
                    estimate_forcings[nn,ii,jj] = 0.
                else:
                    # 3b. If above the threshold,
                    # 4.  Determine the difference between this pixel's AI and the clear-sky
                    #     climatology.
                    delta_ai = OMI_data['AI'][month_idx::6,ii,jj][nn] - \
                               clear_sky_AI[nn % 6, ii,jj]

                    # 6.  Extract the MODIS MYD08 value and weigh the "clear" and "cloud"
                    #     forcing values according to that cloud fraction.
                    cld_frac = MYD08_data['cld_frac_mean'][month_idx::6,ii,jj][nn]

                    # 5.  Extract the NSIDC surface type to figure out which forcing value
                    #     to use.
                    if( not ((pole_mask[nn,ii,jj]) | (unused_mask[nn,ii,jj]))):
                        # Bad grid points. NO forcing
                        calc_forcing =  0.

                    elif( not ((land_mask[nn,ii,jj]) | (coast_mask[nn,ii,jj]))):
                        # Use land forcing value
                        calc_forcing = cld_frac * cloud_dict['land_forcing'][ii] + \
                                       (1 - cld_frac) * clear_dict['land_forcing'][ii]

                    elif( not ice_mask[nn,ii,jj]):
                        # Use either ice, mix, or ocean forcing value

                        if( (NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj][nn] < 20) ):
                            # Use ocean forcing
                            calc_forcing = cld_frac * cloud_dict['ocean_forcing'][ii] + \
                                           (1 - cld_frac) * clear_dict['ocean_forcing'][ii]
                            

                        elif( (NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj][nn] >= 20)  & \
                              (NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj][nn] < 80)):
                            # Use mix forcing
                            calc_forcing = cld_frac * cloud_dict['mix_forcing'][ii] + \
                                           (1 - cld_frac) * clear_dict['mix_forcing'][ii]

                        elif( (NSIDC_data['grid_ice_conc'][month_idx::6,ii,jj][nn] > 80) ):
                            # Use land forcing
                            calc_forcing = cld_frac * cloud_dict['ice_forcing'][ii] + \
                                           (1 - cld_frac) * clear_dict['ice_forcing'][ii]

                        else:
                            # SHOULD NOT GET HERE
                            calc_forcing = 0.
 
                    # 8.  Multiply the ΔAI by the associated forcing value.
                    # NOTE: Modifying to actually calculate it as a true forcing,
                    #       which is "clear-sky - aerosol-sky". Assuming that
                    #       "clear-sky" forcing would give 0 W/m2 of aerosol
                    #       forcing, calculate as 0 - calc_forcing * delta_ai
                    estimate_forcings[nn,ii,jj] = 0 - \
                        calc_forcing * delta_ai

    estimate_forcings = np.ma.masked_invalid(estimate_forcings) 

    if(filter_bad_vals):
        mean_val = np.nanmean(estimate_forcings)
        std_val  = np.nanstd(estimate_forcings)

        estimate_forcings = np.ma.masked_where(estimate_forcings > \
            (mean_val + 8.0 * std_val), estimate_forcings)

    return estimate_forcings


# This verison is set up to use daily data, but still needs the monthly
# OMI data to determine the clear-sky climatology.
# ---------------------------------------------------------------------
def calculate_type_forcing_v3(OMI_daily_data, OMI_monthly_data, coloc_dict, \
        date_str, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw',\
        filter_bad_vals = True, return_modis_nsidc = True, debug = False):
    
    # Necessary pieces:
    # 2. Monthly AI values
    # 4. NSIDC surface type classification for each individual month
    # 5. MODIS monthly cloud fraction 
    #    (/home/bsorenson/data/MODIS/Aqua/MYD08/modis_MYD08_subset_YYYYMM.nc)
    # 6. Arctic comp colocation output dictionary
    # 1. Baseline OMI AI clear-sky climatology (for each month (April, May, etc.))
    # 3. Threshold monthly AI value below which is clear-sky (no forcing) and
    #       above which is aerosol-sky 

    # Load in the daily MODIS and NSIDC data
    # --------------------------------------
    file_strs = np.array([str(tval) for tval in OMI_daily_data['day_values']])
    match_idx = np.where(date_str == file_strs)[0][0]
    local_OMI_daily = np.ma.masked_where(\
        OMI_daily_data['count_AI'][match_idx,:,:] == 0, \
        OMI_daily_data['grid_AI'][match_idx,:,:])
    
    # tidx is the "month_idx"
    # -----------------------    
    tidx = int(date_str[4:6]) - 4

    MYD08_data = read_MODIS_MYD08_single(date_str, minlat = minlat, \
        maxlat = maxlat)
   
    print("HERE:", np.max(MYD08_data['cld_frac_mean']))
 
    # Load in the single-day NSIDC ice concentration
    NSIDC_data =  readNSIDC_daily(date_str, maxlat = maxlat)
    NSIDC_data = grid_data_conc(NSIDC_data, minlat = minlat, maxlat = maxlat)
    NSIDC_data['grid_ice_conc'] = np.ma.masked_where((NSIDC_data['grid_ice_conc'] < 0) | \
        (NSIDC_data['grid_ice_conc'] > 100), NSIDC_data['grid_ice_conc']).squeeze()

    clear_sky_AI = np.array([np.nanmean(\
        np.ma.masked_where(OMI_monthly_data['AI'][tidx::6,:,:] > ai_thresh,\
        OMI_monthly_data['AI'][tidx::6,:,:]), axis = 0) for tidx in range(6)])
    clear_sky_AI = clear_sky_AI[tidx,:,:]
    #clear_sky_AI = np.array([np.nanmean(\
    #    np.ma.masked_where(OMI_data['AI'][tidx::6,:,:] > ai_thresh,\
    #    OMI_data['AI'][tidx::6,:,:]), axis = 0) for tidx in range(6)])
    
    clear_dict = calculate_interp_forcings(coloc_dict, tidx, \
        minlat, maxlat, 'clear', cld_idx = cld_idx, maxerr = maxerr, \
        min_cloud = min_cloud, data_type = data_type)
    cloud_dict = calculate_interp_forcings(coloc_dict, tidx, \
        minlat, maxlat, 'cloud', cld_idx = cld_idx, maxerr = maxerr, \
        min_cloud = min_cloud, data_type = data_type)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # NOTE: the 'ii' latitude subscript on the cloud and clear forcing
    # efficiency values in the loops below is correct: in the 
    # 'calculate_interp_forcings' function, the forcing efficiency 
    # values are calculated as functions of 1-degree latitude bins,
    # which the OMI daily data here are gridded into as well.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    land_mask   = NSIDC_data['grid_land'][:,:].mask
    coast_mask  = NSIDC_data['grid_coastline'][:,:].mask
    pole_mask   = NSIDC_data['grid_pole_hole'][:,:].mask
    unused_mask = NSIDC_data['grid_unused'][:,:].mask
    ice_mask    = NSIDC_data['grid_ice_conc'][:,:].mask

    estimate_forcings = np.full(local_OMI_daily.shape, np.nan)

    if(debug):
        print("OMI SHAPE = ", local_OMI_daily.shape)
        print("MOD SHAPE = ", MYD08_data['day_cld_frac_mean'].shape)
        print("NSI SHAPE = ", NSIDC_data['grid_ice_conc'].shape)

    # 2.  Loop over each individual month
    # -----------------------------------
    for ii in range(local_OMI_daily.shape[0]):
        for jj in range(local_OMI_daily.shape[1]):
    
            # 3.  Determine if a single gridbox has AI above the threshold
            if(local_OMI_daily[ii,jj] < ai_thresh): 
                # 3a. If not above the threshold, set forcing to 0 and continue
                estimate_forcings[ii,jj] = 0.
            else:
                # 3b. If above the threshold,
                # 4.  Determine the difference between this pixel's AI and the clear-sky
                #     climatology.
                delta_ai = local_OMI_daily[ii,jj] - \
                           clear_sky_AI[ii,jj]

                # 6.  Extract the MODIS MYD08 value and weigh the "clear" and "cloud"
                #     forcing values according to that cloud fraction.
                cld_frac = MYD08_data['day_cld_frac_mean'][ii,jj]

                # 5.  Extract the NSIDC surface type to figure out which forcing value
                #     to use.
                if( (pole_mask[ii,jj] == True) & (unused_mask[ii,jj] == True) & \
                    (land_mask[ii,jj] == True) & (coast_mask[ii,jj] == True) & \
                    (ice_mask[ii,jj] == True)):
                    if(debug):
                        print(pole_mask[ii,jj], unused_mask[ii,jj], \
                            land_mask[ii,jj], coast_mask[ii,jj], \
                            ice_mask[ii,jj], 'ALL MASKED')
                    calc_forcing = 0.
                else:
                    if((pole_mask[ii,jj] == False) | \
                       (unused_mask[ii,jj] == False)):
                        # Bad grid points. NO forcing
                        calc_forcing =  0.
                        if(debug):
                            print(pole_mask[ii,jj], unused_mask[ii,jj], \
                                land_mask[ii,jj], coast_mask[ii,jj], \
                                ice_mask[ii,jj], 'POLE/UNUSED')

                    elif((land_mask[ii,jj] == False) | \
                         (coast_mask[ii,jj] == False)):
                        # Use land forcing value
                        if(debug):
                            print(pole_mask[ii,jj], unused_mask[ii,jj], \
                                land_mask[ii,jj], coast_mask[ii,jj], \
                                ice_mask[ii,jj], 'LAND/COAST')

                        # 2023/10/05 - added check to see if the smoke
                        #   is above permanent ice (or above dry snow).
                        #   If it is, use the ice forcing values. 
                        #   Although, is it fair to do this since dry
                        #   snow was removed from the initial forcing
                        #   efficiency calculations...
                        #
                        # Check if the current grid box is permanent ice. 
                        # In that case, use the ice forcing eff. values
                        # -----------------------------------------------
                        if(OMI_daily_data['grid_GPQF'][match_idx,ii,jj] == 3):
                            calc_forcing = cld_frac * \
                                cloud_dict['ice_forcing'][ii] + \
                                (1 - cld_frac) * clear_dict['ice_forcing'][ii]
                        else:
                            calc_forcing = cld_frac * \
                                cloud_dict['land_forcing'][ii] + \
                                (1 - cld_frac) * clear_dict['land_forcing'][ii]
                        estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 

                    elif(ice_mask[ii,jj] == False):
                        # Use either ice, mix, or ocean forcing value
                        if(debug):
                            print(pole_mask[ii,jj], unused_mask[ii,jj], \
                                land_mask[ii,jj], coast_mask[ii,jj], \
                                ice_mask[ii,jj], 'ICE/MIX/OCEAN')

                        if( (NSIDC_data['grid_ice_conc'][ii,jj] < 20) ):
                            # Use ocean forcing
                            calc_forcing = cld_frac * cloud_dict['ocean_forcing'][ii] + \
                                           (1 - cld_frac) * clear_dict['ocean_forcing'][ii]
                            estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 
                            

                        elif( (NSIDC_data['grid_ice_conc'][ii,jj] >= 20)  & \
                              (NSIDC_data['grid_ice_conc'][ii,jj] < 80)):
                            # Use mix forcing
                            calc_forcing = cld_frac * cloud_dict['mix_forcing'][ii] + \
                                           (1 - cld_frac) * clear_dict['mix_forcing'][ii]
                            estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 

                        elif( (NSIDC_data['grid_ice_conc'][ii,jj] > 80) ):
                            # Use land forcing
                            calc_forcing = cld_frac * cloud_dict['ice_forcing'][ii] + \
                                           (1 - cld_frac) * clear_dict['ice_forcing'][ii]
                            estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 

                        else:
                            if(debug):
                                print("FAILED ICE")
                            calc_forcing = 0.

                    else:
                        if(debug):
                            print("FAILED EVERYTHING")
                        calc_forcing = 0.
 
                # 8.  Multiply the ΔAI by the associated forcing value.
                # NOTE: Modifying to actually calculate it as a true forcing,
                #       which is "clear-sky - aerosol-sky". Assuming that
                #       "clear-sky" forcing would give 0 W/m2 of aerosol
                #       forcing, calculate as 0 - calc_forcing * delta_ai
                #estimate_forcings[ii,jj] = 0 - \
                #    calc_forcing * delta_ai

    estimate_forcings = np.ma.masked_invalid(estimate_forcings) 

    if(filter_bad_vals):
        mean_val = np.nanmean(estimate_forcings)
        std_val  = np.nanstd(estimate_forcings)

        estimate_forcings = np.ma.masked_where(estimate_forcings > \
            (mean_val + 8.0 * std_val), estimate_forcings)

    if(return_modis_nsidc):
        return estimate_forcings, MYD08_data, NSIDC_data
    else:
        return estimate_forcings

# This verison is set up to use daily data, but still needs the monthly
# OMI data to determine the clear-sky climatology.
# ---------------------------------------------------------------------
def calculate_type_forcing_v3_monthly(OMI_daily_data, OMI_monthly_data, \
        coloc_dict, month_idx, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw',\
        filter_bad_vals = True, return_modis_nsidc = True, debug = False):

    # Set up arrays to hold the monthly-averaged forcing calculations
    # ---------------------------------------------------------------
    xdim = OMI_daily_data['grid_AI'].shape[1]
    ydim = OMI_daily_data['grid_AI'].shape[2]
    daily_force_vals = np.full( (31, xdim, ydim), np.nan)

    month_force_vals = np.full(\
        (OMI_monthly_data['DATES'][::6].shape[0], xdim, ydim), \
        np.nan)

    begin_date_str = datetime.strptime(\
        OMI_monthly_data['DATES'][month_idx], '%Y%m')
    end_date_str = datetime.strptime(\
        OMI_monthly_data['DATES'][month_idx::6][-1], '%Y%m') + \
        relativedelta(months = 1) - timedelta(days = 1)

    local_date_str = begin_date_str

    day_count = 0
    month_count = 0
    while(local_date_str <= end_date_str):

        date_str = local_date_str.strftime('%Y%m%d')

        print(date_str, day_count, month_count)

        # Calculate the forcing value for this current day
        # ------------------------------------------------
        ##!#estimate_forcing = \
        ##!#    calculate_type_forcing_v3(OMI_daily_data, OMI_monthly_data, \
        ##!#        coloc_dict, date_str, minlat = minlat, maxlat = maxlat, \
        ##!#        ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr,\
        ##!#        min_cloud = min_cloud, data_type = data_type,\
        ##!#        filter_bad_vals = filter_bad_vals, \
        ##!#        return_modis_nsidc = False)

        ##!## Insert the values into the daily holding array
        ##!## ----------------------------------------------
        ##!#daily_force_vals[day_count,:,:] = estimate_forcing[:,:] 

        # Increment working date
        # ----------------------
        new_work_date = local_date_str + timedelta(days = 1)

        # If the new working date has a different month than
        # the previous working date, then average the daily values
        # and move the working date to the next desired month
        # --------------------------------------------------------
        if(new_work_date.month != local_date_str.month):
            ##!#month_force_vals[month_count,:,:] = np.nanmean(\
            ##!#    daily_force_vals[:,:,:], axis = 0)
            day_count = 0
            month_count += 1            

            new_work_date = new_work_date + relativedelta(months = 11)
            #new_work_date = datetime.strptime(
            #    OMI_monthly_data['DATES'][month_idx::6][month_count], \
            #    '%Y%m')

            daily_force_vals[:,:,:] = np.nan 
 
        # If the new working date has the same month as the 
        # previous working date, just increment the day counter
        # and the date string and continue
        # -----------------------------------------------------
        else:
            day_count += 1

        local_date_str = new_work_date

    return month_force_vals

def calc_forcing_grid_trend(forcing_data, trend_type):

    index_jumper = 6

    # Make copy of NSIDC_data array
    local_data   = np.copy(forcing_data)
    local_data = np.ma.masked_invalid(local_data)
    #local_mask = np.ma.masked_where((local_data < 0., local_data)
    local_mask = local_data
    forcing_trends = np.full(local_data.shape[1:], np.nan)
    forcing_pvals  = np.full(local_data.shape[1:], np.nan)
    forcing_uncert = np.full(local_data.shape[1:], np.nan)

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    for i in range(local_data.shape[1]):
        for j in range(local_data.shape[2]):
            # Check the current max and min
            #print(local_mask[:,i,j])
            work_mask = local_mask[:,i,j]
            #work_mask = local_mask[:,i,j][~local_mask[:,i,j].mask][0]
            if(len(work_mask.compressed()) > 1):
                x_vals = np.arange(0,len(work_mask.compressed()))
                # Find the slope of the line of best fit for the time series of
                # average data
                if(trend_type=='standard'): 
                    result = stats.linregress(x_vals, work_mask.compressed())
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,work_mask.compressed())
                    forcing_trends[i,j] = result.slope * len(x_vals)
                    forcing_pvals[i,j]  = result.pvalue
                    forcing_uncert[i,j] = result.stderr * len(x_vals)
                else:
                    res = stats.theilslopes(work_mask.compressed(), x_vals, 0.90)
                    forcing_trends[i,j] = res[0]*len(x_vals)
            #else:
            #    print('no data')

    #nsidc_trends = np.ma.masked_where(((NSIDC_data['grid_lat'] < mingrid_lat) | \
    #    (nsidc_trends == -999.)), nsidc_trends)
    #nsidc_trends = np.ma.masked_where(NSIDC_data['grid_lat'] < mingrid_lat, nsidc_trends)
    #nsidc_pvals  = np.ma.masked_where(NSIDC_data['grid_lat'] < mingrid_lat, nsidc_pvals)
    #nsidc_uncert = np.ma.masked_where(NSIDC_data['grid_lat'] < mingrid_lat, nsidc_uncert)

    return forcing_trends, forcing_pvals, forcing_uncert


def plot_test_forcing_v2(OMI_data, NSIDC_data, MYD08_data, coloc_dict, \
        tidx, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
        save = False):

    estimate_forcing = \
        calculate_type_forcing_v2(OMI_data, NSIDC_data, MYD08_data, \
        coloc_dict, tidx % 6, minlat = minlat, maxlat = maxlat, \
        ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr, \
        min_cloud = min_cloud, data_type = data_type)

    clear_sky_AI = np.array([np.nanmean(\
        np.ma.masked_where(OMI_data['AI'][midx::6,:,:] > ai_thresh,\
        OMI_data['AI'][midx::6,:,:]), axis = 0) for midx in range(6)])

    local_departure = np.where(OMI_data['AI'][tidx,:,:] >= ai_thresh, \
        OMI_data['AI'][tidx,:,:] - clear_sky_AI[tidx % 6,:,:], \
        np.nan)
    local_departure = np.ma.masked_invalid(local_departure)
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 5))
    ax1 = fig.add_subplot(2,3,1, projection = mapcrs)  # Map of original AI
    ax2 = fig.add_subplot(2,3,2, projection = mapcrs)  # Map of climatology departure
    ax3 = fig.add_subplot(2,3,3, projection = mapcrs)  # Map of departure
    ax4 = fig.add_subplot(2,3,4, projection = mapcrs)  # Map of cloud fraction
    ax5 = fig.add_subplot(2,3,5, projection = mapcrs)  # Map of forcing value
    ax6 = fig.add_subplot(2,3,6)  # histogram

    mesh = ax1.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
        OMI_data['AI'][tidx,:,:], shading = 'auto', transform = datacrs, \
        cmap = 'jet', vmin = 0, vmax = 0.5)
    cbar = fig.colorbar(mesh, ax = ax1, label =  'Monthly AI')
    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.coastlines()

    mesh = ax2.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
        clear_sky_AI[tidx % 6,:,:], shading = 'auto', transform = datacrs, \
        cmap = 'jet', vmin = 0, vmax = 0.5)
    cbar = fig.colorbar(mesh, ax = ax2, label = 'AI Climatology')
    ax2.set_extent([-180,180,minlat,90], datacrs)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.coastlines()

    mesh = ax3.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
        local_departure[:,:], shading = 'auto', transform = datacrs, \
        cmap = 'jet', vmin = 0.0, vmax = 0.5)
    cbar = fig.colorbar(mesh, ax = ax3, label = 'ΔAI (AI$_{i}$ - AI$_{clim}$)')
    ax3.set_extent([-180,180,minlat,90], datacrs)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.coastlines()

    mesh = ax4.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
        MYD08_data['cld_frac_mean'][tidx,:,:], shading = 'auto', transform = datacrs, \
        cmap = 'viridis', vmin = 0.0, vmax = 1.0)
    cbar = fig.colorbar(mesh, ax = ax4, label = 'MODIS Cloud Frac')
    ax4.set_extent([-180,180,minlat,90], datacrs)
    ax4.set_boundary(circle, transform=ax4.transAxes)
    ax4.coastlines()

    work_idx = np.where(\
        np.arange(tidx % 6, OMI_data['AI'].shape[0], 6) == tidx)[0][0]
    mesh = ax5.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
        estimate_forcing[work_idx,:,:], shading = 'auto', transform = datacrs, \
        cmap = 'viridis')
    cbar = fig.colorbar(mesh, ax = ax5, label = 'Estimated Forcing [W/m2]')
    ax5.set_extent([-180,180,minlat,90], datacrs)
    ax5.set_boundary(circle, transform=ax5.transAxes)
    ax5.coastlines()

    ax6.hist(np.ma.masked_where(estimate_forcing[work_idx,:,:] == 0, \
        estimate_forcing[work_idx,:,:]).compressed(), bins = 'auto')
    #ax6.set_yscale('log')
    ax6.set_xlabel('Estimated Forcing [W/m2]')
    ax6.set_ylabel('Counts')

    plt.suptitle(OMI_data['DATES'][tidx])

    fig.tight_layout()
    if(save):
        outname = 'test_calc_forcing_v2_' + OMI_data['DATES'][tidx] + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()
 
def plot_test_forcing_v3(OMI_daily_data, OMI_month_data, date_str, \
        coloc_dict, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
        save = False, filter_bad_vals = False):

    estimate_forcing, MYD08_data, NSIDC_data = \
        calculate_type_forcing_v3(OMI_daily_data, OMI_month_data, \
            coloc_dict, date_str, minlat = minlat, maxlat = maxlat, \
            ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr,\
            min_cloud = min_cloud, data_type = data_type,\
            filter_bad_vals = filter_bad_vals, return_modis_nsidc = True)

    file_strs = np.array([str(tval) for tval in OMI_daily_data['day_values']])
    match_idx = np.where(date_str == file_strs)[0][0]
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 5))
    ax1 = fig.add_subplot(2,3,1, projection = mapcrs)  # Map of original AI
    ax4 = fig.add_subplot(2,3,2, projection = mapcrs)  # Map of cloud fraction
    ax5 = fig.add_subplot(2,3,3, projection = mapcrs)  # Map of ice concentration
    ax2 = fig.add_subplot(2,3,4, projection = mapcrs)  # Map of forcing values
    ax6 = fig.add_subplot(2,3,5, projection = mapcrs)  # Map of GPQF values
    ax3 = fig.add_subplot(2,3,6)  # histogram

    mesh = ax1.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], \
        OMI_daily_data['grid_AI'][match_idx,:,:], shading = 'auto', \
        transform = datacrs, cmap = 'jet', vmin = 0, vmax = 4.0)
    cbar = fig.colorbar(mesh, ax = ax1, label =  'Daily AI')
    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.coastlines()

    mesh = ax4.pcolormesh(MYD08_data['lon'], \
        MYD08_data['lat'], MYD08_data['cld_frac_mean'][:,:], shading = 'auto',\
        transform = datacrs, \
        cmap = 'viridis')
    cbar = fig.colorbar(mesh, ax = ax4, label = 'Daily Cloud Fraction')
    ax4.set_extent([-180,180,minlat,90], datacrs)
    ax4.set_boundary(circle, transform=ax4.transAxes)
    ax4.coastlines()

    mesh = ax5.pcolormesh(NSIDC_data['grid_lon'], \
        NSIDC_data['grid_lat'], NSIDC_data['grid_ice_conc'][:,:], shading = 'auto',\
        transform = datacrs, \
        cmap = 'ocean', vmin = 0, vmax = 100)
    cbar = fig.colorbar(mesh, ax = ax5, label = 'Sea Ice Conc')
    ax5.set_extent([-180,180,minlat,90], datacrs)
    ax5.set_boundary(circle, transform=ax5.transAxes)
    ax5.coastlines()

    min_force = np.nanmin(estimate_forcing[:,:])
    max_force = np.nanmax(estimate_forcing[:,:])
    if(abs(max_force) > abs(min_force)):
        lims = [-abs(max_force), abs(max_force)]
    else:
        lims = [-abs(min_force), abs(min_force)]
    mesh = ax2.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], estimate_forcing[:,:], shading = 'auto',\
        transform = datacrs, \
        cmap = 'bwr', vmin = lims[0], vmax = lims[1])
    cbar = fig.colorbar(mesh, ax = ax2, label = 'Estimated Forcing [W/m2]')
    ax2.set_extent([-180,180,minlat,90], datacrs)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.coastlines()

    mesh = ax6.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], \
        OMI_daily_data['grid_GPQF'][match_idx,:,:], shading = 'auto', \
        transform = datacrs, cmap = 'jet', vmin = None, vmax = None)
    cbar = fig.colorbar(mesh, ax = ax6, label =  'Daily GPQF')
    ax6.set_extent([-180,180,minlat,90], datacrs)
    ax6.set_boundary(circle, transform=ax6.transAxes)
    ax6.coastlines()

    ax3.hist(np.ma.masked_where(estimate_forcing[:,:] == 0, \
        estimate_forcing[:,:]).compressed(), bins = 'auto')
    #ax6.set_yscale('log')
    ax3.set_xlabel('Estimated Forcing [W/m2]')
    ax3.set_ylabel('Counts')

    plt.suptitle(date_str)

    fig.tight_layout()
    if(save):
        outname = 'test_calc_forcing_v3_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()
def calculate_type_forcing_old(month_idx, trend_type = 'linear', minlat = 65.):

    # Calculate gridded OMI trends
    OMI_data   = readOMI_NCDF(infile = \
        '/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', \
        minlat = minlat)

    ai_trends, ai_pvals = calcOMI_grid_trend(OMI_data, month_idx, trend_type, \
        minlat)

    ##!#local_data3_Apr  = np.copy(OMI_data3['MONTH_CLIMO'][0,:,:])
    ##!##local_data3_May  = np.copy(OMI_data3['MONTH_CLIMO'][1,:,:])
    ##!##local_data3_Jun  = np.copy(OMI_data3['MONTH_CLIMO'][2,:,:])
    ##!##local_data3_Jul  = np.copy(OMI_data3['MONTH_CLIMO'][3,:,:])
    ##!##local_data3_Aug  = np.copy(OMI_data3['MONTH_CLIMO'][4,:,:])
    ##!##local_data3_Sep  = np.copy(OMI_data3['MONTH_CLIMO'][5,:,:])

    ##!#mask_AI3_Apr = np.ma.masked_where(local_data3_Apr == -999.9, local_data3_Apr)
    ##!##mask_AI3_May = np.ma.masked_where(local_data3_May == -999.9, local_data3_May)
    ##!##mask_AI3_Jun = np.ma.masked_where(local_data3_Jun == -999.9, local_data3_Jun)
    ##!##mask_AI3_Jul = np.ma.masked_where(local_data3_Jul == -999.9, local_data3_Jul)
    ##!##mask_AI3_Aug = np.ma.masked_where(local_data3_Aug == -999.9, local_data3_Aug)
    ##!##mask_AI3_Sep = np.ma.masked_where(local_data3_Sep == -999.9, local_data3_Sep)

    # Calculate the SWF/AI relationships
    relat_dict = plot_compare_all_slopes(date_strs = None, save = True, \
        return_slope_dict = False, return_slope_means = True)

    print(relat_dict)
    print(np.nanmax(ai_trends))

    # Multiply the AI trends by the type trends
    test_trend = 0.5
    for tkey in relat_dict.keys():
        print(tkey, test_trend * relat_dict[tkey]) 

def select_combined_scatter_data(combined_data, xval, yval, 
        swf_min = None, swf_max = None, \
        #cld_min = None, cld_max = None, \
        ch7_min = None, ch7_max = None, \
        sza_min = None, sza_max = None, \
        ice_min = None, ice_max = None, \
        omi_min = None, omi_max = None, \
        ai_diff = None, sza_diff = None, \
        ice_diff = None, ch7_diff = None,\
        cld_diff = None,\
        include_cloud = True, \
        cloud_var = 'ch7'):

    local_xdata = np.copy(combined_data[xval])
    local_ydata = np.copy(combined_data[yval])

    def check_range_vals(local_xdata, local_ydata, lvar,\
            var_min, var_max, var_diff):

        if(var_min is not None):
            if(var_diff is not None):
                checker = var_min - var_diff
            else:
                checker = var_min
            local_xdata = np.ma.masked_where(\
                combined_data[lvar] < checker,local_xdata)
            local_ydata = np.ma.masked_where(\
                combined_data[lvar] < checker,local_ydata)
        if(var_max is not None):
            if(var_diff is not None):
                checker = var_max + var_diff
            else:
                checker = var_max
            local_xdata = np.ma.masked_where(\
                combined_data[lvar] > checker,local_xdata)
            local_ydata = np.ma.masked_where(\
                combined_data[lvar] > checker,local_ydata)

        return local_xdata, local_ydata
        
    #local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
    #    'ceres_cld', cld_min, cld_max, cld_diff)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'modis_' + cloud_var, ch7_min, ch7_max, ch7_diff)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'nsidc_ice', ice_min, ice_max, ice_diff)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'omi_sza', sza_min, sza_max, sza_diff)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'omi_uvai_raw', omi_min, omi_max, ai_diff)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'ceres_swf', swf_min, swf_max, None)

    if(include_cloud):
        local_cloud = np.copy(combined_data['modis_cld'])
        local_cloud = local_cloud[~local_xdata.mask]
    else:
        local_cloud = None   
 
    local_sza = np.copy(combined_data['omi_sza'])
    local_sza = local_sza[~local_xdata.mask]

    local_cld = np.copy(combined_data['modis_' + cloud_var])
    local_cld = local_cld[~local_xdata.mask]

    local_xdata = local_xdata.compressed()
    local_ydata = local_ydata.compressed()

    # Now, remove data that have missing COD/CH7 values
    # -------------------------------------------------
    local_cld = np.ma.masked_invalid(local_cld)
    local_xdata = np.ma.masked_where(local_cld.mask, local_xdata)
    local_ydata = np.ma.masked_where(local_cld.mask, local_ydata)
    local_sza   = np.ma.masked_where(local_cld.mask, local_sza)

    local_xdata = local_xdata[~local_cld.mask]
    local_ydata = local_ydata[~local_cld.mask]
    local_sza   = local_sza[~local_cld.mask]
    local_cld   = local_cld[~local_cld.mask]

    ##!#print("In selecter,")
    ##!#if((ice_max is not None) and (ice_min is not None)):
    ##!#    print("    ice:",ice_min - ice_diff, ice_max + ice_diff, \
    ##!#              "ch7:",ch7_min - ch7_diff, ch7_max + ch7_diff, \
    ##!#              "sza",sza_min - sza_diff, sza_max + sza_diff)
    ##!#elif((ice_max is None) and (ice_min is not None)):
    ##!#    print("    ice:",ice_min - ice_diff, 'None', \
    ##!#              "ch7:",ch7_min - ch7_diff, ch7_max + ch7_diff, \
    ##!#              "sza",sza_min - sza_diff, sza_max + sza_diff)
    ##!#elif((ice_max is not None) and (ice_min is None)):
    ##!#    print("    ice:",'None',ice_max + ice_diff,  \
    ##!#              "ch7:",ch7_min - ch7_diff, ch7_max + ch7_diff, \
    ##!#              "sza",sza_min - sza_diff, sza_max + sza_diff)
    ##!#print("   xdata_size",local_xdata.shape[0])

    out_dict = {}
    out_dict['local_xdata'] = np.round(local_xdata, 4)
    out_dict['local_ydata'] = np.round(local_ydata, 4)
    out_dict['local_sza']   = np.round(local_sza, 4)
    out_dict['local_cld']   = np.round(local_cld, 4)
    out_dict['local_cloud'] = np.round(local_cloud, 4)

    return out_dict
    #return local_xdata, local_ydata, local_cloud

def plot_combined_scatter(combined_data, xval = 'omi_uvai_raw', \
        yval = 'ceres_swf', ax = None, \
        swf_min = None, swf_max = None, \
        #cld_min = None, cld_max = None, \
        cod_min = None, cod_max = None, \
        ch7_min = None, ch7_max = None, \
        sza_min = None, sza_max = None, \
        ice_min = None, ice_max = None, \
        omi_min = None, omi_max = None, \
        ai_diff = None, sza_diff = None, \
        ice_diff = None, ch7_diff = None,\
        cld_diff = None, \
        trend_type = 'theil-sen', \
        show_trend = False, shade_density = False,\
        cloud_var = 'ch7'):

    ##!#local_xdata = np.copy(combined_data[xval])
    ##!#local_ydata = np.copy(combined_data[yval])

    ##!#def check_range_vals(local_xdata, local_ydata, lvar,\
    ##!#        var_min, var_max, var_diff):

    ##!#    if(var_min is not None):
    ##!#        if(var_diff is not None):
    ##!#            checker = var_min - var_diff
    ##!#        else:
    ##!#            checker = var_min
    ##!#        local_xdata = np.ma.masked_where(\
    ##!#            combined_data[lvar] < checker,local_xdata)
    ##!#        local_ydata = np.ma.masked_where(\
    ##!#            combined_data[lvar] < checker,local_ydata)
    ##!#    if(var_max is not None):
    ##!#        if(var_diff is not None):
    ##!#            checker = var_max + var_diff
    ##!#        else:
    ##!#            checker = var_max
    ##!#        local_xdata = np.ma.masked_where(\
    ##!#            combined_data[lvar] > checker,local_xdata)
    ##!#        local_ydata = np.ma.masked_where(\
    ##!#            combined_data[lvar] > checker,local_ydata)

    ##!#    return local_xdata, local_ydata
    ##!#    
    ##!##local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
    ##!##    'ceres_cld', cld_min, cld_max, cld_diff)
    ##!#local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
    ##!#    'modis_ch7', ch7_min, ch7_max, ch7_diff)
    ##!#local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
    ##!#    'nsidc_ice', ice_min, ice_max, ice_diff)
    ##!#local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
    ##!#    'omi_sza', sza_min, sza_max, sza_diff)
    ##!#local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
    ##!#    'omi_uvai_raw', omi_min, omi_max, ai_diff)
    ##!#local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
    ##!#    'ceres_swf', swf_min, swf_max, None)

    ##!#local_xdata = local_xdata.compressed()
    ##!#local_ydata = local_ydata.compressed()

    #local_xdata, local_ydata, local_cloud  = \
    return_dict  = \
        select_combined_scatter_data(combined_data, xval, yval, 
            swf_min = swf_min, swf_max = swf_max, \
            #cld_min = None, cld_max = None, \
            ch7_min = ch7_min, ch7_max = ch7_max, \
            sza_min = sza_min, sza_max = sza_max, \
            ice_min = ice_min, ice_max = ice_max, \
            omi_min = omi_min, omi_max = omi_max, \
            ai_diff = ai_diff, sza_diff = sza_diff, \
            ice_diff = ice_diff, ch7_diff = ch7_diff, \
            #cld_diff = cld_diff, \
            cloud_var = cloud_var 
            )

    local_xdata = return_dict['local_xdata']
    local_ydata = return_dict['local_ydata']
    local_cloud = return_dict['local_cloud']

    if(trend_type == 'theil-sen'):
        if((len(local_xdata) > 20)):
            res = stats.theilslopes(local_ydata, local_xdata, 0.90)
            print("IN HERE: Theil slope:",res[0])

    elif(trend_type == 'linregress'):
        if((len(local_xdata) > 20)):
            result1 = stats.linregress(local_xdata,local_ydata)
            print("IN HERE: Lin slope:",result1.slope)


    print(local_xdata.shape, local_ydata.shape)

    z = 'k'
    if(shade_density):
        if(len(local_xdata) > 10000):
            print("ERROR: Density shading with too many points")
            print("       Preventing this for the sake of computer safety")
        else:
            xy = np.vstack([local_xdata, local_ydata])
            z = stats.gaussian_kde(xy)(xy)       

    in_ax = True
    if(ax is None):
        plt.close('all')
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    
    ax.scatter(local_xdata, local_ydata, c = z, s = 6)
    if(show_trend):
        if(len(local_xdata) > 10000):
            print("ERROR: Plotting trend line with too many points")
            print("       Preventing this for the sake of computer safety")
        else:
            plot_trend_line(ax, local_xdata, local_ydata, color='tab:red', \
                linestyle = '-',  slope = trend_type)

    ax.set_xlabel('OMI UVAI')
    ax.set_ylabel('CERES SWF')

    if(not in_ax):
        plt.show()

    return return_dict

def select_comp_grid_scatter_data(comp_grid_data, xval = 'ai', \
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None):
   
    if('cod_bins' in comp_grid_data.keys()): 
        cloud_var = 'cod'
        ch7_lims  = [np.min(comp_grid_data['cod_bins']), \
            np.max(comp_grid_data['cod_bins'])]
    else:
        cloud_var = 'ch7'
        ch7_lims  = [np.min(comp_grid_data['ch7_bins']), \
            np.max(comp_grid_data['ch7_bins'])]
    #ice_lims  = [101, 105]
    ice_lims  = [np.min(comp_grid_data['ice_bins']), \
        np.max(comp_grid_data['ice_bins'])]
    #sza_lims  = [65, 70]
    sza_lims  = [np.min(comp_grid_data['sza_bins']), \
        np.max(comp_grid_data['sza_bins'])]
    #ai_lims   = [2, 12]
    ai_lims   = [np.min(comp_grid_data['ai_bins']), \
        np.max(comp_grid_data['ai_bins'])]
   
    #if(cld_min is not None):
    #    cld_lims[0] = cld_min
    #if(cld_max is not None):
    #    cld_lims[1] = cld_max
    if(ch7_min is not None):
        ch7_lims[0] = ch7_min
    if(ch7_max is not None):
        ch7_lims[1] = ch7_max
    if(ice_min is not None):
        ice_lims[0] = ice_min
    if(ice_max is not None):
        ice_lims[1] = ice_max
    if(sza_min is not None):
        sza_lims[0] = sza_min
    if(sza_max is not None):
        sza_lims[1] = sza_max
    if(ai_min is not None):
        ai_lims[0] = ai_min
    if(ai_max is not None):
        ai_lims[1] = ai_max
 
    #keep_ai_lims  = np.arange(len(comp_grid_data['ai_bins']))
    keep_ai_lims = np.where((comp_grid_data['ai_bins'] >= ai_lims[0]) & \
                             (comp_grid_data['ai_bins'] <= ai_lims[1]))[0]
    keep_sza_lims = np.where((comp_grid_data['sza_bins'] >= sza_lims[0]) & \
                             (comp_grid_data['sza_bins'] <= sza_lims[1]))[0]
    keep_ice_lims = np.where((comp_grid_data['ice_bins'] >= ice_lims[0]) & \
                             (comp_grid_data['ice_bins'] <= ice_lims[1]))[0]
    keep_ch7_lims = np.where((comp_grid_data[cloud_var + '_bins'] >= ch7_lims[0]) & \
                             (comp_grid_data[cloud_var + '_bins'] <= ch7_lims[1]))[0]
    #keep_cld_lims = np.where((comp_grid_data['cld_bins'] >= cld_lims[0]) & \
    #                         (comp_grid_data['cld_bins'] <= cld_lims[1]))[0]
    #keep_ice_lims  = np.arange(len(comp_grid_data['ice_bins']))
    #keep_ch7_lims  = np.arange(len(comp_grid_data['ch7_bins']))
    
    keep_swf = comp_grid_data['swf_climo'][\
        #keep_cld_lims[0]:keep_cld_lims[-1]+1,\
        keep_ch7_lims[0]:keep_ch7_lims[-1]+1,\
        keep_ice_lims[0]:keep_ice_lims[-1]+1,\
        keep_sza_lims[0]:keep_sza_lims[-1]+1,\
        keep_ai_lims[0]:keep_ai_lims[-1]+1]

    #keep_cld = np.array([[[[[comp_grid_data['cld_bins'][mm] \
    #    for nn in keep_ai_lims] \
    #    for kk in keep_sza_lims] \
    #    for jj in keep_ice_lims] \
    #    for ii in keep_ch7_lims] \
    #    for mm in keep_cld_lims])

    keep_ch7 = np.array([[[[comp_grid_data[cloud_var + '_bins'][ii] \
        for nn in keep_ai_lims] \
        for kk in keep_sza_lims] \
        for jj in keep_ice_lims] \
        for ii in keep_ch7_lims] \
        )
    
    keep_ice = np.array([[[[comp_grid_data['ice_bins'][jj] \
        for nn in keep_ai_lims] \
        for kk in keep_sza_lims] \
        for jj in keep_ice_lims] \
        for ii in keep_ch7_lims] \
        )
    
    keep_sza = np.array([[[[comp_grid_data['sza_bins'][kk] \
        for nn in keep_ai_lims] \
        for kk in keep_sza_lims] \
        for jj in keep_ice_lims] \
        for ii in keep_ch7_lims] \
        )
    
    keep_ai  = np.array([[[[comp_grid_data['ai_bins'][nn] \
        for nn in keep_ai_lims] \
        for kk in keep_sza_lims] \
        for jj in keep_ice_lims] \
        for ii in keep_ch7_lims] \
        )
    


    if(xval == 'ai'): 
        #ax.scatter(final_ai, final_swf, s = 6, color = 'k')
        xlabel = 'OMI UVAI Raw'
        local_xdata = keep_ai.flatten()[~keep_swf.flatten().mask]
    elif(xval == 'sza'):
        #ax.scatter(final_sza, final_swf, s = 6, color = 'k')
        local_xdata = keep_sza.flatten()[~keep_swf.flatten().mask]
        xlabel = 'OMI SZA'
    elif(xval == 'ice'):
        #ax.scatter(final_ice, final_swf, s = 6, color = 'k')
        local_xdata = keep_ice.flatten()[~keep_swf.flatten().mask]
        xlabel = 'NSIDC ICE CONC.'
    elif((xval == 'ch7') | (xval == 'cod')):
        #ax.scatter(final_ch7, final_swf, s = 6, color = 'k')
        local_xdata = keep_ch7.flatten()[~keep_swf.flatten().mask]
        if(xval == 'cod'):
            xlabel = 'MODIS COD'
        else:
            xlabel = 'MODIS CH7 REFL.'
    #elif(xval == 'cld'):
    #    #ax.scatter(final_ch7, final_swf, s = 6, color = 'k')
    #    local_xdata = final_cld
    #    xlabel = 'CERES CLD FRAC.'
    local_ydata = keep_swf.compressed()

    return local_xdata, local_ydata, xlabel

def plot_comp_grid_scatter(comp_grid_data, xval = 'ai', \
        #cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None,\
        ax = None, show_trend = False , \
        trend_type = 'theil-sen', \
        save = False):

    ##!##cld_lims  = [0., 20.]
    ##!##cld_lims  = [np.min(comp_grid_data['cld_bins']), \
    ##!##    np.max(comp_grid_data['cld_bins'])]
    ##!##ch7_lims  = [0., 0.05]
    ##!#ch7_lims  = [np.min(comp_grid_data['ch7_bins']), \
    ##!#    np.max(comp_grid_data['ch7_bins'])]
    ##!##ice_lims  = [101, 105]
    ##!#ice_lims  = [np.min(comp_grid_data['ice_bins']), \
    ##!#    np.max(comp_grid_data['ice_bins'])]
    ##!##sza_lims  = [65, 70]
    ##!#sza_lims  = [np.min(comp_grid_data['sza_bins']), \
    ##!#    np.max(comp_grid_data['sza_bins'])]
    ##!##ai_lims   = [2, 12]
    ##!#ai_lims   = [np.min(comp_grid_data['ai_bins']), \
    ##!#    np.max(comp_grid_data['ai_bins'])]
   
    ##!##if(cld_min is not None):
    ##!##    cld_lims[0] = cld_min
    ##!##if(cld_max is not None):
    ##!##    cld_lims[1] = cld_max
    ##!#if(ch7_min is not None):
    ##!#    ch7_lims[0] = ch7_min
    ##!#if(ch7_max is not None):
    ##!#    ch7_lims[1] = ch7_max
    ##!#if(ice_min is not None):
    ##!#    ice_lims[0] = ice_min
    ##!#if(ice_max is not None):
    ##!#    ice_lims[1] = ice_max
    ##!#if(sza_min is not None):
    ##!#    sza_lims[0] = sza_min
    ##!#if(sza_max is not None):
    ##!#    sza_lims[1] = sza_max
    ##!#if(ai_min is not None):
    ##!#    ai_lims[0] = ai_min
    ##!#if(ai_max is not None):
    ##!#    ai_lims[1] = ai_max
 
    ##!##keep_ai_lims  = np.arange(len(comp_grid_data['ai_bins']))
    ##!#keep_ai_lims = np.where((comp_grid_data['ai_bins'] >= ai_lims[0]) & \
    ##!#                         (comp_grid_data['ai_bins'] <= ai_lims[1]))[0]
    ##!#keep_sza_lims = np.where((comp_grid_data['sza_bins'] >= sza_lims[0]) & \
    ##!#                         (comp_grid_data['sza_bins'] <= sza_lims[1]))[0]
    ##!#keep_ice_lims = np.where((comp_grid_data['ice_bins'] >= ice_lims[0]) & \
    ##!#                         (comp_grid_data['ice_bins'] <= ice_lims[1]))[0]
    ##!#keep_ch7_lims = np.where((comp_grid_data['ch7_bins'] >= ch7_lims[0]) & \
    ##!#                         (comp_grid_data['ch7_bins'] <= ch7_lims[1]))[0]
    ##!##keep_cld_lims = np.where((comp_grid_data['cld_bins'] >= cld_lims[0]) & \
    ##!##                         (comp_grid_data['cld_bins'] <= cld_lims[1]))[0]
    ##!##keep_ice_lims  = np.arange(len(comp_grid_data['ice_bins']))
    ##!##keep_ch7_lims  = np.arange(len(comp_grid_data['ch7_bins']))
    ##!#
    ##!#keep_swf = comp_grid_data['swf_climo'][\
    ##!#    #keep_cld_lims[0]:keep_cld_lims[-1]+1,\
    ##!#    keep_ch7_lims[0]:keep_ch7_lims[-1]+1,\
    ##!#    keep_ice_lims[0]:keep_ice_lims[-1]+1,\
    ##!#    keep_sza_lims[0]:keep_sza_lims[-1]+1,\
    ##!#    keep_ai_lims[0]:keep_ai_lims[-1]+1]

    ##!##keep_cld = np.array([[[[[comp_grid_data['cld_bins'][mm] \
    ##!##    for nn in keep_ai_lims] \
    ##!##    for kk in keep_sza_lims] \
    ##!##    for jj in keep_ice_lims] \
    ##!##    for ii in keep_ch7_lims] \
    ##!##    for mm in keep_cld_lims])

    ##!#keep_ch7 = np.array([[[[comp_grid_data['ch7_bins'][ii] \
    ##!#    for nn in keep_ai_lims] \
    ##!#    for kk in keep_sza_lims] \
    ##!#    for jj in keep_ice_lims] \
    ##!#    for ii in keep_ch7_lims] \
    ##!#    )
    ##!#
    ##!#keep_ice = np.array([[[[comp_grid_data['ice_bins'][jj] \
    ##!#    for nn in keep_ai_lims] \
    ##!#    for kk in keep_sza_lims] \
    ##!#    for jj in keep_ice_lims] \
    ##!#    for ii in keep_ch7_lims] \
    ##!#    )
    ##!#
    ##!#keep_sza = np.array([[[[comp_grid_data['sza_bins'][kk] \
    ##!#    for nn in keep_ai_lims] \
    ##!#    for kk in keep_sza_lims] \
    ##!#    for jj in keep_ice_lims] \
    ##!#    for ii in keep_ch7_lims] \
    ##!#    )
    ##!#
    ##!#keep_ai  = np.array([[[[comp_grid_data['ai_bins'][nn] \
    ##!#    for nn in keep_ai_lims] \
    ##!#    for kk in keep_sza_lims] \
    ##!#    for jj in keep_ice_lims] \
    ##!#    for ii in keep_ch7_lims] \
    ##!#    )
    ##!#

    ##!##final_cld = keep_cld.flatten()[~keep_swf.flatten().mask]
    ##!#final_ch7 = keep_ch7.flatten()[~keep_swf.flatten().mask]
    ##!#final_ice = keep_ice.flatten()[~keep_swf.flatten().mask]
    ##!#final_sza = keep_sza.flatten()[~keep_swf.flatten().mask]
    ##!#final_ai  = keep_ai.flatten()[~keep_swf.flatten().mask]
    ##!#final_swf = keep_swf.compressed()

    ##!#if(xval == 'ai'): 
    ##!#    #ax.scatter(final_ai, final_swf, s = 6, color = 'k')
    ##!#    xlabel = 'OMI UVAI Raw'
    ##!#    local_xdata = final_ai
    ##!#elif(xval == 'sza'):
    ##!#    #ax.scatter(final_sza, final_swf, s = 6, color = 'k')
    ##!#    local_xdata = final_sza
    ##!#    xlabel = 'OMI SZA'
    ##!#elif(xval == 'ice'):
    ##!#    #ax.scatter(final_ice, final_swf, s = 6, color = 'k')
    ##!#    local_xdata = final_ice
    ##!#    xlabel = 'NSIDC ICE CONC.'
    ##!#elif(xval == 'ch7'):
    ##!#    #ax.scatter(final_ch7, final_swf, s = 6, color = 'k')
    ##!#    local_xdata = final_ch7
    ##!#    xlabel = 'MODIS CH7 REFL.'
    ##!##elif(xval == 'cld'):
    ##!##    #ax.scatter(final_ch7, final_swf, s = 6, color = 'k')
    ##!##    local_xdata = final_cld
    ##!##    xlabel = 'CERES CLD FRAC.'
    ##!#local_ydata = final_swf

  
    local_xdata, local_ydata, xlabel = \
            select_comp_grid_scatter_data(comp_grid_data, xval = xval, \
                ch7_min = ch7_min, ch7_max = ch7_max,\
                ice_min = ice_min, ice_max = ice_max,\
                sza_min = sza_min, sza_max = sza_max,\
                ai_min = ai_min,  ai_max = ai_max)

    in_ax = True
    if(ax is None):
        plt.close('all')
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
       
    print('compressed size:',local_ydata.shape)
 
    ax.scatter(local_xdata, local_ydata, s = 6, color = 'k')
    ax.set_xlabel(xlabel)

    if(show_trend):
        if(len(local_xdata) > 10000):
            print("ERROR: Plotting trend line with too many points")
            print("       Preventing this for the sake of computer safety")
        elif(len(local_xdata) < 2):
            print("ERROR: No valid data points for these grid bin ranges")
        else:
            plot_trend_line(ax, local_xdata, local_ydata, color='tab:red', \
                linestyle = '-',  slope = trend_type)

    if(not in_ax):
        plt.show()

def plot_dual_combined_grid_climo(comp_grid_data, combined_data, \
        xval = 'ai', \
        #cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None,\
        save = False, show_trend = False, shade_density = False, \
        trend_type = 'theil-sen'):

    plt.close('all')
    fig = plt.figure(figsize = (9, 4))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    ai_diff = comp_grid_data['ai_bins'][2] - \
        comp_grid_data['ai_edges'][2]
    sza_diff = comp_grid_data['sza_bins'][2] - \
        comp_grid_data['sza_edges'][2]
    ice_diff = comp_grid_data['ice_bins'][2] - \
        comp_grid_data['ice_edges'][2]
    if('cod_bins' in comp_grid_data.keys()):
        ch7_diff = comp_grid_data['cod_bins'][2] - \
            comp_grid_data['cod_edges'][2]
    else:
        ch7_diff = comp_grid_data['ch7_bins'][2] - \
            comp_grid_data['ch7_edges'][2]
    #cld_diff = comp_grid_data['cld_bins'][2] - \
    #    comp_grid_data['cld_edges'][2]

    try:
        plot_combined_scatter(combined_data, ax = ax1, \
            omi_min = ai_min,  omi_max  = ai_max, \
            sza_min = sza_min, sza_max = sza_max, \
            ice_min = ice_min, ice_max = ice_max, \
            ch7_min = ch7_min, ch7_max = ch7_max, \
            #cld_min = cld_min, cld_max = cld_max, \
            ai_diff = ai_diff, sza_diff = sza_diff, \
            ice_diff = ice_diff, ch7_diff = ch7_diff,\
            #cld_diff = cld_diff, \
            trend_type = trend_type, show_trend = show_trend, \
            shade_density = shade_density)
    
    except:
        print("ERROR in plot_combined_scatter")
        return
   
    try: 
        plot_comp_grid_scatter(comp_grid_data, ax = ax2, xval = 'ai', \
            ai_min = ai_min,  ai_max  = ai_max, \
            sza_min = sza_min, sza_max = sza_max, \
            ice_min = ice_min, ice_max = ice_max, \
            ch7_min = ch7_min, ch7_max = ch7_max, \
            #cld_min = cld_min, cld_max = cld_max, \
            show_trend = show_trend, trend_type = trend_type)

    except:
        print("ERROR in plot_comp_grid_scatter")
        return

    ax1.set_title('Raw Colocated L2')
    ax2.set_title('Binned Colocated L2')
    

    title = set_comp_title(\
                           #cld_min, cld_max, \
                           ch7_min, ch7_max, \
                           ice_min, ice_max, \
                           sza_min, sza_max,\
                           ai_min,  ai_max)
    filer = set_file_title(\
                           #cld_min, cld_max, \
                           ch7_min, ch7_max, \
                           ice_min, ice_max, \
                           sza_min, sza_max,\
                           ai_min,  ai_max)
    plt.suptitle(title)

    fig.tight_layout()

    if(save):
        outname = 'comp_dual_' + filer + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image",outname)
    else:
        plt.show()
   
def set_comp_title(\
        #cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None):

    ch7_min = np.round(ch7_min, 2)
    ch7_max = np.round(ch7_max, 2)

    title = '' 
    #if(cld_min is not None):
    #    if(cld_max is not None):
    #        title = title + str(cld_min) + ' < cld < ' + str(cld_max) + '\n'
    #    else:
    #        title = title + 'cld > ' + str(cld_min) + '\n'
    #else:
    #    if(cld_max is not None):
    #        title = title + 'cld < ' + str(cld_max) + '\n'

    if(ch7_min is not None):
        if(ch7_max is not None):
            title = title + str(ch7_min) + ' < ch7 < ' + str(ch7_max) + '\n'
        else:
            title = title + 'ch7 > ' + str(ch7_min) + '\n'
    else:
        if(ch7_max is not None):
            title = title + 'ch7 < ' + str(ch7_max) + '\n'

    if(ice_min is not None):
        if(ice_max is not None):
            title = title + str(ice_min) + ' < ice < ' + str(ice_max) + '\n'
        else:
            title = title + 'ice > ' + str(ice_min) + '\n'
    else:
        if(ice_max is not None):
            title = title + 'ice < ' + str(ice_max) + '\n'

    if(sza_min is not None):
        if(sza_max is not None):
            title = title + str(sza_min) + ' < sza < ' + str(sza_max) + '\n'
        else:
            title = title + 'sza > ' + str(sza_min) + '\n'
    else:
        if(sza_max is not None):
            title = title + 'sza < ' + str(sza_max) + '\n'

    if(ai_min is not None):
        if(ai_max is not None):
            title = title + str(ai_min) + ' < ai < ' + str(ai_max)
        else:
            title = title + 'ai > ' + str(ai_min)
    else:
        if(ai_max is not None):
            title = title + 'ai < ' + str(ai_max)

    return title 

def set_file_title(\
        #cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None):

    ch7_min = np.round(ch7_min, 2)
    ch7_max = np.round(ch7_max, 2)


    title = '' 
    #if(cld_min is not None):
    #    cld_min = int(cld_min)
    #    if(cld_max is not None):
    #        cld_max = int(cld_max)
    #        title = title + str(cld_min) + 'cld' + str(cld_max) + '_'
    #    else:
    #        title = title + str(cld_min) + 'cld_'
    #else:
    #    if(cld_max is not None):
    #        cld_max = int(cld_max)
    #        title = title + 'cld' + str(cld_max) + '_'

    if(ch7_min is not None):
        if(ch7_max is not None):
            title = title + str(ch7_min) + 'ch7' + str(ch7_max) + '_'
        else:
            title = title + str(ch7_min) + 'ch7_'
    else:
        if(ch7_max is not None):
            title = title + 'ch7' + str(ch7_max) + '_'

    if(ice_min is not None):
        ice_min = int(ice_min)
        if(ice_max is not None):
            ice_max = int(ice_max)
            title = title + str(ice_min) + 'ice' + str(ice_max) + '_'
        else:
            title = title + str(ice_min) + 'ice_'
    else:
        if(ice_max is not None):
            ice_max = int(ice_max)
            title = title + 'ice' + str(ice_max) + '_'

    if(sza_min is not None):
        sza_min = int(sza_min)
        if(sza_max is not None):
            sza_max = int(sza_max)
            title = title + str(sza_min) + 'sza' + str(sza_max) + '_'
        else:
            title = title + str(sza_min) + 'sza_'
    else:
        if(sza_max is not None):
            sza_max = int(sza_max)
            title = title + 'sza' + str(sza_max) + '_'

    if(ai_min is not None):
        ai_min  = int(ai_min)
        if(ai_max is not None):
            ai_max  = int(ai_max)
            title = title + str(ai_min) + 'ai' + str(ai_max)
        else:
            title = title + str(ai_min) + 'ai'
    else:
        if(ai_max is not None):
            ai_max  = int(ai_max)
            title = title + 'ai' + str(ai_max)

    print("HERE",title)

    #if(ch7_min is not None):
    #if(ch7_max is not None):

    #if(ice_min is not None):
    #if(ice_max is not None):

    #if(sza_min is not None):
    #if(sza_max is not None):

    #if(ai_min is not None):
    #if(ai_max is not None):

    return title 

# smoother: either 'none','smooth','smoother'
# sizer   : either 1, 2, or 3
def calc_raw_grid_slopes(combined_data, comp_grid_data, \
        ai_min = None, ai_max = None, \
        trend_type = 'theil-sen', \
        smoother = 'None', sizer = 1):

    ice_mins = np.array([0,  80, 105, 21])
    ice_maxs = np.array([20,100,None, 79])
    #ice_mins = np.array([0,80,105])
    #ice_maxs = np.array([20,100,None])

    if(smoother == 'None'):    
        if('cod_bins' in comp_grid_data.keys()):
            cloud_var = 'cod'
            ch7_mins = comp_grid_data['cod_bins'][0::sizer]
            ch7_maxs = comp_grid_data['cod_bins'][0::sizer]
        else:
            cloud_var = 'ch7'
            ch7_mins = comp_grid_data['ch7_bins'][0::sizer]
            ch7_maxs = comp_grid_data['ch7_bins'][0::sizer]
        
        sza_mins = comp_grid_data['sza_bins'][4::sizer]
        sza_maxs = comp_grid_data['sza_bins'][4::sizer]
    elif(smoother == 'smooth'): 
        if('cod_bins' in comp_grid_data.keys()):
            cloud_var = 'cod'
            ch7_mins = comp_grid_data['cod_bins'][0:-1:sizer]
            ch7_maxs = comp_grid_data['cod_bins'][1::sizer]
        else:
            cloud_var = 'ch7'
            ch7_mins = comp_grid_data['ch7_bins'][0:-1:sizer]
            ch7_maxs = comp_grid_data['ch7_bins'][1::sizer]
        
        sza_mins = comp_grid_data['sza_bins'][4:-1:sizer]
        sza_maxs = comp_grid_data['sza_bins'][5::sizer]
    elif(smoother == 'smoother'):
        if('cod_bins' in comp_grid_data.keys()):
            cloud_var = 'cod'
            ch7_mins = comp_grid_data['cod_bins'][0:-2:sizer]
            ch7_maxs = comp_grid_data['cod_bins'][2::sizer]
        else:
            cloud_var = 'ch7'
            ch7_mins = comp_grid_data['ch7_bins'][0:-2:sizer]
            ch7_maxs = comp_grid_data['ch7_bins'][2::sizer]
        
        sza_mins = comp_grid_data['sza_bins'][4:-2:sizer]
        sza_maxs = comp_grid_data['sza_bins'][6::sizer]
    
    raw_slopes   = np.full((ice_mins.size, ch7_mins.size, sza_mins.size), \
        np.nan)
    grid_slopes  = np.full((ice_mins.size, ch7_mins.size, sza_mins.size), \
        np.nan)
    raw_stderr   = np.full((ice_mins.size, ch7_mins.size, sza_mins.size), \
        np.nan)
    grid_stderr  = np.full((ice_mins.size, ch7_mins.size, sza_mins.size), \
        np.nan)
    raw_pvals    = np.full((ice_mins.size, ch7_mins.size, sza_mins.size), \
        np.nan)
    grid_pvals   = np.full((ice_mins.size, ch7_mins.size, sza_mins.size), \
        np.nan)
    raw_counts   = np.full((ice_mins.size, ch7_mins.size, sza_mins.size), \
        np.nan)
    grid_counts  = np.full((ice_mins.size, ch7_mins.size, sza_mins.size), \
        np.nan)

    # 0: "cloudy"
    # 1: "probably cloudy"
    # 2: "probably clear"
    # 3: "clear"
    # 4: "cloudy" + "probably cloudy"
    # 5: "clear"  + "probably clear"
    raw_cldvals  = np.full((ice_mins.size, ch7_mins.size, sza_mins.size, 6), \
        np.nan)
    cldval_names = ['cloudy','probably cloudy','probably clear','clear',\
        'cloudy + probably cloudy','clear + probably clear']
    
    ai_diff  = comp_grid_data['ai_bins'][2] - comp_grid_data['ai_edges'][2]
    sza_diff = comp_grid_data['sza_bins'][2] - comp_grid_data['sza_edges'][2]
    ice_diff = comp_grid_data['ice_bins'][2] - comp_grid_data['ice_edges'][2]
    ch7_diff = comp_grid_data[cloud_var + '_bins'][2] - \
        comp_grid_data[cloud_var + '_edges'][2]
    
    for ii in range(ice_mins.size):
        for jj in range(ch7_mins.size):
            for kk in range(sza_mins.size):
                print(ii,jj,kk)
   
                 
                #raw_xdata, raw_ydata, raw_cloud = \
                return_dict = \
                    select_combined_scatter_data(combined_data, \
                        xval = 'omi_uvai_raw', yval = 'ceres_swf', 
                        #cld_min = None, cld_max = None, \
                        ch7_min = ch7_mins[jj], ch7_max = ch7_maxs[jj], \
                        sza_min = sza_mins[kk], sza_max = sza_maxs[kk], \
                        ice_min = ice_mins[ii], ice_max = ice_maxs[ii], \
                        omi_min = ai_min, omi_max = ai_max, \
                        ai_diff = ai_diff, sza_diff = sza_diff, \
                        ice_diff = ice_diff, ch7_diff = ch7_diff, \
                        include_cloud = True,
                        #cld_diff = cld_diff, \
                        cloud_var = cloud_var 
                        )
   
                raw_xdata = return_dict['local_xdata']
                raw_ydata = return_dict['local_ydata']
                raw_cloud = return_dict['local_cloud']
 
                grid_xdata, grid_ydata, xlabel = \
                        select_comp_grid_scatter_data(comp_grid_data, xval = 'ai', \
                            ch7_min = ch7_mins[jj], ch7_max = ch7_maxs[jj],\
                            ice_min = ice_mins[ii], ice_max = ice_maxs[ii],\
                            sza_min = sza_mins[kk], sza_max = sza_maxs[kk],\
                            ai_min = ai_min,  ai_max = ai_max)
   
                # Handle the cloud fraction stuff
                # -------------------------------
                totalsize = raw_cloud.size
                if(totalsize > 0):
                    for nn in range(4):
                        raw_cldvals[ii,jj,kk,nn] = \
                            np.where(raw_cloud == nn)[0].size / totalsize
                    raw_cldvals[ii,jj,kk,4] = \
                        np.where((raw_cloud == 0) | (raw_cloud == 1))[0].size / totalsize
                    raw_cldvals[ii,jj,kk,5] = \
                        np.where((raw_cloud == 2) | (raw_cloud == 3))[0].size / totalsize
 
                # Calculate slopes 
                # ----------------
                #if((len(grid_xdata) > 20) & (len(raw_xdata) > 20)):
                if(trend_type == 'theil-sen'):
                    if((len(raw_xdata) >= 20)):
                        res = stats.theilslopes(raw_ydata, raw_xdata, 0.90)
                        raw_slopes[ii,jj,kk]   = res[0]

                    if((len(grid_xdata) > 20)):
                        res = stats.theilslopes(grid_ydata, grid_xdata, 0.90)
                        grid_slopes[ii,jj,kk]   = res[0]
                elif(trend_type == 'linregress'):
                    if((len(raw_xdata) >= 20)):
                        result1 = stats.linregress(raw_xdata,raw_ydata)
                        raw_slopes[ii,jj,kk] = result1.slope 
                        raw_stderr[ii,jj,kk] = result1.stderr
                        raw_pvals[ii,jj,kk]  = result1.pvalue

                    if((len(grid_xdata) >= 20)):
                        result2 = stats.linregress(grid_xdata,grid_ydata)
                        grid_slopes[ii,jj,kk] = result2.slope 
                        grid_stderr[ii,jj,kk] = result2.stderr
                        grid_pvals[ii,jj,kk]  = result2.pvalue
                else:
                    print("WARNING: invalid trend_type specified. Using theil-sen")
                    if((len(raw_xdata) >= 20)):
                        res = stats.theilslopes(raw_ydata, raw_xdata, 0.90)
                        raw_slopes[ii,jj,kk]   = res[0]

                    if((len(grid_xdata) >= 20)):
                        res = stats.theilslopes(grid_ydata, grid_xdata, 0.90)
                        grid_slopes[ii,jj,kk]   = res[0]
   
                raw_counts[ii,jj,kk]  = raw_xdata.shape[0]
                grid_counts[ii,jj,kk] = grid_xdata.shape[0]
                print(raw_xdata.shape, grid_xdata.shape, raw_slopes[ii,jj,kk], grid_slopes[ii,jj,kk])
    
                #raw_slopes[ii,jj,kk]   = 
    
    raw_slopes  = np.ma.masked_invalid(raw_slopes)
    raw_cldvals = np.ma.masked_invalid(raw_cldvals)
    grid_slopes = np.ma.masked_invalid(grid_slopes)
    raw_stderr  = np.ma.masked_invalid(raw_stderr)
    grid_stderr = np.ma.masked_invalid(grid_stderr)
    raw_pvals   = np.ma.masked_invalid(raw_pvals)
    grid_pvals  = np.ma.masked_invalid(grid_pvals)
    raw_counts  = np.ma.masked_invalid(raw_counts)
    grid_counts = np.ma.masked_invalid(grid_counts)

    out_dict = {}
    out_dict['raw_slopes'] = raw_slopes
    out_dict['grid_slopes'] = grid_slopes
    out_dict['raw_cldvals'] = raw_cldvals
    out_dict['cldval_names'] = cldval_names
    out_dict['raw_stderr'] = raw_stderr
    out_dict['grid_stderr'] = grid_stderr
    out_dict['raw_pvals'] = raw_pvals
    out_dict['grid_pvals'] = grid_pvals
    out_dict['raw_counts'] = raw_counts
    out_dict['grid_counts'] = grid_counts
    out_dict['sza_mins'] = sza_mins
    out_dict['sza_maxs'] = sza_maxs
    out_dict[cloud_var + '_mins'] = ch7_mins
    out_dict[cloud_var + '_maxs'] = ch7_maxs
    out_dict['ice_mins'] = ice_mins
    out_dict['ice_maxs'] = ice_maxs
    out_dict['smoother'] = smoother
    out_dict['sizer'] = sizer
    out_dict['trend_type'] = trend_type
    
    return out_dict

def plot_raw_grid_slopes(out_dict, vmin = -10, vmax = 10, \
        save = False, hatch_stderr = False, hatch_style = '.', \
        maxerr = 2):

    plt.close('all')
    fig = plt.figure(figsize = (10, 6))
    ax1 = fig.add_subplot(2,3,1)
    ax2 = fig.add_subplot(2,3,2)
    ax3 = fig.add_subplot(2,3,3)
    
    ax4 = fig.add_subplot(2,3,4)
    ax5 = fig.add_subplot(2,3,5)
    ax6 = fig.add_subplot(2,3,6)
    
    cmap = 'bwr'
   
    if('cod_mins' in out_dict.keys()):
        cloud_var = 'cod'
    else:
        cloud_var = 'ch7'
     
    shrk = 1.0
    mesh = ax1.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        out_dict['raw_slopes'][0,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    #cbar = plt.colorbar(mesh,\
    #    ax = ax1, orientation='vertical',shrink = shrk, extend = 'both')
    ax1.set_xlabel('SZA')
    ax1.set_ylabel(cloud_var.upper())
    ax1.set_title('Raw Ocean')
    if(hatch_stderr):
        hasher = np.ma.masked_where(out_dict['raw_stderr'][0,:,:] > maxerr, \
            out_dict['raw_stderr'][0,:,:])
        ax1.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')
    
    mesh = ax2.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        out_dict['raw_slopes'][1,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    #cbar = plt.colorbar(mesh,\
    #    ax = ax2, orientation='vertical',shrink = shrk, extend = 'both')
    ax2.set_xlabel('SZA')
    ax2.set_ylabel(cloud_var.upper())
    ax2.set_title('Raw Ice')
    if(hatch_stderr):
        hasher = np.ma.masked_where(out_dict['raw_stderr'][1,:,:] > maxerr, \
            out_dict['grid_stderr'][2,:,:])
        ax2.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')
    
    mesh =ax3.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        out_dict['raw_slopes'][2,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    cbar = plt.colorbar(mesh,\
        ax = ax3, orientation='vertical',shrink = shrk, extend = 'both')
    cbar.set_label('AI - SWF slope [Wm$^{-2}$/AI]')
    ax3.set_xlabel('SZA')
    ax3.set_ylabel(cloud_var.upper())
    ax3.set_title('Raw Land')
    if(hatch_stderr):
        hasher = np.ma.masked_where(out_dict['raw_stderr'][2,:,:] > maxerr, \
            out_dict['raw_stderr'][2,:,:])
        ax3.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')
    
    mesh = ax4.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        out_dict['grid_slopes'][0,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    #cbar = plt.colorbar(mesh,\
    #    ax = ax4, orientation='vertical',shrink = shrk, extend = 'both')
    ax4.set_xlabel('SZA')
    ax4.set_ylabel(cloud_var.upper())
    ax4.set_title('Grid Ocean')
    if(hatch_stderr):
        hasher = np.ma.masked_where(out_dict['grid_stderr'][0,:,:] > maxerr, \
            out_dict['grid_stderr'][0,:,:])
        ax4.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')
    
    mesh = ax5.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'],\
        out_dict['grid_slopes'][1,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    #cbar = plt.colorbar(mesh,\
    #    ax = ax5, orientation='vertical',shrink = shrk, extend = 'both')
    ax5.set_xlabel('SZA')
    ax5.set_ylabel(cloud_var.upper())
    ax5.set_title('Grid Ice')
    if(hatch_stderr):
        hasher = np.ma.masked_where(out_dict['grid_stderr'][1,:,:] > maxerr, \
            out_dict['grid_stderr'][1,:,:])
        ax5.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')
    
    mesh = ax6.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        out_dict['grid_slopes'][2,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    cbar = plt.colorbar(mesh,\
        ax = ax6, orientation='vertical',shrink = shrk, extend = 'both')
    cbar.set_label('AI - SWF slope [Wm$^{-2}$/AI]')
    ax6.set_xlabel('SZA')
    ax6.set_ylabel(cloud_var.upper())
    ax6.set_title('Grid Land')
    if(hatch_stderr):
        hasher = np.ma.masked_where(out_dict['grid_stderr'][2,:,:] > maxerr, \
            out_dict['grid_stderr'][2,:,:])
        ax6.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')
  
    title_str = 'MODIS CH7 Bin Size = ' + str(out_dict[cloud_var + '_mins'][1] - \
                out_dict[cloud_var + '_mins'][0]) + '\n' + \
                 'OMI SZA Bin Size = ' + str(out_dict['sza_mins'][1] - \
                out_dict['sza_mins'][0])
    if(out_dict['smoother'] == 'smooth'):
        title_str = title_str + '\n1-Bin Smoothed'
        out_adder = '_smooth1bin'
    elif(out_dict['smoother'] == 'smoother'):
        title_str = title_str + '\n2-Bin Smoothed'
        out_adder = '_smooth2bin'
    else:
        out_adder = ''
    if(hatch_stderr):
        title_str = title_str + '\nHatching: Slope Standard Error < ' + \
            str(maxerr)
        hatch_adder = '_hatch' + str(maxerr)
    else:
        hatch_adder = ''
    plt.suptitle(title_str)
    
    fig.tight_layout()
   
    if(save): 
        if(out_dict['trend_type'] == 'linregress'):
            type_adder = '_lin'
        elif(out_dict['trend_type'] == 'theil-sen'):
            type_adder = '_thl'

        outname = 'comp_grid_ai_swf_slopes_highres' + type_adder + \
            out_adder + hatch_adder + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved",outname)
    else: 
        plt.show()

# dtype: 'raw','grid'
# cld_idx:
#   - 0: "cloudy"
#   - 1: "probably cloudy"
#   - 2: "probably clear"
#   - 3: "clear"
def plot_slopes_cloud_types(out_dict, vmin = -10, vmax = 10, \
        save = False, hatch_stderr = False, hatch_style = '.', \
        cld_idx = 0, maxerr = 2, remove_high_error = True, \
        hatch_cloud = False, cloud_style = '//', \
        min_cloud = 0.95):

    plt.close('all')
    fig = plt.figure(figsize = (10, 6))
    # Trends
    ax1 = fig.add_subplot(2,3,1)
    ax2 = fig.add_subplot(2,3,2)
    ax3 = fig.add_subplot(2,3,3)
   
    # Cloud fracs 
    ax4 = fig.add_subplot(2,3,4)
    ax5 = fig.add_subplot(2,3,5)
    ax6 = fig.add_subplot(2,3,6)
    
    cmap = 'bwr'
   
    if('cod_mins' in out_dict.keys()):
        cloud_var = 'cod'
    else:
        cloud_var = 'ch7'
 
    shrk = 1.0
    if(remove_high_error):
        plotter = np.ma.masked_where(out_dict['raw_stderr'][0,:,:] > maxerr, \
            out_dict['raw_slopes'][0,:,:])
    else:
        plotter = out_dict['raw_slopes'][0,:,:]
    mesh = ax1.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        plotter, \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    #cbar = plt.colorbar(mesh,\
    #    ax = ax1, orientation='vertical',shrink = shrk, extend = 'both')
    ax1.set_xlabel('SZA')
    ax1.set_ylabel(cloud_var.upper())
    ax1.set_title('Raw Ocean')
    if(hatch_stderr):
        hasher = np.ma.masked_where(out_dict['raw_stderr'][0,:,:] > maxerr, \
            out_dict['raw_stderr'][0,:,:])
        ax1.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')
    
    if(remove_high_error):
        plotter = np.ma.masked_where(out_dict['raw_stderr'][1,:,:] > maxerr, \
            out_dict['raw_slopes'][1,:,:])
    else:
        plotter = out_dict['raw_slopes'][1,:,:]
    mesh = ax2.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        plotter, \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    #cbar = plt.colorbar(mesh,\
    #    ax = ax2, orientation='vertical',shrink = shrk, extend = 'both')
    ax2.set_xlabel('SZA')
    ax2.set_ylabel(cloud_var.upper())
    ax2.set_title('Raw Ice')
    if(hatch_stderr):
        hasher = np.ma.masked_where(out_dict['raw_stderr'][1,:,:] > maxerr, \
            out_dict['grid_stderr'][2,:,:])
        ax2.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')
    
    if(remove_high_error):
        plotter = np.ma.masked_where(out_dict['raw_stderr'][2,:,:] > maxerr, \
            out_dict['raw_slopes'][2,:,:])
    else:
        plotter = out_dict['raw_slopes'][2,:,:]
    mesh =ax3.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        plotter, \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    cbar = plt.colorbar(mesh,\
        ax = ax3, orientation='vertical',shrink = shrk, extend = 'both')
    cbar.set_label('AI - SWF slope [Wm$^{-2}$/AI]')
    ax3.set_xlabel('SZA')
    ax3.set_ylabel(cloud_var.upper())
    ax3.set_title('Raw Land')
    if(hatch_stderr):
        hasher = np.ma.masked_where(out_dict['raw_stderr'][2,:,:] > maxerr, \
            out_dict['raw_stderr'][2,:,:])
        ax3.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')
    
    mask_cloud = np.ma.masked_where(\
        out_dict['raw_cldvals'][0,:,:,cld_idx] < -9, \
        out_dict['raw_cldvals'][0,:,:,cld_idx])
    mesh = ax4.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        mask_cloud, cmap = 'plasma', shading = 'auto', vmin = 0, vmax = 1)
    #cbar = plt.colorbar(mesh,\
    #    ax = ax4, orientation='vertical',shrink = shrk, extend = 'both')
    ax4.set_xlabel('SZA')
    ax4.set_ylabel(cloud_var.upper())
    ax4.set_title('Ocean Cloud')
    if(hatch_cloud):
        if((cld_idx <= 1) | (cld_idx == 4)):
            hasher = np.ma.masked_where(mask_cloud < min_cloud, \
                mask_cloud)
        else:
            hasher = np.ma.masked_where(mask_cloud > min_cloud, \
                mask_cloud)
        ax1.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = cloud_style, alpha=0., shading = 'auto')
        ax4.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = cloud_style, alpha=0., shading = 'auto')
  
    
    mask_cloud = np.ma.masked_where(\
        out_dict['raw_cldvals'][1,:,:,cld_idx] < -9, \
        out_dict['raw_cldvals'][1,:,:,cld_idx])
    mesh = ax5.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        mask_cloud, cmap = 'plasma', shading = 'auto', vmin = 0, vmax = 1)
    #cbar = plt.colorbar(mesh,\
    #    ax = ax5, orientation='vertical',shrink = shrk, extend = 'both')
    ax5.set_xlabel('SZA')
    ax5.set_ylabel(cloud_var.upper())
    ax5.set_title('Ice Cloud')
    if(hatch_cloud):
        if((cld_idx <= 1) | (cld_idx == 4)):
            hasher = np.ma.masked_where(mask_cloud < min_cloud, \
                mask_cloud)
        else:
            hasher = np.ma.masked_where(mask_cloud > min_cloud, \
                mask_cloud)
        ax2.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = cloud_style, alpha=0., shading = 'auto')
        ax5.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = cloud_style, alpha=0., shading = 'auto')
   
    mask_cloud = np.ma.masked_where(\
        out_dict['raw_cldvals'][2,:,:,cld_idx] < -9, \
        out_dict['raw_cldvals'][2,:,:,cld_idx])
    mesh = ax6.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        mask_cloud, cmap = 'plasma', shading = 'auto', vmin = 0, vmax = 1)
    cbar = plt.colorbar(mesh,\
        ax = ax6, orientation='vertical',shrink = shrk, extend = 'both')
    cbar.set_label('Percent of\n\"' + out_dict['cldval_names'][cld_idx] + '\"')
    ax6.set_xlabel('SZA')
    ax6.set_ylabel(cloud_var.upper())
    ax6.set_title('Land Cloud')
    if(hatch_cloud):
        if((cld_idx <= 1) | (cld_idx == 4)):
            hasher = np.ma.masked_where(mask_cloud < min_cloud, \
                mask_cloud)
        else:
            hasher = np.ma.masked_where(mask_cloud > min_cloud, \
                mask_cloud)
        ax3.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = cloud_style, alpha=0., shading = 'auto')
        ax6.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            hasher, hatch = cloud_style, alpha=0., shading = 'auto')
  
    title_str = 'MODIS ' + cloud_var.upper() + ' Bin Size = ' + \
                str(out_dict[cloud_var + '_mins'][1] - \
                out_dict[cloud_var + '_mins'][0]) + '\n' + \
                 'OMI SZA Bin Size = ' + str(out_dict['sza_mins'][1] - \
                out_dict['sza_mins'][0])
    if(out_dict['smoother'] == 'smooth'):
        title_str = title_str + '\n1-Bin Smoothed'
        out_adder = '_smooth1bin'
    elif(out_dict['smoother'] == 'smoother'):
        title_str = title_str + '\n2-Bin Smoothed'
        out_adder = '_smooth2bin'
    else:
        out_adder = ''
    if((not hatch_stderr) and (not hatch_cloud)):
        hatch_adder = ''
    else:
        if(hatch_stderr):
            title_str = title_str + '\nDotted: Slope Standard Error < ' + \
                str(maxerr)
            hatch_adder = '_hatcherr' + str(maxerr)
        if(hatch_cloud):
            title_str = title_str + '\nSlashed: Pcnt of ' + \
                out_dict['cldval_names'][cld_idx] + ' pixels > ' + \
                str(min_cloud * 100) + '%'
            hatch_adder = '_hatchcld'
    if(remove_high_error):
        error_adder = '_nohigherror'
    else:
        error_adder = ''

    plt.suptitle(title_str)
    
    fig.tight_layout()
   
    if(save): 
        if(out_dict['trend_type'] == 'linregress'):
            type_adder = '_lin'
        elif(out_dict['trend_type'] == 'theil-sen'):
            type_adder = '_thl'

        outname = 'comp_grid_ai_swf_slopecloud_highres' + type_adder + \
            out_adder + hatch_adder + error_adder + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved",outname)
    else: 
        plt.show()


def calc_slope_clear_clean_sfctype(out_dict, sfc_type_idx, \
    cld_idx, min_cloud = 0.95, maxerr = 2, data_type = 'raw'):

    mask_cloud = np.ma.masked_where(\
        out_dict['raw_cldvals'][sfc_type_idx,:,:,cld_idx] < -9, \
        out_dict['raw_cldvals'][sfc_type_idx,:,:,cld_idx])
    
    cloud_slopes = np.ma.masked_where(mask_cloud < min_cloud, \
        out_dict[data_type + '_slopes'][sfc_type_idx,:,:])
    clear_slopes = np.ma.masked_where(mask_cloud >= min_cloud, \
        out_dict[data_type + '_slopes'][sfc_type_idx,:,:])
    
    cloud_slopes = np.ma.masked_where(\
        out_dict[data_type + '_stderr'][sfc_type_idx,:,:] > maxerr, \
        cloud_slopes)
    clear_slopes = np.ma.masked_where(\
        out_dict[data_type + '_stderr'][sfc_type_idx,:,:] > maxerr, \
        clear_slopes)
    
    sza_cloud_means = np.nanmean(cloud_slopes, axis = 0)
    sza_clear_means = np.nanmean(clear_slopes, axis = 0)
    sza_cloud_std   = np.std(cloud_slopes, axis = 0)
    sza_clear_std   = np.std(clear_slopes, axis = 0)
    ch7_cloud_means = np.nanmean(cloud_slopes, axis = 1)
    ch7_clear_means = np.nanmean(clear_slopes, axis = 1)
    ch7_cloud_std   = np.std(cloud_slopes, axis = 1)
    ch7_clear_std   = np.std(clear_slopes, axis = 1)

    return_dict = {}
    if('cod_mins' in out_dict.keys()):
        out_var = 'cod'
    else:
        out_var = 'ch7'

    return_dict['sza_cloud_means'] = sza_cloud_means
    return_dict['sza_cloud_std']   = sza_cloud_std
    return_dict['sza_clear_means'] = sza_clear_means
    return_dict['sza_clear_std']   = sza_clear_std
    return_dict[out_var + '_cloud_means'] = ch7_cloud_means
    return_dict[out_var + '_cloud_std']   = ch7_cloud_std
    return_dict[out_var + '_clear_means'] = ch7_clear_means
    return_dict[out_var + '_clear_std']   = ch7_clear_std

    return return_dict   
 
    
# dtype: 'raw','grid'
# cld_idx:
#   - 0: "cloudy"
#   - 1: "probably cloudy"
#   - 2: "probably clear"
#   - 3: "clear"
# xvar: 'sza', or 'ch7'
def plot_slopes_cloud_types_szamean(out_dict,cld_idx = 0, maxerr = 2, \
        data_type = 'raw', xvar = 'sza', remove_high_error = True, \
        hatch_cloud = False, \
        min_cloud = 0.95, save = False):

    ocean_slopes = calc_slope_clear_clean_sfctype(out_dict, 0, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
    ice_slopes   = calc_slope_clear_clean_sfctype(out_dict, 1, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
    land_slopes  = calc_slope_clear_clean_sfctype(out_dict, 2, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
    if(out_dict['raw_slopes'].shape[0] == 4):
        include_mix = True
        mix_slopes  = calc_slope_clear_clean_sfctype(out_dict, 3, cld_idx, \
            maxerr = maxerr, min_cloud = min_cloud, data_type = data_type)
        figsize = (6, 6)
    else:
        include_mix = False
        figsize = (9, 3)
        

    plt.close('all')
    fig = plt.figure(figsize = figsize)
    # Trends
    if(include_mix):
        axs = fig.subplots(nrows = 2, ncols = 2, sharey = True)
        ax1 = axs[0,0]
        ax2 = axs[0,1]
        ax3 = axs[1,0]
        ax4 = axs[1,1]
    else:
        axs = fig.subplots(nrows = 1, ncols = 3, sharey = True)
        ax1 = axs[0]
        ax2 = axs[1]

    ax1.axhline(0, color = 'k', linestyle = '--')
    ax2.axhline(0, color = 'k', linestyle = '--')
    ax3.axhline(0, color = 'k', linestyle = '--')
    if(include_mix): ax4.axhline(0, color = 'k', linestyle = '--')
 
    # Ocean stuff 
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_cloud_means'], \
        color = 'tab:blue', label = 'Cloud')
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_cloud_means'] - \
        ocean_slopes[xvar + '_cloud_std'], ':', \
        color = 'tab:blue')
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_cloud_means'] + \
        ocean_slopes[xvar + '_cloud_std'], ':', \
        color = 'tab:blue')

    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_clear_means'], \
        color = 'tab:orange', label = 'Clear')
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_clear_means'] - \
        ocean_slopes[xvar + '_clear_std'], ':', \
        color = 'tab:orange')
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_clear_means'] + \
        ocean_slopes[xvar + '_clear_std'], ':', \
        color = 'tab:orange')

    # Ice stuff 
    ax2.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_cloud_means'], \
        color = 'tab:blue', label = 'Cloud')
    ax2.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_cloud_means'] - \
        ice_slopes[xvar + '_cloud_std'], ':', \
        color = 'tab:blue')
    ax2.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_cloud_means'] + \
        ice_slopes[xvar + '_cloud_std'], ':', \
        color = 'tab:blue')

    ax2.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_clear_means'], \
        color = 'tab:orange', label = 'Clear')
    ax2.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_clear_means'] - \
        ice_slopes[xvar + '_clear_std'], ':', \
        color = 'tab:orange')
    ax2.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_clear_means'] + \
        ice_slopes[xvar + '_clear_std'], ':', \
        color = 'tab:orange')

    # Land stuff 
    ax3.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_cloud_means'], \
        color = 'tab:blue', label = 'Cloud')
    ax3.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_cloud_means'] - \
        land_slopes[xvar + '_cloud_std'], ':', \
        color = 'tab:blue')
    ax3.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_cloud_means'] + \
        land_slopes[xvar + '_cloud_std'], ':', \
        color = 'tab:blue')

    ax3.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_clear_means'], \
        color = 'tab:orange', label = 'Clear')
    ax3.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_clear_means'] - \
        land_slopes[xvar + '_clear_std'], ':', \
        color = 'tab:orange')
    ax3.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_clear_means'] + \
        land_slopes[xvar + '_clear_std'], ':', \
        color = 'tab:orange')

    ax1_mins = ax1.get_ylim()[0]
    ax2_mins = ax2.get_ylim()[0]
    ax3_mins = ax3.get_ylim()[0]
    ax1_maxs = ax1.get_ylim()[1]
    ax2_maxs = ax2.get_ylim()[1]
    ax3_maxs = ax3.get_ylim()[1]

    if(include_mix):
        ax4.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_cloud_means'], \
            color = 'tab:blue', label = 'Cloud')
        ax4.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_cloud_means'] - \
            mix_slopes[xvar + '_cloud_std'], ':', \
            color = 'tab:blue')
        ax4.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_cloud_means'] + \
            mix_slopes[xvar + '_cloud_std'], ':', \
            color = 'tab:blue')

        ax4.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_clear_means'], \
            color = 'tab:orange', label = 'Clear')
        ax4.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_clear_means'] - \
            mix_slopes[xvar + '_clear_std'], ':', \
            color = 'tab:orange')
        ax4.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_clear_means'] + \
            mix_slopes[xvar + '_clear_std'], ':', \
            color = 'tab:orange')
        ax4.legend()

        ax4_mins = ax4.get_ylim()[0]
        ax4_maxs = ax4.get_ylim()[1]
    else:
        ax3.legend()

    all_mins = [ax1_mins, ax2_mins, ax3_mins]
    all_maxs = [ax1_maxs, ax2_maxs, ax3_maxs]
    if(include_mix):
        all_mins += [ax4_mins]
        all_maxs += [ax4_maxs]

    ax1.set_ylim(np.min(all_mins), np.max(all_maxs))
    ax2.set_ylim(np.min(all_mins), np.max(all_maxs))
    ax3.set_ylim(np.min(all_mins), np.max(all_maxs))
    if(include_mix): 
        ax4.set_ylim(np.min(all_mins), np.max(all_maxs))

    ax1.grid(alpha = 0.5)
    ax2.grid(alpha = 0.5)
    ax3.grid(alpha = 0.5)
    if(include_mix): 
        ax4.grid(alpha = 0.5)

    ax1.set_title('Ocean')
    ax2.set_title('Ice')
    ax3.set_title('Land')
    if(include_mix):  
        ax4.set_title('Mix')
 
    if(xvar == 'sza'): labeltxt = 'OMI SZA'
    elif(xvar == 'ch7'): labeltxt = 'MODIS CH7'
    elif(xvar == 'cod'): labeltxt = 'MODIS COD'
    ax1.set_xlabel(labeltxt)
    ax2.set_xlabel(labeltxt)
    ax3.set_xlabel(labeltxt)
    if(include_mix):
        ax4.set_xlabel(labeltxt)
 
    ax1.set_ylabel('AI-SWF Slope [Wm$^{-2}$AI$^{-1}$]')
    if(include_mix):
        ax3.set_ylabel('AI-SWF Slope [Wm$^{-2}$AI$^{-1}$]')
    #ax2.set_ylabel('AI-SWF Slope [Wm$^{-2}$]')
    #ax3.set_ylabel('AI-SWF Slope [Wm$^{-2}$]')

    fig.tight_layout()
   
    if(save):

        if(out_dict['smoother'] == 'None'):
            smth_adder = '_raw'
        elif(out_dict['smoother'] == 'smooth'):
            smth_adder = '_smth'
        elif(out_dict['smoother'] == 'smoother'):
            smth_adder = '_smth2'
 
        if(out_dict['trend_type'] == 'linregress'):
            type_adder = '_lin'
        elif(out_dict['trend_type'] == 'theil-sen'):
            type_adder = '_thl'

        if(include_mix):
            mix_add = '_mix'
        else:
            mix_add = ''
        outname = 'comp_grid_ai_swf_slopecloud_' + xvar + 'means' + type_adder + \
            smth_adder + mix_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved",outname)
    else: 
        plt.show()


# dtype: 'raw','grid'
def plot_compare_grid_climo_stderr(comp_grid_data,  \
        dtype, ice_idx, save = False):

    vmin = -15
    vmax = 15
   
    label_dict = {
        0: 'Ocean',
        1: 'Ice',\
        2: 'Land',\
        3: 'Mix',
    }

    if('cod_mins' in comp_grid_data.keys()):
        cloud_var = 'cod'
    else:
        cloud_var = 'ch7'
    plt.close('all') 
    fig = plt.figure(figsize = (7.5, 3.5))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
   
    cmap = 'bwr' 
    shrk = 1.0
    cmap2 = plt.get_cmap('jet')
    colorvals = np.arange(0, 6, 1)
    norm = cm.BoundaryNorm(colorvals, cmap2.N, extend = 'upper')
        
    mesh = ax1.pcolormesh(comp_grid_data['sza_mins'], comp_grid_data[cloud_var + '_mins'], \
        comp_grid_data[dtype + '_slopes'][ice_idx,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    cbar = plt.colorbar(mesh,\
        ax = ax1, orientation='vertical',shrink = shrk, extend = 'both')
    cbar.set_label('AI - SWF slope [Wm$^{-2}$/AI]')
    ax1.set_xlabel('SZA')
    ax1.set_ylabel(cloud_var.upper())
    ax1.set_title(dtype.title() + ' ' + label_dict[ice_idx] + ' Forcing')


    mesh = ax2.pcolormesh(comp_grid_data['sza_mins'], comp_grid_data[cloud_var + '_mins'], \
        comp_grid_data[dtype + '_stderr'][ice_idx,:,:], \
        norm = norm, cmap = 'jet', shading = 'auto')
    cbar = plt.colorbar(ScalarMappable(norm = norm, cmap = cmap2),\
        ax = ax2, orientation='vertical',shrink = shrk, extend = 'both')
    cbar.set_label('AI - SWF Slope Standard Error')
    ax2.set_xlabel('SZA')
    ax2.set_ylabel(cloud_var.upper())
    ax2.set_title(dtype.title() + ' ' + label_dict[ice_idx] + ' Slope Error')

    
    if(comp_grid_data['smoother'] == 'None'):
        title_str = 'Unsmoothed'
        out_adder = ''
    elif(comp_grid_data['smoother'] == 'smooth'):
        title_str = '1-Bin Smoothed'
        out_adder = '_smooth1bin'
    elif(comp_grid_data['smoother'] == 'smoother'):
        title_str = '2-Bin Smoothed'
        out_adder = '_smooth2bin'

    plt.suptitle(title_str)

    fig.tight_layout()
    if(save):
        outname = 'comp_grid_climo_stderr_' + dtype + out_adder + '_' + label_dict[ice_idx] + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

# dtype: 'raw','grid'
def plot_compare_grid_climo_counts(comp_grid_data,  \
        dtype, ice_idx, ax1 = None, ax2 = None, save = False, \
        hatch_counts = False, hatch_style = '///'):

    vmin = -15
    vmax = 15
   
    label_dict = {
        0: 'Ocean',
        1: 'Ice',\
        2: 'Land',\
        3: 'Mix',
    }

    if('cod_mins' in comp_grid_data.keys()):
        cloud_var = 'cod'
    else:
        cloud_var = 'ch7'

    in_ax = True
    if((ax1 is None) and (ax2 is None)):
        in_ax = False
        plt.close('all') 
        fig = plt.figure(figsize = (7.5, 3.5))
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)
   
    cmap = 'bwr' 
    shrk = 1.0
    cmap2 = plt.get_cmap('viridis')
        
    mesh = ax1.pcolormesh(comp_grid_data['sza_mins'], \
        comp_grid_data[cloud_var + '_mins'], \
        comp_grid_data[dtype + '_slopes'][ice_idx,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    cbar = plt.colorbar(mesh,\
        ax = ax1, orientation='vertical',shrink = shrk, extend = 'both')
    cbar.set_label('AI - SWF slope [Wm$^{-2}$/AI]')
    ax1.set_xlabel('SZA')
    ax1.set_ylabel(cloud_var.upper())
    ax1.set_title(dtype.title() + ' ' + label_dict[ice_idx] + ' Forcing')


    mesh = ax2.pcolormesh(comp_grid_data['sza_mins'], \
        comp_grid_data[cloud_var + '_mins'], \
        comp_grid_data[dtype + '_counts'][ice_idx,:,:], \
        cmap = 'viridis', shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax2, orientation='vertical', \
        shrink = shrk, extend = 'both')
    cbar.set_label('Bin Counts')
    ax2.set_xlabel('SZA')
    ax2.set_ylabel(cloud_var.upper())
    ax2.set_title(dtype.title() + ' ' + label_dict[ice_idx] + ' Bin Counts')

    if(hatch_counts):
        hasher = np.ma.masked_where(comp_grid_data[dtype + '_counts'][ice_idx,:,:] > 20, \
            comp_grid_data[dtype + '_counts'][ice_idx,:,:])
        ax2.pcolor(comp_grid_data['sza_mins'], comp_grid_data[cloud_var + '_mins'], \
            hasher, hatch = hatch_style, alpha=0., shading = 'auto')

    
    if(comp_grid_data['smoother'] == 'None'):
        title_str = 'Unsmoothed'
        out_adder = ''
    elif(comp_grid_data['smoother'] == 'smooth'):
        title_str = '1-Bin Smoothed'
        out_adder = '_smooth1bin'
    elif(comp_grid_data['smoother'] == 'smoother'):
        title_str = '2-Bin Smoothed'
        out_adder = '_smooth2bin'

    plt.suptitle(title_str)

    if(not in_ax):
        fig.tight_layout()
        if(save):
            outname = 'comp_grid_climo_counts_' + dtype + out_adder + '_' + label_dict[ice_idx] + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show()

def plot_compare_grid_climo_counts_3version(comp_data1, comp_data2, \
        comp_data3, \
        dtype, ice_idx, save = False, \
        hatch_counts = False, hatch_style = '///'):

    label_dict = {
        0: 'Ocean',
        1: 'Ice',\
        2: 'Land',\
        3: 'Mix',
    }

    plt.close('all')
    fig = plt.figure(figsize = (7.5, 9))
    ax1 = fig.add_subplot(3,2,1)
    ax2 = fig.add_subplot(3,2,2)
    ax3 = fig.add_subplot(3,2,3)
    ax4 = fig.add_subplot(3,2,4)
    ax5 = fig.add_subplot(3,2,5)
    ax6 = fig.add_subplot(3,2,6)

    plot_compare_grid_climo_counts(comp_data1,  \
        dtype, ice_idx, ax1 = ax1, ax2 = ax2, save = False, \
        hatch_counts = hatch_counts, hatch_style = hatch_style)
    plot_compare_grid_climo_counts(comp_data2,  \
        dtype, ice_idx, ax1 = ax3, ax2 = ax4, save = False, \
        hatch_counts = hatch_counts, hatch_style = hatch_style)
    plot_compare_grid_climo_counts(comp_data3,  \
        dtype, ice_idx, ax1 = ax5, ax2 = ax6, save = False, \
        hatch_counts = hatch_counts, hatch_style = hatch_style)

    if(comp_data1['smoother'] == 'None'):
        title_str = 'Unsmoothed'
        out_adder = ''
    elif(comp_data1['smoother'] == 'smooth'):
        title_str = '1-Bin Smoothed'
        out_adder = '_smooth1bin'
    elif(comp_data1['smoother'] == 'smoother'):
        title_str = '2-Bin Smoothed'
        out_adder = '_smooth2bin'


    fig.tight_layout()
    if(save):
        outname = 'comp_grid_climo_counts_3version_' + dtype + out_adder + \
            '_' + label_dict[ice_idx].lower() + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

# Plots the calculated slopes in the 2-d pcolormesh on one subplot.
# Then, on the other subplot, shows the scatter data for a desired
# bin value from the first subplot
def plot_compare_slopes_scatter(out_dict, combined_data, comp_grid_data, \
        x_idx, y_idx, dtype = 'raw', ice_idx = 0, ai_min = 2, \
        ai_max = None, show_trend = False, save = False):

    if('cod_mins' in out_dict.keys()):
        cloud_var = 'cod'
    else:
        cloud_var = 'ch7'

    #ice_mins = np.array([0,80,105])
    #ice_maxs = np.array([20,100,None])
    ice_mins = np.array([0,  80, 105, 21])
    ice_maxs = np.array([20,100,None, 79])

    sizer = out_dict['sizer']
    if(out_dict['smoother'] == 'None'):    
        if('cod_bins' in comp_grid_data.keys()):
            cloud_var = 'cod'
            ch7_mins = comp_grid_data['cod_bins'][0::sizer]
            ch7_maxs = comp_grid_data['cod_bins'][0::sizer]
        else:
            cloud_var = 'ch7'
            ch7_mins = comp_grid_data['ch7_bins'][0::sizer]
            ch7_maxs = comp_grid_data['ch7_bins'][0::sizer]
        
        sza_mins = comp_grid_data['sza_bins'][4::sizer]
        sza_maxs = comp_grid_data['sza_bins'][4::sizer]
    elif(out_dict['smoother'] == 'smooth'): 
        if('cod_bins' in comp_grid_data.keys()):
            cloud_var = 'cod'
            ch7_mins = comp_grid_data['cod_bins'][0:-1:sizer]
            ch7_maxs = comp_grid_data['cod_bins'][1::sizer]
        else:
            cloud_var = 'ch7'
            ch7_mins = comp_grid_data['ch7_bins'][0:-1:sizer]
            ch7_maxs = comp_grid_data['ch7_bins'][1::sizer]
        
        sza_mins = comp_grid_data['sza_bins'][4:-1:sizer]
        sza_maxs = comp_grid_data['sza_bins'][5::sizer]
    elif(out_dict['smoother'] == 'smoother'):
        if('cod_bins' in comp_grid_data.keys()):
            cloud_var = 'cod'
            ch7_mins = comp_grid_data['cod_bins'][0:-2:sizer]
            ch7_maxs = comp_grid_data['cod_bins'][2::sizer]
        else:
            cloud_var = 'ch7'
            ch7_mins = comp_grid_data['ch7_bins'][0:-2:sizer]
            ch7_maxs = comp_grid_data['ch7_bins'][2::sizer]
        
        sza_mins = comp_grid_data['sza_bins'][4:-2:sizer]
        sza_maxs = comp_grid_data['sza_bins'][6::sizer]
    
    ai_diff  = comp_grid_data['ai_bins'][2] - comp_grid_data['ai_edges'][2]
    sza_diff = comp_grid_data['sza_bins'][2] - comp_grid_data['sza_edges'][2]
    ice_diff = comp_grid_data['ice_bins'][2] - comp_grid_data['ice_edges'][2]
    ch7_diff = comp_grid_data[cloud_var + '_bins'][2] - \
        comp_grid_data[cloud_var + '_edges'][2]
    
    ## Calculate slopes 
    ## ----------------
    ##if((len(grid_xdata) > 20) & (len(raw_xdata) > 20)):
    #if(trend_type == 'theil-sen'):
    #    if((len(raw_xdata) > 20)):
    #        res = stats.theilslopes(raw_ydata, raw_xdata, 0.90)
    #        raw_slopes[ii,jj,kk]   = res[0]

    #    if((len(grid_xdata) > 20)):
    #        res = stats.theilslopes(grid_ydata, grid_xdata, 0.90)
    #        grid_slopes[ii,jj,kk]   = res[0]
    #elif(trend_type == 'linregress'):
    #    if((len(raw_xdata) > 20)):
    #        result1 = stats.linregress(raw_xdata,raw_ydata)
    #        raw_slopes[ii,jj,kk] = result1.slope 
    #        raw_stderr[ii,jj,kk] = result1.stderr
    #        raw_pvals[ii,jj,kk]  = result1.pvalue

    #    if((len(grid_xdata) > 20)):
    #        result2 = stats.linregress(grid_xdata,grid_ydata)
    #        grid_slopes[ii,jj,kk] = result2.slope 
    #        grid_stderr[ii,jj,kk] = result2.stderr
    #        grid_pvals[ii,jj,kk]  = result2.pvalue
    #else:
    #    print("WARNING: invalid trend_type specified. Using theil-sen")
    #    if((len(raw_xdata) > 20)):
    #        res = stats.theilslopes(raw_ydata, raw_xdata, 0.90)
    #        raw_slopes[ii,jj,kk]   = res[0]

    #    if((len(grid_xdata) > 20)):
    #        res = stats.theilslopes(grid_ydata, grid_xdata, 0.90)
    #        grid_slopes[ii,jj,kk]   = res[0]
   

    titles = ['Ocean','Ice','Land','Mix']
   
    plt.close('all') 
    fig = plt.figure(figsize = (7,6))
    ax1 = fig.add_subplot(2,2,1) # Raw pcolor
    ax2 = fig.add_subplot(2,2,2) # Raw scatter
    ax3 = fig.add_subplot(2,2,3) # Grid pcolor
    ax4 = fig.add_subplot(2,2,4) # Grid scatter

    vmin = -15
    vmax = 15 
    shrk = 1.0
    cmap = 'bwr'
    mesh = ax1.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        out_dict['raw_slopes'][ice_idx,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    cbar = plt.colorbar(mesh,\
        ax = ax1, orientation='vertical',shrink = shrk, extend = 'both')
    ax1.set_xlabel('SZA')
    ax1.set_ylabel(cloud_var.upper())
    ax1.set_title('Raw ' + titles[ice_idx])
    #if(hatch_stderr):
    #    hasher = np.ma.masked_where(out_dict['raw_stderr'][0,:,:] > maxerr, \
    #        out_dict['raw_stderr'][0,:,:])
    #    ax1.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
    #        hasher, hatch = hatch_style, alpha=0., shading = 'auto')

    print("SZA Mins:")
    print("    ",sza_mins[x_idx], sza_maxs[x_idx])
    print("CH7 Ranges:")
    print("    ",ch7_mins[y_idx],ch7_maxs[y_idx])
    ax1.plot(sza_mins[x_idx], ch7_mins[y_idx], 'o', color = 'tab:purple')


    return_dict = \
        plot_combined_scatter(combined_data, ax = ax2, \
        #cld_min = None, cld_max = None, \
        ch7_min = ch7_mins[y_idx], ch7_max = ch7_maxs[y_idx], \
        sza_min = sza_mins[x_idx], sza_max = sza_maxs[x_idx], \
        ice_min = ice_mins[ice_idx], ice_max = ice_maxs[ice_idx], \
        omi_min = ai_min, omi_max = ai_max, \
        ai_diff = ai_diff, sza_diff = sza_diff, \
        ice_diff = ice_diff, ch7_diff = ch7_diff,\
        cld_diff = None, \
        trend_type = out_dict['trend_type'], \
        show_trend = show_trend, shade_density = False,\
        cloud_var = cloud_var)

    mesh = ax3.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        out_dict['grid_slopes'][ice_idx,:,:], \
        cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
    cbar = plt.colorbar(mesh,\
        ax = ax3, orientation='vertical',shrink = shrk, extend = 'both')
    ax3.set_xlabel('SZA')
    ax3.set_ylabel(cloud_var.upper())
    ax3.set_title('Grid ' + titles[ice_idx])
    
    ax3.plot(sza_mins[x_idx], ch7_mins[y_idx], 'o', color = 'tab:purple')


    plot_comp_grid_scatter(comp_grid_data, ax = ax4, xval = 'ai', \
        ai_min = ai_min,  ai_max  = ai_max, \
        ch7_min = ch7_mins[y_idx], ch7_max = ch7_mins[y_idx], \
        sza_min = sza_mins[x_idx], sza_max = sza_maxs[x_idx], \
        ice_min = ice_mins[ice_idx], ice_max = ice_maxs[ice_idx], \
        #cld_min = cld_min, cld_max = cld_max, \
        show_trend = show_trend, trend_type = out_dict['trend_type'])


    fig.tight_layout()

    if(save):
        outname = 'comp_grid_scatter_compare_' + titles[ice_idx] + '_' + \
            str(x_idx) + 'x' + str(y_idx) + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)

    plt.show()

    return return_dict

def calc_pcnt_aerosol_over_type(date_list, min_AI, ax = None): 

    count_data  = np.full(len(date_list), np.nan)
    total_count = np.full(len(date_list), np.nan)
    total2_count = np.full(len(date_list), np.nan)
    pcnt_ice    = np.full(len(date_list), np.nan)
    pcnt_mix    = np.full(len(date_list), np.nan)
    pcnt_ocn    = np.full(len(date_list), np.nan)
    pcnt_lnd    = np.full(len(date_list), np.nan)
    pcnt_oth    = np.full(len(date_list), np.nan)
    pcnt_ice_cld    = np.full(len(date_list), np.nan)
    pcnt_mix_cld    = np.full(len(date_list), np.nan)
    pcnt_ocn_cld    = np.full(len(date_list), np.nan)
    pcnt_lnd_cld    = np.full(len(date_list), np.nan)
    pcnt_oth_cld    = np.full(len(date_list), np.nan)
    pcnt_ice_clr    = np.full(len(date_list), np.nan)
    pcnt_mix_clr    = np.full(len(date_list), np.nan)
    pcnt_ocn_clr    = np.full(len(date_list), np.nan)
    pcnt_lnd_clr    = np.full(len(date_list), np.nan)
    pcnt_oth_clr    = np.full(len(date_list), np.nan)
    
    for ii, date in enumerate(date_list):
        coloc_data = read_colocated(date)
    
        # Test finding the indices with high aerosol
        # ------------------------------------------
        
        mask_data = np.ma.masked_where(coloc_data['OMI_RAW'] < min_AI, coloc_data['OMI_RAW'])
        
        mask_ice = np.ma.masked_where((coloc_data['NSIDC_ICE'].mask == True), \
            mask_data)
        mask_io_mix = np.ma.masked_where((coloc_data['NSIDC_IO_MIX'].mask == True), \
            mask_data)
        mask_ocean = np.ma.masked_where((coloc_data['NSIDC_OCEAN'].mask == True), \
            mask_data)
        mask_land = np.ma.masked_where((coloc_data['NSIDC_LAND'].mask == True), \
            mask_data)
        mask_othr = np.ma.masked_where((coloc_data['NSIDC_OTHER'].mask == True), \
            mask_data)
       
        mask_ice_clr = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_ice)
        mask_ice_cld = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_ice)
         
        mask_mix_clr = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_io_mix)
        mask_mix_cld = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_io_mix)
         
        mask_ocn_clr = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_ocean)
        mask_ocn_cld = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_ocean)
         
        mask_lnd_clr = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_land)
        mask_lnd_cld = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_land)
     
        mask_oth_clr = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_othr)
        mask_oth_cld = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_othr)
     
        count_data_val = mask_data.compressed().size
        count_ice      = mask_ice.compressed().size
        count_mix      = mask_io_mix.compressed().size
        count_ocean    = mask_ocean.compressed().size
        count_land     = mask_land.compressed().size
        count_othr     = mask_othr.compressed().size

        count_cld_ice      = mask_ice_cld.compressed().size
        count_cld_mix      = mask_mix_cld.compressed().size
        count_cld_ocn      = mask_ocn_cld.compressed().size
        count_cld_lnd      = mask_lnd_cld.compressed().size
        count_cld_oth      = mask_oth_cld.compressed().size

        count_clr_ice      = mask_ice_clr.compressed().size
        count_clr_mix      = mask_mix_clr.compressed().size
        count_clr_ocn      = mask_ocn_clr.compressed().size
        count_clr_lnd      = mask_lnd_clr.compressed().size
        count_clr_oth      = mask_oth_clr.compressed().size

        totals         = count_ice + count_mix + count_ocean + count_land + \
            count_othr
        total2         = count_cld_ice + count_clr_ice + \
                         count_cld_mix + count_clr_mix + \
                         count_cld_ocn + count_clr_ocn + \
                         count_cld_lnd + count_clr_lnd + \
                         count_cld_oth + count_clr_oth
    
        pcnt_ice_val = np.round((count_ice / count_data_val) * 100., 3)
        pcnt_mix_val = np.round((count_mix / count_data_val) * 100., 3)
        pcnt_ocn_val = np.round((count_ocean / count_data_val) * 100., 3)
        pcnt_lnd_val = np.round((count_land / count_data_val) * 100., 3)
        pcnt_oth_val = np.round((count_othr / count_data_val) * 100., 3)
    
        count_data[ii]  = count_data_val
        pcnt_ice[ii]    = pcnt_ice_val 
        pcnt_mix[ii]    = pcnt_mix_val 
        pcnt_ocn[ii]    = pcnt_ocn_val 
        pcnt_lnd[ii]    = pcnt_lnd_val 
        pcnt_oth[ii]    = pcnt_oth_val 
        total_count[ii] = totals
        total2_count[ii] = total2

    print('    date      cnt    ice   mix    ocn    lnd    oth     tot')
    for ii, date in enumerate(date_list):
        print("%s %5d %6.3f %6.3f %6.3f %6.3f %6.3f %5d %5d" % \
            (date, count_data[ii], pcnt_ice[ii], pcnt_mix[ii], pcnt_ocn[ii], \
            pcnt_lnd[ii], pcnt_oth[ii], total_count[ii], total2_count[ii]))

    in_ax = True 
    if(ax is None): 
        plt.close('all')
        in_ax = False
        fig = plt.figure(figsize = (10, 4))
        ax = fig.add_subplot(1,1,1)

    xvals = np.arange(len(date_list))
    
    ax.bar(xvals, pcnt_ice, label = 'Ice')
    ax.bar(xvals, pcnt_mix, bottom = pcnt_ice, label = 'Mix')
    ax.bar(xvals, pcnt_ocn, bottom = pcnt_ice + pcnt_mix, label = 'Ocean')
    ax.bar(xvals, pcnt_lnd, bottom = pcnt_ice + pcnt_mix + pcnt_ocn, label = 'Land')
    ax.bar(xvals, pcnt_oth, bottom = pcnt_ice + pcnt_mix + pcnt_ocn + pcnt_lnd, label = 'Other')
    ax.legend()


    if(not in_ax):
        fig.tight_layout()
        plt.show() 

def calc_pcnt_aerosol_over_type_dayavgs(day_data, min_AI, ax = None, \
        area_calc = False, hatch_cloud = False, plot_map = False): 

    count_data   = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    total_count  = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    total2_count = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_ice     = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_mix     = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_ocn     = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_lnd     = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_oth     = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_ice_cld = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_mix_cld = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_ocn_cld = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_lnd_cld = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_oth_cld = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_ice_clr = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_mix_clr = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_ocn_clr = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_lnd_clr = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
    pcnt_oth_clr = np.full(day_data['omi_ai_raw'].shape[2], np.nan)
   
    base_ones   = np.full(day_data['omi_ai_raw'][:,:,0].shape, 1)
    base_twos   = np.full(day_data['omi_ai_raw'][:,:,0].shape, 2)
    base_threes = np.full(day_data['omi_ai_raw'][:,:,0].shape, 3)
    base_fours  = np.full(day_data['omi_ai_raw'][:,:,0].shape, 4)
    base_fives  = np.full(day_data['omi_ai_raw'][:,:,0].shape, 5)
    
    for ii in range(day_data['omi_ai_raw'].shape[2]):
    
        # Test finding the indices with high aerosol
        # ------------------------------------------
        mask_data = np.ma.masked_where(day_data['omi_ai_raw'][:,:,ii] < min_AI, \
            day_data['omi_ai_raw'][:,:,ii])
      
        mask_ice = np.ma.masked_where((\
            #data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.), \
            day_data['nsidc_ice'][:,:,ii] == -999.) | \
            (day_data['nsidc_ice'][:,:,ii] > 100.) | \
            (day_data['nsidc_ice'][:,:,ii] < 80.), \
            mask_data)
        #coloc_data['NSIDC_ICE'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        #    coloc_data['NSIDC_ICE'])
        mask_io_mix = np.ma.masked_where((\
            #data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.), \
            day_data['nsidc_ice'][:,:,ii] == -999.) | \
            (day_data['nsidc_ice'][:,:,ii] >= 80.) | \
            (day_data['nsidc_ice'][:,:,ii] <= 0.), \
            mask_data)
        #coloc_data['NSIDC_IOMIX'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        #    coloc_data['NSIDC_IO_MIX'])
        mask_ocean = np.ma.masked_where((\
            day_data['nsidc_ice'][:,:,ii] == -999.) | \
            (day_data['nsidc_ice'][:,:,ii] > 0.), \
            mask_data)
        #coloc_data['NSIDC_OCEAN'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        #    coloc_data['NSIDC_OCEAN'])
        mask_land = np.ma.masked_where((\
            day_data['nsidc_ice'][:,:,ii] == -999.) | \
            (day_data['nsidc_ice'][:,:,ii] != 254.), \
            mask_data)
        mask_othr = np.ma.masked_where((\
            day_data['nsidc_ice'][:,:,ii] == -999.) | \
            (day_data['nsidc_ice'][:,:,ii] < 251) | \
            (day_data['nsidc_ice'][:,:,ii] > 253), \
            mask_data) 
        #coloc_data['NSIDC_LAND'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        #    coloc_data['NSIDC_LAND'])

        if(plot_map):
            #plt.close('all')
            fig = plt.figure(figsize = (9, 4)) 
            #fig = plt.figure(figsize = (6, 6)) 
            #ax1 = fig.add_subplot(1,1,1, projection = mapcrs)
            ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
            ax2 = fig.add_subplot(1,2,2, projection = mapcrs)

            print(mask_ice.compressed().shape[0], mask_io_mix.compressed().shape[0], \
                mask_ocean.compressed().shape[0], mask_land.compressed().shape[0], \
                mask_othr.compressed().shape[0])

            if(mask_ice.compressed().shape[0] > 0):
                ones = np.ma.masked_where(mask_ice.mask == True, base_ones).T
                ax1.pcolormesh(day_data['longitude'][:], day_data['latitude'][:], \
                    ones, transform = datacrs, shading = 'auto', \
                    vmin = 1, vmax = 10, cmap = 'tab10')
            if(mask_io_mix.compressed().shape[0] > 0):
                twos = np.ma.masked_where(mask_io_mix.mask == True, base_twos).T
                ax1.pcolormesh(day_data['longitude'][:], day_data['latitude'][:], \
                    twos, transform = datacrs, shading = 'auto', \
                    vmin = 1, vmax = 10, cmap = 'tab10')
            if(mask_ocean.compressed().shape[0] > 0):
                threes = np.ma.masked_where(mask_ocean.mask == True, base_threes).T
                ax1.pcolormesh(day_data['longitude'][:], day_data['latitude'][:], \
                    threes, transform = datacrs, shading = 'auto', \
                    vmin = 1, vmax = 10, cmap = 'tab10')
            if(mask_land.compressed().shape[0] > 0):
                fours = np.ma.masked_where(mask_land.mask == True, base_fours).T
                ax1.pcolormesh(day_data['longitude'][:], day_data['latitude'][:], \
                    fours, transform = datacrs, shading = 'auto', \
                    vmin = 1, vmax = 10, cmap = 'tab10')
            if(mask_othr.compressed().shape[0] > 0):
                fives = np.ma.masked_where(mask_othr.mask == True, base_fives).T
                ax1.pcolormesh(day_data['longitude'][:], day_data['latitude'][:], \
                    fives, transform = datacrs, shading = 'auto', \
                    vmin = 1, vmax = 10, cmap = 'tab10')

            ax1.coastlines()
            ax1.set_extent([-180, 180, 65, 90], datacrs)

            masker = np.ma.masked_where(day_data['omi_ai_raw'][:,:,ii] < -20, \
                day_data['omi_ai_raw'][:,:,ii]).T
            ax2.pcolormesh(day_data['longitude'][:], day_data['latitude'][:], \
                masker, transform = datacrs, shading = 'auto', \
                cmap = 'jet')
            ax2.coastlines()
            ax2.set_extent([-180, 180, 65, 90], datacrs)
       
            plt.suptitle(str(day_data['dates'][ii]) + ' ' + str(ii))
 
            fig.tight_layout()

            plt.show()
    
        mask_ice_clr = np.ma.masked_where(day_data['modis_cld'][:,:,ii] <= 2, mask_ice)
        mask_ice_cld = np.ma.masked_where(day_data['modis_cld'][:,:,ii] > 2, mask_ice)
         
        mask_mix_clr = np.ma.masked_where(day_data['modis_cld'][:,:,ii] <= 2, mask_io_mix)
        mask_mix_cld = np.ma.masked_where(day_data['modis_cld'][:,:,ii] > 2, mask_io_mix)
         
        mask_ocn_clr = np.ma.masked_where(day_data['modis_cld'][:,:,ii] <= 2, mask_ocean)
        mask_ocn_cld = np.ma.masked_where(day_data['modis_cld'][:,:,ii] > 2, mask_ocean)
         
        mask_lnd_clr = np.ma.masked_where(day_data['modis_cld'][:,:,ii] <= 2, mask_land)
        mask_lnd_cld = np.ma.masked_where(day_data['modis_cld'][:,:,ii] > 2, mask_land)
     
        mask_oth_clr = np.ma.masked_where(day_data['modis_cld'][:,:,ii] <= 2, mask_othr)
        mask_oth_cld = np.ma.masked_where(day_data['modis_cld'][:,:,ii] > 2, mask_othr)
    
        if(area_calc):
            mask_data   = np.ma.masked_where(mask_data.mask == True, day_data['grid_areas'][:,:])
            mask_ice    = np.ma.masked_where(mask_ice.mask == True, day_data['grid_areas'][:,:])
            mask_io_mix = np.ma.masked_where(mask_io_mix.mask == True, day_data['grid_areas'][:,:])
            mask_ocean  = np.ma.masked_where(mask_ocean.mask == True, day_data['grid_areas'][:,:])
            mask_land   = np.ma.masked_where(mask_land.mask == True, day_data['grid_areas'][:,:])
            mask_othr   = np.ma.masked_where(mask_othr.mask == True, day_data['grid_areas'][:,:])

            mask_ice_clr = np.ma.masked_where(mask_ice_clr.mask == True, day_data['grid_areas'][:,:])
            mask_ice_cld = np.ma.masked_where(mask_ice_cld.mask == True, day_data['grid_areas'][:,:])

            mask_mix_clr = np.ma.masked_where(mask_mix_clr.mask == True, day_data['grid_areas'][:,:])
            mask_mix_cld = np.ma.masked_where(mask_mix_cld.mask == True, day_data['grid_areas'][:,:])

            mask_ocn_clr = np.ma.masked_where(mask_ocn_clr.mask == True, day_data['grid_areas'][:,:])
            mask_ocn_cld = np.ma.masked_where(mask_ocn_cld.mask == True, day_data['grid_areas'][:,:])

            mask_lnd_clr = np.ma.masked_where(mask_lnd_clr.mask == True, day_data['grid_areas'][:,:])
            mask_lnd_cld = np.ma.masked_where(mask_lnd_cld.mask == True, day_data['grid_areas'][:,:])

            mask_oth_clr = np.ma.masked_where(mask_oth_clr.mask == True, day_data['grid_areas'][:,:])
            mask_oth_cld = np.ma.masked_where(mask_oth_cld.mask == True, day_data['grid_areas'][:,:])

            count_data_val = sum(mask_data.compressed())
            count_ice      = sum(mask_ice.compressed())
            count_mix      = sum(mask_io_mix.compressed())
            count_ocean    = sum(mask_ocean.compressed())
            count_land     = sum(mask_land.compressed())
            count_othr     = sum(mask_othr.compressed())

            count_cld_ice  = sum(mask_ice_cld.compressed())
            count_cld_mix  = sum(mask_mix_cld.compressed())
            count_cld_ocn  = sum(mask_ocn_cld.compressed())
            count_cld_lnd  = sum(mask_lnd_cld.compressed())
            count_cld_oth  = sum(mask_oth_cld.compressed())

            count_clr_ice  = sum(mask_ice_clr.compressed())
            count_clr_mix  = sum(mask_mix_clr.compressed())
            count_clr_ocn  = sum(mask_ocn_clr.compressed())
            count_clr_lnd  = sum(mask_lnd_clr.compressed())
            count_clr_oth  = sum(mask_oth_clr.compressed())

        else:
    
            count_data_val = mask_data.compressed().size
            count_ice      = mask_ice.compressed().size
            count_mix      = mask_io_mix.compressed().size
            count_ocean    = mask_ocean.compressed().size
            count_land     = mask_land.compressed().size
            count_othr     = mask_othr.compressed().size

            count_cld_ice  = mask_ice_cld.compressed().size
            count_cld_mix  = mask_mix_cld.compressed().size
            count_cld_ocn  = mask_ocn_cld.compressed().size
            count_cld_lnd  = mask_lnd_cld.compressed().size
            count_cld_oth  = mask_oth_cld.compressed().size

            count_clr_ice  = mask_ice_clr.compressed().size
            count_clr_mix  = mask_mix_clr.compressed().size
            count_clr_ocn  = mask_ocn_clr.compressed().size
            count_clr_lnd  = mask_lnd_clr.compressed().size
            count_clr_oth  = mask_oth_clr.compressed().size

        totals         = count_ice + count_mix + count_ocean + count_land + \
            count_othr
        total2         = count_cld_ice + count_clr_ice + \
                         count_cld_mix + count_clr_mix + \
                         count_cld_ocn + count_clr_ocn + \
                         count_cld_lnd + count_clr_lnd + \
                         count_cld_oth + count_clr_oth
    
        pcnt_ice_val = np.round((count_ice / count_data_val) * 100., 3)
        pcnt_mix_val = np.round((count_mix / count_data_val) * 100., 3)
        pcnt_ocn_val = np.round((count_ocean / count_data_val) * 100., 3)
        pcnt_lnd_val = np.round((count_land / count_data_val) * 100., 3)
        pcnt_oth_val = np.round((count_othr / count_data_val) * 100., 3)


        pcnt_ice_clr_val = np.round((count_clr_ice / count_data_val) * 100., 3)
        pcnt_ice_cld_val = np.round((count_cld_ice / count_data_val) * 100., 3)
        pcnt_mix_clr_val = np.round((count_clr_mix / count_data_val) * 100., 3)
        pcnt_mix_cld_val = np.round((count_cld_mix / count_data_val) * 100., 3)
        pcnt_ocn_clr_val = np.round((count_clr_ocn / count_data_val) * 100., 3)
        pcnt_ocn_cld_val = np.round((count_cld_ocn / count_data_val) * 100., 3)
        pcnt_lnd_clr_val = np.round((count_clr_lnd / count_data_val) * 100., 3)
        pcnt_lnd_cld_val = np.round((count_cld_lnd / count_data_val) * 100., 3)
        pcnt_oth_clr_val = np.round((count_clr_oth / count_data_val) * 100., 3)
        pcnt_oth_cld_val = np.round((count_cld_oth / count_data_val) * 100., 3)
    
        count_data[ii]  = count_data_val
        pcnt_ice[ii]    = pcnt_ice_val 
        pcnt_mix[ii]    = pcnt_mix_val 
        pcnt_ocn[ii]    = pcnt_ocn_val 
        pcnt_lnd[ii]    = pcnt_lnd_val 
        pcnt_oth[ii]    = pcnt_oth_val 

        pcnt_ice_clr[ii]    = pcnt_ice_clr_val 
        pcnt_ice_cld[ii]    = pcnt_ice_cld_val 
        pcnt_mix_clr[ii]    = pcnt_mix_clr_val 
        pcnt_mix_cld[ii]    = pcnt_mix_cld_val 
        pcnt_ocn_clr[ii]    = pcnt_ocn_clr_val 
        pcnt_ocn_cld[ii]    = pcnt_ocn_cld_val 
        pcnt_lnd_clr[ii]    = pcnt_lnd_clr_val 
        pcnt_lnd_cld[ii]    = pcnt_lnd_cld_val 
        pcnt_oth_clr[ii]    = pcnt_oth_clr_val 
        pcnt_oth_cld[ii]    = pcnt_oth_cld_val 

        total_count[ii] = totals
        total2_count[ii] = total2

    if(area_calc):
        print_fmt = "%s %8d %6.3f %6.3f %6.3f %6.3f %6.3f %8d %8d"
    else:
        print_fmt = "%s %5d %6.3f %6.3f %6.3f %6.3f %6.3f %5d %5d"
    print('  date    cnt    ice   mix    ocn    lnd    oth    tot')
    for ii in range(day_data['omi_ai_raw'].shape[2]):
    
        print(print_fmt % \
            (day_data['dates'][ii], count_data[ii], pcnt_ice[ii], pcnt_mix[ii], pcnt_ocn[ii], \
            pcnt_lnd[ii], pcnt_oth[ii], total_count[ii], total2_count[ii]))

    in_ax = True 
    if(ax is None): 
        plt.close('all')
        in_ax = False
        fig = plt.figure(figsize = (10, 4))
        ax = fig.add_subplot(1,1,1)

    xvals = np.arange(day_data['omi_ai_raw'].shape[2])
   
    if(not hatch_cloud): 
        ax.bar(xvals, pcnt_ice, label = 'Ice')
        ax.bar(xvals, pcnt_mix, bottom = pcnt_ice, label = 'Mix')
        ax.bar(xvals, pcnt_ocn, bottom = pcnt_ice + pcnt_mix, label = 'Ocean')
        ax.bar(xvals, pcnt_lnd, bottom = pcnt_ice + pcnt_mix + pcnt_ocn, label = 'Land')
        ax.bar(xvals, pcnt_oth, bottom = pcnt_ice + pcnt_mix + pcnt_ocn + pcnt_lnd, \
            label = 'Other')
    else:
        ax.bar(xvals, pcnt_ice_clr, label = 'Ice')
        ax.bar(xvals, pcnt_ice_cld, bottom = pcnt_ice_clr, color = 'tab:blue', hatch = '///')
        ax.bar(xvals, pcnt_mix_clr, bottom = pcnt_ice_clr + pcnt_ice_cld, \
            label = 'Mix')
        ax.bar(xvals, pcnt_mix_cld, bottom = pcnt_ice_clr + pcnt_ice_cld + \
            pcnt_mix_clr, color = 'tab:orange', hatch = '///')
        ax.bar(xvals, pcnt_ocn_clr, bottom = pcnt_ice_clr + pcnt_ice_cld + \
            pcnt_mix_clr + pcnt_mix_cld, \
            label = 'Ocean')
        ax.bar(xvals, pcnt_ocn_cld, bottom = pcnt_ice_clr + pcnt_ice_cld + \
            pcnt_mix_clr + pcnt_mix_cld + pcnt_ocn_clr, \
            color = 'tab:green', hatch = '///')
        ax.bar(xvals, pcnt_lnd_clr, bottom = pcnt_ice_clr + pcnt_ice_cld + \
            pcnt_mix_clr + pcnt_mix_cld + pcnt_ocn_clr + pcnt_ocn_cld, \
            label = 'Land')
        ax.bar(xvals, pcnt_lnd_cld, bottom = pcnt_ice_clr + pcnt_ice_cld + \
            pcnt_mix_clr + pcnt_mix_cld + pcnt_ocn_clr + pcnt_ocn_cld + \
            pcnt_lnd_clr, \
            color = 'tab:red', hatch = '///')
        ax.bar(xvals, pcnt_oth_clr, bottom = pcnt_ice_clr + pcnt_ice_cld + \
            pcnt_mix_clr + pcnt_mix_cld + pcnt_ocn_clr + pcnt_ocn_cld + \
            pcnt_lnd_clr + pcnt_lnd_cld, \
            label = 'Other')
        ax.bar(xvals, pcnt_oth_cld, bottom = pcnt_ice_clr + pcnt_ice_cld + \
            pcnt_mix_clr + pcnt_mix_cld + pcnt_ocn_clr + pcnt_ocn_cld + \
            pcnt_lnd_clr + pcnt_lnd_cld + pcnt_oth_clr, \
            color = 'tab:purple', hatch = '///')
    ax.legend()

    if(not in_ax):
        fig.tight_layout()
        plt.show() 


def plot_aerosol_over_types(coloc_data, min_AI = 1.0, minlat = 70, \
        ai_val = 'OMI_RAW', save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data)
    else:
        dt_date_str = datetime.strptime(coloc_data['date_str'], '%Y%m%d%H%M')

    mask_data = np.ma.masked_where(coloc_data[ai_val] < min_AI, \
        coloc_data[ai_val])
    
    mask_ice = np.ma.masked_where((coloc_data['NSIDC_ICE'].mask == True), \
        mask_data)
    mask_mix = np.ma.masked_where((coloc_data['NSIDC_IO_MIX'].mask == True), \
        mask_data)
    mask_ocn = np.ma.masked_where((coloc_data['NSIDC_OCEAN'].mask == True), \
        mask_data)
    mask_lnd = np.ma.masked_where((coloc_data['NSIDC_LAND'].mask == True), \
        mask_data)

    mask_data_clear = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_data)
    mask_data_cloud = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_data)

    mask_ice_clear = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_ice)
    mask_ice_cloud = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_ice)
     
    mask_mix_clear = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_mix)
    mask_mix_cloud = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_mix)
     
    mask_ocn_clear = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_ocn)
    mask_ocn_cloud = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_ocn)
     
    mask_lnd_clear = np.ma.masked_where(coloc_data['MODIS_CLD'] <= 2, mask_lnd)
    mask_lnd_cloud = np.ma.masked_where(coloc_data['MODIS_CLD'] > 2, mask_lnd)

    fig = plt.figure(figsize = (9, 12))
    ax1  = fig.add_subplot(5,3,1 , projection = mapcrs)
    ax2  = fig.add_subplot(5,3,2 , projection = mapcrs)
    ax3  = fig.add_subplot(5,3,3 , projection = mapcrs)
    ax4  = fig.add_subplot(5,3,4 , projection = mapcrs)
    ax5  = fig.add_subplot(5,3,5 , projection = mapcrs)
    ax6  = fig.add_subplot(5,3,6 , projection = mapcrs)
    ax7  = fig.add_subplot(5,3,7 , projection = mapcrs)
    ax8  = fig.add_subplot(5,3,8 , projection = mapcrs)
    ax9  = fig.add_subplot(5,3,9 , projection = mapcrs)
    ax10 = fig.add_subplot(5,3,10, projection = mapcrs)
    ax11 = fig.add_subplot(5,3,11, projection = mapcrs)
    ax12 = fig.add_subplot(5,3,12, projection = mapcrs)
    ax13 = fig.add_subplot(5,3,13, projection = mapcrs)
    ax14 = fig.add_subplot(5,3,14, projection = mapcrs)
    ax15 = fig.add_subplot(5,3,15, projection = mapcrs)

    plt.suptitle(dt_date_str.strftime('%Y%m%d%H%M'))

    ax1.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_data, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax2.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_data_clear, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax3.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_data_cloud, \
        transform = datacrs, shading = 'auto', cmap = 'jet')

    ax4.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_ice, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax5.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_ice_clear, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax6.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_ice_cloud, \
        transform = datacrs, shading = 'auto', cmap = 'jet')

    ax7.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_mix, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax8.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_mix_clear, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax9.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_mix_cloud, \
        transform = datacrs, shading = 'auto', cmap = 'jet')

    ax10.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_ocn, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax11.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_ocn_clear, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax12.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_ocn_cloud, \
        transform = datacrs, shading = 'auto', cmap = 'jet')

    ax13.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_lnd, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax14.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_lnd_clear, \
        transform = datacrs, shading = 'auto', cmap = 'jet')
    ax15.pcolormesh(coloc_data['LON'], coloc_data['LAT'], mask_lnd_cloud, \
        transform = datacrs, shading = 'auto', cmap = 'jet')

    ax1.set_title('Total AI')
    ax2.set_title('Total Clear AI')
    ax3.set_title('Total Cloud AI')

    ax4.set_title('Ice AI')
    ax5.set_title('Ice Clear AI')
    ax6.set_title('Ice Cloud AI')

    ax7.set_title('Mix AI')
    ax8.set_title('Mix Clear AI')
    ax9.set_title('Mix Cloud AI')

    ax10.set_title('Ocean AI')
    ax11.set_title('Ocean Clear AI')
    ax12.set_title('Ocean Cloud AI')

    ax13.set_title('Land AI')
    ax14.set_title('Land Clear AI')
    ax15.set_title('Land Cloud AI')

    ax1.coastlines()
    ax2.coastlines()
    ax3.coastlines()
    ax4.coastlines()
    ax5.coastlines()
    ax6.coastlines()
    ax7.coastlines()
    ax8.coastlines()
    ax9.coastlines()
    ax10.coastlines()
    ax11.coastlines()
    ax12.coastlines()
    ax13.coastlines()
    ax14.coastlines()
    ax15.coastlines()

    ax1.set_extent( [-180, 180, minlat, 90], datacrs)
    ax2.set_extent( [-180, 180, minlat, 90], datacrs)
    ax3.set_extent( [-180, 180, minlat, 90], datacrs)
    ax4.set_extent( [-180, 180, minlat, 90], datacrs)
    ax5.set_extent( [-180, 180, minlat, 90], datacrs)
    ax6.set_extent( [-180, 180, minlat, 90], datacrs)
    ax7.set_extent( [-180, 180, minlat, 90], datacrs)
    ax8.set_extent( [-180, 180, minlat, 90], datacrs)
    ax9.set_extent( [-180, 180, minlat, 90], datacrs)
    ax10.set_extent([-180, 180, minlat, 90], datacrs)
    ax11.set_extent([-180, 180, minlat, 90], datacrs)
    ax12.set_extent([-180, 180, minlat, 90], datacrs)
    ax13.set_extent([-180, 180, minlat, 90], datacrs)
    ax14.set_extent([-180, 180, minlat, 90], datacrs)
    ax15.set_extent([-180, 180, minlat, 90], datacrs)

    fig.tight_layout()

    date_str = dt_date_str.strftime('%Y%m%d%H%M')

    fig2 = plt.figure()
    ax21 = fig2.add_subplot(2,2,1, projection = mapcrs)
    ax22 = fig2.add_subplot(2,2,2, projection = mapcrs)
    ax23 = fig2.add_subplot(2,2,3, projection = mapcrs)
    ax24 = fig2.add_subplot(2,2,4, projection = mapcrs)

    plot_spatial(ax21, coloc_data['LON'], coloc_data['LAT'], coloc_data[ai_val], \
        date_str, cmap = 'jet', zoom = False, minlat = minlat)
    plot_spatial(ax22, coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH1'], \
        date_str, cmap = 'Greys_r', zoom = False, minlat = minlat)
    plot_spatial(ax23, coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CLD'], \
        date_str, cmap = 'jet', zoom = False, minlat = minlat)
    plot_spatial(ax24, coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH7'], \
        date_str, cmap = 'Greys_r', zoom = False, minlat = minlat)

    fig2.tight_layout()

    plt.show()

def plot_aerosol_over_type_combined(data, dates, min_ai = 1.5, save = False, plot_map = False):

    fig = plt.figure(figsize = (9, 6))
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    #ax3 = fig.add_subplot(3,1,3)
    calc_pcnt_aerosol_over_type(dates, min_ai, ax = ax1)
    calc_pcnt_aerosol_over_type_dayavgs(data, min_ai, ax = ax2, area_calc = True, hatch_cloud = True, \
        plot_map = plot_map)
    #calc_pcnt_aerosol_over_type_dayavgs(data, 1.5, ax = ax3, area_calc = True, hatch_cloud = True)
    
    ax1.set_ylabel('Pcnt of OMI Pixels')
    #ax2.set_ylabel('Pcnt of 0.5 deg. Grid Boxes')
    ax2.set_ylabel('Pcnt oF Aerosol Area')
   
    plt.suptitle('Minimum AI = ' + str(min_ai))
 
    fig.tight_layout()

    if(save):
        outname = 'aerosol_over_type_minAI_' + str(int(min_ai*10)) + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

#def calc_pcnt_aerosol_mid_swath(date_list, min_AI):
#    
#    for ii, date in enumerate(date_list):
#        coloc_data = read_colocated(date)
#

def match_aeronet_to_grid_AI(data, aeronet_file = 'aeronet_site_info.txt', \
    min_ai = 1.5):

    tester = pd.read_csv(aeronet_file)

    mask_all_data = np.ma.masked_where((data['omi_ai_raw_count'][:,:,:] == 0) | \
        (data['omi_ai_raw'][:,:,:] < min_ai), data['omi_ai_raw'][:,:,:])

    grid_lat, grid_lon = np.meshgrid(data['latitude'][:], data['longitude'])

    # Figure out the grid indices for each AERONET site
    # ------------------------------------------------
    idx_dict = {}
    for ii in range(len(tester)):
        print(tester['Site_name'][ii], tester['Lat'][ii], tester['Lon'][ii])
        idxs = nearest_gridpoint(float(tester['Lat'][ii]), \
                                 float(tester['Lon'][ii]),\
                                 grid_lat, grid_lon)

        # Make sure the lats and lons are actually within the grid
        # --------------------------------------------------------
        if(    (abs(float(tester['Lat'][ii]) - grid_lat[idxs]) < 1.) & \
               (abs(float(tester['Lon'][ii]) - grid_lon[idxs]) < 1.)):
            print("This site is in the grid")
            idx_dict[tester['Site_name'][ii]] = idxs
        else:
            print("This site is NOT in the grid")

    # Loop over the dates
    # -------------------
    for ii in range(mask_all_data.shape[2]):
        print(data['dates'][ii])
        # Loop over each site 
        # -------------------
        for key in idx_dict.keys():
            # Check the mask of the grid value for this site
            # ----------------------------------------------
            #print("HERE:",mask_all_data[:,:,ii][idx_dict[key]])
            mask_val = mask_all_data[:,:,ii][idx_dict[key]]
            if(not isinstance(mask_val, list)):
                mask_val = [mask_val] 

            mask_val = np.ma.masked_array(mask_all_data[:,:,ii][idx_dict[key]]).mask
                
            if(mask_val == False):
                print('Aeronet data for site',key)


    return idx_dict, grid_lat, grid_lon

# dtype: 'cloud', or 'clear', or 'average'
def plot_type_forcing_all_months(OMI_data, NSIDC_data, dtype, minlat = 70., \
        maxlat = 87., use_szas = False, save = False, coloc_dict = None, \
        debug = False):

    full_forcings = np.full((6, OMI_data['LAT'].shape[0], \
            OMI_data['LAT'].shape[1]), np.nan)
    if(dtype == 'average'):
        full_cloud = np.full((6, OMI_data['LAT'].shape[0], \
                OMI_data['LAT'].shape[1]), np.nan)
        full_clear = np.full((6, OMI_data['LAT'].shape[0], \
                OMI_data['LAT'].shape[1]), np.nan)
    
    for month_idx in range(6):
        if(dtype == 'average'):
            full_cloud[month_idx,:,:] = \
                calculate_type_forcing(OMI_data, NSIDC_data, month_idx, \
                    'cloud', minlat = minlat, maxlat = maxlat, \
                     coloc_dict = coloc_dict, debug = debug)
            full_clear[month_idx,:,:] = \
                calculate_type_forcing(OMI_data, NSIDC_data, month_idx, \
                    'clear', minlat = minlat, maxlat = maxlat, \
                    coloc_dict = coloc_dict, debug = debug)
        else:
            full_forcings[month_idx,:,:] = \
                calculate_type_forcing(OMI_data, NSIDC_data, month_idx, \
                    dtype, minlat = minlat, maxlat = maxlat, \
                    coloc_dict = coloc_dict, debug = debug)

    if(dtype == 'average'):
        full_forcings = \
            np.nanmean(np.ma.masked_invalid(\
            np.array([full_cloud, full_clear])), axis = 0)
    else:
        full_forcings = np.ma.masked_invalid(full_forcings)
 
    fig = plt.figure(figsize = (9, 12))
    
    ax1 = fig.add_subplot(4,3,1, projection = mapcrs)
    ax2 = fig.add_subplot(4,3,2, projection = mapcrs)
    ax3 = fig.add_subplot(4,3,3, projection = mapcrs)
    ax4 = fig.add_subplot(4,3,4, projection = mapcrs)
    ax5 = fig.add_subplot(4,3,5, projection = mapcrs)
    ax6 = fig.add_subplot(4,3,6, projection = mapcrs)

    ax7  = fig.add_subplot(4,3,7)
    ax8  = fig.add_subplot(4,3,8)
    ax9  = fig.add_subplot(4,3,9)
    ax10 = fig.add_subplot(4,3,10)
    ax11 = fig.add_subplot(4,3,11)
    ax12 = fig.add_subplot(4,3,12)
    
    mesh = ax1.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            full_forcings[0,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax1)
    ax1.coastlines()
    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_title('April')   
 
    mesh = ax2.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            full_forcings[1,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax2)
    ax2.coastlines()
    ax2.set_extent([-180,180,minlat,90], datacrs)
    ax2.set_title('May')   
    
    mesh = ax3.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            full_forcings[2,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax3)
    ax3.coastlines()
    ax3.set_extent([-180,180,minlat,90], datacrs)
    ax3.set_title('June')   
    
    mesh = ax4.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            full_forcings[3,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax4)
    ax4.coastlines()
    ax4.set_extent([-180,180,minlat,90], datacrs)
    ax4.set_title('July')   
    
    mesh = ax5.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            full_forcings[4,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax5)
    ax5.coastlines()
    ax5.set_extent([-180,180,minlat,90], datacrs)
    ax5.set_title('August')  
    
    mesh = ax6.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            full_forcings[5,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax6)
    ax6.coastlines()
    ax6.set_extent([-180,180,minlat,90], datacrs)
    ax6.set_title('September')
  
    # Plot line-graph zonal averages
    # ------------------------------
 
    # Calculate the average of each of the daily standard deviations,
    # following the methodology I found on :
    # https://www.statology.org/averaging-standard-deviations/
    # -----------------------------------------------------------------
    #out_dict['cld_frac_std'] = \
    #    np.sqrt((np.sum((all_counts_data - 1) * (all_std_data**2.), \
    #    axis = 0)) / \
    #    (np.sum(all_counts_data, axis = 0) - all_counts_data.shape[0]))
    zonal_avgs = np.nanmean(full_forcings, axis = 2)
    if(dtype == 'average'):
        #zonal_cloud_avgs = np.nanmean(full_cloud, axis = 2)
        #zonal_clear_avgs = np.nanmean(full_clear, axis = 2)
        zonal_cloud_stds = np.nanstd(full_cloud, axis = 2)
        zonal_clear_stds = np.nanstd(full_clear, axis = 2)

        # Average the cloud and clear data
        # --------------------------------
        #zonal_avgs = np.concatenate([zonal_cloud_avgs, zonal_clear_avgs])
        #zonal_avgs = np.nanmean(zonal_avgs, axis = 0)

        # Calculate the average of the standard deviations
        # ------------------------------------------------
        zonal_stds = np.sqrt( \
            (zonal_cloud_stds**2. + zonal_clear_stds**2.) / 2)
        #np.sqrt((np.sum((all_counts_data - 1) * (all_std_data**2.), axis = 0)) / \
        #        (np.sum(all_counts_data, axis = 0) - all_counts_data.shape[0]))
        

    else:
        zonal_stds = np.nanstd(full_forcings, axis = 2)

    ax7.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[0])
    ax7.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[0] + zonal_stds[0], linestyle = '--', color = 'tab:blue')
    ax7.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[0] - zonal_stds[0], linestyle = '--', color = 'tab:blue')
    ax7.axhline(0, linestyle = ':', color = 'k')
    ax7.set_title('April')
    ax7.set_xlabel('Latitude')    
    ax7.set_ylabel('Zonal Avg. Forcing')    
 
    ax8.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[1])
    ax8.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[1] + zonal_stds[1], linestyle = '--', color = 'tab:blue')
    ax8.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[1] - zonal_stds[1], linestyle = '--', color = 'tab:blue')
    ax8.axhline(0, linestyle = ':', color = 'k')
    ax8.set_title('May')
    ax8.set_xlabel('Latitude')    
    ax8.set_ylabel('Zonal Avg. Forcing')    
 
    ax9.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[2])
    ax9.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[2] + zonal_stds[2], linestyle = '--', color = 'tab:blue')
    ax9.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[2] - zonal_stds[2], linestyle = '--', color = 'tab:blue')
    ax9.axhline(0, linestyle = ':', color = 'k')
    ax9.set_title('June')
    ax9.set_xlabel('Latitude')    
    ax9.set_ylabel('Zonal Avg. Forcing')    
 
    ax10.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[3])
    ax10.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[3] + zonal_stds[3], linestyle = '--', color = 'tab:blue')
    ax10.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[3] - zonal_stds[3], linestyle = '--', color = 'tab:blue')
    ax10.axhline(0, linestyle = ':', color = 'k')
    ax10.set_title('July')
    ax10.set_xlabel('Latitude')    
    ax10.set_ylabel('Zonal Avg. Forcing')    
 
    ax11.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[4])
    ax11.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[4] + zonal_stds[4], linestyle = '--', color = 'tab:blue')
    ax11.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[4] - zonal_stds[4], linestyle = '--', color = 'tab:blue')
    ax11.axhline(0, linestyle = ':', color = 'k')
    ax11.set_title('August')
    ax11.set_xlabel('Latitude')    
    ax11.set_ylabel('Zonal Avg. Forcing')    
 
    ax12.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[5])
    ax12.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[5] + zonal_stds[5], linestyle = '--', color = 'tab:blue')
    ax12.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[5] - zonal_stds[5], linestyle = '--', color = 'tab:blue')
    ax12.axhline(0, linestyle = ':', color = 'k')
    ax12.set_title('September')
    ax12.set_xlabel('Latitude')    
    ax12.set_ylabel('Zonal Avg. Forcing')    

    min_val = np.min([np.min(ax7.get_ylim()),  np.min(ax8.get_ylim()), \
              np.min(ax9.get_ylim()),  np.min(ax10.get_ylim()),
              np.min(ax11.get_ylim()), np.min(ax12.get_ylim())])
    max_val = np.max([np.max(ax7.get_ylim()),  np.max(ax8.get_ylim()), \
            np.max(ax9.get_ylim()),  np.max(ax10.get_ylim()),
            np.max(ax11.get_ylim()), np.max(ax12.get_ylim())])

    rangers = [min_val, max_val]
    ax7.set_ylim(rangers)
    ax8.set_ylim(rangers)
    ax9.set_ylim(rangers)
    ax10.set_ylim(rangers)
    ax11.set_ylim(rangers)
    ax12.set_ylim(rangers)
 
    plt.suptitle(dtype)
    fig.tight_layout()

    if(save):
        outname = 'calc_arctic_forcing_' + dtype + '_mix.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else: 
        plt.show()


def plot_type_forcing_v2_all_months(OMI_data, NSIDC_data, MYD08_data, \
        coloc_dict, minlat = 70., maxlat = 87., use_szas = False, \
        ai_thresh = -0.15, cld_idx = 0, maxerr = 2, \
        min_cloud = 0.95, data_type = 'raw', save = False,debug = False):

    forcing_trends = np.full((6, OMI_data['LAT'].shape[0], \
            OMI_data['LAT'].shape[1]), np.nan)
    forcing_pval   = np.full((6, OMI_data['LAT'].shape[0], \
            OMI_data['LAT'].shape[1]), np.nan)
    forcing_uncert = np.full((6, OMI_data['LAT'].shape[0], \
            OMI_data['LAT'].shape[1]), np.nan)
  
    for ii in range(6): 
        # Find the estimated forcing values for the current month series
        # --------------------------------------------------------------
        estimated_forcing = \
            calculate_type_forcing_v2(OMI_data, NSIDC_data, MYD08_data, \
            coloc_dict, ii, minlat = minlat, maxlat = maxlat, \
            ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr, \
            min_cloud = min_cloud, data_type = data_type)

        # Calculate the trend in the forcings
        # -----------------------------------
        forcing_trends[ii,:,:], forcing_pval[ii,:,:], \
            forcing_uncert[ii,:,:] = calc_forcing_grid_trend(\
            estimated_forcing, 'standard')

    fig = plt.figure(figsize = (9, 12))
    
    ax1 = fig.add_subplot(4,3,1, projection = mapcrs)
    ax2 = fig.add_subplot(4,3,2, projection = mapcrs)
    ax3 = fig.add_subplot(4,3,3, projection = mapcrs)
    ax4 = fig.add_subplot(4,3,4, projection = mapcrs)
    ax5 = fig.add_subplot(4,3,5, projection = mapcrs)
    ax6 = fig.add_subplot(4,3,6, projection = mapcrs)

    ax7  = fig.add_subplot(4,3,7)
    ax8  = fig.add_subplot(4,3,8)
    ax9  = fig.add_subplot(4,3,9)
    ax10 = fig.add_subplot(4,3,10)
    ax11 = fig.add_subplot(4,3,11)
    ax12 = fig.add_subplot(4,3,12)
    
    mesh = ax1.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            forcing_trends[0,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax1, label = 'Forcing [W/m2]')
    ax1.coastlines()
    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_title('April')   
 
    mesh = ax2.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            forcing_trends[1,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax2, label = 'Forcing [W/m2]')
    ax2.coastlines()
    ax2.set_extent([-180,180,minlat,90], datacrs)
    ax2.set_title('May')   
    
    mesh = ax3.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            forcing_trends[2,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax3, label = 'Forcing [W/m2]')
    ax3.coastlines()
    ax3.set_extent([-180,180,minlat,90], datacrs)
    ax3.set_title('June')   
    
    mesh = ax4.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            forcing_trends[3,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax4, label = 'Forcing [W/m2]')
    ax4.coastlines()
    ax4.set_extent([-180,180,minlat,90], datacrs)
    ax4.set_title('July')   
    
    mesh = ax5.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            forcing_trends[4,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax5, label = 'Forcing [W/m2]')
    ax5.coastlines()
    ax5.set_extent([-180,180,minlat,90], datacrs)
    ax5.set_title('August')  
    
    mesh = ax6.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
            forcing_trends[5,:,:], transform = datacrs, shading = 'auto', \
            vmin = -4, vmax = 4, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax6, label = 'Forcing [W/m2]')
    ax6.coastlines()
    ax6.set_extent([-180,180,minlat,90], datacrs)
    ax6.set_title('September')
  
    # Plot line-graph zonal averages
    # ------------------------------
 
    # Calculate the average of each of the daily standard deviations,
    # following the methodology I found on :
    # https://www.statology.org/averaging-standard-deviations/
    # -----------------------------------------------------------------
    #out_dict['cld_frac_std'] = \
    #    np.sqrt((np.sum((all_counts_data - 1) * (all_std_data**2.), \
    #    axis = 0)) / \
    #    (np.sum(all_counts_data, axis = 0) - all_counts_data.shape[0]))
    zonal_avgs = np.nanmean(forcing_trends, axis = 2)
    zonal_stds = np.nanstd(forcing_trends, axis = 2)

    ax7.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[0])
    ax7.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[0] + zonal_stds[0], linestyle = '--', color = 'tab:blue')
    ax7.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[0] - zonal_stds[0], linestyle = '--', color = 'tab:blue')
    ax7.axhline(0, linestyle = ':', color = 'k')
    ax7.set_title('April')
    ax7.set_xlabel('Latitude')    
    ax7.set_ylabel('Zonal Avg. Forcing')    
 
    ax8.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[1])
    ax8.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[1] + zonal_stds[1], linestyle = '--', color = 'tab:blue')
    ax8.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[1] - zonal_stds[1], linestyle = '--', color = 'tab:blue')
    ax8.axhline(0, linestyle = ':', color = 'k')
    ax8.set_title('May')
    ax8.set_xlabel('Latitude')    
    ax8.set_ylabel('Zonal Avg. Forcing')    
 
    ax9.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[2])
    ax9.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[2] + zonal_stds[2], linestyle = '--', color = 'tab:blue')
    ax9.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[2] - zonal_stds[2], linestyle = '--', color = 'tab:blue')
    ax9.axhline(0, linestyle = ':', color = 'k')
    ax9.set_title('June')
    ax9.set_xlabel('Latitude')    
    ax9.set_ylabel('Zonal Avg. Forcing')    
 
    ax10.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[3])
    ax10.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[3] + zonal_stds[3], linestyle = '--', color = 'tab:blue')
    ax10.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[3] - zonal_stds[3], linestyle = '--', color = 'tab:blue')
    ax10.axhline(0, linestyle = ':', color = 'k')
    ax10.set_title('July')
    ax10.set_xlabel('Latitude')    
    ax10.set_ylabel('Zonal Avg. Forcing')    
 
    ax11.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[4])
    ax11.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[4] + zonal_stds[4], linestyle = '--', color = 'tab:blue')
    ax11.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[4] - zonal_stds[4], linestyle = '--', color = 'tab:blue')
    ax11.axhline(0, linestyle = ':', color = 'k')
    ax11.set_title('August')
    ax11.set_xlabel('Latitude')    
    ax11.set_ylabel('Zonal Avg. Forcing')    
 
    ax12.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[5])
    ax12.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[5] + zonal_stds[5], linestyle = '--', color = 'tab:blue')
    ax12.plot(NSIDC_data['grid_lat'][:,0], zonal_avgs[5] - zonal_stds[5], linestyle = '--', color = 'tab:blue')
    ax12.axhline(0, linestyle = ':', color = 'k')
    ax12.set_title('September')
    ax12.set_xlabel('Latitude')    
    ax12.set_ylabel('Zonal Avg. Forcing')    

    min_val = np.min([np.min(ax7.get_ylim()),  np.min(ax8.get_ylim()), \
              np.min(ax9.get_ylim()),  np.min(ax10.get_ylim()),
              np.min(ax11.get_ylim()), np.min(ax12.get_ylim())])
    max_val = np.max([np.max(ax7.get_ylim()),  np.max(ax8.get_ylim()), \
            np.max(ax9.get_ylim()),  np.max(ax10.get_ylim()),
            np.max(ax11.get_ylim()), np.max(ax12.get_ylim())])

    rangers = [min_val, max_val]
    ax7.set_ylim(rangers)
    ax8.set_ylim(rangers)
    ax9.set_ylim(rangers)
    ax10.set_ylim(rangers)
    ax11.set_ylim(rangers)
    ax12.set_ylim(rangers)

    plt.suptitle('Forcing Estimate 2: Single Month Based', weight = 'bold', fontsize = 12) 
    fig.tight_layout()

    if(save):
        outname = 'calc_arctic_forcing_v2.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else: 
        plt.show()
