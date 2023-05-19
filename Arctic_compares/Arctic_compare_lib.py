"""
  NAME:

  PURPOSE:


"""
import os
home_dir = os.environ['HOME']
import sys
import json
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
        minlat = 65.):
    if(isinstance(date_str, str)):
        date_str = [date_str]

    if(download):
        # If desired, download all data and/or match up the files in the JSON
        # -------------------------------------------------------------------
        auto_all_download(date_str, download = download, rewrite_json = True)

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

            if(include_tropomi & (dt_date_str.year > 2017)):                    
                generate_TROPOMI_prep_data(dstr, copy_to_raindrop = \
                    copy_to_raindrop, minlat = minlat)

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
            write_MODIS_to_HDF5(MODIS_date, channel = 2, swath = True, \
                save_path = save_dir2, minlat = minlat, \
                remove_empty_scans = remove_empty_scans)
            write_MODIS_to_HDF5(MODIS_date, channel = 7, swath = True, \
                save_path = save_dir2, minlat = minlat, \
                remove_empty_scans = remove_empty_scans)

            NSIDC_date = file_date_dict[dstr]['NSIDC'][:8]
            print(NSIDC_date)
            writeNSIDC_to_HDF5(NSIDC_date, save_path = save_dir2, \
                minlat = minlat, remove_empty_scans = remove_empty_scans)

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

def single_wrap_function(date_str, minlat, min_AI, out_time_dict, download, \
        images, process, include_tropomi, new_only, copy_to_raindrop, \
        remove_empty_scans):

    print(date_str)

    good_list = download_identify_OMI_swaths(date_str, \
        minlat = minlat, min_AI = min_AI, remove_bad = True)

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
                remove_empty_scans = remove_empty_scans)

        else:
            print(ttime + " in json database. Reprocessing")
            automate_all_preprocess(ttime, download = download, \
                images = images, process = process, \
                omi_dtype = 'ltc3', include_tropomi = include_tropomi, \
                copy_to_raindrop = copy_to_raindrop, \
                minlat = minlat, \
                remove_empty_scans = remove_empty_scans)

def entire_wrapper(min_AI = 1.0, minlat = 70., new_only = True, \
        download = True, images = False, process = False, run_list = None, \
        include_tropomi = True, copy_to_raindrop = True, \
        remove_empty_scans = True):

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
                copy_to_raindrop, remove_empty_scans)

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
                    copy_to_raindrop, remove_empty_scans)
                

        #if((date_str == '20140812') | \
        #   (date_str == '20140814') | \
        #   (date_str == '20140816')):

        ##!#    print(date_str, ai_val)

        ##!#    good_list = download_identify_OMI_swaths(date_str, \
        ##!#        minlat = minlat, remove_bad = True)

        ##!#    final_good_list = final_good_list + good_list

        ##!#    for ttime in good_list:
        ##!#        if(new_only and (ttime not in out_time_dict.keys())):
        ##!#        #if(ttime not in out_time_dict.keys()):
        ##!#            print("NEW TIME: " + ttime + " not in json database. Run preprocessor")
        ##!#            automate_all_preprocess(ttime, download = download, \
        ##!#                images = images, process = process, \
        ##!#                omi_dtype = 'ltc3', include_tropomi = include_tropomi)

        ##!#        else:
        ##!#            print(ttime + " in json database. Reprocessing")
        ##!#            automate_all_preprocess(ttime, download = download, \
        ##!#                images = images, process = process, \
        ##!#                omi_dtype = 'ltc3', include_tropomi = include_tropomi)

    #return final_good_list

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
        omi_dtype = 'ltc3'):
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
                modis_date_list, file_list = download_MODIS_swath(CERES_date_str, \
                        dest_dir = home_dir + '/data/MODIS/Aqua/', download = download)
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
    
            ##!#print(date_str)
            ##!#print('    OMI - ', dstr)
            ##!#print('  CERES - ', CERES_date_str)
            ##!#print('  MODIS - ', *modis_date_list)
            ##!#print('  NSIDC - ', CERES_date_str[:8])

            ##!#print('   OMI file  - ', omi_shawn)
            ##!#print(' CERES file  - ', ceres_file)
            ##!#print(' MODIS files - ', file_list)
            ##!#print('  NSIDC file - ', nsidc_file)

            out_time_dict[dstr] = {}
            out_time_dict[dstr]['CERES'] = CERES_date_str
            out_time_dict[dstr]['MODIS'] = modis_date_list
            out_time_dict[dstr]['NSIDC'] = CERES_date_str[:8]

            ##out_file_dict[dstr] = {}
            ##out_file_dict[dstr]['OMI'] = omi_shawn
            ##out_file_dict[dstr]['CERES'] = ceres_file
            ##out_file_dict[dstr]['MODIS'] = file_list
            ##out_file_dict[dstr]['NSIDC'] = nsidc_file

            dt_date_str = datetime.strptime(modis_date_list[0], '%Y%m%d%H%M')
            local_date = dt_date_str.strftime('%Y-%m-%d')
            if(local_date not in out_file_dict.keys()):
                out_file_dict[local_date] = {}
            if('Lat' not in out_file_dict[local_date].keys()):
                out_file_dict[local_date]['Lat'] = [60., 90.]
            if('Lon' not in out_file_dict[local_date].keys()):
                out_file_dict[local_date]['Lon'] = [-180., 180.]
       
            for ttime, tfile in zip(modis_date_list, file_list):
                dt_date_str = datetime.strptime(ttime, '%Y%m%d%H%M')

                local_time = dt_date_str.strftime('%H%M')

                if(local_time not in out_file_dict[local_date].keys()):
                    out_file_dict[local_date][local_time] = {}

                # Make sure the MODIS files have the path attached
                # ------------------------------------------------
                if(len(tfile.split('/')) == 1):
                    tfile = '/home/bsorenson/data/MODIS/Aqua/' + tfile

                out_file_dict[local_date][local_time]['omi']   = omi_file
                out_file_dict[local_date][local_time]['ceres'] = ceres_file
                out_file_dict[local_date][local_time]['modis'] = tfile
                out_file_dict[local_date][local_time]['ceres_time'] = CERES_date_str[8:10]
                out_file_dict[local_date][local_time]['swath']      = modis_date_list
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
    coloc_data['MODIS_CH2'] = np.ma.masked_where((\
        data['modis_ch2'][:,:] == -999.) | (data['modis_ch2'][:,:] > 1.0), \
        data['modis_ch2'][:,:])
    coloc_data['MODIS_CH2'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['MODIS_CH2'])
    coloc_data['MODIS_CH7'] = np.ma.masked_where((\
        data['modis_ch7'][:,:] == -999.) | (data['modis_ch7'][:,:] > 1.0), \
        data['modis_ch7'])
    coloc_data['MODIS_CH7'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['MODIS_CH7'])
    coloc_data['NSIDC_ICE'] = np.ma.masked_where((\
        #data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.), \
        data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.) | \
        (data['nsidc_ice'][:,:] < 80.), \
        data['nsidc_ice'])
    coloc_data['NSIDC_ICE'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['NSIDC_ICE'])
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
    total_data['MODIS_CH2'] = np.ma.masked_where((\
        temp_data['MODIS_CH2'][:,:] == -999.) | (temp_data['MODIS_CH2'][:,:] > 1.0), \
        temp_data['MODIS_CH2'][:,:])
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
    comp_grid_data['sza_bins'] = data['omi_sza_bins'][:]
    comp_grid_data['ice_bins'] = data['nsidc_ice_bins'][:]
    comp_grid_data['ch7_bins'] = data['modis_ch7_bins'][:]
    comp_grid_data['cld_bins'] = data['ceres_cld_bins'][:]
    comp_grid_data['ai_edges'] = data['omi_ai_edges'][:]
    comp_grid_data['sza_edges'] = data['omi_sza_edges'][:]
    comp_grid_data['ice_edges'] = data['nsidc_ice_edges'][:]
    comp_grid_data['ch7_edges'] = data['modis_ch7_edges'][:]
    comp_grid_data['cld_edges'] = data['ceres_cld_edges'][:]

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
    #plot_MODIS_channel(modis_date_str, ch2, swath = True, \
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
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_spatial(ax, lon, lat, data, date_str, cmap = 'Greys_r', \
        plabel = '', colorbar = True, zoom = False, xmin = None):
  
    if(xmin is not None):
        data = np.ma.masked_where(data < xmin, data) 
     
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
        ax.set_extent([-180, 180, 60, 90], datacrs)

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
    plot_spatial(ax1, coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH2'], \
        pdate, cmap = 'Greys_r', zoom = zoom)
    ##!#ax1.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH2'], \
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
    plot_spatial(ax4, coloc_data['LON'], coloc_data['LAT'], coloc_data['OMI'], \
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
        fig1.savefig(outname, dpi = 300)
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
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname) 
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
        select_category(coloc_data, 'MODIS_CH2', cat), \
        #np.ma.masked_where(coloc_data[plot_var].mask, coloc_data['MODIS_CH2']), \
        pdate, cmap = 'Greys_r', zoom = zoom)
    ##!#ax1.pcolormesh(coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH2'], \
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
        fig1.savefig(outname, dpi = 300)
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
            fig.savefig(outname, dpi = 300)
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
        fig1.savefig(outname, dpi = 300)
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
    plot_spatial(ax1, coloc_data['LON'], coloc_data['LAT'], coloc_data['MODIS_CH2'], \
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
        fig1.savefig(outname, dpi = 300)
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
        fig.savefig(outname, dpi = 300)
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
        fig.savefig(outname, dpi = 300)
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

     
def calculate_type_forcing(month_idx, trend_type = 'linear', minlat = 65.):

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


def plot_combined_scatter(combined_data, xval = 'omi_uvai_raw', \
        yval = 'ceres_swf', ax = None, \
        swf_min = None, swf_max = None, \
        cld_min = None, cld_max = None, \
        ch7_min = None, ch7_max = None, \
        sza_min = None, sza_max = None, \
        ice_min = None, ice_max = None, \
        omi_min = None, omi_max = None, \
        trend_type = 'theil-sen', \
        show_trend = False, shade_density = False):

    local_xdata = np.copy(combined_data[xval])
    local_ydata = np.copy(combined_data[yval])

    def check_range_vals(local_xdata, local_ydata, lvar,\
            var_min, var_max):

        if(var_min is not None):
            local_xdata = np.ma.masked_where(\
                combined_data[lvar] < var_min,local_xdata)
            local_ydata = np.ma.masked_where(\
                combined_data[lvar] < var_min,local_ydata)
        if(var_max is not None):
            local_xdata = np.ma.masked_where(\
                combined_data[lvar] > var_max,local_xdata)
            local_ydata = np.ma.masked_where(\
                combined_data[lvar] > var_max,local_ydata)

        return local_xdata, local_ydata
        
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'ceres_cld', cld_min, cld_max)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'ceres_swf', swf_min, swf_max)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'modis_ch7', ch7_min, ch7_max)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'omi_sza', sza_min, sza_max)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'nsidc_ice', ice_min, ice_max)
    local_xdata, local_ydata = check_range_vals(local_xdata, local_ydata,\
        'omi_uvai_raw', omi_min, omi_max)

    local_xdata = local_xdata.compressed()
    local_ydata = local_ydata.compressed()

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


def plot_comp_grid_scatter(comp_grid_data, xval = 'ai', \
        cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None,\
        ax = None, show_trend = False , \
        trend_type = 'theil-sen', \
        save = False):

    #cld_lims  = [0., 20.]
    cld_lims  = [np.min(comp_grid_data['cld_bins']), np.max(comp_grid_data['cld_bins'])]
    #ch7_lims  = [0., 0.05]
    ch7_lims  = [np.min(comp_grid_data['ch7_bins']), np.max(comp_grid_data['ch7_bins'])]
    #ice_lims  = [101, 105]
    ice_lims  = [np.min(comp_grid_data['ice_bins']), np.max(comp_grid_data['ice_bins'])]
    #sza_lims  = [65, 70]
    sza_lims  = [np.min(comp_grid_data['sza_bins']), np.max(comp_grid_data['sza_bins'])]
    #ai_lims   = [2, 12]
    ai_lims   = [np.min(comp_grid_data['ai_bins']), np.max(comp_grid_data['ai_bins'])]
   
    if(cld_min is not None):
        cld_lims[0] = cld_min
    if(cld_max is not None):
        cld_lims[1] = cld_max
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
    keep_ch7_lims = np.where((comp_grid_data['ch7_bins'] >= ch7_lims[0]) & \
                             (comp_grid_data['ch7_bins'] <= ch7_lims[1]))[0]
    keep_cld_lims = np.where((comp_grid_data['cld_bins'] >= cld_lims[0]) & \
                             (comp_grid_data['cld_bins'] <= cld_lims[1]))[0]
    #keep_ice_lims  = np.arange(len(comp_grid_data['ice_bins']))
    #keep_ch7_lims  = np.arange(len(comp_grid_data['ch7_bins']))
    
    keep_swf = comp_grid_data['swf_climo'][keep_cld_lims[0]:keep_cld_lims[-1]+1,\
                                           keep_ch7_lims[0]:keep_ch7_lims[-1]+1,\
                                           keep_ice_lims[0]:keep_ice_lims[-1]+1,\
                                           keep_sza_lims[0]:keep_sza_lims[-1]+1,\
                                           keep_ai_lims[0]:keep_ai_lims[-1]+1]
    
    #keep_swf_new = np.array([[[[comp_grid_data['swf_climo'][ii,jj,kk,nn] \
    #    for nn in keep_ai_lims] \
    #    for kk in keep_sza_lims] \
    #    for jj in keep_ice_lims] \
    #    for ii in keep_ch7_lims])

    keep_cld = np.array([[[[[comp_grid_data['cld_bins'][mm] \
        for nn in keep_ai_lims] \
        for kk in keep_sza_lims] \
        for jj in keep_ice_lims] \
        for ii in keep_ch7_lims] \
        for mm in keep_cld_lims])

    keep_ch7 = np.array([[[[[comp_grid_data['ch7_bins'][ii] \
        for nn in keep_ai_lims] \
        for kk in keep_sza_lims] \
        for jj in keep_ice_lims] \
        for ii in keep_ch7_lims] \
        for mm in keep_cld_lims])
    
    keep_ice = np.array([[[[[comp_grid_data['ice_bins'][jj] \
        for nn in keep_ai_lims] \
        for kk in keep_sza_lims] \
        for jj in keep_ice_lims] \
        for ii in keep_ch7_lims] \
        for mm in keep_cld_lims])
    
    keep_sza = np.array([[[[[comp_grid_data['sza_bins'][kk] \
        for nn in keep_ai_lims] \
        for kk in keep_sza_lims] \
        for jj in keep_ice_lims] \
        for ii in keep_ch7_lims] \
        for mm in keep_cld_lims])
    
    keep_ai  = np.array([[[[[comp_grid_data['ai_bins'][nn] \
        for nn in keep_ai_lims] \
        for kk in keep_sza_lims] \
        for jj in keep_ice_lims] \
        for ii in keep_ch7_lims] \
        for mm in keep_cld_lims])
    

    final_cld = keep_cld.flatten()[~keep_swf.flatten().mask]
    final_ch7 = keep_ch7.flatten()[~keep_swf.flatten().mask]
    final_ice = keep_ice.flatten()[~keep_swf.flatten().mask]
    final_sza = keep_sza.flatten()[~keep_swf.flatten().mask]
    final_ai  = keep_ai.flatten()[~keep_swf.flatten().mask]
    final_swf = keep_swf.compressed()
  
    in_ax = True
    if(ax is None):
        plt.close('all')
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
       
    print('compressed size:',final_swf.shape)
 
    if(xval == 'ai'): 
        #ax.scatter(final_ai, final_swf, s = 6, color = 'k')
        xlabel = 'OMI UVAI Raw'
        local_xdata = final_ai
    elif(xval == 'sza'):
        #ax.scatter(final_sza, final_swf, s = 6, color = 'k')
        local_xdata = final_sza
        xlabel = 'OMI SZA'
    elif(xval == 'ice'):
        #ax.scatter(final_ice, final_swf, s = 6, color = 'k')
        local_xdata = final_ice
        xlabel = 'NSIDC ICE CONC.'
    elif(xval == 'ch7'):
        #ax.scatter(final_ch7, final_swf, s = 6, color = 'k')
        local_xdata = final_ch7
        xlabel = 'MODIS CH7 REFL.'
    elif(xval == 'cld'):
        #ax.scatter(final_ch7, final_swf, s = 6, color = 'k')
        local_xdata = final_cld
        xlabel = 'CERES CLD FRAC.'
    local_ydata = final_swf

    ax.scatter(local_xdata, local_ydata, s = 6, color = 'k')
    ax.set_xlabel(xlabel)

    if(show_trend):
        if(len(local_xdata) > 10000):
            print("ERROR: Plotting trend line with too many points")
            print("       Preventing this for the sake of computer safety")
        else:
            plot_trend_line(ax, local_xdata, local_ydata, color='tab:red', \
                linestyle = '-',  slope = trend_type)

    if(not in_ax):
        plt.show()

def plot_dual_combined_grid_climo(comp_grid_data, combined_data, \
        xval = 'ai', \
        cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None,\
        save = False, show_trend = False, trend_type = 'theil-sen'):

    plt.close('all')
    fig = plt.figure(figsize = (9, 4))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    try:
        plot_combined_scatter(combined_data, ax = ax1, \
            omi_min = ai_min,  omi_max  = ai_max, \
            sza_min = sza_min, sza_max = sza_max, \
            ice_min = ice_min, ice_max = ice_max, \
            ch7_min = ch7_min, ch7_max = ch7_max, \
            cld_min = cld_min, cld_max = cld_max, \
            trend_type = 'theil-sen', show_trend = show_trend)
        
        plot_comp_grid_scatter(comp_grid_data, ax = ax2, xval = 'ai', \
            ai_min = ai_min,  ai_max  = ai_max, \
            sza_min = sza_min, sza_max = sza_max, \
            ice_min = ice_min, ice_max = ice_max, \
            ch7_min = ch7_min, ch7_max = ch7_max, \
            cld_min = cld_min, cld_max = cld_max, \
            show_trend = show_trend)

        ax1.set_title('Raw Colocated L2')
        ax2.set_title('Binned Colocated L2')
        

        title = set_comp_title(cld_min, cld_max, \
                               ch7_min, ch7_max, \
                               ice_min, ice_max, \
                               sza_min, sza_max,\
                               ai_min,  ai_max)
        filer = set_file_title(cld_min, cld_max, \
                               ch7_min, ch7_max, \
                               ice_min, ice_max, \
                               sza_min, sza_max,\
                               ai_min,  ai_max)
        plt.suptitle(title)

        fig.tight_layout()

        if(save):
            outname = 'comp_dual_' + filer + '.png'
            fig.savefig(outname, dpi = 300)
            print("Saved image",outname)
        else:
            plt.show()
    except TypeError:
        print("ERROR: no data for these criteria")
   
def set_comp_title(\
        cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None):

    ch7_min = np.round(ch7_min, 2)
    ch7_max = np.round(ch7_max, 2)

    title = '' 
    if(cld_min is not None):
        if(cld_max is not None):
            title = title + str(cld_min) + ' < cld < ' + str(cld_max) + '\n'
        else:
            title = title + 'cld > ' + str(cld_min) + '\n'
    else:
        if(cld_max is not None):
            title = title + 'cld < ' + str(cld_max) + '\n'

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
        cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None):

    ch7_min = np.round(ch7_min, 2)
    ch7_max = np.round(ch7_max, 2)


    title = '' 
    if(cld_min is not None):
        cld_min = int(cld_min)
        if(cld_max is not None):
            cld_max = int(cld_max)
            title = title + str(cld_min) + 'cld' + str(cld_max) + '_'
        else:
            title = title + str(cld_min) + 'cld_'
    else:
        if(cld_max is not None):
            cld_max = int(cld_max)
            title = title + 'cld' + str(cld_max) + '_'

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
