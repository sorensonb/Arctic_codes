"""
  NAME:

  PURPOSE:


"""
import os
home_dir = os.environ['HOME']
import sys
import json
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import ImageGrid
sys.path.append(home_dir)
import python_lib
import importlib
from matplotlib.cm import turbo
import matplotlib.colors as mcolors
#from scipy.stats import sem
from scipy.stats import sem, norm as statnorm
from astropy.modeling import models,fitting
from sklearn.metrics import r2_score, mean_squared_error
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
from matplotlib.cm import ScalarMappable
import matplotlib.cm as cm
import matplotlib as mpl

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
        reprocess_only_omi = False, \
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


            if(not reprocess_only_omi):
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
        remove_empty_scans, remove_ch2_file, reprocess_only_omi, skiprows = [52]):

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
                reprocess_only_omi = reprocess_only_omi, \
                remove_ch2_file = remove_ch2_file)

        else:
            print(ttime + " in json database. Reprocessing")
            automate_all_preprocess(ttime, download = download, \
                images = images, process = process, \
                omi_dtype = 'ltc3', include_tropomi = include_tropomi, \
                copy_to_raindrop = copy_to_raindrop, \
                minlat = minlat, \
                remove_empty_scans = remove_empty_scans, \
                reprocess_only_omi = reprocess_only_omi, \
                remove_ch2_file = remove_ch2_file)

def entire_wrapper(min_AI = 1.0, minlat = 70., new_only = True, \
        download = True, images = False, process = False, run_list = None, \
        include_tropomi = True, copy_to_raindrop = True, \
        remove_empty_scans = True, remove_ch2_file = False, \
        reprocess_only_omi = False, \
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
                reprocess_only_omi, skiprows = skiprows)

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
                    reprocess_only_omi, skiprows = skiprows)
                
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
        compare_tropomi = False, local_dir = None):
  
    if(local_dir is None):
        local_dir = data_dir

    filename =  local_dir + 'colocated_subset_' + date_str + '.hdf5'
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
            filename =  local_dir + 'colocated_subset_' + date_str2 + '.hdf5'
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
    coloc_data['OMI_PERT'] = np.ma.masked_where(coloc_data['OMI_PERT'] < -15, \
        coloc_data['OMI_PERT'])
    coloc_data['OMI_PERT'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['OMI_PERT'])
    coloc_data['OMI_PERT'] = np.ma.masked_where(data['nsidc_ice'][:,:] == -999. , \
        coloc_data['OMI_PERT'])
    coloc_data['OMI_RAW']  = np.ma.masked_invalid(data['omi_uvai_raw'][:,:])
    #coloc_data['OMI_RAW'][:,23:53] = -9e9
    coloc_data['OMI_RAW'] = np.ma.masked_where(coloc_data['OMI_RAW'] < -15, \
        coloc_data['OMI_RAW'])
    coloc_data['OMI_RAW'] = np.ma.masked_where(data['nsidc_ice'][:,:] == -999. , \
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
        (data['nsidc_ice'][:,:] <= 20.), \
        data['nsidc_ice'])
    coloc_data['NSIDC_IOMIX'] = np.ma.masked_where(coloc_data['LAT'] < minlat, \
        coloc_data['NSIDC_IO_MIX'])
    coloc_data['NSIDC_OCEAN'] = np.ma.masked_where((\
        data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 20.), \
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


def read_daily_month_force_HDF5(infile):

    data = h5py.File(infile)

    out_dict = {}
    out_dict['TYPE_ADDER'] = infile.strip().split('/')[-1].split('.')[0].split('_')[4:]
    out_dict['LAT'] = data['latitude'][:]
    out_dict['LON'] = data['longitude'][:]
    out_dict['DATES'] = data['dates'][:]
    out_dict['FORCE_EST'] = data['force_estimate'][:,:,:]

    data.close()
    
    return out_dict

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Writing functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Writes a single Shawn file to HDF5. 
def write_daily_month_force_to_HDF5(all_month_values, OMI_monthly_data, \
        maxerr = None, ai_thresh = None, minlat = None, 
        dtype = None, 
        save_path = './', name_add = ''):

    # Create a new HDF5 dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'arctic_month_est_forcing' + name_add + '.hdf5'
    dset = h5py.File(outfile,'w')

    local_dates = np.array([int(tdate) for tdate in OMI_monthly_data['DATES']])
 
    cdt = dset.create_dataset('latitude',  data = OMI_monthly_data['LAT'][:,0].squeeze())
    cdt = dset.create_dataset('longitude', data = OMI_monthly_data['LON'][0,:].squeeze())
    cdt = dset.create_dataset('dates', data = local_dates)
    cdt = dset.create_dataset('force_estimate', data = all_month_values)
    if(maxerr is not None):
        cdt.attrs['maxerr'] = str(maxerr)
    if(ai_thresh is not None):
        cdt.attrs['ai_thresh'] = str(ai_thresh)
    if(minlat is not None):
        cdt.attrs['minlat'] = str(minlat)
    if(dtype is not None):
        cdt.attrs['omi_data_type'] = dtype

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  

# Read the HDF5 file generated by write_daily_month_force_L2L3_error_to_HDF5(
def read_daily_month_force_L2L3_error_from_HDF5(filename):

    data = h5py.File(filename, 'r')

    out_dict = {}

    out_dict['dates']          = data['dates'][:]
    out_dict['latitude']       = data['latitude'][:]
    out_dict['longitude']      = data['longitude'][:]
    if('force_estimate_orig' in data.keys()):
        out_dict['force_estimate_orig'] = data['force_estimate_orig'][:,:]
    if('force_estimate' in data.keys()):
        out_dict['force_estimate'] = data['force_estimate'][:]
    out_dict['attributes'] = {}

    if ('force_estimate' in data.keys()): 
        if( len(data['force_estimate'].attrs.keys()) != 0):
            for attr in data['force_estimate'].attrs.keys():
                out_dict['attributes'][attr] = data['force_estimate'].attrs[attr]
            """
            if('description' in data['force_estimate'].attrs.keys()):
                out_dict['description']     = data['force_estimate'].attrs['description']
            if('ai_thresh' in data['force_estimate'].attrs.keys()):
                out_dict['ai_thresh']     = float(data['force_estimate'].attrs['ai_thresh'])
            if('maxerr' in data['force_estimate'].attrs.keys()):
                out_dict['maxerr']        = float(data['force_estimate'].attrs['maxerr'])
            if('minlat' in data['force_estimate'].attrs.keys()):
                out_dict['minlat']        = float(data['force_estimate'].attrs['minlat'])
            if('L2L3_err_mean' in data['force_estimate'].attrs.keys()):
                out_dict['L2L3_err_mean'] = float(data['force_estimate'].attrs['L2L3_err_mean'])
            if('L2L3_err_std' in data['force_estimate'].attrs.keys()):
                out_dict['L2L3_err_std']  = float(data['force_estimate'].attrs['L2L3_err_std'])
            if('ice_err_mean' in data['force_estimate'].attrs.keys()):
                out_dict['ice_err_mean'] = float(data['force_estimate'].attrs['ice_err_mean'])
            if('ice_err_std' in data['force_estimate'].attrs.keys()):
                out_dict['ice_err_std']  = float(data['force_estimate'].attrs['ice_err_std'])
            """

    if ('force_estimate_orig' in data.keys()): 
        if( len(data['force_estimate_orig'].attrs.keys()) != 0):
            for attr in data['force_estimate_orig'].attrs.keys():
                out_dict['attributes'][attr] = data['force_estimate_orig'].attrs[attr]
            """
            if('description' in data['force_estimate_orig'].attrs.keys()):
                out_dict['orig_description']     = data['force_estimate_orig'].attrs['description']
            if('ai_thresh' in data['force_estimate_orig'].attrs.keys()):
                out_dict['ai_thresh']     = float(data['force_estimate_orig'].attrs['ai_thresh'])
            if('maxerr' in data['force_estimate_orig'].attrs.keys()):
                out_dict['maxerr']        = float(data['force_estimate_orig'].attrs['maxerr'])
            if('minlat' in data['force_estimate_orig'].attrs.keys()):
                out_dict['minlat']        = float(data['force_estimate_orig'].attrs['minlat'])
            if('L2L3_err_mean' in data['force_estimate_orig'].attrs.keys()):
                out_dict['L2L3_err_mean'] = float(data['force_estimate_orig'].attrs['L2L3_err_mean'])
            if('L2L3_err_std' in data['force_estimate_orig'].attrs.keys()):
                out_dict['L2L3_err_std']  = float(data['force_estimate_orig'].attrs['L2L3_err_std'])
            if('ice_err_mean' in data['force_estimate_orig'].attrs.keys()):
                out_dict['ice_err_mean'] = float(data['force_estimate_orig'].attrs['ice_err_mean'])
            if('ice3_err_std' in data['force_estimate_orig'].attrs.keys()):
                out_dict['ice_err_std']  = float(data['force_estimate_orig'].attrs['ice_err_std'])
            """

    if('bin_dict' in data.keys()):
        out_dict['bin_dict'] = {}
        out_dict['bin_dict']['cod_bin_edges'] = data['bin_dict/cod_bin_edges'][:]
        out_dict['bin_dict']['cod_bin_means'] = data['bin_dict/cod_bin_means'][:]
        out_dict['bin_dict']['sza_bin_edges'] = data['bin_dict/sza_bin_edges'][:]
        out_dict['bin_dict']['sza_bin_means'] = data['bin_dict/sza_bin_means'][:]
        out_dict['bin_dict']['ice_bin_edges'] = data['bin_dict/ice_bin_edges'][:]
        out_dict['bin_dict']['ice_bin_means'] = data['bin_dict/ice_bin_means'][:]

        if( len(data['bin_dict/ice_bin_means'].attrs.keys()) != 0):
            for attr in data['bin_dict/ice_min_means'].attrs.keys():
                out_dict['attributes'][attr] = data['ice_bin_means'].attrs[attr]
            
            """
            out_dict['ai_thresh']     = float(data['bin_dict/ice_bin_means'].attrs['ai_thresh'])
            out_dict['maxerr']        = float(data['bin_dict/ice_bin_means'].attrs['maxerr'])
            out_dict['minlat']        = float(data['bin_dict/ice_bin_means'].attrs['minlat'])
            out_dict['L2L3_err_mean'] = float(data['bin_dict/ice_bin_means'].attrs['L2L3_err_mean'])
            out_dict['L2L3_err_std']  = float(data['bin_dict/ice_bin_means'].attrs['L2L3_err_std'])
            """

    if('slope_dict' in data.keys()):
        out_dict['slope_dict'] = {}
        #out_dict['slope_dict']['ai_min'] = data['slope_dict/ai_min'] 
        #out_dict['slope_dict']['min_ob'] = data['slope_dict/min_ob'] 
        out_dict['slope_dict']['counts']        = data['slope_dict/counts'][:,:,:]
        out_dict['slope_dict']['slopess']       = data['slope_dict/intercepts'][:,:,:] 
        out_dict['slope_dict']['intercepts']    = data['slope_dict/intercepts'][:,:,:] 
        out_dict['slope_dict']['intcpt_stderr'] = data['slope_dict/intcpt_stderr'][:,:,:]
        out_dict['slope_dict']['slope_stderr']  = data['slope_dict/slope_stderr'][:,:,:] 
        out_dict['slope_dict']['slope_pvals']   = data['slope_dict/slope_pvals'][:,:,:] 

        if( len(data['slope_dict'].attrs.keys()) > 0):
            if('trend_type' in data['slope_dict'].attrs.keys()):
                out_dict['slope_dict']['trend_type'] = data['slope_dict'].attrs['trend_type']
            if('sim_name' in data['slope_dict'].attrs.keys()):
                out_dict['slope_dict']['sim_name'] = data['slope_dict'].attrs['sim_name']
            if('ai_min' in data['slope_dict'].attrs.keys()):
                out_dict['slope_dict']['ai_min'] = data['slope_dict'].attrs['ai_min']
            if('min_ob' in data['slope_dict'].attrs.keys()):
                out_dict['slope_dict']['min_ob'] = data['slope_dict'].attrs['min_ob']

    data.close()

    out_dict['filename'] = filename

    print_L2L3_sim_info(out_dict)

    return out_dict

def print_L2L3_sim_info(sim_dict):

    if('filename' in sim_dict.keys()):
        print('File: ', sim_dict['filename'])
    print('\nGlobal Attributes')
    for key in sim_dict['attributes'].keys():
        print('    ',key,': ',sim_dict['attributes'][key])
    if('maxerr' in sim_dict.keys()):
        print('    maxerr: ',  sim_dict['maxerr'])
    print('\nBin Dict Info:')
    print('    cod_bins: ', sim_dict['bin_dict']['cod_bin_edges'])
    print('    sza_bins: ', sim_dict['bin_dict']['sza_bin_edges'])
    print('    ice_bins: ', sim_dict['bin_dict']['ice_bin_edges'])
    print('\nSlope Dict Info:')
    if('ai_min' in sim_dict['slope_dict'].keys()):
        print('    ai_min: ',  sim_dict['slope_dict']['ai_min'])
    if('min_ob' in sim_dict['slope_dict'].keys()):
        print('    min_ob: ',  sim_dict['slope_dict']['min_ob'])
    if('trend_type' in sim_dict['slope_dict'].keys()):
        print('    trend_type: ',  sim_dict['slope_dict']['trend_type'])
    if('sim_name' in sim_dict['slope_dict'].keys()):
        print('    sim_name: ',  sim_dict['slope_dict']['sim_name'])

    # maxerr = 
    # ai_thresh = 
    # L2L3_err_mean  = 
    # L2L3_err_std   = 
    # use_intercepts = 
    # mod_slopes     = 
    # mod_intcpts    = 
    # mod_cod        = 
    # mod_ice        = 
    # cod_bins = [....]
    # sza_bins = [....]
    # ice_bins = [....]
    # slope_dict
    # - ai_min = 
    # - min_ob = 
    # - trend_type = 
    # - sim_name = 

# Writes a single Shawn file to HDF5. 
def write_daily_month_force_L2L3_error_to_HDF5(\
        all_month_values, OMI_monthly_data, \
        slope_dict, bin_dict, control_run_values, \
        minlat = None, maxlat = None, \
        maxerr = None, ai_thresh = None, \
        reference_ice = None, reference_cld = None, \
        L2L3_err_mean = None, L2L3_err_std = None, \
        ice_err_mean = None, ice_err_std = None, \
        cod_err_mean = None, cod_err_std = None, \
        calc_from_bckgd = None, \
        dtype = None, \
        overwrite_old_file = False, \
        write_daily_values = False, \
        OMI_daily_data = None, \
        save_path = './', name_add = ''):

    if(write_daily_values):
        outfile = save_path + 'arctic_daily_est_forcing' + name_add + '.hdf5'
        if(os.path.exists(outfile)):
            ii = 1
            while(os.path.exists(outfile)):
                outfile = save_path + 'arctic_daily_est_forcing' + name_add + '_v' + str(ii) + '.hdf5'
                ii += 1
    else:
        # Create a new HDF5 dataset to write to the file
        # ------------------------------------------------
        num_runs = int(all_month_values.shape[0])
        outfile = save_path + 'arctic_month_est_forcing_L2L3err_count' + str(num_runs) + name_add + '.hdf5'
        if(os.path.exists(outfile)):
            ii = 1
            while(os.path.exists(outfile)):
                outfile = save_path + 'arctic_month_est_forcing_L2L3err_count' + \
                    str(num_runs) + name_add + '_v' + str(ii) + '.hdf5'
                ii += 1

    print(outfile)

    dset = h5py.File(outfile,'w')

    if(not write_daily_values):
        local_dates = np.array([int(tdate) for tdate in OMI_monthly_data['DATES']])
    else:
        local_dates = OMI_daily_data['day_values']
 
    cdt = dset.create_dataset('latitude',  data = OMI_monthly_data['LAT'][:,0].squeeze())
    cdt = dset.create_dataset('longitude', data = OMI_monthly_data['LON'][0,:].squeeze())
    cdt = dset.create_dataset('dates', data = local_dates)

    if(write_daily_values):
        cdt = dset.create_dataset('force_estimate_orig', data = all_month_values)
        if( (reference_cld is not None) and (reference_ice is not None)):
                cdt.attrs['description'] = 'Daily forcing values with clouds set to ' + \
                    reference_cld + ' values and ice set to ' + reference_ice + ' values.'
        else:
            if(L2L3_err_mean is not None):
                cdt.attrs['description'] = 'Daily forcing values with L2L3 err modifications'
            elif(ice_err_mean is not None):
                cdt.attrs['description'] = 'Daily forcing values with ice err modifications'
            elif(cod_err_mean is not None):
                cdt.attrs['description'] = 'Daily forcing values with COD err modifications'
            elif(reference_cld is not None):
                cdt.attrs['description'] = 'Daily forcing values with clouds set to ' + \
                    reference_cld + ' values'
            elif(reference_ice is not None):
                cdt.attrs['description'] = 'Daily forcing values with ice set to ' + \
                    reference_ice + ' values'
            else:
                cdt.attrs['description'] = 'Daily forcing values with no modifications'
    else:

        cdt = dset.create_dataset('force_estimate_orig', data = control_run_values)
        cdt.attrs['description'] = 'Monthly forcing values (average of daily) with no modifications'
        cdt = dset.create_dataset('force_estimate', data = all_month_values)
        cdt.attrs['description'] = 'Monthly forcing values (average of daily) with L2L3 error mods'

    if(maxerr is not None):
        cdt.attrs['maxerr'] = str(maxerr)
    if(ai_thresh is not None):
        cdt.attrs['ai_thresh'] = str(ai_thresh)
    if(minlat is not None):
        cdt.attrs['minlat'] = str(minlat)
    if(dtype is not None):
        cdt.attrs['omi_data_type'] = dtype
    if(L2L3_err_mean is not None):
        cdt.attrs['L2L3_err_mean'] = L2L3_err_mean
    if(L2L3_err_std is not None):
        cdt.attrs['L2L3_err_std'] = L2L3_err_std
    if(ice_err_mean is not None):
        cdt.attrs['ice_err_mean'] = ice_err_mean
    if(ice_err_std is not None):
        cdt.attrs['ice_err_std'] = ice_err_std
    if(cod_err_mean is not None):
        cdt.attrs['cod_err_mean'] = cod_err_mean
    if(cod_err_std is not None):
        cdt.attrs['cod_err_std'] = cod_err_std
    if(reference_cld is not None):
        cdt.attrs['ref_cld'] = reference_cld
    if(reference_ice is not None):
        cdt.attrs['ref_ice'] = reference_ice
    if(calc_from_bckgd is not None):
        cdt.attrs['calc_from_bckgd'] = calc_from_bckgd


    grp = dset.create_group('slope_dict')
    cdt = grp.create_dataset('counts', data = slope_dict['counts'])
    cdt = grp.create_dataset('slopes', data = slope_dict['slopes'])
    cdt = grp.create_dataset('intercepts', data = slope_dict['intercepts'])
    cdt = grp.create_dataset('min_ob', data = slope_dict['min_ob'])
    cdt = grp.create_dataset('ai_min', data = slope_dict['ai_min'])
    cdt = grp.create_dataset('intcpt_stderr', data = slope_dict['intcpt_stderr'])
    cdt = grp.create_dataset('slope_stderr', data = slope_dict['slope_stderr'])
    cdt = grp.create_dataset('slope_pvals', data = slope_dict['slope_pvals'])
    grp.attrs['trend_type'] = slope_dict['trend_type']
    grp.attrs['sim_name'] = slope_dict['sim_name']
    grp.attrs['ai_min'] = slope_dict['ai_min']
    grp.attrs['min_ob'] = slope_dict['min_ob']

    grp = dset.create_group('bin_dict')
    cdt = grp.create_dataset('cod_bin_edges', data = bin_dict['cod_bin_edges'])
    cdt = grp.create_dataset('cod_bin_means', data = bin_dict['cod_bin_means'])
    cdt = grp.create_dataset('sza_bin_edges', data = bin_dict['sza_bin_edges'])
    cdt = grp.create_dataset('sza_bin_means', data = bin_dict['sza_bin_means'])
    cdt = grp.create_dataset('ice_bin_edges', data = bin_dict['ice_bin_edges'])
    cdt = grp.create_dataset('ice_bin_means', data = bin_dict['ice_bin_means'])

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  


# Writes a single Shawn file to HDF5. 
def write_monthly_force_trend_sims_to_HDF5(daily_dict, forcing_trends, \
        err_mean, err_std, save_path = './', name_add = ''):

    outfile = save_path + 'arctic_monthly_force_trends_count' + \
        str(int(forcing_trends.shape[0])) + name_add + '.hdf5'
    if(os.path.exists(outfile)):
        ii = 1
        while(os.path.exists(outfile)):
            outfile = save_path + 'arctic_monthly_force_trends_count' + \
                str(int(forcing_trends.shape[0])) + name_add + '_v' + str(ii) + '.hdf5'
            ii += 1

    print(outfile)

    dset = h5py.File(outfile,'w')

    cdt = dset.create_dataset('latitude',  data = daily_dict['latitude'][:])
    cdt = dset.create_dataset('longitude', data = daily_dict['longitude'][:])

    cdt = dset.create_dataset('force_trends', data = forcing_trends)
    cdt.attrs['description'] = 'Monthly forcing trends calculated using Monte Carlo error methods'

    cdt.attrs['err_mean'] = str(err_mean)
    cdt.attrs['err_std'] = str(err_std)
    cdt.attrs['num_sims'] = str(forcing_trends.shape[0])

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  

# Writes a single Shawn file to HDF5. 
def write_monthly_force_vals_sims_to_HDF5(daily_dict, forcing_values, \
        err_mean, err_std, save_path = './', name_add = ''):

    outfile = save_path + 'arctic_monthly_force_values_count' + \
        str(int(forcing_values.shape[0])) + name_add + '.hdf5'
    if(os.path.exists(outfile)):
        ii = 1
        while(os.path.exists(outfile)):
            outfile = save_path + 'arctic_monthly_force_values_count' + \
                str(int(forcing_values.shape[0])) + name_add + '_v' + str(ii) + '.hdf5'
            ii += 1

    print(outfile)

    dset = h5py.File(outfile,'w')

    cdt = dset.create_dataset('latitude',  data = daily_dict['latitude'][:])
    cdt = dset.create_dataset('longitude', data = daily_dict['longitude'][:])

    cdt = dset.create_dataset('monthly_force_vals', data = forcing_values)
    cdt.attrs['description'] = 'Monthly forcing trends calculated using Monte Carlo error methods'

    cdt.attrs['err_mean'] = err_mean
    cdt.attrs['err_std'] = err_std
    cdt.attrs['num_sims'] = forcing_values.shape[0]

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  





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
        #size = (14, 5)
        size = (7.5, 8.5)
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
    ax1 = fig1.add_subplot(3,2,1, projection = mapcrs)
    ax2 = fig1.add_subplot(3,2,2, projection = mapcrs)
    ax3 = fig1.add_subplot(3,2,3, projection = mapcrs)
    ax4 = fig1.add_subplot(3,2,4, projection = mapcrs)
    ax5 = fig1.add_subplot(3,2,5, projection = mapcrs)
    ax6 = fig1.add_subplot(3,2,6, projection = mapcrs)
    #ax1 = fig1.add_subplot(2,3,1, projection = mapcrs)
    #ax2 = fig1.add_subplot(2,3,2, projection = mapcrs)
    #ax3 = fig1.add_subplot(2,3,3, projection = mapcrs)
    #ax4 = fig1.add_subplot(2,3,4, projection = mapcrs)
    #ax5 = fig1.add_subplot(2,3,5, projection = mapcrs)
    #ax6 = fig1.add_subplot(2,3,6, projection = mapcrs)

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
    ax1.set_title('Aqua MODIS True Color')
    plot_MODIS_channel(modis_date, 'cloud_mask', swath = True, \
        zoom = zoom, ax = ax2, vmax = None)
    ax2.set_title('Aqua MODIS Cloud Mask')
    plot_MODIS_channel(modis_date, ch1, swath = True, \
        zoom = zoom, ax = ax3, vmax = 0.4)
    ax3.set_title('Aqua MODIS Ch7\n2.105 μm - 2.155 μm')
    #plot_MODIS_channel(modis_date_str, ch1, swath = True, \
    #    zoom = zoom, ax = ax3)

    # Plot the NSIDC data
    # -------------------
    plotNSIDC_daily_figure(nsidc_date, minlat = minlat, \
        zoom = zoom, ax = ax4, gridlines = False, \
        title = 'SSMI/S Sea Ice Concentration', save = False)

    # Plot the OMI data
    # -----------------
    plotOMI_single_swath_figure(omi_date, \
            dtype = omi_dtype, only_sea_ice = False, minlat = minlat, \
            ax = ax5, skiprows = [52], lat_circles = None, save = False, \
            zoom = zoom, shawn_path = shawn_path)
    ax5.set_title('OMI UVAI')   
 
    # Plot the CERES data
    # -------------------
    plotCERES_hrly_figure(ceres_date, 'SWF',  \
        minlat = minlat, lat_circles = None, ax = ax6, title = 'SWF',\
        grid_data = True, zoom = zoom, vmax = 450, vmin = None, save = False)
    ax6.set_title('Aqua CERES TOA SWF')   
    #plotCERES_hrly_figure(ceres_date, 'ALB',  \
    #    minlat = minlat, lat_circles = None, ax = ax6, title = 'ALB',\
    #    grid_data = True, zoom = zoom, vmax = None, vmin = None, save = False)

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))

    ax1.set_extent([ 155,200, 64, 80], datacrs)
    ax2.set_extent([ 155,200, 64, 80], datacrs)
    ax3.set_extent([ 155,200, 64, 80], datacrs)
    ax4.set_extent([ 155,200, 64, 80], datacrs)
    ax5.set_extent([ 155,200, 64, 80], datacrs)
    ax6.set_extent([ 155,200, 64, 80], datacrs)

    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax3, 'c)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax4, 'd)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax5, 'e)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax6, 'f)', fontsize = 11, backgroundcolor = 'white')

    fig1.tight_layout()

    if(save):
        outname = 'omi_ceres_modis_nsidc_compare_' + omi_date + '_v2.png'
        fig1.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_compare_OMI_MODIS_NSIDC_v2(date_str, ch1, \
        omi_dtype = 'shawn', minlat = 65., zoom = False, save = False):

    if(home_dir + '/Research/OMI/' not in sys.path):
        sys.path.append(home_dir + '/Research/OMI/')
    if(home_dir + '/Research/MODIS/obs_smoke_forcing/' not in sys.path):
        sys.path.append(home_dir + '/Research/MODIS/obs_smoke_forcing/')
    if(home_dir + '/Research/NSIDC/' not in sys.path):
        sys.path.append(home_dir + '/Research/NSIDC/')

    from OMILib import plotOMI_single_swath_figure
    from MODISLib import plot_MODIS_channel

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
        #size = (14, 5)
        size = (9, 4)
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0)
    elif(date_str[:8] == '20190811'\
        ):
        size = (12, 6.5)
    elif(date_str[:7] == '2006072'\
        ):
        size = (12, 8)
    elif(date_str[:8] == '20140811'\
        ):
        size = (9, 4)
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 135)
        #mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = ccrs.LambertConformal(central_longitude = 135, central_latitude = 75)
    else:
        size = (10, 9)
        mapcrs = mapcrs

    # Plot the MODIS true-color and channel data
    # ------------------------------------------
    modis_date = file_date_dict[date_str]['MODIS'][0]
    nsidc_date = file_date_dict[date_str]['NSIDC'][:8]
    omi_date   = date_str

    # Read in data for value printing
    # -------------------------------
    MODIS_ch7 = read_MODIS_channel(modis_date, 7, swath = True)

    # Make the overall figure
    # -----------------------
    plt.close('all')
    #fig1 = plt.figure(figsize = file_date_dict[modis_date_str]['size'])
    fig1 = plt.figure(figsize = size)
    ax1 = fig1.add_subplot(1,2,1, projection = mapcrs)
    ax2 = fig1.add_subplot(1,2,2, projection = mapcrs)
    #ax3 = fig1.add_subplot(2,2,3, projection = mapcrs)
    #ax4 = fig1.add_subplot(2,2,4, projection = mapcrs)
    #ax5 = fig1.add_subplot(2,3,5, projection = mapcrs)
    #ax6 = fig1.add_subplot(2,3,6, projection = mapcrs)
   
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = zoom, ax = ax1)

    # Plot the OMI data
    # -----------------
    plotOMI_single_swath_figure(omi_date, \
            dtype = omi_dtype, only_sea_ice = False, minlat = minlat, \
            ax = ax2, skiprows = [52], lat_circles = None, save = False, \
            zoom = zoom, shawn_path = shawn_path, colorbar = True)
    
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    if(date_str[:8] == '20140811'):
        ax1.set_title('Aqua MODIS True Color')
        ax2.set_title('OMI UVAI Aerosol Index')
        ax1.set_extent([ 125,165, 70, 85], datacrs)
        ax2.set_extent([ 125,165, 70, 85], datacrs)
        plt.suptitle(dt_date_str.strftime('Smoke over Ocean\n%H:%M UTC %d %B %Y'), weight = 'bold')
    elif(date_str[:8] == '20180705'):
        ax1.set_title('Aqua MODIS True Color')
        ax2.set_title('OMI UVAI Aerosol Index')
        ax1.set_extent([ 145,200, 64, 80], datacrs)
        ax2.set_extent([ 145,200, 64, 80], datacrs)
        plt.suptitle(dt_date_str.strftime('Smoke over Ice\n%H:%M UTC %d %B %Y'), weight = 'bold')
    else:
        plt.suptitle(dt_date_str.strftime('%H:%M UTC %d %B %Y'))

    
    fig1.tight_layout()

    if(save):
        outname = 'omi_modis_nsidc_compare_v2_' + omi_date + '.png'
        fig1.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_compare_OMI_MODIS_v2(date_str, ch1, \
        omi_dtype = 'shawn', minlat = 65., zoom = False, save = False):

    if(home_dir + '/Research/OMI/' not in sys.path):
        sys.path.append(home_dir + '/Research/OMI/')
    if(home_dir + '/Research/MODIS/obs_smoke_forcing/' not in sys.path):
        sys.path.append(home_dir + '/Research/MODIS/obs_smoke_forcing/')

    from OMILib import plotOMI_single_swath_figure
    from MODISLib import plot_MODIS_channel

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
        #size = (14, 5)
        size = (9, 4)
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0)
    elif(date_str[:8] == '20190811'\
        ):
        size = (12, 6.5)
    elif(date_str[:7] == '2006072'\
        ):
        size = (12, 8)
    elif(date_str[:8] == '20140811'\
        ):
        size = (9, 4)
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 135)
        #mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = ccrs.LambertConformal(central_longitude = 135, central_latitude = 75)
    else:
        size = (10, 9)
        mapcrs = mapcrs

    # Plot the MODIS true-color and channel data
    # ------------------------------------------
    modis_date = file_date_dict[date_str]['MODIS'][0]
    omi_date   = date_str

    # Read in data for value printing
    # -------------------------------
    MODIS_ch7 = read_MODIS_channel(modis_date, 7, swath = True, \
        include_cloud_mask = True)

    # Make the overall figure
    # -----------------------
    plt.close('all')
    #fig1 = plt.figure(figsize = file_date_dict[modis_date_str]['size'])
    fig1 = plt.figure(figsize = (7, 5))
    ax1 = fig1.add_subplot(2,2,1, projection = mapcrs)
    ax2 = fig1.add_subplot(2,2,2, projection = mapcrs)
    ax3 = fig1.add_subplot(2,2,3, projection = mapcrs)
    ax4 = fig1.add_subplot(2,2,4, projection = mapcrs)

    # Plot the OMI data
    # -----------------
    plotOMI_single_swath_figure(omi_date, \
            dtype = omi_dtype, only_sea_ice = False, minlat = minlat, \
            ax = ax1, skiprows = [52], lat_circles = None, save = False, \
            zoom = zoom, shawn_path = shawn_path, colorbar = True)
    ax1.set_title('OMI UVAI')
   
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = zoom, ax = ax2)
    ax2.set_title('Aqua MODIS True Color')

    plot_MODIS_channel(modis_date, ch1, swath = True, \
        zoom = zoom, ax = ax3, vmax = 0.4)
    ax3.set_title('Aqua MODIS Ch7\n2.105 μm - 2.155 μm')

    # Plot the CH7 data with cloud mask overlaid
    # ------------------------------------------
    plot_MODIS_channel(modis_date, ch1, swath = True, \
        zoom = zoom, ax = ax4, vmax = 0.4, plot_cbar = False)
    col_dict = {0: 'tab:blue',1: 'tab:orange',2: 'tab:green',3: 'tab:red'}
    labels = np.array(['Cloudy','Prob.\nCloudy','Prob.\nClear','Clear'])
    ccmm = mcolors.ListedColormap([col_dict[x] for x in col_dict.keys()])
    bounds = np.array([-0.5,0.5,1.5,2.5,3.5])   
    norm = mcolors.BoundaryNorm(bounds, ccmm.N)
    cld_mask = ax4.pcolormesh(MODIS_ch7['lon'], MODIS_ch7['lat'], \
        MODIS_ch7['cloud_mask'], cmap = ccmm, norm = norm, \
        transform = ccrs.PlateCarree(),
        shading = 'auto', alpha = 0.3) 
    cbar = plt.colorbar(ScalarMappable(norm = norm, cmap = ccmm), \
        ticks = [0,1,2,3], ax = ax4, orientation='vertical', pad = 0.03, \
        fraction = 0.052)
    cbar.ax.set_yticklabels(labels)
    #plot_MODIS_channel(modis_date, 'cloud_mask', swath = True, \
    #    zoom = zoom, ax = ax2, vmax = None)
    ax4.set_title('Aqua MODIS Ch7\nCloud Mask Overlay')

    x1 = 162.8
    y1 = 65.6
    x2 = 172.7
    y2 = 68.3
    dx = x2 - x1
    dy = y2 - y1
    ax4.arrow(x1, y1, dx, dy, facecolor = 'red', width = 0.2, head_width = 0.2, edgecolor = 'red', transform = ccrs.PlateCarree())
    x1 = 155.2
    y1 = 64.8
    plot_figure_text(ax4, 'Misclassified\nSmoke', xval = x1, yval = y1, \
        transform = ccrs.PlateCarree(), \
        color = 'black', fontsize = 8, backgroundcolor = None,\
        halign = 'center', location = 'lower_right', weight = None)

    x1 = -159.9
    y1 = 79.3
    x2 = -173.5
    y2 = 71.6
    dx = x2 - x1
    dy = y2 - y1
    #ax4.annotate('', xy = (x2, y2), xytext = (x1, y1), arrowprops = dict(headwidth = 1.0, headlength = 1.0, width = 0.2), transform = ccrs.PlateCarree())
    ax4.arrow(x1, y1, dx, dy, facecolor = 'red', width = 0.2, head_width = 0.2, edgecolor = 'red', transform = ccrs.PlateCarree())
    #plot_figure_text(ax4, 'Misclass.\nClear/Ice', xval = 148.1, yval = 67.6, \
    plot_figure_text(ax4, 'Misclassified\nSmoke', xval = x1, yval = y1, \
        transform = ccrs.PlateCarree(), \
        color = 'black', fontsize = 8, backgroundcolor = 'white',\
        halign = 'center', location = 'lower_right', weight = None)
    
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    if(date_str[:8] == '20140811'):
        #ax1.set_title('OMI UVAI')
        #ax2.set_title('Aqua MODIS True Color')
        #ax3.set_title('Aqua MODIS 2.1 μm Refl.')
        #ax4.set_title('Aqua MODIS Cloud Mask')
        ax1.set_extent([ 125,165, 70, 85], datacrs)
        ax2.set_extent([ 125,165, 70, 85], datacrs)
        ax3.set_extent([ 125,165, 70, 85], datacrs)
        ax4.set_extent([ 125,165, 70, 85], datacrs)
        plt.suptitle(dt_date_str.strftime('%H:%M UTC %d %B %Y'))
    elif(date_str[:8] == '20180705'):
        #ax1.set_title('OMI UVAI')
        #ax2.set_title('Aqua MODIS True Color')
        #ax3.set_title('Aqua MODIS 2.1 μm Refl.')
        #ax4.set_title('Aqua MODIS Cloud Mask')
        ax1.set_extent([ 145,200, 64, 80], datacrs)
        ax2.set_extent([ 145,200, 64, 80], datacrs)
        ax3.set_extent([ 145,200, 64, 80], datacrs)
        ax4.set_extent([ 145,200, 64, 80], datacrs)
        plt.suptitle(dt_date_str.strftime('%H:%M UTC %d %B %Y'))
    else:
        plt.suptitle(dt_date_str.strftime('%H:%M UTC %d %B %Y'))

    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax3, 'c)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax4, 'd)', fontsize = 11, backgroundcolor = 'white')
    
    fig1.tight_layout()

    if(save):
        outname = 'omi_modis_compare_v2_' + omi_date + '.png'
        fig1.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


def plot_compare_OMI_MODIS_NSIDC_v2_combined(date_str1, date_str2, ch1, \
        omi_dtype = 'shawn', minlat = 65., zoom = False, save = False):

    if(home_dir + '/Research/OMI/' not in sys.path):
        sys.path.append(home_dir + '/Research/OMI/')
    if(home_dir + '/Research/MODIS/obs_smoke_forcing/' not in sys.path):
        sys.path.append(home_dir + '/Research/MODIS/obs_smoke_forcing/')
    if(home_dir + '/Research/NSIDC/' not in sys.path):
        sys.path.append(home_dir + '/Research/NSIDC/')

    from OMILib import plotOMI_single_swath_figure
    from MODISLib import plot_MODIS_channel

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

    if((date_str1[:8] == '20080422') | \
       (date_str1[:8] == '20180705')\
        ):
        #size = (14, 5)
        mapcrs1 = ccrs.NorthPolarStereo(central_longitude = 0)
    #elif(date_str[:8] == '20190811'\
    #    ):
    #    size = (12, 6.5)
    #elif(date_str[:7] == '2006072'\
    #    ):
    #    size = (12, 8)
    elif(date_str1[:8] == '20140811'\
        ):
        mapcrs1 = ccrs.NorthPolarStereo(central_longitude = 135)
        #mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = ccrs.LambertConformal(central_longitude = 135, central_latitude = 75)
    else:
        mapcrs1 = mapcrs

    if((date_str2[:8] == '20080422') | \
       (date_str2[:8] == '20180705')\
        ):
        #size = (14, 5)
        mapcrs2 = ccrs.NorthPolarStereo(central_longitude = 0)
    #elif(date_str[:8] == '20190811'\
    #    ):
    #    size = (12, 6.5)
    #elif(date_str[:7] == '2006072'\
    #    ):
    #    size = (12, 8)
    elif(date_str2[:8] == '20140811'\
        ):
        mapcrs2 = ccrs.NorthPolarStereo(central_longitude = 135)
        #mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = ccrs.LambertConformal(central_longitude = 135, central_latitude = 75)
    else:
        mapcrs2 = mapcrs


    # Plot the MODIS true-color and channel data
    # ------------------------------------------
    modis_date1 = file_date_dict[date_str1]['MODIS'][0]
    nsidc_date1 = file_date_dict[date_str1]['NSIDC'][:8]
    ceres_date1 = file_date_dict[date_str1]['CERES']
    omi_date1   = date_str1
    modis_date2 = file_date_dict[date_str2]['MODIS'][0]
    nsidc_date2 = file_date_dict[date_str2]['NSIDC'][:8]
    ceres_date2 = file_date_dict[date_str2]['CERES']
    omi_date2   = date_str2

    # Read in data for value printing
    # -------------------------------
    MODIS_ch7_1 = read_MODIS_channel(modis_date1, 7, swath = True)
    MODIS_ch7_2 = read_MODIS_channel(modis_date2, 7, swath = True)

    # Make the overall figure
    # -----------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (11, 5.5))
    ax1 = fig1.add_subplot(2,3,1, projection = mapcrs1)
    ax2 = fig1.add_subplot(2,3,2, projection = mapcrs1)
    ax3 = fig1.add_subplot(2,3,4, projection = mapcrs2)
    ax4 = fig1.add_subplot(2,3,5, projection = mapcrs2)
    ax5 = fig1.add_subplot(2,3,3, projection = mapcrs1)
    ax6 = fig1.add_subplot(2,3,6, projection = mapcrs2)
   
    plot_MODIS_channel(modis_date1, 'true_color', swath = True, \
        zoom = zoom, ax = ax1)
    plot_MODIS_channel(modis_date2, 'true_color', swath = True, \
        zoom = zoom, ax = ax3, ptitle = None)

    # Plot the OMI data
    # -----------------
    plotOMI_single_swath_figure(omi_date1, \
            dtype = omi_dtype, only_sea_ice = False, minlat = minlat, \
            ax = ax2, skiprows = [52], lat_circles = None, save = False, 
            title = None, \
            zoom = zoom, shawn_path = shawn_path, colorbar = True)
    plotOMI_single_swath_figure(omi_date2, \
            dtype = omi_dtype, only_sea_ice = False, minlat = minlat, \
            ax = ax4, skiprows = [52], lat_circles = None, save = False, \
            title = None, \
            zoom = zoom, shawn_path = shawn_path, colorbar = True)

    # Plot the CERES data
    # -------------------
    plotCERES_hrly_figure(ceres_date1, 'SWF',  \
        minlat = minlat, lat_circles = None, ax = ax5, title = 'SWF',\
        grid_data = True, zoom = zoom, vmax = 350, vmin = None, save = False)
    plotCERES_hrly_figure(ceres_date2, 'SWF',  \
        minlat = minlat, lat_circles = None, ax = ax6, title = None,\
        grid_data = True, zoom = zoom, vmax = 450, vmin = None, save = False)
    
    dt_date_str1 = datetime.strptime(date_str1, '%Y%m%d%H%M')
    dt_date_str2 = datetime.strptime(date_str2, '%Y%m%d%H%M')
    ax1.set_title('Aqua MODIS True Color')
    ax2.set_title('OMI UV Aerosol Index')
    ax5.set_title('Aqua CERES TOA SWF')
    row_label_size = 12
    fig1.text(0.025, 0.65, dt_date_str1.strftime('%Y-%m-%d %H:%M'), ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 0)
    fig1.text(0.025, 0.22, dt_date_str2.strftime('%Y-%m-%d %H:%M'), ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 0)
    #ax3.set_ylabel('SMOKE OVER ICE', weight = 'bold')
    #plt.suptitle(dt_date_str.strftime('Smoke over Ocean\n%H:%M UTC %d %B %Y'), weight = 'bold')
    #ax3.set_title('Aqua MODIS True Color')
    #ax4.set_title('OMI UVAI Aerosol Index')
    #plt.suptitle(dt_date_str.strftime('Smoke over Ice\n%H:%M UTC %d %B %Y'), weight = 'bold')

    if(date_str1[:8] == '20140811'):
        ax1.set_extent([ 117,170, 70, 85], datacrs)
        ax2.set_extent([ 117,170, 70, 85], datacrs)
        ax5.set_extent([ 117,170, 70, 85], datacrs)
    elif(date_str1[:8] == '20180705'):
        #ax1.set_extent([ 145,200, 64, 80], datacrs)
        #ax2.set_extent([ 145,200, 64, 80], datacrs)
        ax1.set_extent([ 155,200, 64, 80], datacrs)
        ax2.set_extent([ 155,200, 64, 80], datacrs)
        ax5.set_extent([ 155,200, 64, 80], datacrs)

    if(date_str2[:8] == '20140811'):
        ax3.set_extent([ 117,170, 70, 85], datacrs)
        ax4.set_extent([ 117,170, 70, 85], datacrs)
        ax6.set_extent([ 117,170, 70, 85], datacrs)
    elif(date_str2[:8] == '20180705'):
        ax3.set_extent([ 155,200, 64, 80], datacrs)
        ax4.set_extent([ 155,200, 64, 80], datacrs)
        ax6.set_extent([ 155,200, 64, 80], datacrs)

    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax3, 'd)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax4, 'e)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax5, 'c)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax6, 'f)', fontsize = 11, backgroundcolor = 'white')
  
    #title_str =  'Arctic Biomass Burning Smoke\nOver Ocean: ' + \
    #    dt_date_str1.strftime('%Y-%m-%d %H:%M') + \
    #    '\nOver Ice:      ' + dt_date_str2.strftime('%Y-%m-%d %H:%M')
    title_str =  'Arctic Biomass Burning Smoke\nTop: Over Ocean\n' + \
        'Bottom: Over Ice'
    plt.suptitle(title_str)
 
    fig1.tight_layout()

    if(save):
        #outname = 'omi_modis_nsidc_compare_v2_combined_' + omi_date1 + \
        #outname = 'omi_modis_ceres_compare_v2_combined_' + omi_date1 + \
        #    '_' + omi_date2 + '.png'
        outname = 'omi_modis_ceres_compare_v3_combined_' + omi_date1 + \
            '_' + omi_date2 + '.png'
        fig1.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


def plot_compare_OMI_MODIS_CERES_v3_combined(calc_data1, calc_data2, \
        omi_dtype = 'shawn', minlat = 65., zoom = False, save = False):
   
    in_calc1 = h5py.File(calc_data1)
    mask_AI  = np.ma.masked_invalid(in_calc1['omi_uvai_pert'])

    # Figure out where the pixels containing AI are located
    # -----------------------------------------------------
    high_AI_idxs = np.where(mask_AI > 1.5)

    # Switch the longitudes to 0 - 360 rather than -180 - 180
    # -------------------------------------------------------
    work_lons = np.copy(in_calc1['omi_lon'][:,:])
    work_lons[work_lons < 0] = 360 + work_lons[work_lons < 0]
    

    # Figure out the minimum/maximum latitude/longitudes in the AI ranges
    # -------------------------------------------------------------------
    min_lat1 = np.min(in_calc1['omi_lat'][:,:][high_AI_idxs]) - 0
    max_lat1 = np.max(in_calc1['omi_lat'][:,:][high_AI_idxs]) + 2
    
    if(max_lat1 > 90):
        max_lat1 = 90

    # Calculate the 5th and 95th percentiles of longitude to prevent
    # outliers from ruining the mean or max/min values
    # --------------------------------------------------------------
    min_lon1 = np.percentile(work_lons[high_AI_idxs], 5) - 2
    max_lon1 = np.percentile(work_lons[high_AI_idxs], 95) + 2

    center_lon1 = (min_lon1 + max_lon1) / 2.
    workcrs1 = ccrs.NorthPolarStereo(central_longitude = center_lon1)

    # Now set up the crs for plot 2. 
    # ------------------------------
    workcrs2 = ccrs.NorthPolarStereo(central_longitude = 0)

    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure(figsize = (12.8, 5.5))
    ax1 = fig.add_subplot(2,4,1, projection = workcrs1)
    ax2 = fig.add_subplot(2,4,2, projection = workcrs1)
    ax3 = fig.add_subplot(2,4,3, projection = workcrs1)
    ax4 = fig.add_subplot(2,4,4, projection = workcrs1)
    ax5 = fig.add_subplot(2,4,5, projection = workcrs2)
    ax6 = fig.add_subplot(2,4,6, projection = workcrs2)
    ax7 = fig.add_subplot(2,4,7, projection = workcrs2)
    ax8 = fig.add_subplot(2,4,8, projection = workcrs2)

    axs1 = [ax1, ax2, ax3, ax4]
    axs2 = [ax5, ax6, ax7, ax8]

    plot_compare_OMI_MODIS_CERES(calc_data1, axs1, vmin = 75, \
        vmax = 250, auto_zoom = True, zoom = False, plot_titles = True, \
        label_xloc = 115.7, label_yloc = 78.6, labels = ['a)','b)','c)','d)'])
    plot_compare_OMI_MODIS_CERES(calc_data2, axs2, vmin = 75, \
        vmax = 450, auto_zoom = False, zoom = True, \
        label_xloc = -166, label_yloc = 65.4, labels = ['e)','f)','g)','h)'])

    #mask_AI  = np.ma.masked_invalid(in_calc1['omi_uvai_pert'])
    #hasher = np.ma.masked_invalid(mask_AI)
    #hasher = np.ma.masked_where( (hasher < 1.5), \
    #    hasher)

    #print("MIN MASK AI:", np.min(hasher), 'MAX MASK AI:', np.nanmax(hasher))


    #mask_AI = np.ma.masked_where(mask_AI < 1.5, mask_AI)
    #ax3.pcolormesh(in_calc1['omi_lon'][:,:], in_calc1['omi_lat'][:,:], mask_AI, \
    #    transform = ccrs.PlateCarree(), shading = 'auto', vmin = -1, vmax = 5, cmap = 'jet', alpha = 0.5)



    #ax3.pcolor(in_calc1['omi_lon'][:,:], in_calc1['omi_lat'][:,:],\
    #    hasher, hatch = '...', alpha=0.0, transform = datacrs, color = 'r', shading = 'auto')
    #ax1.pcolor(in_calc1['omi_lon'][:,:], in_calc1['omi_lat'][:,:], \
    #    hasher, alpha=0., shading = 'auto', \
    #    transform = ccrs.PlateCarree())

    dt_date_str1 = datetime.strptime(calc_data1.strip().split('/')[-1].split('_')[-1][:12], '%Y%m%d%H%M')
    dt_date_str2 = datetime.strptime(calc_data2.strip().split('/')[-1].split('_')[-1][:12], '%Y%m%d%H%M')

    plt.suptitle('Top: Biomass Burning Smoke over Ocean\nBottom: Biomass Burning Smoke over Ice and Land')

    fig.tight_layout(rect = [0.04,0.00,1.00,1.00])
    row_label_size = 10
    plotloc = axs1[0].get_position() 
    fig.text(plotloc.xmin - 0.01, (plotloc.ymax + plotloc.ymin) / 2., \
        dt_date_str1.strftime('%H:%M UTC %d-%b-%Y'), ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    plotloc = axs2[0].get_position() 
    fig.text(plotloc.xmin - 0.01, (plotloc.ymax + plotloc.ymin) / 2., \
        dt_date_str2.strftime('%H:%M UTC %d-%b-%Y'), ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)


    in_calc1.close()
   
    if(save):
        outname = 'omi_modis_ceres_compare_v4_combined_' + dt_date_str1.strftime('%Y%m%d%H%M') + \
            '_' + dt_date_str2.strftime('%Y%m%d%H%M') + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else: 
        plt.show()
 
def plot_compare_OMI_MODIS_CERES(calc_data, axs = None, \
        label_xloc = None, label_yloc = None, \
        labels = None, 
        omi_dtype = 'shawn', minlat = 65., auto_zoom = True, \
        vmin = None, vmax = None, \
        plot_titles = False, \
        zoom = False, save = False):

    # Load in the JSON file string relations
    # --------------------------------------
    with open(json_time_database, 'r') as fin:
        file_date_dict = json.load(fin)

    file_name = calc_data.strip().split('/')[-1] 
    date_str = file_name.split('.')[0].split('_')[-1]
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    
    sim_name = file_name.split('_')[3]
    print(sim_name, date_str)
    
    in_calc = h5py.File(calc_data)
    
    mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                   (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
    mask_lwf  = np.ma.masked_where((in_calc['ceres_lwf'][:,:] == -999.) | \
                                   (in_calc['ceres_lwf'][:,:] > 3000), in_calc['ceres_lwf'][:,:])
    mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
    mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)
    mask_AI   = np.ma.masked_invalid(in_calc['omi_uvai_pert'])
    
    # Plot the MODIS true-color and channel data
    # ------------------------------------------
    modis_date = file_date_dict[date_str]['MODIS'][0]
    nsidc_date = file_date_dict[date_str]['NSIDC'][:8]
    ceres_date = file_date_dict[date_str]['CERES']
    omi_date   = date_str

   
    if(auto_zoom):

        # Figure out where the pixels containing AI are located
        # -----------------------------------------------------
        high_AI_idxs = np.where(mask_AI > 1.5)

        # Switch the longitudes to 0 - 360 rather than -180 - 180
        # -------------------------------------------------------
        work_lons = np.copy(in_calc['omi_lon'][:,:])
        work_lons[work_lons < 0] = 360 + work_lons[work_lons < 0]
        

        # Figure out the minimum/maximum latitude/longitudes in the AI ranges
        # -------------------------------------------------------------------
        min_lat = np.min(in_calc['omi_lat'][:,:][high_AI_idxs]) - 0
        max_lat = np.max(in_calc['omi_lat'][:,:][high_AI_idxs]) + 2
        
        if(max_lat > 90):
            max_lat = 90

        # Calculate the 5th and 95th percentiles of longitude to prevent
        # outliers from ruining the mean or max/min values
        # --------------------------------------------------------------
        min_lon = np.percentile(work_lons[high_AI_idxs], 5) - 10
        max_lon = np.percentile(work_lons[high_AI_idxs], 95) + 2

        center_lon = (min_lon + max_lon) / 2.
        workcrs = ccrs.NorthPolarStereo(central_longitude = center_lon)

        #fig = plt.figure(figsize = (10.5, 6))
        figsize = (10, 4)
        #fig = plt.figure(figsize = (10, 4))
        #ax1 = fig.add_subplot(1,3,1, projection = workcrs) # OMI AI
        #ax2 = fig.add_subplot(1,3,2, projection = workcrs) # True color
        #ax3 = fig.add_subplot(1,3,3, projection = workcrs) # True color
        
    else:
        if((date_str[:8] == '20080422') | \
           (date_str[:8] == '20180705')\
            ):
            workcrs = ccrs.NorthPolarStereo(central_longitude = 0)
            min_lat = 64
            max_lat = 80
            min_lon = 155
            max_lon = 195

        elif(date_str[:8] == '20140811'\
            ):
            workcrs = ccrs.NorthPolarStereo(central_longitude = 135)
            min_lat = 70
            max_lat = 85
            min_lon = 117
            max_lon = 170
        else:
            #workcrs = mapcrs
            workcrs = ccrs.NorthPolarStereo(central_longitude = 0)

        figsize = (9, 4)


    in_ax = True
    if(axs is None):
        in_ax = False
        fig = plt.figure(figsize = figsize)
        ax1 = fig.add_subplot(1,4,1, projection = workcrs)
        ax2 = fig.add_subplot(1,4,2, projection = workcrs)
        ax3 = fig.add_subplot(1,4,3, projection = workcrs)
        ax4 = fig.add_subplot(1,4,4, projection = workcrs)
    else:
        ax1 = axs[0]
        ax2 = axs[1]    
        ax3 = axs[2]
        ax4 = axs[3]
 
    CERES_data_hrly = readgridCERES_hrly_grid(ceres_date, 'swf', \
        satellite = 'Aqua', minlat = 65.5)
    
    # Plot the MODIS data
    # ------------------- 
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = True, ax = ax1)
    #mesh = ax2.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_calc, \
    #    transform = ccrs.PlateCarree(), shading = 'auto')
    #cbar = fig.colorbar(mesh, ax = ax2, label = 'SWF [Wm$^{-2}$]')
    #ax2.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    ax1.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    #if(auto_zoom):
    #    ax2.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    #else:
    #    ax2.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
    #    ax2.set_boundary(circle, transform=ax2.transAxes)
    if(plot_titles):
        ax1.set_title('Aqua MODIS')
    else:
        ax1.set_title(None)
    ax1.coastlines()

    # Plot the OMI data
    # ----------------- 
    mesh = ax2.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = -2, vmax = 5, cmap = 'jet')
    cbar = plt.colorbar(mesh, ax = ax2, fraction = 0.046, label = 'OMI UVAI')
    #ax1.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    #if((auto_zoom) | (zoom)):
    #if(auto_zoom):
    ax2.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    #else:
    #    ax1.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
    #    ax1.set_boundary(circle, transform=ax1.transAxes)
    if(plot_titles):
        ax2.set_title('OMI UVAI')
    ax2.coastlines()
   
  
    # Plot the CERES data
    # -------------------
    plot_lon  = CERES_data_hrly['lon']
    plot_lat  = CERES_data_hrly['lat']
    mask_flux = CERES_data_hrly['swf']

    mesh = ax3.pcolormesh(plot_lon, plot_lat,mask_flux,transform = datacrs,\
        cmap = 'viridis', vmin = vmin, vmax = vmax, shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax3, fraction = 0.046, label = 'SWF [Wm$^{-2}$]')
    ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    #if(auto_zoom):
    #    ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    #else:
    #    ax3.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
    #    ax3.set_boundary(circle, transform=ax3.transAxes)
    if(plot_titles):
        ax3.set_title('Aqua CERES TOA SWF')
    ax3.coastlines()

    mesh = ax4.pcolormesh(plot_lon, plot_lat,mask_flux,transform = datacrs,\
        cmap = 'viridis', vmin = vmin, vmax = vmax, shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax4, fraction = 0.046, label = 'SWF [Wm$^{-2}$]')
    ax4.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    #if(auto_zoom):
    #    ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    #else:
    #    ax3.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
    #    ax3.set_boundary(circle, transform=ax3.transAxes)
    if(plot_titles):
        ax4.set_title('Aqua CERES TOA SWF\nOverlay: UVAI > 1.0')
    ax4.coastlines()


    mask_AI   = np.ma.masked_invalid(in_calc['omi_uvai_pert'])
    mask_AI = np.ma.masked_where(mask_AI < 1.0, mask_AI)
    ax4.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = -2, vmax = 5, cmap = 'jet', alpha = 0.5)


    ##!#plotCERES_hrly_figure(ceres_date, 'SWF',  \
    ##!#    minlat = minlat, lat_circles = None, ax = ax3, title = 'SWF',\
    ##!#    grid_data = True, zoom = zoom, vmax = 350, vmin = None, save = False)
    ##!#if(auto_zoom):
    ##!#    ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
 

    font_size = 10 
    if((label_xloc is not None) & (label_yloc is not None)):
        if(labels is not None):
            label_text1 = labels[0]
            label_text2 = labels[1]
            label_text3 = labels[2]
            label_text4 = labels[3]
        else:
            label_text1 = 'a)'
            label_text2 = 'b)'
            label_text3 = 'c)'
            label_text4 = 'd)'
        plot_figure_text(ax1, label_text1, \
            xval = label_xloc, yval = label_yloc, transform = ccrs.PlateCarree(), \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax2, label_text2, \
            xval = label_xloc, yval = label_yloc, transform = ccrs.PlateCarree(), \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax3, label_text3, \
            xval = label_xloc, yval = label_yloc, transform = ccrs.PlateCarree(), \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax4, label_text4, \
            xval = label_xloc, yval = label_yloc, transform = ccrs.PlateCarree(), \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
    else:
        plot_figure_text(ax1, 'a)', \
            xval = None, yval = None, transform = None, \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax2, 'b)', \
            xval = None, yval = None, transform = None, \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax3, 'c)', \
            xval = None, yval = None, transform = None, \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax4, 'd)', \
            xval = None, yval = None, transform = None, \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
    
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    
    #in_base.close()
    in_calc.close()
   
    if(not in_ax):  
        fig.tight_layout()

        if(save):    
            if(auto_zoom):
                zoom_add = '_zoom'
            else:
                zoom_add = ''
            outname = 'omi_modis_ceres_nn_compare_'+ date_str + '_' + zoom_add + '_v3.png'
            fig.savefig(outname, dpi = 200)
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
# mod_slopes: either "add" or "subtract", either adds or subtracts
 #              the calculated standard error of the SZA-mean slopes
def calculate_interp_forcings(coloc_dict, month_idx, minlat, maxlat, \
        dtype, cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
        mod_slopes = None):

    ocean_slopes = calc_slope_clear_clean_sfctype(coloc_dict, 0, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type, \
        mod_slopes = mod_slopes)
    ice_slopes   = calc_slope_clear_clean_sfctype(coloc_dict, 1, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type, \
        mod_slopes = mod_slopes)
    land_slopes  = calc_slope_clear_clean_sfctype(coloc_dict, 2, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type, \
        mod_slopes = mod_slopes)
    if(coloc_dict['raw_slopes'].shape[0] == 4):
        include_mix = True
        mix_slopes  = calc_slope_clear_clean_sfctype(coloc_dict, 3, cld_idx, \
            maxerr = maxerr, min_cloud = min_cloud, data_type = data_type, \
            mod_slopes = mod_slopes)
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
# 
# reference_ice: if set to a string (e.g. '2005'), the ice cocentration
# from that day in 2005 will be used to determine which forcing
# efficiency value is used. 
# reference_cld: same as reference_ice but for the MODIS daily cloud
# fraction.
# ---------------------------------------------------------------------
def calculate_type_forcing_v3(OMI_daily_data, OMI_monthly_data, coloc_dict, \
        date_str, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw',\
        reference_ice = None, reference_cld = None, mod_slopes = None, \
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
    if(not (date_str in file_strs)):
        # Return a nan array
        # ------------------
        print("WARNING: Date", date_str, "not in daily data. Returning nans")
        return np.full(OMI_daily_data['grid_AI'][10,:,:].shape, np.nan)
         
    match_idx = np.where(date_str == file_strs)[0][0]
    local_OMI_daily = np.ma.masked_where(\
        OMI_daily_data['count_AI'][match_idx,:,:] == 0, \
        OMI_daily_data['grid_AI'][match_idx,:,:])
    
    # tidx is the "month_idx"
    # -----------------------    
    tidx = int(date_str[4:6]) - 4

    if(reference_cld is None):
        cld_str = date_str
    else:
        cld_str = reference_cld + date_str[4:]

    MYD08_data = read_MODIS_MYD08_single(cld_str, minlat = minlat, \
        maxlat = maxlat)
   
    # Load in the single-day NSIDC ice concentration
    #
    # If the user wants to use a different year's ice values as a reference,
    # modify the date string here accordingly
    # ------------------------------------------------------------------------
    if(reference_ice is None):
        ice_str = date_str
    else:
        ice_str = reference_ice + date_str[4:]

    NSIDC_data =  readNSIDC_daily(ice_str, maxlat = maxlat)
    NSIDC_data = grid_data_conc(NSIDC_data, minlat = minlat, maxlat = maxlat)
    NSIDC_data['grid_ice_conc'] = np.ma.masked_where((NSIDC_data['grid_ice_conc'] < 0) | \
        (NSIDC_data['grid_ice_conc'] > 100), NSIDC_data['grid_ice_conc']).squeeze()

    clear_sky_AI = np.array([np.nanmean(\
        np.ma.masked_where(OMI_monthly_data['AI'][midx::6,:,:] > ai_thresh,\
        OMI_monthly_data['AI'][midx::6,:,:]), axis = 0) for midx in range(6)])
    clear_sky_AI = clear_sky_AI[tidx,:,:]

    #clear_sky_AI = np.array([np.nanmean(\
    #    np.ma.masked_where(OMI_data['AI'][tidx::6,:,:] > ai_thresh,\
    #    OMI_data['AI'][tidx::6,:,:]), axis = 0) for tidx in range(6)])
    
    clear_dict = calculate_interp_forcings(coloc_dict, tidx, \
        minlat, maxlat, 'clear', cld_idx = cld_idx, maxerr = maxerr, \
        min_cloud = min_cloud, data_type = data_type, mod_slopes = mod_slopes)
    cloud_dict = calculate_interp_forcings(coloc_dict, tidx, \
        minlat, maxlat, 'cloud', cld_idx = cld_idx, maxerr = maxerr, \
        min_cloud = min_cloud, data_type = data_type, mod_slopes = mod_slopes)

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
        return estimate_forcings, MYD08_data, NSIDC_data, clear_sky_AI
    else:
        return estimate_forcings

# This verison is set up to use daily data, but still needs the monthly
# OMI data to determine the clear-sky climatology.
# ---------------------------------------------------------------------
def calculate_type_forcing_v3_monthly(OMI_daily_data, OMI_monthly_data, \
        coloc_dict, month_idx, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw',\
        reference_ice = None, reference_cld = None, mod_slopes = None, \
        filter_bad_vals = True, return_modis_nsidc = True, debug = False):


    # Set up arrays to hold the monthly-averaged forcing calculations
    # ---------------------------------------------------------------
    xdim = OMI_daily_data['grid_AI'].shape[1]
    ydim = OMI_daily_data['grid_AI'].shape[2]
    daily_force_vals = np.full( (31, xdim, ydim), np.nan)

    if(str(month_idx) == 'all'):
        l_all_months = True
        month_force_vals = np.full(\
            (OMI_monthly_data['DATES'].shape[0], xdim, ydim), \
            np.nan)

        begin_date_str = datetime.strptime(\
            OMI_monthly_data['DATES'][0], '%Y%m')
        end_date_str = datetime.strptime(\
            OMI_monthly_data['DATES'][-1], '%Y%m') + \
            relativedelta(months = 1) - timedelta(days = 1)
    else:
        l_all_months = False
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
        estimate_forcing = \
            calculate_type_forcing_v3(OMI_daily_data, OMI_monthly_data, \
                coloc_dict, date_str, minlat = minlat, maxlat = maxlat, \
                ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr,\
                min_cloud = min_cloud, data_type = data_type,\
                filter_bad_vals = filter_bad_vals, \
                reference_ice = reference_ice, \
                reference_cld = reference_cld, \
                mod_slopes = mod_slopes, \
                return_modis_nsidc = False)

        # Insert the values into the daily holding array
        # ----------------------------------------------
        daily_force_vals[day_count,:,:] = estimate_forcing[:,:] 

        # Increment working date
        # ----------------------
        new_work_date = local_date_str + timedelta(days = 1)

        # If the new working date has a different month than
        # the previous working date, then average the daily values
        # and move the working date to the next desired month
        # --------------------------------------------------------
        if(new_work_date.month != local_date_str.month):
            month_force_vals[month_count,:,:] = np.nanmean(\
                daily_force_vals[:,:,:], axis = 0)
            day_count = 0
            month_count += 1            

            if(l_all_months):
                if(new_work_date.month == 10):
                    new_work_date = new_work_date + relativedelta(months = 6)
            else:
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
                if((trend_type == 'standard') | (trend_type == 'linregress')): 
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
        coloc_dict, minlat = 65., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
        save = False, filter_bad_vals = False):

    print("HERE1", date_str, minlat, maxlat, ai_thresh, cld_idx, maxerr, min_cloud, data_type)
    estimate_forcing, MYD08_data, NSIDC_data, clear_sky_AI = \
        calculate_type_forcing_v3(OMI_daily_data, OMI_month_data, \
            coloc_dict, date_str, minlat = minlat, maxlat = maxlat, \
            ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr,\
            min_cloud = min_cloud, data_type = data_type,\
            filter_bad_vals = filter_bad_vals, return_modis_nsidc = True)

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')
    
    file_strs = np.array([str(tval) for tval in OMI_daily_data['day_values']])
    match_idx = np.where(date_str == file_strs)[0][0]
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 5))
    ax1 = fig.add_subplot(2,3,1, projection = mapcrs)  # Map of original AI
    ax2 = fig.add_subplot(2,3,2, projection = mapcrs)  # Map of background clear AI
    ax3 = fig.add_subplot(2,3,3, projection = mapcrs)  # Map of AI above thresh
    ax4 = fig.add_subplot(2,3,4, projection = mapcrs)  # Map of forcing values
    ax5 = fig.add_subplot(2,3,5, projection = mapcrs)  # Map of cloud fraction
    ax6 = fig.add_subplot(2,3,6, projection = mapcrs)  # Map of ice concentration
    #ax6 = fig.add_subplot(2,3,2, projection = mapcrs)  # Map of GPQF values
    #ax3 = fig.add_subplot(2,3,6)  # histogram

    above_threshold = np.ma.masked_where(\
        OMI_daily_data['grid_AI'][match_idx,:,:] < ai_thresh, \
        OMI_daily_data['grid_AI'][match_idx,:,:])

    mesh = ax1.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], \
        OMI_daily_data['grid_AI'][match_idx,:,:], shading = 'auto', \
        transform = datacrs, cmap = 'jet', vmin = 0, vmax = 4.0)
    cbar = fig.colorbar(mesh, ax = ax1, label =  'UVAI')
    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.coastlines()
    ax1.set_title('Daily OMI UVAI')

    mesh = ax2.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], \
        #OMI_daily_data['grid_GPQF'][match_idx,:,:], shading = 'auto', \
        clear_sky_AI, shading = 'auto', \
        transform = datacrs, cmap = 'jet', vmin = None, vmax = None)
    cbar = fig.colorbar(mesh, ax = ax2, label =  'UVAI')
    ax2.set_extent([-180,180,minlat,90], datacrs)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.coastlines()
    ax2.set_title('Daily OMI UVAI\nClear-sky Background')

    mesh = ax3.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], \
        #OMI_daily_data['grid_GPQF'][match_idx,:,:], shading = 'auto', \
        above_threshold, shading = 'auto', \
        transform = datacrs, cmap = 'jet', vmin = 0, vmax = 4.0)
    cbar = fig.colorbar(mesh, ax = ax3, label =  'UVAI')
    ax3.set_extent([-180,180,minlat,90], datacrs)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.coastlines()
    ax3.set_title('Daily OMI UVAI\nUVAI > ' + str(ai_thresh))

    mesh = ax4.pcolormesh(MYD08_data['lon'], \
        MYD08_data['lat'], MYD08_data['cld_frac_mean'][:,:], shading = 'auto',\
        transform = datacrs, \
        cmap = 'viridis')
    cbar = fig.colorbar(mesh, ax = ax4, label = 'Cloud Fraction')
    ax4.set_extent([-180,180,minlat,90], datacrs)
    ax4.set_boundary(circle, transform=ax4.transAxes)
    ax4.coastlines()
    ax4.set_title('Aqua MODIS\nCloud Fraction')

    mask_ice = np.ma.masked_where(NSIDC_data['grid_ice_conc'][:,:] == 0, \
        NSIDC_data['grid_ice_conc'][:,:])
    mesh = ax5.pcolormesh(NSIDC_data['grid_lon'], \
        NSIDC_data['grid_lat'], mask_ice, shading = 'auto',\
        transform = datacrs, \
        cmap = 'ocean', vmin = 1, vmax = 100)
    cbar = fig.colorbar(mesh, ax = ax5, label = 'Ice Concentration [%]')
    ax5.set_extent([-180,180,minlat,90], datacrs)
    ax5.set_boundary(circle, transform=ax5.transAxes)
    ax5.coastlines()
    ax5.set_title('SSMI/S Sea\nIce Concentration')

    min_force = np.nanmin(estimate_forcing[:,:])
    max_force = np.nanmax(estimate_forcing[:,:])
    if(abs(max_force) > abs(min_force)):
        lims = [-abs(max_force), abs(max_force)]
    else:
        lims = [-abs(min_force), abs(min_force)]
    mesh = ax6.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], estimate_forcing[:,:], shading = 'auto',\
        transform = datacrs, \
        cmap = 'bwr', vmin = lims[0], vmax = lims[1])
    cbar = fig.colorbar(mesh, ax = ax6, label = 'Aerosol Forcing [W/m2]')
    ax6.set_extent([-180,180,minlat,90], datacrs)
    ax6.set_boundary(circle, transform=ax6.transAxes)
    ax6.coastlines()
    ax6.set_title('Estimated\nAerosol Forcing')

    ###ax3.hist(np.ma.masked_where(estimate_forcing[:,:] == 0, \
    ###    estimate_forcing[:,:]).compressed(), bins = 'auto')
    ####ax6.set_yscale('log')
    ###ax3.set_xlabel('Estimated Forcing [W/m2]')
    ###ax3.set_ylabel('Counts')

    plt.suptitle(dt_date_str.strftime('%Y-%m-%d'))

    plot_subplot_label(ax1, 'a)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax4, 'd)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax5, 'e)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax6, 'f)', fontsize = 10, backgroundcolor = None)

    fig.tight_layout()
    if(save):
        #outname = 'test_calc_forcing_v3_' + date_str + '.png'
        dtype_add = ''
        if('data_type' in coloc_dict.keys()):
            if(coloc_dict['data_type'] == 'omi_uvai_pert'):
                dtype_add = '_pert'
        minlat_add = ''
        if('minlat' in coloc_dict.keys()):
            minlat_add = '_minlat' + str(int(coloc_dict['minlat']))
                
        outname = 'test_calc_forcing_v3_' + date_str + dtype_add + \
            minlat_add + '_v2.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


# all_month_values: returned from calculate_type_forcing_v3_monthly with
# 'all' used in place of the month_idx argument
# ----------------------------------------------------------------------
def plot_test_forcing_v3_monthly(all_month_values, OMI_monthly_data, date_str, \
        minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
        save = False, filter_bad_vals = False):

    file_strs = np.array([str(tval) for tval in OMI_monthly_data['DATES']])
    match_idx = np.where(date_str == file_strs)[0][0]
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 4))
    ax1 = fig.add_subplot(1,2,1, projection = mapcrs)  # Map of monthy force est
    ax2 = fig.add_subplot(1,2,2)  # histogram

    max_data = np.max(all_month_values[match_idx,:,:])
    mesh = ax1.pcolormesh(OMI_monthly_data['LON'], \
        OMI_monthly_data['LAT'], \
        all_month_values[match_idx,:,:], shading = 'auto', \
        transform = datacrs, cmap = 'bwr', vmin = -max_data, vmax = max_data)
    cbar = fig.colorbar(mesh, ax = ax1, label =  'Monthly Forcing [W/m2]')
    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.coastlines()

    ax2.hist(np.ma.masked_where(all_month_values[match_idx,:,:] == 0, \
        all_month_values[match_idx,:,:]).compressed(), bins = 'auto')
    #ax6.set_yscale('log')
    ax2.set_xlabel('Estimated Forcing [W/m2]')
    ax2.set_ylabel('Counts')

    plt.suptitle(date_str)

    fig.tight_layout()
    if(save):
        outname = 'test_calc_forcing_v3_monthly_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


# all_month_values: returned from calculate_type_forcing_v3_monthly with
# 'all' used in place of the month_idx argument
# 
# This function plots the estimated forcing values for a single, 
# provided date and for its associated month
# ----------------------------------------------------------------------
def plot_test_forcing_v3_daily_monthly(date_str, all_month_values, \
        OMI_daily_data, OMI_monthly_data, coloc_dict,\
        minlat = 65., maxlat = 87., ai_thresh = 0.7, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
        save = False, filter_bad_vals = False):

    #estimate_forcing, MYD08_data, NSIDC_data = \
    # Grab the daily forcing values here
    # ----------------------------------
    print("HERE2", date_str, minlat, maxlat, ai_thresh, cld_idx, maxerr, min_cloud, data_type)
    estimate_forcing = \
        calculate_type_forcing_v3(OMI_daily_data, OMI_monthly_data, \
            coloc_dict, date_str, minlat = minlat, maxlat = maxlat, \
            ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr,\
            min_cloud = min_cloud, data_type = data_type,\
            filter_bad_vals = filter_bad_vals, return_modis_nsidc = False)

    file_strs = np.array([str(tval) for tval in OMI_daily_data['day_values']])
    match_idx = np.where(date_str == file_strs)[0][0]
    
    plt.close('all')
    fig = plt.figure(figsize = (10, 3))
    ax1 = fig.add_subplot(1,3,1, projection = mapcrs)  # Map of daily OMI AI
    ax2 = fig.add_subplot(1,3,2, projection = mapcrs)  # Map of daily forcing
    ax3 = fig.add_subplot(1,3,3, projection = mapcrs)  # map of monthly forcing estimate
   
    # Plot daily OMI AI
    # ----------------- 
    mesh = ax1.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], \
        OMI_daily_data['grid_AI'][match_idx,:,:], shading = 'auto', \
        transform = datacrs, cmap = 'jet', vmin = 0, vmax = 4.0)
    cbar = fig.colorbar(mesh, ax = ax1, label =  'Daily AI')
    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_title('Gridded Perturbed OMI AI\n' + date_str)
    ax1.coastlines()

    # Plot daily forcing estimate
    # ---------------------------
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
    cbar = fig.colorbar(mesh, ax = ax2, label = 'Est. Aerosol Forcing [W/m2]')
    ax2.set_extent([-180,180,minlat,90], datacrs)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.set_title('Estimated Daily Forcing\n' + date_str)
    ax2.coastlines()
    
    # Plot monthly forcing estimate
    # -----------------------------
    file_strs = np.array([str(tval) for tval in OMI_monthly_data['DATES']])
    match_idx = np.where(date_str[:6] == file_strs)[0][0]
    max_data = np.max(all_month_values[match_idx,:,:])
    mesh = ax3.pcolormesh(OMI_monthly_data['LON'], \
        OMI_monthly_data['LAT'], \
        all_month_values[match_idx,:,:], shading = 'auto', \
        transform = datacrs, cmap = 'bwr', vmin = -max_data, vmax = max_data)
    cbar = fig.colorbar(mesh, ax = ax3, label =  'Est. Aerosol Forcing [W/m2]')
    ax3.set_extent([-180,180,minlat,90], datacrs)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.set_title('Estimated Monthly Forcing\n' + date_str[:6])
    ax3.coastlines()

    plot_subplot_label(ax1, 'a)', fontsize = 14, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 14, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 14, backgroundcolor = None)

    fig.tight_layout()
    if(save):
        outname = 'test_calc_forcing_v3_daily_monthly_' + date_str + '.png'
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
    ax.set_ylabel('CERES SWF [W/m2]')

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
  
# Allows a user to compare raw colocated data for two sfc types 
def plot_dual_combined_multi_type(comp_grid_data, combined_data, \
        xval = 'omi_uvai_raw', \
        #cld_min = None, cld_max = None,\
        ch7_min1 = None, ch7_max1 = None,\
        ch7_min2 = None, ch7_max2 = None,\
        ice_min1 = None, ice_max1 = None,\
        ice_min2 = None, ice_max2 = None,\
        sza_min1 = None, sza_max1 = None,\
        sza_min2 = None, sza_max2 = None,\
        ai_min = 2,  ai_max = None,\
        save = False, show_trend = False, shade_density = False, \
        minlat = None, \
        trend_type = 'theil-sen'):

    if(ch7_min2 == None): ch7_min2 = ch7_min1
    if(ch7_max2 == None): ch7_max2 = ch7_max1
    if(sza_min2 == None): sza_min2 = sza_min1
    if(sza_max2 == None): sza_max2 = sza_max1

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
            xval = xval, \
            omi_min = ai_min,  omi_max  = ai_max, \
            sza_min = sza_min1, sza_max = sza_max1, \
            ice_min = ice_min1, ice_max = ice_max1, \
            ch7_min = ch7_min1, ch7_max = ch7_max1, \
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
        plot_combined_scatter(combined_data, ax = ax2, \
            xval = xval, \
            omi_min = ai_min,  omi_max  = ai_max, \
            sza_min = sza_min2, sza_max = sza_max2, \
            ice_min = ice_min2, ice_max = ice_max2, \
            ch7_min = ch7_min2, ch7_max = ch7_max2, \
            #cld_min = cld_min, cld_max = cld_max, \
            ai_diff = ai_diff, sza_diff = sza_diff, \
            ice_diff = ice_diff, ch7_diff = ch7_diff,\
            #cld_diff = cld_diff, \
            trend_type = trend_type, show_trend = show_trend, \
            shade_density = shade_density)

    except:
        print("ERROR in plot_combined_scatter for sfc type 2")
        return

    if(ice_min1 == 0):
        title1 = 'Ocean\nIce Conc. < 20%'
    elif(ice_min1 == 21):
        title1 = 'Mix\n21% < Ice Conc. < 79%'
    elif(ice_min1 == 80):
        title1 = 'Ice\nIce Conc. > 80%'
    elif(ice_min1 == 105):
        title1 = 'Land'

    if(ice_min2 == 0):
        title2 = 'Ocean\nIce Conc. <= 20%'
    elif(ice_min2 == 21):
        title2 = 'Mix\n20% < Ice Conc. < 80%'
    elif(ice_min2 == 80):
        title2 = 'Ice\nIce Conc. >= 80%'
    elif(ice_min2 == 105):
        title2 = 'Land'
    ax1.set_title(title1)
    ax2.set_title(title2)
    

    title = set_comp_title(\
                           #cld_min, cld_max, \
                           ch7_min1, ch7_max1, \
                           ice_min1, ice_max1, \
                           sza_min1, sza_max1,\
                           ai_min,  ai_max, version2 = True)
    filer = set_file_title(\
                           #cld_min, cld_max, \
                           ch7_min1, ch7_max1, \
                           ice_min1, ice_max1, \
                           sza_min1, sza_max1,\
                           ai_min,  ai_max, \
                           minlat, xval)
    plt.suptitle(title)

    fig.tight_layout()

    if(save):
        outname = 'comp_dual_rawL2_' + filer + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image",outname)
    else:
        plt.show()
   
def set_comp_title(\
        #cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None, \
        minlat = None, 
        version2 = False):

    ch7_min = np.round(ch7_min, 2)
    ch7_max = np.round(ch7_max, 2)

    if(version2):
        title = 'Raw Colocated L2\n'
    else:
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
            title = title + str(ch7_min) + ' < 2.1 μm Refl. < ' + str(ch7_max) + '\n'
        else:
            title = title + '2.1 μm Refl. > ' + str(ch7_min) + '\n'
    else:
        if(ch7_max is not None):
            title = title + '2.1 μm Refl. < ' + str(ch7_max) + '\n'

    if(not version2):
        if(ice_min is not None):
            if(ice_max is not None):
                title = title + str(ice_min) + ' < ice < ' + str(ice_max) + '\n'
            else:
                title = title + 'ICE % > ' + str(ice_min) + '\n'
        else:
            if(ice_max is not None):
                title = title + 'ICE % < ' + str(ice_max) + '\n'

    if(sza_min is not None):
        if(sza_max is not None):
            title = title + str(sza_min) + ' < SZA < ' + str(sza_max)
        else:
            title = title + 'SZA > ' + str(sza_min)
    else:
        if(sza_max is not None):
            title = title + 'SZA < ' + str(sza_max)

    if(not version2):
        if(ai_min is not None):
            if(ai_max is not None):
                title = '\n' + title + str(ai_min) + ' < AI < ' + str(ai_max)
            else:
                title = '\n' + title + 'AI > ' + str(ai_min)
        else:
            if(ai_max is not None):
                title = '\n' + title + 'AI < ' + str(ai_max)

    return title 

def set_file_title(\
        #cld_min = None, cld_max = None,\
        ch7_min = None, ch7_max = None,\
        ice_min = None, ice_max = None,\
        sza_min = None, sza_max = None,\
        ai_min = None,  ai_max = None, \
        minlat = None, xval = None):

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

    if(minlat is not None):
        title = title + '_minlat' + str(int(minlat))

    title = title + '_' + xval

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
        smoother = 'None', sizer = 1, \
        xval = 'omi_uvai_raw'):

    ice_mins = np.array([0,  80, 105, 20])
    ice_maxs = np.array([20,100,None, 80])
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
                        xval = xval, yval = 'ceres_swf', 
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
    out_dict['data_type']  = xval   
 
    return out_dict

# smoother: either 'none','smooth','smoother'
# sizer   : either 1, 2, or 3
# Made on 2024/02/29 to remove the CH7 bin stuff and 
def calc_raw_grid_slopes_v2(combined_data, comp_grid_data, \
        ai_min = None, ai_max = None, \
        trend_type = 'theil-sen', \
        smoother = 'None', sizer = 1, \
        xval = 'omi_uvai_raw'):

    ice_mins = np.array([0,  80, 105, 20])
    ice_maxs = np.array([20,100,None, 80])
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
                        xval = xval, yval = 'ceres_swf', 
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
    out_dict['data_type']  = xval   
 
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
        min_cloud = 0.95, dtype = 'raw', \
        plot_counts = False):

    plt.close('all')
    #fig = plt.figure(figsize = (10, 6))
    if(plot_counts):
        fig = plt.figure(figsize = (11, 9))
        # Trends
        axs = fig.subplots(nrows = 3, ncols = 4, sharex = True, sharey = True)
        ax1 = axs[0,0]  # ocean force
        ax7 = axs[0,1]  # mix force
        ax2 = axs[0,2]  # ice force
        ax3 = axs[0,3]  # land force
        ax4 = axs[1,0]  # ocean cloud
        ax8 = axs[1,1]  # mix cloud
        ax5 = axs[1,2]  # ice cloud
        ax6 = axs[1,3]  # land cloud
        ax9 = axs[2,0]  # ocean count
        ax10= axs[2,1]  # mix count
        ax11= axs[2,2]  # ice count
        ax12= axs[2,3]  # land count
    else:
        fig = plt.figure(figsize = (11, 6))
        # Trends
        axs = fig.subplots(nrows = 2, ncols = 4, sharex = True, sharey = True)
        ax1 = axs[0,0]  # ocean force
        ax7 = axs[0,1]  # mix force
        ax2 = axs[0,2]  # ice force
        ax3 = axs[0,3]  # land force
        ax4 = axs[1,0]  # ocean cloud
        ax8 = axs[1,1]  # mix cloud
        ax5 = axs[1,2]  # ice cloud
        ax6 = axs[1,3]  # land cloud
    #ax1 = fig.add_subplot(2,4,1)
    #ax7 = fig.add_subplot(2,4,2)
    #ax2 = fig.add_subplot(2,4,3)
    #ax3 = fig.add_subplot(2,4,4)
   
    ## Cloud fracs 
    #ax4 = fig.add_subplot(2,4,5)
    #ax8 = fig.add_subplot(2,4,6)
    #ax5 = fig.add_subplot(2,4,7)
    #ax6 = fig.add_subplot(2,4,8)
    
    cmap = 'bwr'
  
    xlabel = 'OMI Solar Zenith Angle [$^{o}$]' 
    if('cod_mins' in out_dict.keys()):
        cloud_var = 'cod'
        ylabel = 'MODIS COD'
    else:
        cloud_var = 'ch7'
        ylabel = 'MODIS 2.1 μm Refl.'
 
    shrk = 1.0

    # ice_idx: 0 = ocean, 1 = ice, 2 = land, 3 = mix
    ice_idx = 0

    def plot_force_efficiency(pax, out_dict, dtype, ice_idx, maxerr, vmin,\
            vmax, cmap, cloud_var, show_cbar = False, hatch_stderr = False):

        # Plot the first type: Ocean
        # --------------------------
        if(remove_high_error):
            plotter = np.ma.masked_where(out_dict[dtype + '_stderr'][ice_idx,:,:] > maxerr, \
                out_dict[dtype + '_slopes'][ice_idx,:,:])
        else:
            plotter = out_dict[dtype + '_slopes'][ice_idx,:,:]
        mesh = pax.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            plotter, \
            cmap = cmap, shading = 'auto', vmin = vmin, vmax = vmax)
        if(show_cbar):
            cbar = plt.colorbar(mesh,\
                ax = pax, orientation='vertical',shrink = shrk, extend = 'both')
            cbar.set_label('Forcing Efficiency\n[Wm$^{-2}$ AI$^{-1}$]')
        if(hatch_stderr):
            hasher = np.ma.masked_where(out_dict[dtype + '_stderr'][ice_idx,:,:] > maxerr, \
                out_dict[dtype + '_stderr'][ice_idx,:,:])
            pax.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
                hasher, hatch = hatch_style, alpha=0., shading = 'auto')

    def plot_cloud_percent(pax, fax, out_dict, dtype, ice_idx, cld_idx, min_cloud,\
            cloud_var, show_cbar = False, cloud_style = '//'):

        # Plot the ocean cloud
        # --------------------
        mask_cloud = np.ma.masked_where(\
            out_dict['raw_cldvals'][ice_idx,:,:,cld_idx] < -9, \
            out_dict['raw_cldvals'][ice_idx,:,:,cld_idx])
        mesh = pax.pcolormesh(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
            mask_cloud, cmap = 'plasma', shading = 'auto', vmin = 0, vmax = 1)
        if(show_cbar):
            cbar = plt.colorbar(mesh,\
                ax = pax, orientation='vertical',shrink = shrk, extend = 'both')
            cbar.set_label('Percent of\n\"' + out_dict['cldval_names'][cld_idx] + '\"')
        if(hatch_cloud):
            if((cld_idx <= 1) | (cld_idx == 4)):
                hasher = np.ma.masked_where(mask_cloud < min_cloud, \
                    mask_cloud)
            else:
                hasher = np.ma.masked_where(mask_cloud > min_cloud, \
                    mask_cloud)
            fax.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
                hasher, hatch = cloud_style, alpha=0., shading = 'auto')
            pax.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
                hasher, hatch = cloud_style, alpha=0., shading = 'auto')

    plot_force_efficiency(ax1, out_dict, dtype, 0, maxerr, vmin, vmax, \
        cmap, cloud_var, show_cbar = False, hatch_stderr = hatch_stderr)
    plot_force_efficiency(ax7, out_dict, dtype, 3, maxerr, vmin, vmax, \
        cmap, cloud_var, show_cbar = False, hatch_stderr = hatch_stderr)
    plot_force_efficiency(ax2, out_dict, dtype, 1, maxerr, vmin, vmax, \
        cmap, cloud_var, show_cbar = False, hatch_stderr = hatch_stderr)
    plot_force_efficiency(ax3, out_dict, dtype, 2, maxerr, vmin, vmax, \
        cmap, cloud_var, show_cbar = True, hatch_stderr = hatch_stderr)

    ax1.set_ylabel(ylabel)
    ax1.set_title('Ocean\n(Ice Conc. <= ' + str(out_dict['ice_maxs'][0]) + '%)')
    ax7.set_title('Mix\n(' + str(out_dict['ice_mins'][3]) + \
        '% < Ice Conc. < ' + str(out_dict['ice_maxs'][3]) + '%)')
    ax2.set_title('Ice\n(Ice Conc. >= ' + str(out_dict['ice_mins'][1]) + '%)')
    ax3.set_title('Land')

    plot_cloud_percent(ax4, ax1, out_dict, dtype, 0, cld_idx, min_cloud, \
        cloud_var, show_cbar = False, cloud_style = cloud_style)
    plot_cloud_percent(ax8, ax7, out_dict, dtype, 3, cld_idx, min_cloud, \
        cloud_var, show_cbar = False, cloud_style = cloud_style)
    plot_cloud_percent(ax5, ax2, out_dict, dtype, 1, cld_idx, min_cloud, \
        cloud_var, show_cbar = False, cloud_style = cloud_style)
    plot_cloud_percent(ax6, ax3, out_dict, dtype, 2, cld_idx, min_cloud, \
        cloud_var, show_cbar = True, cloud_style = cloud_style)

    ax4.set_ylabel(ylabel)
    ax4.set_xlabel(xlabel)
    ax8.set_xlabel(xlabel)
    ax5.set_xlabel(xlabel)
    ax6.set_xlabel(xlabel)
  
    plot_subplot_label(ax1, 'a)', fontsize = 12, yval = 0.27, xval = 75.,  backgroundcolor = 'white')
    plot_subplot_label(ax7, 'b)', fontsize = 12, yval = 0.27, xval = 75.,  backgroundcolor = 'white')
    plot_subplot_label(ax2, 'c)', fontsize = 12, yval = 0.27, xval = 75.,  backgroundcolor = 'white')
    plot_subplot_label(ax3, 'd)', fontsize = 12, yval = 0.27, xval = 75.,  backgroundcolor = 'white')
    plot_subplot_label(ax4, 'e)', fontsize = 12, yval = 0.27, xval = 75.,  backgroundcolor = 'white')
    plot_subplot_label(ax8, 'f)', fontsize = 12, yval = 0.27, xval = 75.,  backgroundcolor = 'white')
    plot_subplot_label(ax5, 'g)', fontsize = 12, yval = 0.27, xval = 75.,  backgroundcolor = 'white')
    plot_subplot_label(ax6, 'h)', fontsize = 12, yval = 0.27, xval = 75.,  backgroundcolor = 'white')

  
    #title_str = 'MODIS ' + cloud_var.upper() + ' Bin Size = ' + \
    #            str(out_dict[cloud_var + '_mins'][1] - \
    #            out_dict[cloud_var + '_mins'][0]) + '\n' + \
    #             'OMI SZA Bin Size = ' + str(out_dict['sza_mins'][1] - \
    #            out_dict['sza_mins'][0])
    title_str = 'Observation-derived Arctic Aerosol Forcing Efficiency'
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

        if(out_dict['data_type'] == 'omi_uvai_pert'):
            dtype_adder = '_pert'
        else:
            dtype_adder = ''

        if('minlat' in out_dict.keys()):
            minlat_adder = '_minlat' + str(int(out_dict['minlat']))
        else:
            minlat_adder = ''

        outname = 'comp_grid_ai_swf_slopecloud_highres' + type_adder + \
            out_adder + hatch_adder + error_adder + dtype_adder + \
            minlat_adder + '_v2.png'
        fig.savefig(outname, dpi = 200)
        print("Saved",outname)
    else: 
        plt.show()


def identify_slope_clear_clean_sfctype(out_dict, sfc_type_idx, \
        cld_idx, min_cloud = 0.95, maxerr = 2, data_type = 'raw', \
        use_error = False):

    mask_cloud = np.ma.masked_where(\
        out_dict['raw_cldvals'][sfc_type_idx,:,:,cld_idx] < -9, \
        out_dict['raw_cldvals'][sfc_type_idx,:,:,cld_idx])
   
    if(use_error): 
        slope_val = data_type + '_stderr'
    else:
        slope_val = data_type + '_slopes'

    cloud_slopes = np.ma.masked_where(mask_cloud < min_cloud, \
        out_dict[slope_val][sfc_type_idx,:,:])
    clear_slopes = np.ma.masked_where(mask_cloud >= min_cloud, \
        out_dict[slope_val][sfc_type_idx,:,:])
    
    cloud_slopes = np.ma.masked_where(\
        out_dict[data_type + '_stderr'][sfc_type_idx,:,:] > maxerr, \
        cloud_slopes)
    clear_slopes = np.ma.masked_where(\
        out_dict[data_type + '_stderr'][sfc_type_idx,:,:] > maxerr, \
        clear_slopes)
   
    return clear_slopes, cloud_slopes

def calc_slope_error(clear_slopes, cloud_slopes):

    clear_mean = np.nanmean(clear_slopes, axis = 0)    
    cloud_mean = np.nanmean(cloud_slopes, axis = 0)    
    clear_std  = np.nanstd(clear_slopes, axis = 0)    
    cloud_std  = np.nanstd(cloud_slopes, axis = 0)    

    test_slopes = np.copy(clear_slopes)
    test_slopes[clear_slopes.mask] = np.nan
    clear_count = np.count_nonzero(~np.isnan(test_slopes), axis = 0)
    test_slopes = np.copy(cloud_slopes)
    test_slopes[cloud_slopes.mask] = np.nan
    cloud_count = np.count_nonzero(~np.isnan(test_slopes), axis = 0)

    clear_error = clear_std / np.sqrt(clear_count)
    cloud_error = cloud_std / np.sqrt(cloud_count)

    return clear_error, cloud_error

def calc_slope_clear_clean_sfctype(out_dict, sfc_type_idx, \
    cld_idx, min_cloud = 0.95, maxerr = 2, data_type = 'raw',\
    use_error = False, mod_slopes = None):

    clear_slopes, cloud_slopes = identify_slope_clear_clean_sfctype(\
        out_dict, sfc_type_idx, cld_idx, min_cloud = min_cloud, \
        maxerr = maxerr, data_type = data_type, use_error = use_error)

    clear_error, cloud_error = calc_slope_error(clear_slopes, cloud_slopes)

    sza_cloud_means = np.nanmean(cloud_slopes, axis = 0)
    sza_clear_means = np.nanmean(clear_slopes, axis = 0)
    sza_cloud_std   = np.std(cloud_slopes, axis = 0)
    sza_clear_std   = np.std(clear_slopes, axis = 0)
    ch7_cloud_means = np.nanmean(cloud_slopes, axis = 1)
    ch7_clear_means = np.nanmean(clear_slopes, axis = 1)
    ch7_cloud_std   = np.std(cloud_slopes, axis = 1)
    ch7_clear_std   = np.std(clear_slopes, axis = 1)

    if(mod_slopes is not None):
        if(mod_slopes == 'add'):
            #print('ADDED    ERROR: mean(CLEAR) = ', np.nanmean(clear_error))
            #print('                mean(CLOUD) = ', np.nanmean(cloud_error))
            sza_clear_means = sza_clear_means + clear_error
            sza_cloud_means = sza_cloud_means + cloud_error
        elif(mod_slopes == 'subtract'):
            #print('SUBTRACT ERROR: mean(CLEAR) = ', np.nanmean(clear_error))
            #print('                mean(CLOUD) = ', np.nanmean(cloud_error))
            sza_clear_means = sza_clear_means - clear_error
            sza_cloud_means = sza_cloud_means - cloud_error
        else:
            print("INVALID MOD SLOPE VALUE. MUST BE \"add\" OR " + \
                "\subtract\"")  
            print("CONTINUING WITH UNALTERED SLOPES")

    return_dict = {}
    if('cod_mins' in out_dict.keys()):
        out_var = 'cod'
    else:
        out_var = 'ch7'

    return_dict['sza_cloud_means'] = sza_cloud_means
    return_dict['sza_cloud_std']   = sza_cloud_std
    return_dict['sza_cloud_err']   = cloud_error
    return_dict['sza_clear_means'] = sza_clear_means
    return_dict['sza_clear_std']   = sza_clear_std
    return_dict['sza_clear_err']   = clear_error
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
# mod_var: 'std' or 'err'
# reverse_sign: Switches the sign of the forcing efficiency to reflect
#               the forcing definition used in the paper (positive = 
#               more energy into surface: darkening) (negative = 
#               less energy into the surface: brightening)                
def plot_slopes_cloud_types_szamean(out_dict,cld_idx = 0, maxerr = 2, \
        data_type = 'raw', xvar = 'sza', remove_high_error = True, \
        hatch_cloud = False, mod_var = 'err', \
        min_cloud = 0.95, use_error = False, reverse_sign = False, \
        save = False):

    ocean_slopes = calc_slope_clear_clean_sfctype(out_dict, 0, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type, \
        use_error = use_error)
    ice_slopes   = calc_slope_clear_clean_sfctype(out_dict, 1, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type, \
        use_error = use_error)
    land_slopes  = calc_slope_clear_clean_sfctype(out_dict, 2, cld_idx, \
        maxerr = maxerr, min_cloud = min_cloud, data_type = data_type, \
        use_error = use_error)

    if(reverse_sign):
        ocean_slopes[xvar + '_cloud_means'] *= -1
        ice_slopes[xvar + '_cloud_means']   *= -1
        land_slopes[xvar + '_cloud_means']  *= -1
        ocean_slopes[xvar + '_clear_means'] *= -1
        ice_slopes[xvar + '_clear_means']   *= -1
        land_slopes[xvar + '_clear_means']  *= -1

    if(out_dict['raw_slopes'].shape[0] == 4):
        include_mix = True
        mix_slopes  = calc_slope_clear_clean_sfctype(out_dict, 3, cld_idx, \
            maxerr = maxerr, min_cloud = min_cloud, data_type = data_type, \
            use_error = use_error)
        if(reverse_sign):
            mix_slopes[xvar + '_cloud_means'] *= -1
            mix_slopes[xvar + '_clear_means'] *= -1
        figsize = (10, 3)
    else:
        include_mix = False
        figsize = (9, 3)
        

    plt.close('all')
    fig = plt.figure(figsize = figsize)
    # Trends
    if(include_mix):
        axs = fig.subplots(nrows = 1, ncols = 4, sharey = True)
        ax1 = axs[0]
        ax2 = axs[1]
        ax3 = axs[2]
        ax4 = axs[3]
        #ax1 = axs[0,0]
        #ax2 = axs[0,1]
        #ax3 = axs[1,0]
        #ax4 = axs[1,1]
    else:
        axs = fig.subplots(nrows = 1, ncols = 3, sharey = True)
        ax1 = axs[0]
        ax3 = axs[1]
        ax4 = axs[2]

    ax1.axhline(0, color = 'k', linestyle = '--')
    ax3.axhline(0, color = 'k', linestyle = '--')
    ax4.axhline(0, color = 'k', linestyle = '--')
    if(include_mix): ax2.axhline(0, color = 'k', linestyle = '--')

 
    # Ocean stuff 
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_cloud_means'], \
        color = 'tab:blue', label = 'Cloud')
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_cloud_means'] - \
        ocean_slopes[xvar + '_cloud_' + mod_var], ':', \
        color = 'tab:blue')
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_cloud_means'] + \
        ocean_slopes[xvar + '_cloud_' + mod_var], ':', \
        color = 'tab:blue')

    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_clear_means'], \
        color = 'tab:orange', label = 'Clear')
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_clear_means'] - \
        ocean_slopes[xvar + '_clear_' + mod_var], ':', \
        color = 'tab:orange')
    ax1.plot(out_dict[xvar + '_mins'], ocean_slopes[xvar + '_clear_means'] + \
        ocean_slopes[xvar + '_clear_' + mod_var], ':', \
        color = 'tab:orange')

    # Ice stuff 
    ax3.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_cloud_means'], \
        color = 'tab:blue', label = 'Cloud')
    ax3.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_cloud_means'] - \
        ice_slopes[xvar + '_cloud_' + mod_var], ':', \
        color = 'tab:blue')
    ax3.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_cloud_means'] + \
        ice_slopes[xvar + '_cloud_' + mod_var], ':', \
        color = 'tab:blue')

    ax3.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_clear_means'], \
        color = 'tab:orange', label = 'Clear')
    ax3.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_clear_means'] - \
        ice_slopes[xvar + '_clear_' + mod_var], ':', \
        color = 'tab:orange')
    ax3.plot(out_dict[xvar + '_mins'], ice_slopes[xvar + '_clear_means'] + \
        ice_slopes[xvar + '_clear_' + mod_var], ':', \
        color = 'tab:orange')

    # Land stuff 
    ax4.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_cloud_means'], \
        color = 'tab:blue', label = 'Cloud')
    ax4.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_cloud_means'] - \
        land_slopes[xvar + '_cloud_' + mod_var], ':', \
        color = 'tab:blue')
    ax4.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_cloud_means'] + \
        land_slopes[xvar + '_cloud_' + mod_var], ':', \
        color = 'tab:blue')

    ax4.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_clear_means'], \
        color = 'tab:orange', label = 'Clear')
    ax4.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_clear_means'] - \
        land_slopes[xvar + '_clear_' + mod_var], ':', \
        color = 'tab:orange')
    ax4.plot(out_dict[xvar + '_mins'], land_slopes[xvar + '_clear_means'] + \
        land_slopes[xvar + '_clear_' + mod_var], ':', \
        color = 'tab:orange')

    ax1_mins = ax1.get_ylim()[0]
    ax3_mins = ax3.get_ylim()[0]
    ax4_mins = ax4.get_ylim()[0]
    ax1_maxs = ax1.get_ylim()[1]
    ax3_maxs = ax3.get_ylim()[1]
    ax4_maxs = ax4.get_ylim()[1]

    if(include_mix):
        ax2.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_cloud_means'], \
            color = 'tab:blue', label = 'Cloud')
        ax2.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_cloud_means'] - \
            mix_slopes[xvar + '_cloud_std'], ':', \
            color = 'tab:blue')
        ax2.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_cloud_means'] + \
            mix_slopes[xvar + '_cloud_std'], ':', \
            color = 'tab:blue')

        ax2.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_clear_means'], \
            color = 'tab:orange', label = 'Clear')
        ax2.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_clear_means'] - \
            mix_slopes[xvar + '_clear_std'], ':', \
            color = 'tab:orange')
        ax2.plot(out_dict[xvar + '_mins'], mix_slopes[xvar + '_clear_means'] + \
            mix_slopes[xvar + '_clear_std'], ':', \
            color = 'tab:orange')

        ax2_mins = ax2.get_ylim()[0]
        ax2_maxs = ax2.get_ylim()[1]

    if(reverse_sign):
        ax1.legend(loc = 2)
    else:
        ax4.legend()

    all_mins = [ax1_mins, ax3_mins, ax4_mins]
    all_maxs = [ax1_maxs, ax3_maxs, ax4_maxs]
    if(include_mix):
        all_mins += [ax2_mins]
        all_maxs += [ax2_maxs]

    ax1.set_ylim(np.min(all_mins), np.max(all_maxs))
    ax3.set_ylim(np.min(all_mins), np.max(all_maxs))
    ax4.set_ylim(np.min(all_mins), np.max(all_maxs))
    if(include_mix): 
        ax2.set_ylim(np.min(all_mins), np.max(all_maxs))

    ax1.grid(alpha = 0.5)
    ax3.grid(alpha = 0.5)
    ax4.grid(alpha = 0.5)
    if(include_mix): 
        ax2.grid(alpha = 0.5)

    ax1.set_title('Ocean\n(Ice Conc. <= ' + str(out_dict['ice_maxs'][0]) + '%)')
    ax3.set_title('Ice\n(Ice Conc. >= ' + str(out_dict['ice_mins'][1]) + '%)')
    ax4.set_title('Land')
    if(include_mix):  
        ax2.set_title('Mix\n(' + str(out_dict['ice_mins'][3]) + '% < Ice Conc. < ' + str(out_dict['ice_maxs'][3]) + '%)')
 
    if(xvar == 'sza'): labeltxt = 'OMI Solar Zenith Angle [$^{o}$]' 
    elif(xvar == 'ch7'): labeltxt = 'MODIS CH7'
    elif(xvar == 'cod'): labeltxt = 'MODIS COD'
    ax1.set_xlabel(labeltxt)
    ax3.set_xlabel(labeltxt)
    ax4.set_xlabel(labeltxt)
    if(include_mix):
        ax2.set_xlabel(labeltxt)

    if(use_error):
        labeltext =  'AI-SWF Error [Wm$^{-2}$AI$^{-1}$]'
    else:
        labeltext =  'Forcing Efficiency\n[Wm$^{-2}$AI$^{-1}$]'
    ax1.set_ylabel(labeltext)
    #if(include_mix):
    #    ax3.set_ylabel(labeltext)
    #ax2.set_ylabel('AI-SWF Slope [Wm$^{-2}$]')
    #ax3.set_ylabel('AI-SWF Slope [Wm$^{-2}$]')

    if(reverse_sign):
        ymax = 0.80 * np.max(all_maxs)
    else:
        ymax = 0.60 * np.max(all_maxs)
    plot_subplot_label(ax1, 'a)', fontsize = 12, xval = 76., yval = ymax)
    plot_subplot_label(ax2, 'b)', fontsize = 12, xval = 76., yval = ymax)
    plot_subplot_label(ax3, 'c)', fontsize = 12, xval = 76., yval = ymax)
    plot_subplot_label(ax4, 'd)', fontsize = 12, xval = 76., yval = ymax)
    #plot_subplot_label(ax1, 'a)', fontsize = 12, xval = 76., yval = 7.8, backgroundcolor = 'white')
    #plot_subplot_label(ax2, 'b)', fontsize = 12, xval = 76., yval = 7.8, backgroundcolor = 'white')
    #plot_subplot_label(ax3, 'c)', fontsize = 12, xval = 76., yval = 7.8, backgroundcolor = 'white')
    #plot_subplot_label(ax4, 'd)', fontsize = 12, xval = 76., yval = 7.8, backgroundcolor = 'white')

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

        if(out_dict['data_type'] == 'omi_uvai_pert'):
            dtype_adder = '_pert'
        else:
            dtype_adder = ''

        if('minlat' in out_dict.keys()):
            minlat_adder = '_minlat' + str(int(out_dict['minlat']))
        else:
            minlat_adder = ''


        if(include_mix):
            mix_add = '_mix_v2'
            #mix_add = '_mix'
        else:
            mix_add = ''

        if(reverse_sign):
            rvs_adder = '_revsign'
        else:
            rvs_adder = ''

        if(use_error):
            outname = 'comp_grid_ai_swf_errorcloud_' + xvar + 'means' + type_adder + \
                smth_adder + mix_add + dtype_adder + minlat_adder + '.png'
        else:
            outname = 'comp_grid_ai_swf_slopecloud_' + xvar + 'means' + type_adder + \
                smth_adder + mix_add + dtype_adder + minlat_adder + rvs_adder + '.png'
                #smth_adder + mix_add + '.png'
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
    norm = mpl.colors.BoundaryNorm(colorvals, cmap2.N, extend = 'upper')
        
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

# dtype: 'RAW' or 'PERT'
def calc_pcnt_aerosol_over_type(date_list, min_AI, ax = None, \
        minlat = 70., dtype = 'RAW', save = False, local_dir = None): 

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
        coloc_data = read_colocated(date, minlat = minlat, local_dir = local_dir)
    
        # Test finding the indices with high aerosol
        # ------------------------------------------
        
        mask_data = np.ma.masked_where(coloc_data['OMI_' + dtype] < min_AI, \
            coloc_data['OMI_' + dtype])
        
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
    
        pcnt_ice_val = (count_ice / count_data_val) * 100.
        pcnt_mix_val = (count_mix / count_data_val) * 100.
        pcnt_ocn_val = (count_ocean / count_data_val) * 100.
        pcnt_lnd_val = (count_land / count_data_val) * 100.
        pcnt_oth_val = (count_othr / count_data_val) * 100.
   
        if((pcnt_ice_val + pcnt_mix_val + pcnt_ocn_val + pcnt_lnd_val + pcnt_oth_val) != 100):
            print("BAD SUM: ", date)
 
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
        print("%s %5d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %5d" % \
            (date, count_data[ii], pcnt_ice[ii], pcnt_mix[ii], pcnt_ocn[ii], \
            pcnt_lnd[ii], pcnt_oth[ii], (pcnt_ice[ii]+ pcnt_mix[ii]+ pcnt_ocn[ii]+ \
            pcnt_lnd[ii]+ pcnt_oth[ii]), total_count[ii]))
            #pcnt_lnd[ii], pcnt_oth[ii], total_count[ii], total2_count[ii]))

    print("TOTAL ACROSS SWATHS:", np.sum(total_count))

    in_ax = True 
    if(ax is None): 
        plt.close('all')
        in_ax = False
        fig = plt.figure(figsize = (9, 3.5))
        ax = fig.add_subplot(1,1,1)

    int_dates = np.array([int(dl) for dl in date_list])
    #print(int_dates)
    xvals = np.arange(len(date_list))
    
    ax.bar(xvals, pcnt_ice, label = 'Ice')
    ax.bar(xvals, pcnt_mix, bottom = pcnt_ice, label = 'Mix')
    ax.bar(xvals, pcnt_ocn, bottom = pcnt_ice + pcnt_mix, label = 'Ocean')
    ax.bar(xvals, pcnt_lnd, bottom = pcnt_ice + pcnt_mix + pcnt_ocn, label = 'Land')
    ax.bar(xvals, pcnt_oth, bottom = pcnt_ice + pcnt_mix + pcnt_ocn + pcnt_lnd, label = 'Other')
    ax.set_xlim([np.min(xvals), np.max(xvals)])
    #lines, labels = ax.get_legend_handles_labels()
    #lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    #ax.legend()
    #plt.legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0, 0.01, 1, 1),\
    ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5),\
        ncol=1)


    if(not in_ax):
        ax.set_title('SSMI/S Surface Types of L2 Arctic OMI Pixels' + \
            ' with AI > ' + str(min_AI))
        ax.set_xlabel('Swath Number')
        ax.set_ylabel('Percent of OMI Aerosol Pixels [%]')
        ax.set_ylim([0, 100])
        fig.tight_layout()
        if(save):
            outname = 'aerosol_over_type_swath_minAI' + dtype + '_'+ str(int(min_AI*10)) + '_minlat' + str(int(minlat)) + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show() 

# = = = = = = =
#
# These taken from: https://stackoverflow.com/questions/18386210/annotating-ranges-of-data
#
# = = = = = = =

##def rotate_point(x, y, angle_rad):
##    cos,sin = np.cos(angle_rad),np.sin(angle_rad)
##    return cos*x-sin*y,sin*x+cos*y
##
##def draw_brace(ax, span, position, text, text_pos, brace_scale=1.0, beta_scale=300., rotate=False, rotate_text=False):
##    '''
##        all positions and sizes are in axes units
##        span: size of the curl
##        position: placement of the tip of the curl
##        text: label to place somewhere
##        text_pos: position for the label
##        beta_scale: scaling for the curl, higher makes a smaller radius
##        rotate: true rotates to place the curl vertically
##        rotate_text: true rotates the text vertically        
##    '''
##    # get the total width to help scale the figure
##    ax_xmin, ax_xmax = ax.get_xlim()
##    xax_span = ax_xmax - ax_xmin
##    resolution = int(span/xax_span*100)*2+1 # guaranteed uneven
##    beta = beta_scale/xax_span # the higher this is, the smaller the radius
##    # center the shape at (0, 0)
##    x = np.linspace(-span/2., span/2., resolution)
##    # calculate the shape
##    x_half = x[:int(resolution/2)+1]
##    y_half_brace = (1/(1.+np.exp(-beta*(x_half-x_half[0])))
##                + 1/(1.+np.exp(-beta*(x_half-x_half[-1]))))
##    y = np.concatenate((y_half_brace, y_half_brace[-2::-1]))
##    # put the tip of the curl at (0, 0)
##    max_y = np.max(y)    
##    min_y = np.min(y)
##    y /= (max_y-min_y)
##    y *= brace_scale
##    y -= max_y
##    # rotate the trace before shifting
##    if rotate:
##        x,y = rotate_point(x, y, np.pi/2)
##    # shift to the user's spot   
##    x += position[0]        
##    y += position[1]
##    ax.autoscale(False)
##    ax.plot(x, y, color='black', lw=1, clip_on=False)
##    # put the text
##    ax.text(text_pos[0], text_pos[1], text, ha='center', va='bottom', rotation=90 if rotate_text else 0)


def draw_brace(ax, xspan, yy, text):
    """Draws an annotated brace outside the axes."""
    xmin, xmax = xspan
    xspan = xmax - xmin
    ax_xmin, ax_xmax = ax.get_xlim()
    xax_span = ax_xmax - ax_xmin

    ymin, ymax = ax.get_ylim()
    yspan = ymax - ymin
    resolution = int(xspan/xax_span*100)*2+1 # guaranteed uneven
    beta = 300./xax_span # the higher this is, the smaller the radius

    x = np.linspace(xmin, xmax, resolution)
    x_half = x[:int(resolution/2)+1]
    y_half_brace = (1/(1.+np.exp(-beta*(x_half-x_half[0])))
                + 1/(1.+np.exp(-beta*(x_half-x_half[-1]))))
    y = np.concatenate((y_half_brace, y_half_brace[-2::-1]))
    y = yy + (.05*y - .01)*yspan # adjust vertical position

    ax.autoscale(False)
    ax.plot(x, -y, color='black', lw=1, clip_on=False)

    ax.text((xmax+xmin)/2., -yy-.17*yspan, text, ha='center', va='bottom', fontsize = 8)



# dtype: 'RAW' or 'PERT'
def calc_pcnt_aerosol_over_type_v2(sim_name, min_AI, ax = None, \
        minlat = 70., dtype = 'RAW', save = False): 

    file_list = glob('neuralnet_output/test_calc_out_' + sim_name + '*.hdf5')

    count_data   = np.full(len(file_list), np.nan)
    total_count  = np.full(len(file_list), np.nan)
    total2_count = np.full(len(file_list), np.nan)
    pcnt_ice     = np.full(len(file_list), np.nan)
    pcnt_mix     = np.full(len(file_list), np.nan)
    pcnt_ocn     = np.full(len(file_list), np.nan)
    pcnt_lnd     = np.full(len(file_list), np.nan)
    pcnt_oth     = np.full(len(file_list), np.nan)
    pcnt_ice_cld = np.full(len(file_list), np.nan)
    pcnt_mix_cld = np.full(len(file_list), np.nan)
    pcnt_ocn_cld = np.full(len(file_list), np.nan)
    pcnt_lnd_cld = np.full(len(file_list), np.nan)
    pcnt_oth_cld = np.full(len(file_list), np.nan)
    pcnt_ice_clr = np.full(len(file_list), np.nan)
    pcnt_mix_clr = np.full(len(file_list), np.nan)
    pcnt_ocn_clr = np.full(len(file_list), np.nan)
    pcnt_lnd_clr = np.full(len(file_list), np.nan)
    pcnt_oth_clr = np.full(len(file_list), np.nan)

    for ii, ff in enumerate(file_list):
    
        # Read data from the file
        # -----------------------
        data = h5py.File(ff,'r')
        print(ff)

        # Test finding the indices with high aerosol
        # ------------------------------------------
        mask_data = np.ma.masked_invalid(data['omi_uvai_pert'][:,:])
        mask_data = np.ma.masked_where(data['omi_uvai_pert'][:,:] < min_AI, \
            mask_data)
        mask_data = np.ma.masked_where(data['omi_lat'][:,:] < minlat, mask_data)
        mask_data = np.ma.masked_where(data['nsidc_ice'][:,:] == -999., mask_data)

        mask_ice = np.ma.masked_where((\
            #data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.), \
             data['nsidc_ice'][:,:] == -999.) | \
            (data['nsidc_ice'][:,:] > 100.) | \
            (data['nsidc_ice'][:,:] < 80.), \
            mask_data)
        mask_io_mix = np.ma.masked_where((\
            #data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.), \
             data['nsidc_ice'][:,:] == -999.) | \
            (data['nsidc_ice'][:,:] >= 80.) | \
            (data['nsidc_ice'][:,:] < 20.), \
            mask_data)
        mask_ocean = np.ma.masked_where((\
             data['nsidc_ice'][:,:] == -999.) | \
            (data['nsidc_ice'][:,:] >= 20.) | \
            (data['nsidc_ice'][:,:] < 0.), \
            mask_data)
        mask_land = np.ma.masked_where((\
             data['nsidc_ice'][:,:] == -999.) | \
            (data['nsidc_ice'][:,:] != 254.), \
            mask_data)
        mask_othr = np.ma.masked_where((\
             data['nsidc_ice'][:,:] == -999.) | \
            (data['nsidc_ice'][:,:] < 251) | \
            (data['nsidc_ice'][:,:] > 253), \
            mask_data) 

        count_data_val = mask_data.compressed().size
        count_ice      = mask_ice.compressed().size
        count_mix      = mask_io_mix.compressed().size
        count_ocean    = mask_ocean.compressed().size
        count_land     = mask_land.compressed().size
        count_othr     = mask_othr.compressed().size

        totals         = count_ice + count_mix + count_ocean + count_land + \
            count_othr
    
        pcnt_ice_val = (count_ice / count_data_val) * 100.
        pcnt_mix_val = (count_mix / count_data_val) * 100.
        pcnt_ocn_val = (count_ocean / count_data_val) * 100.
        pcnt_lnd_val = (count_land / count_data_val) * 100.
        pcnt_oth_val = (count_othr / count_data_val) * 100.
  
        calc_sum = pcnt_ice_val + pcnt_mix_val + pcnt_ocn_val + pcnt_lnd_val + pcnt_oth_val  
        if((calc_sum) != 100):
            print("BAD SUM: ", ff, calc_sum)
 
        count_data[ii]  = count_data_val
        pcnt_ice[ii]    = pcnt_ice_val 
        pcnt_mix[ii]    = pcnt_mix_val 
        pcnt_ocn[ii]    = pcnt_ocn_val 
        pcnt_lnd[ii]    = pcnt_lnd_val 
        pcnt_oth[ii]    = pcnt_oth_val 
        total_count[ii] = totals

    print('    date      cnt    ice   mix    ocn    lnd    oth     tot')
    for ii, date in enumerate(file_list):
        print("%s %5d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %5d" % \
            (date, count_data[ii], pcnt_ice[ii], pcnt_mix[ii], pcnt_ocn[ii], \
            pcnt_lnd[ii], pcnt_oth[ii], (pcnt_ice[ii]+ pcnt_mix[ii]+ pcnt_ocn[ii]+ \
            pcnt_lnd[ii]+ pcnt_oth[ii]), total_count[ii]))
            #pcnt_lnd[ii], pcnt_oth[ii], total_count[ii], total2_count[ii]))

    print("TOTAL ACROSS SWATHS:", np.sum(total_count))

    in_ax = True 
    if(ax is None): 
        plt.close('all')
        in_ax = False
        fig = plt.figure(figsize = (11, 3.5))
        ax = fig.add_subplot(1,1,1)

    xvals = np.arange(len(file_list))
    
    ax.bar(xvals, pcnt_ice, label = 'Ice')
    ax.bar(xvals, pcnt_mix, bottom = pcnt_ice, label = 'Mix')
    ax.bar(xvals, pcnt_ocn, bottom = pcnt_ice + pcnt_mix, label = 'Ocean')
    ax.bar(xvals, pcnt_lnd, bottom = pcnt_ice + pcnt_mix + pcnt_ocn, label = 'Land')
    ax.bar(xvals, pcnt_oth, bottom = pcnt_ice + pcnt_mix + pcnt_ocn + pcnt_lnd, label = 'Other')
    ax.set_xlim([np.min(xvals), np.max(xvals)])
    ax.xaxis.set_visible(False)
    #lines, labels = ax.get_legend_handles_labels()
    #lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    #ax.legend()
    #plt.legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0, 0.01, 1, 1),\
    ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5),\
        ncol=1)
    
    # Add brackets to the figure
    # --------------------------
    int_dates = np.array([int(ff.strip().split('/')[-1].split('_')[-1][:6]) for ff in file_list])
    unique_dates = np.unique(int_dates)
    print(int_dates)
    for idate in unique_dates:
        # Find indices of each date
        # -------------------------
        date_idxs = np.where(int_dates == idate)
        beg_idx = date_idxs[0][0]
        end_idx = date_idxs[0][-1]
        print(idate, beg_idx, end_idx)
    
        #draw_brace(ax, 20, (20, 10), 'test', (20, -5), brace_scale=1.0, beta_scale=300., rotate=True, rotate_text=False)
        dt_date = datetime.strptime(str(idate), '%Y%m')
        draw_brace(ax, (beg_idx, end_idx), 0.1, dt_date.strftime('%b\n%Y'))

    if(not in_ax):
        ax.set_title('SSMIS Surface Types of L2 Arctic OMI Aerosol Swaths' + \
            '\nPlotted for UVAI > ' + str(min_AI))
        ax.set_xlabel('Swath Number')
        ax.set_ylabel('Percent of OMI Aerosol Pixels [%]')
        ax.set_ylim([0, 100])
        fig.tight_layout()
        if(save):
            outname = 'aerosol_over_type_swath_minAI' + dtype + '_'+ str(int(min_AI*10)) + '_minlat' + str(int(minlat)) + '_v2.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
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
   
    #plt.suptitle('Minimum AI = ' + str(min_ai))
    plt.suptitle('Percent Coverage of OMI  = ' + str(min_ai))
 
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

# omi_data_type = 'thl' or 'lin'
def plot_type_forcing_v3_all_months(all_month_values, OMI_monthly_data, \
        minlat = 70., maxlat = 87., use_szas = False, \
        ai_thresh = 0.7, cld_idx = 0, maxerr = 2, \
        min_cloud = 0.95, data_type = 'raw', trend_type = 'standard', \
        hstyle = '.....', \
        plot_zonal_avgs = False, plot_pvals = True, \
        omi_data_type = 'linregress', \
        version4 = False,
        save = False,debug = False):

    forcing_trends = np.full((6, OMI_monthly_data['LAT'].shape[0], \
            OMI_monthly_data['LAT'].shape[1]), np.nan)
    forcing_pval   = np.full((6, OMI_monthly_data['LAT'].shape[0], \
            OMI_monthly_data['LAT'].shape[1]), np.nan)
    forcing_uncert = np.full((6, OMI_monthly_data['LAT'].shape[0], \
            OMI_monthly_data['LAT'].shape[1]), np.nan)
  
    for ii in range(6): 
        # Calculate the trend in the forcings
        # -----------------------------------
        forcing_trends[ii,:,:], forcing_pval[ii,:,:], \
            forcing_uncert[ii,:,:] = calc_forcing_grid_trend(\
            all_month_values[ii::6,:,:], trend_type)

    mask_pvals = np.ma.masked_where((OMI_monthly_data['LAT'] < minlat) | \
        (forcing_pval > 0.05), forcing_pval)

    plt.close('all')
    if(plot_zonal_avgs):
        figsize = (9, 10)
        ncols = 4
        nrows = 3
    else:
        figsize = (9, 5)
        ncols = 2
        nrows = 3

    fig = plt.figure(figsize = figsize)
    
    ax1 = fig.add_subplot(ncols,nrows,1, projection = mapcrs)
    ax2 = fig.add_subplot(ncols,nrows,2, projection = mapcrs)
    ax3 = fig.add_subplot(ncols,nrows,3, projection = mapcrs)
    ax4 = fig.add_subplot(ncols,nrows,4, projection = mapcrs)
    ax5 = fig.add_subplot(ncols,nrows,5, projection = mapcrs)
    ax6 = fig.add_subplot(ncols,nrows,6, projection = mapcrs)

    if(plot_zonal_avgs):
        ax7  = fig.add_subplot(4,3,7)
        ax8  = fig.add_subplot(4,3,8)
        ax9  = fig.add_subplot(4,3,9)
        ax10 = fig.add_subplot(4,3,10)
        ax11 = fig.add_subplot(4,3,11)
        ax12 = fig.add_subplot(4,3,12)
   
    slope_max = 0.8 
    mesh = ax1.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            forcing_trends[0,:,:], transform = datacrs, shading = 'auto', \
            vmin = -slope_max, vmax = slope_max, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax1, label = 'Forcing Trend\n[W/m2/period]')
    ax1.coastlines()
    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_title('April')   
    ax1.set_boundary(circle, transform=ax1.transAxes)
    if(plot_pvals):
        ax1.pcolor(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            mask_pvals[0,:,:], hatch = hstyle, alpha = 0.0, \
            shading = 'auto', transform = datacrs)
 
    mesh = ax2.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            forcing_trends[1,:,:], transform = datacrs, shading = 'auto', \
            vmin = -slope_max, vmax = slope_max, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax2, label = 'Forcing Trend\n[W/m2/period]')
    ax2.coastlines()
    ax2.set_extent([-180,180,minlat,90], datacrs)
    ax2.set_title('May')   
    ax2.set_boundary(circle, transform=ax2.transAxes)
    if(plot_pvals):
        ax2.pcolor(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            mask_pvals[1,:,:], hatch = hstyle, alpha = 0.0, \
            shading = 'auto', transform = datacrs)
    
    mesh = ax3.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            forcing_trends[2,:,:], transform = datacrs, shading = 'auto', \
            vmin = -slope_max, vmax = slope_max, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax3, label = 'Forcing Trend\n[W/m2/period]')
    ax3.coastlines()
    ax3.set_extent([-180,180,minlat,90], datacrs)
    ax3.set_title('June')   
    ax3.set_boundary(circle, transform=ax3.transAxes)
    if(plot_pvals):
        ax3.pcolor(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            mask_pvals[2,:,:], hatch = hstyle, alpha = 0.0, \
            shading = 'auto', transform = datacrs)
    
    mesh = ax4.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            forcing_trends[3,:,:], transform = datacrs, shading = 'auto', \
            vmin = -slope_max, vmax = slope_max, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax4, label = 'Forcing Trend\n[W/m2/period]')
    ax4.coastlines()
    ax4.set_extent([-180,180,minlat,90], datacrs)
    ax4.set_title('July')   
    ax4.set_boundary(circle, transform=ax4.transAxes)
    if(plot_pvals):
        ax4.pcolor(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            mask_pvals[3,:,:], hatch = hstyle, alpha = 0.0, \
            shading = 'auto', transform = datacrs)
    
    mesh = ax5.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            forcing_trends[4,:,:], transform = datacrs, shading = 'auto', \
            vmin = -slope_max, vmax = slope_max, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax5, label = 'Forcing Trend\n[W/m2/period]')
    ax5.coastlines()
    ax5.set_extent([-180,180,minlat,90], datacrs)
    ax5.set_title('August')  
    ax5.set_boundary(circle, transform=ax5.transAxes)
    if(plot_pvals):
        ax5.pcolor(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            mask_pvals[4,:,:], hatch = hstyle, alpha = 0.0, \
            shading = 'auto', transform = datacrs)
    
    mesh = ax6.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            forcing_trends[5,:,:], transform = datacrs, shading = 'auto', \
            vmin = -slope_max, vmax = slope_max, cmap = 'bwr')
    cbar = plt.colorbar(mesh, ax = ax6, label = 'Forcing Trend\n[W/m2/period]')
    ax6.coastlines()
    ax6.set_extent([-180,180,minlat,90], datacrs)
    ax6.set_title('September')
    ax6.set_boundary(circle, transform=ax6.transAxes)
    if(plot_pvals):
        ax6.pcolor(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
            mask_pvals[5,:,:], hatch = hstyle, alpha = 0.0, \
            shading = 'auto', transform = datacrs)

    if(plot_zonal_avgs):
  
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

        ax7.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[0])
        ax7.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[0] + zonal_stds[0], linestyle = '--', color = 'tab:blue')
        ax7.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[0] - zonal_stds[0], linestyle = '--', color = 'tab:blue')
        ax7.axhline(0, linestyle = ':', color = 'k')
        ax7.set_title('April')
        ax7.set_xlabel('Latitude')    
        ax7.set_ylabel('Zonal Avg. Forcing Trend')    
 
        ax8.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[1])
        ax8.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[1] + zonal_stds[1], linestyle = '--', color = 'tab:blue')
        ax8.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[1] - zonal_stds[1], linestyle = '--', color = 'tab:blue')
        ax8.axhline(0, linestyle = ':', color = 'k')
        ax8.set_title('May')
        ax8.set_xlabel('Latitude')    
        ax8.set_ylabel('Zonal Avg. Forcing Trend')    
 
        ax9.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[2])
        ax9.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[2] + zonal_stds[2], linestyle = '--', color = 'tab:blue')
        ax9.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[2] - zonal_stds[2], linestyle = '--', color = 'tab:blue')
        ax9.axhline(0, linestyle = ':', color = 'k')
        ax9.set_title('June')
        ax9.set_xlabel('Latitude')    
        ax9.set_ylabel('Zonal Avg. Forcing Trend')    
 
        ax10.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[3])
        ax10.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[3] + zonal_stds[3], linestyle = '--', color = 'tab:blue')
        ax10.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[3] - zonal_stds[3], linestyle = '--', color = 'tab:blue')
        ax10.axhline(0, linestyle = ':', color = 'k')
        ax10.set_title('July')
        ax10.set_xlabel('Latitude')    
        ax10.set_ylabel('Zonal Avg. Forcing Trend')    
 
        ax11.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[4])
        ax11.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[4] + zonal_stds[4], linestyle = '--', color = 'tab:blue')
        ax11.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[4] - zonal_stds[4], linestyle = '--', color = 'tab:blue')
        ax11.axhline(0, linestyle = ':', color = 'k')
        ax11.set_title('August')
        ax11.set_xlabel('Latitude')    
        ax11.set_ylabel('Zonal Avg. Forcing Trend')    
 
        ax12.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[5], label = 'Avg. Trend')
        ax12.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[5] + zonal_stds[5], linestyle = '--', color = 'tab:blue')
        ax12.plot(OMI_monthly_data['LAT'][:,0], zonal_avgs[5] - zonal_stds[5], linestyle = '--', color = 'tab:blue', \
            label = '+/- St. Dev. Trend')
        ax12.legend(fontsize = 10)
        ax12.axhline(0, linestyle = ':', color = 'k')
        ax12.set_title('September')
        ax12.set_xlabel('Latitude')    
        ax12.set_ylabel('Zonal Avg. Forcing Trend')    

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

    plot_subplot_label(ax1, 'a)', fontsize = 14, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 14, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 14, backgroundcolor = None)
    plot_subplot_label(ax4, 'd)', fontsize = 14, backgroundcolor = None)
    plot_subplot_label(ax5, 'e)', fontsize = 14, backgroundcolor = None)
    plot_subplot_label(ax6, 'f)', fontsize = 14, backgroundcolor = None)

    #plt.suptitle('Forcing Estimate 3:\nDaily-averaged Single Month Based', weight = 'bold', fontsize = 12) 
    plt.suptitle('Estimated Aerosol Forcing Trend\n2005 - 2020', \
        weight = 'bold', fontsize = 12) 
    fig.tight_layout()

    if(save):
        if(plot_zonal_avgs):
            zonal_add = '_zonal'
        else:
            zonal_add = ''

        if(plot_pvals):
            pval_add = '_pvals'
        else:
            pval_add = ''

        if(version4):
            version_add = 'v4'
        else:
            version_add = 'v3'

        outname = 'calc_arctic_forcing_' + version_add + '_monthly_' + \
            trend_type + zonal_add + pval_add + '_' + omi_data_type + \
            '_minlat' + str(int(minlat)) + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else: 
        plt.show()


def plot_type_forcing_v3_all_months_arctic_avg_combined(all_month_values, \
        OMI_monthly_data, version4 = False, max_pval = 0.05, \
        slope_type = 'lin', 
        horiz_orient = True, save = False):

    # Prep the overall figure
    # -----------------------
    plt.close('all')
    if(horiz_orient):
        fig = plt.figure(figsize = (13, 6.5))
        axs  = fig.subplots(nrows = 3, ncols = 6, sharex = True, sharey = True)
        ax1  = axs[0,0]  # control April 
        ax2  = axs[0,1]  # control May 
        ax3  = axs[0,2]  # control June 
        ax4  = axs[0,3]  # control July 
        ax5  = axs[0,4]  # control August
        ax6  = axs[0,5]  # control September
        ax7  = axs[1,0]  # low     April 
        ax8  = axs[1,1]  # low     May 
        ax9  = axs[1,2]  # low     June 
        ax10 = axs[1,3]  # low     July 
        ax11 = axs[1,4]  # low     August
        ax12 = axs[1,5]  # low     September
        ax13 = axs[2,0]  # high    April 
        ax14 = axs[2,1]  # high    May 
        ax15 = axs[2,2]  # high    June 
        ax16 = axs[2,3]  # high    July 
        ax17 = axs[2,4]  # high    August
        ax18 = axs[2,5]  # high    September
        orient_add = '_horiz'
    else:
        #fig = plt.figure(figsize = (7.5, 11))
        fig = plt.figure(figsize = (8, 9.5))
        axs  = fig.subplots(nrows = 6, ncols = 3, sharex = True, sharey = True)
        ax1  = axs[0,0]  # control April 
        ax2  = axs[1,0]  # control May 
        ax3  = axs[2,0]  # control June 
        ax4  = axs[3,0]  # control July 
        ax5  = axs[4,0]  # control August
        ax6  = axs[5,0]  # control September
        ax7  = axs[0,1]  # low     April 
        ax8  = axs[1,1]  # low     May 
        ax9  = axs[2,1]  # low     June 
        ax10 = axs[3,1]  # low     July 
        ax11 = axs[4,1]  # low     August
        ax12 = axs[5,1]  # low     September
        ax13 = axs[0,2]  # high    April 
        ax14 = axs[1,2]  # high    May 
        ax15 = axs[2,2]  # high    June 
        ax16 = axs[3,2]  # high    July 
        ax17 = axs[4,2]  # high    August
        ax18 = axs[5,2]  # high    September
        orient_add = '_vert'

    # Plot the first row: control (65 - 87)
    # -------------------------------------
    ax_list = [ax1, ax2, ax3, ax4, ax5, ax6]
    first_min, first_max = \
        plot_type_forcing_v3_all_months_arctic_avg(all_month_values, \
        OMI_monthly_data, \
        axs = ax_list, \
        minlat = 65., maxlat = 87., omi_data_type = 'lin', 
        version4 = True, 
        save = False, debug = False,max_pval = 0.05)


    # Plot the second row: low (65 - 75)
    # ----------------------------------
    ax_list = [ax7, ax8, ax9, ax10, ax11, ax12]
    secnd_min, secnd_max = \
        plot_type_forcing_v3_all_months_arctic_avg(all_month_values, \
        OMI_monthly_data, \
        axs = ax_list, \
        minlat = 65., maxlat = 75., omi_data_type = 'lin', 
        version4 = True, 
        save = False, debug = False,max_pval = 0.05)


    # Plot the third row: high (75 - 87)
    # ----------------------------------
    ax_list = [ax13, ax14, ax15, ax16, ax17, ax18]
    third_min, third_max = \
    plot_type_forcing_v3_all_months_arctic_avg(all_month_values, \
        OMI_monthly_data, \
        axs = ax_list, \
        minlat = 75., maxlat = 87., omi_data_type = 'lin', 
        version4 = True, 
        save = False, debug = False,max_pval = 0.05)

    if(horiz_orient):
        ax1.set_ylabel('Forcing [Wm$^{-2}$]')
        ax7.set_ylabel('Forcing [Wm$^{-2}$]')
        ax13.set_ylabel('Forcing [Wm$^{-2}$]')
    else:
        ax1.set_ylabel('Forcing [Wm$^{-2}$]')
        ax2.set_ylabel('Forcing [Wm$^{-2}$]')
        ax3.set_ylabel('Forcing [Wm$^{-2}$]')
        ax4.set_ylabel('Forcing [Wm$^{-2}$]')
        ax5.set_ylabel('Forcing [Wm$^{-2}$]')
        ax6.set_ylabel('Forcing [Wm$^{-2}$]')

        ax1.set_title('65$^{o}$ - 87$^{o}$ N', weight = 'bold')
        ax7.set_title('65$^{o}$ - 75$^{o}$ N', weight = 'bold')
        ax13.set_title('75$^{o}$ - 87$^{o}$ N', weight = 'bold')

        fig.subplots_adjust(left = 0.2)

        row_label_size = 12 
        fig.text(0.07, 0.82, 'April', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.69, 'May', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.555, 'June', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.42, 'July', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.295, 'August', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.16, 'September', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)



    plt.suptitle('Arctic Regional Average Aerosol Forcing\n', 
        weight = 'bold') 

    all_ranges = [first_min, first_max, secnd_min, secnd_max, \
        third_min, third_max]

    for pax in axs.flatten():
        pax.set_ylim(np.min(all_ranges), np.max(all_ranges))

    #fig.tight_layout()

    if(save):
        file_add = ''
        trend_add = ''
        outname = 'calc_arctic_forcing_v3_monthly_arcticavg_aimin0_useintcpt_' + \
            slope_type + orient_add + '_combined.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


# omi_data_type: 'lin' or 'thl'
def plot_type_forcing_v3_all_months_arctic_avg(all_month_values, \
        OMI_monthly_data, \
        axs = None, \
        minlat = 70., maxlat = 87., use_szas = False, \
        ai_thresh = 0.7, cld_idx = 0, maxerr = 2, \
        min_cloud = 0.95, data_type = 'raw', trend_type = 'standard', \
        omi_data_type = 'lin', 
        version4 = False, 
        save = False,debug = False, month_values2 = None, \
        month_values3 = None, labels = None, max_pval = 0.05):

    in_ax = True
    if(axs is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure(figsize = (9, 5))
        axs = fig.subplots(nrows = 2, ncols = 3)
 
        ax1 = axs[0,0]
        ax2 = axs[0,1]
        ax3 = axs[0,2]
        ax4 = axs[1,0]
        ax5 = axs[1,1]
        ax6 = axs[1,2]
    else:
        ax1 = axs[0]   
        ax2 = axs[1]   
        ax3 = axs[2]   
        ax4 = axs[3]   
        ax5 = axs[4]   
        ax6 = axs[5]   
 
    #ax1 = fig.add_subplot(2,3,1)
    #ax2 = fig.add_subplot(2,3,2)
    #ax3 = fig.add_subplot(2,3,3)
    #ax4 = fig.add_subplot(2,3,4)
    #ax5 = fig.add_subplot(2,3,5)
    #ax6 = fig.add_subplot(2,3,6)
  
    xvals = np.arange(2005, 2021)
    lwidth = 0.5
    fsize = 8

    months = ['April','May','June','July','August','September']
    axs = [ax1, ax2, ax3, ax4, ax5, ax6]

    def process_monthly_data(month_values, lats, axs, data_type = None, \
            linestyle = '-', color = 'tab:blue', alpha = 1.0, minlat = 70., \
            maxlat = 87., in_ax = False):

        # Mask the month values that are outside the desired
        # lat bounds
        # --------------------------------------------------
        keep_idxs = np.where((lats[:,0] >= (minlat + 0.5)) & \
            (lats[:,0] < (maxlat + 0.5)))[0]

        print("KEEP LATS = ", lats[keep_idxs,0])

        # Calculate the monthly averages of the estimated forcings
        # over the entire Arctic region
        # --------------------------------------------------------
        #arctic_avgs = np.array([np.nanmean(month_values[idx::6,:,:], \
        arctic_avgs = np.array([np.nanmean(month_values[idx::6,keep_idxs,:], \
            axis = (1,2)) for idx in range(6)])

        plot_max = np.max(arctic_avgs) + 0.05 
        plot_min = np.min(arctic_avgs) - 0.05 
        ylims = [plot_min, plot_max]

        trend_vals = np.full(len(axs), np.nan)
    
        for ii, ax in enumerate(axs):
            print(months[ii])
            ax.axhline(0, linestyle = '--', color = 'k', alpha = 0.5)
            ax.plot(xvals, arctic_avgs[ii], linestyle = linestyle, \
                color = color, alpha = alpha)
            ax.grid(alpha = 0.25, color = 'grey')
            if(not in_ax):
                ax.set_title(months[ii])   
                ax.set_ylabel("Avg. Forcing [W/m2]")
                ax.set_ylim(plot_min, plot_max)
            #if(data_type is None):
            rvals = plot_trend_line(ax, xvals, arctic_avgs[ii], color='black', \
                linestyle = '-', linewidth = lwidth, slope = trend_type)
            calc_trend = rvals.slope * len(xvals)
            ptext = 'ΔF$_{aer}$ = '+str(np.round(calc_trend, 3))
            if(version4):
                text_loc = 'lower_left'
                #ax.set_title(ptext, fontsize = 11)
                if(rvals.pvalue <= max_pval):
                    plot_figure_text(ax, ptext, xval = 2005, yval = -0.35, \
                        fontsize = fsize, color = color, weight = 'bold', backgroundcolor = 'white')
                else:
                    plot_figure_text(ax, ptext, xval = 2005, yval = -0.35, \
                        fontsize = fsize, color = color, backgroundcolor = 'white')
            else:
                text_loc = 'upper_left'
                if(rvals.pvalue <= max_pval):
                    plot_figure_text(ax, ptext, location = text_loc, \
                        fontsize = fsize, color = color, weight = 'bold', backgroundcolor = 'white')
                else:
                    plot_figure_text(ax, ptext, location = text_loc, \
                        fontsize = fsize, color = color, backgroundcolor = 'white')

            trend_vals[ii] = calc_trend
   
        return trend_vals, ylims
 
    file_add = ''
    title_add = ''
    orig_trends, orig_lims = process_monthly_data(all_month_values, \
        OMI_monthly_data['LAT'], axs, minlat = minlat, maxlat = maxlat, \
        in_ax = in_ax)
    total_lims = orig_lims
    if(month_values2 is not None):
        v2_trends, v2_lims = process_monthly_data(month_values2, \
            OMI_monthly_data['LAT'], axs, color = 'tab:orange', \
            minlat = minlat, maxlat = maxlat, \
            data_type = labels[0])
        file_add += '_' + labels[0]
        title_add += '\nOrange - ' + labels[0]
        total_lims = total_lims + v2_lims
    if(month_values3 is not None):
        v3_trends, v3_lims = process_monthly_data(month_values3, \
            OMI_monthly_data['LAT'], axs, color = 'tab:green', \
            minlat = minlat, maxlat = maxlat, \
            data_type = labels[1])
        file_add += labels[1]
        title_add += ', Green - ' + labels[1]
        total_lims = total_lims + v3_lims
    
    print(total_lims)
    total_min = np.min(total_lims)
    total_max = np.max(total_lims)
    for ii, ax in enumerate(axs):
        ax.set_ylim(total_min, total_max) 

    # Calculate and print the percent error of the trends calculated from
    # the adderror and suberror simulations
    # -------------------------------------------------------------------
    if(month_values2 is not None):
        for ii in range(orig_trends.shape[0]):
            pcnt_err1 = ((v2_trends[ii] - orig_trends[ii]) / orig_trends[ii]) * 100.
            pcnt_err2 = ((v3_trends[ii] - orig_trends[ii]) / orig_trends[ii]) * 100.
           
            print('{0:>6.2e} {1:>6.2e} {2:>6.2e} {3:>6.1f} {4:>6.1f}'.format(\
                orig_trends[ii], v2_trends[ii], v3_trends[ii], pcnt_err1, \
                pcnt_err2))
 
    #plt.suptitle('Forcing Estimate 3: Daily-averaged Single Month Based' + \
        #'\nBold - Significant at p = ' + str(max_pval) + title_add, \
    if(not in_ax):
        plt.suptitle('Arctic Regional Average Aerosol Forcing: ' + str(minlat) + \
            '$^{o}$ N - ' + str(maxlat) + '$^{o}$ N\n' + \
            'Bold - Significant at p = ' + str(max_pval) + title_add, \
            weight = 'bold') 

        plot_subplot_label(ax1, 'a)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax2, 'b)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax3, 'c)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax4, 'd)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax5, 'e)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax6, 'f)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)

        fig.tight_layout()

        if(save):
            if(version4):
                version_add = 'v4'
            else:
                version_add = 'v3'
            outname = 'calc_arctic_forcing_' + version_add + \
                '_monthly_arcticavg_' + trend_type + \
                file_add + '_' + omi_data_type + '_minlat' + str(int(minlat)) + \
                '_maxlat' + str(int(maxlat)) + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else: 
            plt.show()

    return total_min, total_max

# omi_data_type: 'raw' or 'pert' NOT NEEDED ANYMORE
# Takes a list of input HDF5 files and plots all the results
# Assumes that the first argument is the base simulation
# -----------------------------------------------------------
def plot_type_forcing_v4_all_months_arctic_avg(all_month_files, \
        minlat = 70., maxlat = 87., axs = None, \
        data_type = 'raw', trend_type = 'standard', \
        save = False,debug = False, labels = None, max_pval = 0.05,\
        linestyle_list = None, linecolor_list = None):

    in_ax = True
    if(axs is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure(figsize = (9, 5))
        axs = fig.subplots(nrows = 2, ncols = 3)
 
        ax1 = axs[0,0]
        ax2 = axs[0,1]
        ax3 = axs[0,2]
        ax4 = axs[1,0]
        ax5 = axs[1,1]
        ax6 = axs[1,2]
    else:
        ax1 = axs[0]   
        ax2 = axs[1]   
        ax3 = axs[2]   
        ax4 = axs[3]   
        ax5 = axs[4]   
        ax6 = axs[5]   
  
    xvals = np.arange(2005, 2021)
    lwidth = 0.5
    fsize = 10

    months = ['April','May','June','July','August','September']
    axs = [ax1, ax2, ax3, ax4, ax5, ax6]

    def process_monthly_data(infile, axs, data_type = None, \
            linestyle = '-', color = 'tab:blue', alpha = 1.0, linewidth = None, \
            minlat = 70., maxlat = 87., in_ax = False):

        data_dict = read_daily_month_force_HDF5(infile)
        month_values = data_dict['FORCE_EST']
        lats = data_dict['LAT']

        # Mask the month values that are outside the desired
        # lat bounds
        # --------------------------------------------------
        keep_idxs = np.where((lats[:] >= (minlat + 0.5)) & \
            (lats[:] < (maxlat + 0.5)))[0]

        #print("KEEP LATS = ", lats[keep_idxs])

        # Calculate the monthly averages of the estimated forcings
        # over the entire Arctic region
        # --------------------------------------------------------
        #arctic_avgs = np.array([np.nanmean(month_values[idx::6,:,:], \
        arctic_avgs = np.array([np.nanmean(month_values[idx::6,keep_idxs,:], \
            axis = (1,2)) for idx in range(6)])


        plot_max = np.max(arctic_avgs) + 0.01 
        plot_min = np.min(arctic_avgs) - 0.01 
        ylims = [plot_min, plot_max]

        trend_vals = np.full(len(axs), np.nan)
    
        for ii, ax in enumerate(axs):
            print(months[ii])
            ax.axhline(0, linestyle = '--', color = 'tab:grey', alpha = 0.5)
            ax.grid(alpha = 0.25, color = 'grey')
            ax.plot(xvals, arctic_avgs[ii], linestyle = linestyle, \
                color = color, alpha = alpha, linewidth = linewidth)
            if(not in_ax):
                ax.set_title(months[ii])   
                ax.set_ylabel("Avg. Forcing [W/m2]")
                ax.set_ylim(plot_min, plot_max)


            # First, calculate the trend
            zdata = stats.linregress(xvals,  arctic_avgs[ii])
            calc_trend = len(xvals) * zdata.slope
            #zdata = np.polyfit(xdata, ydata, 1)

            #if(data_type is None):
            #rvals = plot_trend_line(ax, xvals, arctic_avgs[ii], color='black', \
            #    linestyle = '-', linewidth = linewidth, slope = trend_type)
            #calc_trend = rvals.slope * len(xvals)
            ptext = 'ΔF$_{aer}$ = '+str(np.round(calc_trend, 3))

            #if(version4):
            #    text_loc = 'lower_left'
            #else:
            #    text_loc = 'upper_left'
            #if(rvals.pvalue <= max_pval):
            #    plot_figure_text(ax, ptext, location = text_loc, \
            #        fontsize = fsize, color = color, weight = 'bold', backgroundcolor = 'white')
            #else:
            #    plot_figure_text(ax, ptext, location = text_loc, \
            #        fontsize = fsize, color = color, backgroundcolor = 'white')

            trend_vals[ii] = calc_trend
   
        return trend_vals, ylims

    for ii, infile in enumerate(all_month_files):

        if(ii == 0):
            pcolor = 'k'
            linewidth = 2
            alpha = 1
        else:
            pcolor = None
            linewidth = 1
            alpha = 1.00

        if(linestyle_list is not None):
            linestyle = linestyle_list[ii]
            linecolor = linecolor_list[ii]
            trend_vals, ylims = process_monthly_data(infile, axs, \
                color = linecolor, linestyle = linestyle, \
                linewidth = linewidth, alpha = alpha, minlat = minlat, \
                maxlat = maxlat, in_ax = in_ax)
        else:
            trend_vals, ylims = process_monthly_data(infile, axs, \
                color = pcolor, \
                linewidth = linewidth, alpha = alpha, minlat = minlat, \
                maxlat = maxlat, in_ax = in_ax)

        if(ii == 0):
            total_lims = ylims
        else:
            total_lims = total_lims + ylims
     
    # Replot the first line
    # ---------------------    
    _, _  = process_monthly_data(all_month_files[0], axs, color = 'k', \
        linewidth = 1.5, alpha = 1.0, minlat = minlat, maxlat = maxlat, in_ax = in_ax)
 
    #file_add = ''
    title_add = ''
    #orig_trends, orig_lims = process_monthly_data(all_month_values, axs)
    #total_lims = orig_lims
    #if(month_values2 is not None):
    #    v2_trends, v2_lims = process_monthly_data(month_values2, axs, color = 'tab:orange', data_type = labels[0])
    #    file_add += '_' + labels[0]
    #    title_add += '\nOrange - ' + labels[0]
    #    total_lims = total_lims + v2_lims
    #if(month_values3 is not None):
    #    v3_trends, v3_lims = process_monthly_data(month_values3, axs, color = 'tab:green', data_type = labels[1])
    #    file_add += labels[1]
    #    title_add += ', Green - ' + labels[1]
    #    total_lims = total_lims + v3_lims
    
    print(total_lims)
    total_min = np.min(total_lims)
    total_max = np.max(total_lims)
    for ii, ax in enumerate(axs):
        ax.set_ylim(total_min, total_max) 

    ##!## Calculate and print the percent error of the trends calculated from
    ##!## the adderror and suberror simulations
    ##!## -------------------------------------------------------------------
    ##!#if(month_values2 is not None):
    ##!#    for ii in range(orig_trends.shape[0]):
    ##!#        pcnt_err1 = ((v2_trends[ii] - orig_trends[ii]) / orig_trends[ii]) * 100.
    ##!#        pcnt_err2 = ((v3_trends[ii] - orig_trends[ii]) / orig_trends[ii]) * 100.
    ##!#       
    ##!#        print('{0:>6.2e} {1:>6.2e} {2:>6.2e} {3:>6.1f} {4:>6.1f}'.format(\
    ##!#            orig_trends[ii], v2_trends[ii], v3_trends[ii], pcnt_err1, \
    ##!#            pcnt_err2))
 
    #plt.suptitle('Forcing Estimate 3: Daily-averaged Single Month Based' + \
        #'\nBold - Significant at p = ' + str(max_pval) + title_add, \
    #plt.suptitle('Observational Arctic Aerosol Forcing Estimate\n' + \
    #    'Bold - Significant at p = ' + str(max_pval) + title_add, \
    #    weight = 'bold') 
    if(not in_ax):
        plt.suptitle('Arctic Regional Average Aerosol Forcing: ' + str(minlat) + \
            '$^{o}$ N - ' + str(maxlat) + '$^{o}$ N\n' + \
            'Bold - Significant at p = ' + str(max_pval) + title_add, \
            weight = 'bold') 

        plot_subplot_label(ax1, 'a)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax2, 'b)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax3, 'c)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax4, 'd)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax5, 'e)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax6, 'f)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)

        fig.tight_layout()

        if(save):
            version_add = 'v4'
            outname = 'calc_arctic_forcing_uncert_' + version_add + \
                '_monthly_arcticavg_' + trend_type + \
                file_add + '_' + omi_data_type + str(int(minlat)) + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else: 
            plt.show()

    return total_min, total_max

def plot_type_forcing_v4_all_months_arctic_avg_combined(all_month_files, \
        version4 = False, max_pval = 0.05, \
        slope_type = 'lin', \
        horiz_orient = True, save = False):

    # Prep the overall figure
    # -----------------------
    plt.close('all')
    if(horiz_orient):
        fig = plt.figure(figsize = (13, 6.5))
        axs  = fig.subplots(nrows = 3, ncols = 6, sharex = True, sharey = True)
        ax1  = axs[0,0]  # control April 
        ax2  = axs[0,1]  # control May 
        ax3  = axs[0,2]  # control June 
        ax4  = axs[0,3]  # control July 
        ax5  = axs[0,4]  # control August
        ax6  = axs[0,5]  # control September
        ax7  = axs[1,0]  # low     April 
        ax8  = axs[1,1]  # low     May 
        ax9  = axs[1,2]  # low     June 
        ax10 = axs[1,3]  # low     July 
        ax11 = axs[1,4]  # low     August
        ax12 = axs[1,5]  # low     September
        ax13 = axs[2,0]  # high    April 
        ax14 = axs[2,1]  # high    May 
        ax15 = axs[2,2]  # high    June 
        ax16 = axs[2,3]  # high    July 
        ax17 = axs[2,4]  # high    August
        ax18 = axs[2,5]  # high    September
        orient_add = '_horiz'
    else:
        #fig = plt.figure(figsize = (7.5, 11))
        fig = plt.figure(figsize = (8.5, 9.5))
        axs  = fig.subplots(nrows = 6, ncols = 3, sharex = True, sharey = True)
        ax1  = axs[0,0]  # control April 
        ax2  = axs[1,0]  # control May 
        ax3  = axs[2,0]  # control June 
        ax4  = axs[3,0]  # control July 
        ax5  = axs[4,0]  # control August
        ax6  = axs[5,0]  # control September
        ax7  = axs[0,1]  # low     April 
        ax8  = axs[1,1]  # low     May 
        ax9  = axs[2,1]  # low     June 
        ax10 = axs[3,1]  # low     July 
        ax11 = axs[4,1]  # low     August
        ax12 = axs[5,1]  # low     September
        ax13 = axs[0,2]  # high    April 
        ax14 = axs[1,2]  # high    May 
        ax15 = axs[2,2]  # high    June 
        ax16 = axs[3,2]  # high    July 
        ax17 = axs[4,2]  # high    August
        ax18 = axs[5,2]  # high    September
        orient_add = '_vert'


    linestyle_dict = {
        '': '-', 
        'addslopeerror': '-', 
        'subslopeerror': '--', 
        'addintcpterror': '-', 
        'subintcpterror': '--', 
        'iceplus5': '-', 
        'iceminus5': '--', 
        'iceplus15': '-', 
        'iceminus15': '--', 
        'codplus2': '-', 
        'codminus2': '--', 
        'codplus5': '-', 
        'codminus5': '--', 
        'addL2L3error30': '-', 
        'subL2L3error30': '--', 
    }

    linecolor_dict = {
        '': 'k', 
        'addslopeerror': 'tab:blue', 
        'subslopeerror': 'tab:blue', 
        'addintcpterror': 'tab:orange', 
        'subintcpterror': 'tab:orange', 
        'iceplus5': 'tab:olive', 
        'iceminus5': 'tab:olive', 
        'iceplus15': 'tab:red', 
        'iceminus15': 'tab:red', 
        'codplus2': 'tab:purple', 
        'codminus2': 'tab:purple', 
        'codplus5': 'tab:cyan', 
        'codminus5': 'tab:cyan', 
        'addL2L3error30': 'tab:brown', 
        'subL2L3error30': 'tab:brown', 
    }

    linelabel_dict = {
        '': 'k', 
        'addslopeerror': '+ slope error', 
        'subslopeerror': '-  slope error', 
        'addintcpterror': '+ intcpt error', 
        'subintcpterror': '-  intcpt error', 
        'iceplus5': 'ice conc. + 5%', 
        'iceminus5': 'ice conc.  -  5%', 
        'iceplus15': 'ice conc. + 15%', 
        'iceminus15': 'ice conc.  -  15%', 
        'codplus2': 'COD + 2', 
        'codminus2': 'COD  -  2', 
        'codplus5': 'COD + 5', 
        'codminus5': 'COD  -  5', 
        'addL2L3error30': '+ 30 W/m2', 
        'subL2L3error30': '- 30 W/m2', 
    }

    file_types = [ptitle.strip().split('/')[-1].split('_')[-1][3:].split('.')[0] \
        for ptitle in all_month_files]

    linestyle_list = [linestyle_dict[ftype] for ftype in file_types]
    linecolor_list = [linecolor_dict[ftype] for ftype in file_types]

    print(linestyle_list, linecolor_list)

    # Plot the first row: control (65 - 87)
    # -------------------------------------
    ax_list = [ax1, ax2, ax3, ax4, ax5, ax6]
    ##!#first_min, first_max = \
    ##!#    plot_type_forcing_v3_all_months_arctic_avg(all_month_values, \
    ##!#    OMI_monthly_data, \
    ##!#    axs = ax_list, \
    ##!#    minlat = 65., maxlat = 87., omi_data_type = 'lin', 
    ##!#    version4 = True, 
    ##!#    save = False, debug = False,max_pval = 0.05)

    first_min, first_max = \
        plot_type_forcing_v4_all_months_arctic_avg(all_month_files, \
        axs = ax_list, \
        minlat = 65., maxlat = 87., \
        data_type = 'raw', trend_type = 'standard', \
        save = False,debug = False, labels = None,\
        linestyle_list = linestyle_list, \
        linecolor_list = linecolor_list)


    # Plot the second row: low (65 - 75)
    # ----------------------------------
    ax_list = [ax7, ax8, ax9, ax10, ax11, ax12]
    ##!#secnd_min, secnd_max = \
    ##!#    plot_type_forcing_v3_all_months_arctic_avg(all_month_values, \
    ##!#    OMI_monthly_data, \
    ##!#    axs = ax_list, \
    ##!#    minlat = 65., maxlat = 75., omi_data_type = 'lin', 
    ##!#    version4 = True, 
    ##!#    save = False, debug = False,max_pval = 0.05)

    secnd_min, secnd_max = \
        plot_type_forcing_v4_all_months_arctic_avg(all_month_files, \
        axs = ax_list, \
        minlat = 65., maxlat = 75., \
        data_type = 'raw', trend_type = 'standard', \
        save = False,debug = False, labels = None, \
        linestyle_list = linestyle_list, \
        linecolor_list = linecolor_list)

    # Plot the third row: high (75 - 87)
    # ----------------------------------
    ax_list = [ax13, ax14, ax15, ax16, ax17, ax18]
    ##!#third_min, third_max = \
    ##!#plot_type_forcing_v3_all_months_arctic_avg(all_month_values, \
    ##!#    OMI_monthly_data, \
    ##!#    axs = ax_list, \
    ##!#    minlat = 75., maxlat = 87., omi_data_type = 'lin', 
    ##!#    version4 = True, 
    ##!#    save = False, debug = False,max_pval = 0.05)

    third_min, third_max = \
        plot_type_forcing_v4_all_months_arctic_avg(all_month_files, \
        axs = ax_list, \
        minlat = 75., maxlat = 87., \
        data_type = 'raw', trend_type = 'standard', \
        save = False,debug = False, labels = None, \
        linestyle_list = linestyle_list, \
        linecolor_list = linecolor_list)


    if(horiz_orient):
        ax1.set_ylabel('Forcing [Wm$^{-2}$]')
        ax7.set_ylabel('Forcing [Wm$^{-2}$]')
        ax13.set_ylabel('Forcing [Wm$^{-2}$]')
    else:
        ax1.set_ylabel('Forcing [Wm$^{-2}$]')
        ax2.set_ylabel('Forcing [Wm$^{-2}$]')
        ax3.set_ylabel('Forcing [Wm$^{-2}$]')
        ax4.set_ylabel('Forcing [Wm$^{-2}$]')
        ax5.set_ylabel('Forcing [Wm$^{-2}$]')
        ax6.set_ylabel('Forcing [Wm$^{-2}$]')

        ax1.set_title('65$^{o}$ - 87$^{o}$ N', weight = 'bold')
        ax7.set_title('65$^{o}$ - 75$^{o}$ N', weight = 'bold')
        ax13.set_title('75$^{o}$ - 87$^{o}$ N', weight = 'bold')

        fig.subplots_adjust(left = 0.2)

        row_label_size = 12 
        fig.text(0.08, 0.82, 'April', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.08, 0.69, 'May', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.08, 0.555, 'June', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.08, 0.42, 'July', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.08, 0.295, 'August', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.08, 0.16, 'September', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)

        linestyle_list = [linestyle_dict[ftype] for ftype in file_types]
        linecolor_list = [linecolor_dict[ftype] for ftype in file_types]
        linelabel_list = [linelabel_dict[ftype] for ftype in file_types][1:]

        # Extract all the file types here
        # -------------------------------
        custom_lines = [Line2D([0], [0], color = lcolor, linestyle = lstyle) \
            for lcolor, lstyle in zip(linecolor_list[1:], linestyle_list[1:])]
        
        #custom_lines = [Line2D([0], [0], color='k'),
        #                Line2D([0], [0], color='k', linestyle = '--'),\
        #                Line2D([0], [0], color='k', linestyle = ':')]
        #ax10.legend(custom_lines, ['0.64 μm', '2.25 μm', '10.35 μm'],\
        #    fontsize = 8, loc = 2)

        ##!#norm = mc.Normalize(vmin=2005, vmax = 2020)
        ##!#bounds = np.arange(2005 - 0.5, 2020 + 1.5) 
        ##!#ticks = np.arange(2005, 2020 + 1) 
        ##!#cmap = 'turbo'

        ##!#cbar_ax = fig.add_axes([0.20, 0.07, 0.7, 0.01])
        ##!#if(stype == 'ice'):
        ##!#    cbar_adder = 'sea ice'
        ##!#elif(stype == 'cld'):
        ##!#    cbar_adder = 'cloud optical depth'
        #cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),
        #             cax=cbar_ax, orientation='horizontal', ticks = ticks[::3], \
        #             boundaries = bounds, label='Year of ' + cbar_adder + ' used in forcing calculation')

        fig.legend(custom_lines, linelabel_list, loc = 'lower center', \
            bbox_to_anchor = (0, 0.01, 1, 1),\
            bbox_transform = plt.gcf().transFigure, ncol=6, fontsize = 9)

    all_ranges = [first_min, first_max, secnd_min, secnd_max, \
        third_min, third_max]

    plt.suptitle('Arctic Regional Average Aerosol Forcing\nUncertainty Analysis', \
        weight = 'bold') 

    for pax in axs.flatten():
        pax.set_ylim(np.min(all_ranges), np.max(all_ranges))

    #fig.tight_layout()

    if(save):
        file_add = ''
        trend_add = ''
        #outname = 'calc_arctic_forcing_v4_monthly_arcticavg_aimin0_useintcpt_ref' + stype + \
        #    '_' + trend_type + file_add + trend_add + '_' + slope_type + \
        #    orient_add + '_combined.png'

        outname = 'calc_arctic_forcing_uncert_v4_' + \
            'monthly_arcticavg_' + slope_type + orient_add + '_combined.png'

        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()





# Used for plotting all the different ref-ice simulations
# ptype: 'forcing', 'error', 'pcnt_error'
# vtype: '' or 'v3'
# stype: 'ice' or 'cld'. Allows users to either plot 'refice' simulations
#        or 'refcld' simulations.
# return_slopes: returns the calculated slopes from each of the yearly
#       simulations
# slope_type: 'lin' or 'thl'
def plot_type_forcing_v3_all_months_arctic_avg_manyrefice(all_month_vals, \
        OMI_monthly_data, axs = None, \
        min_year = 2005, max_year = 2020, \
        minlat = 70., maxlat = 87., use_szas = False, \
        ai_thresh = 0.7, cld_idx = 0, maxerr = 2, \
        min_cloud = 0.95, data_type = 'raw', trend_type = 'standard', \
        save = False,debug = False, max_pval = 0.05, \
        ptype = 'forcing', vtype = '', \
        slope_type = 'lin', \
        stype = 'ice', show_trends = False, \
        version4 = False):

    if(stype == 'ice'):
        title_add = 'Sea Ice Retreat'
    elif(stype == 'cld'):
        title_add = 'Arctic Cloud Cover'
    else:
        print("ERROR: INVALID STYPE")
        return 

    in_ax = True
    if(axs is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure(figsize = (9, 5))
        axs = fig.subplots(nrows = 2, ncols = 3)
 
        ax1 = axs[0,0]
        ax2 = axs[0,1]
        ax3 = axs[0,2]
        ax4 = axs[1,0]
        ax5 = axs[1,1]
        ax6 = axs[1,2]
    else:
        ax1 = axs[0]   
        ax2 = axs[1]   
        ax3 = axs[2]   
        ax4 = axs[3]   
        ax5 = axs[4]   
        ax6 = axs[5]   
  


 
    years = np.arange(min_year, max_year + 1)
    lwidth = 0.5
    fsize = 10

    months = ['April','May','June','July','August','September']
    axs = [ax1, ax2, ax3, ax4, ax5, ax6]

    plot_c = turbo((years-np.min(years))/\
                     (np.max(years)-np.min(years)))

    if((ptype == 'error') | (ptype == 'pcnt_error')):
        ax1.axhline(0, linestyle = '--', color = 'black')
        ax2.axhline(0, linestyle = '--', color = 'black')
        ax3.axhline(0, linestyle = '--', color = 'black')
        ax4.axhline(0, linestyle = '--', color = 'black')
        ax5.axhline(0, linestyle = '--', color = 'black')
        ax6.axhline(0, linestyle = '--', color = 'black')
    
    # Mask the month values that are outside the desired
    # lat bounds
    # --------------------------------------------------
    keep_idxs = np.where((OMI_monthly_data['LAT'][:,0] >= (minlat + 0.5)) & \
        (OMI_monthly_data['LAT'][:,0] < (maxlat + 0.5)))[0]

    # Convert the base values to Arctic and monthly averages
    # ------------------------------------------------------
    all_month_vals =  np.array([np.nanmean(all_month_vals[idx::6,keep_idxs,:], \
        axis = (1,2)) for idx in range(6)])

    # Combined all the single-year stuff into one array
    # ------------------------------------------------- 
    combined_vals = np.full( (years.shape[0], \
        6, all_month_vals.shape[1]), np.nan)

    calc_slopes = np.full( (years.shape[0], 6), np.nan)

    if(vtype != ''):
        vtype = vtype + '_'

    # Loop over each reference ice year, reading in the 
    # individual forcing values.
    for ii, year in enumerate(years):
        # Read in the force vals with this ref ice
        # ----------------------------------------
        if(version4):
            if(slope_type == 'thl'):
                filename = home_dir + '/Research/Arctic_compares/' + \
                    'arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_ref' + \
                    stype + str(year) + '.hdf5'
            else:
                filename = home_dir + '/Research/Arctic_compares/' + \
                    'arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_lin_ref' + \
                    stype + str(year) + '.hdf5'
        else:
            filename = home_dir + '/Research/Arctic_compares/' + \
                'arctic_month_est_forcing_dayaithresh07_' + vtype + 'ref' + \
                stype + str(year) + '.hdf5'

        print(filename)
        all_month_dict = read_daily_month_force_HDF5(filename)
        month_values = all_month_dict['FORCE_EST']
        lats = all_month_dict['LAT']

        # Mask the month values that are outside the desired
        # lat bounds
        # --------------------------------------------------
        keep_idxs = np.where((lats[:] >= (minlat + 0.5)) & \
            (lats[:] < (maxlat + 0.5)))[0]

        # Calculate the monthly averages of the estimated forcings
        # over the entire Arctic region
        # --------------------------------------------------------
        #arctic_avgs = np.array([np.nanmean(month_values[idx::6,:,:], \
        combined_vals[ii,:,:] = np.array([np.nanmean(month_values[idx::6,keep_idxs,:], \
            axis = (1,2)) for idx in range(6)])

        if(ptype == 'error'):
            combined_vals[ii,:,:] = (combined_vals[ii,:,:] - \
                all_month_vals[:,:])
        elif(ptype == 'pcnt_error'):
            combined_vals[ii,:,:] = ((combined_vals[ii,:,:] - all_month_vals[:,:]) / \
                all_month_vals[:,:]) * 100.
   
        # Now, loop over each month and plot accordingly
        # ----------------------------------------------
        for jj, ax in enumerate(axs):
            if(show_trends):
                rvals = plot_trend_line(ax, years, combined_vals[ii,jj,:], color=plot_c[ii], \
                    linestyle = '-', linewidth = lwidth, slope = trend_type)
            ax.plot(years, combined_vals[ii,jj,:], linestyle = '-', \
                c = plot_c[ii], alpha = 1.0)



    local_min = np.min(combined_vals)
    local_max = np.max(combined_vals)
    ranges = np.abs(np.array([local_min, local_max]))
    local_min = local_min - 0.1 * np.max(ranges)
    local_max = local_max + 0.1 * np.max(ranges)

    if(ptype == 'error'):
        plabel = 'Avg. Forcing\nError [W/m2]'
        file_add = '_error'
    elif(ptype == 'pcnt_error'):
        plabel = 'Avg. Forcing\nPcnt. Error [%]'
        file_add = '_pcnterror'
    else:
        plabel = 'Avg. Forcing [W/m2]'
        file_add = '_forcing'

    if(show_trends):
        trend_add = '_trends'
    else:
        trend_add = ''

    for jj in range(6):
  
        axs[jj].axhline(0, color = 'black', linestyle = '--', alpha = 0.75)
        axs[jj].grid(alpha = 0.25, color = 'grey')
        if(ptype == 'forcing'):
            axs[jj].plot(years, all_month_vals[jj,:], linestyle = '-', \
                c = 'black', alpha = 1.0)
        if(not in_ax):
            axs[jj].set_ylabel(plabel) 
            axs[jj].set_title(months[jj])   
            axs[jj].set_ylim(local_min, local_max)

    if(not in_ax):
        plot_subplot_label(ax1, 'a)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax2, 'b)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax3, 'c)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax4, 'd)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax5, 'e)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)
        plot_subplot_label(ax6, 'f)', fontsize = 12, location = 'upper_upper_left', backgroundcolor = None)

        plt.suptitle('Arctic Regional Average Aerosol Forcing: ' + str(minlat) + \
            '$^{o}$ N - ' + str(maxlat) + '$^{o}$ N\n' + \
            'Sensitivity to ' + title_add,
            weight = 'bold') 

        fig.tight_layout()

        if(save):
            if(version4):
                outname = 'calc_arctic_forcing_v4_monthly_arcticavg_aimin0_useintcpt_ref' + stype + \
                    '_' + trend_type + file_add + trend_add + '_' + slope_type + \
                '_minlat' + str(int(minlat)) + \
                '_maxlat' + str(int(maxlat)) + '.png'
            else:
                outname = 'calc_arctic_forcing_v3_monthly_arcticavg_ref' + stype + \
                    '_' + vtype + trend_type + file_add + trend_add + '_' + slope_type + \
                    '_minlat' + str(int(minlat)) + \
                '_maxlat' + str(int(maxlat)) + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else: 
            plt.show()

    return local_min, local_max

def plot_type_forcing_v4_all_months_arctic_avg_manyrefice_combined(all_month_vals, \
        OMI_monthly_data, \
        min_year = 2005, max_year = 2020, \
        minlat = 70., maxlat = 87., trend_type = 'standard', \
        save = False,debug = False, max_pval = 0.05, \
        ptype = 'forcing', vtype = '', \
        slope_type = 'lin', horiz_orient = True, \
        stype = 'ice', show_trends = False, \
        version4 = False):

    if(stype == 'ice'):
        title_add = 'Sea Ice Retreat'
    elif(stype == 'cld'):
        title_add = 'Arctic Cloud Optical Depth'
    else:
        print("ERROR: INVALID STYPE")
        return 

    # Prep the overall figure
    # -----------------------
    plt.close('all')
    if(horiz_orient):
        fig = plt.figure(figsize = (12, 6.5))
        #axs  = fig.subplots(nrows = 3, ncols = 6, sharex = True, sharey = True)
        ##ax1  = axs[0,0]  # control April 
        ##ax2  = axs[0,1]  # control May 
        ##ax3  = axs[0,2]  # control June 
        ##ax4  = axs[0,3]  # control July 
        ##ax5  = axs[0,4]  # control August
        ##ax6  = axs[0,5]  # control September
        ##ax7  = axs[1,0]  # low     April 
        ##ax8  = axs[1,1]  # low     May 
        ##ax9  = axs[1,2]  # low     June 
        ##ax10 = axs[1,3]  # low     July 
        ##ax11 = axs[1,4]  # low     August
        ##ax12 = axs[1,5]  # low     September
        ##ax13 = axs[2,0]  # high    April 
        ##ax14 = axs[2,1]  # high    May 
        ##ax15 = axs[2,2]  # high    June 
        ##ax16 = axs[2,3]  # high    July 
        ##ax17 = axs[2,4]  # high    August
        ##ax18 = axs[2,5]  # high    September

        grid = ImageGrid(fig, 111, \
            nrows_ncols = (3, 6), \
            axes_pad = 0.20, \
            share_all = True, \
            cbar_location = 'right', \
            cbar_mode = 'edge', \
            cbar_size = '7%', \
            cbar_pad = 0.15)

        ax1  = grid[0,0]  # control April 
        ax2  = grid[0,1]  # control May 
        ax3  = grid[0,2]  # control June 
        ax4  = grid[0,3]  # control July 
        ax5  = grid[0,4]  # control August
        ax6  = grid[0,5]  # control September
        ax7  = grid[1,0]  # low     April 
        ax8  = grid[1,1]  # low     May 
        ax9  = grid[1,2]  # low     June 
        ax10 = grid[1,3]  # low     July 
        ax11 = grid[1,4]  # low     August
        ax12 = grid[1,5]  # low     September
        ax13 = grid[2,0]  # high    April 
        ax14 = grid[2,1]  # high    May 
        ax15 = grid[2,2]  # high    June 
        ax16 = grid[2,3]  # high    July 
        ax17 = grid[2,4]  # high    August
        ax18 = grid[2,5]  # high    September
        orient_add = '_horiz'
    else:
        fig = plt.figure(figsize = (8.5, 9.5))
        axs  = fig.subplots(nrows = 6, ncols = 3, sharex = True, sharey = True)
        ax1  = axs[0,0]  # control April 
        ax2  = axs[1,0]  # control May 
        ax3  = axs[2,0]  # control June 
        ax4  = axs[3,0]  # control July 
        ax5  = axs[4,0]  # control August
        ax6  = axs[5,0]  # control September
        ax7  = axs[0,1]  # low     April 
        ax8  = axs[1,1]  # low     May 
        ax9  = axs[2,1]  # low     June 
        ax10 = axs[3,1]  # low     July 
        ax11 = axs[4,1]  # low     August
        ax12 = axs[5,1]  # low     September
        ax13 = axs[0,2]  # high    April 
        ax14 = axs[1,2]  # high    May 
        ax15 = axs[2,2]  # high    June 
        ax16 = axs[3,2]  # high    July 
        ax17 = axs[4,2]  # high    August
        ax18 = axs[5,2]  # high    September

        """
        grid = ImageGrid(fig, 111, \
            nrows_ncols = (6, 3), \
            axes_pad = 0.20, \
            share_all = True, \
            cbar_location = 'bottom', \
            cbar_mode = 'edge', \
            cbar_size = '7%', \
            cbar_pad = 0.15)

        ax1  = axs[0,0]  # control April 
        ax2  = axs[1,0]  # control May 
        ax3  = axs[2,0]  # control June 
        ax4  = axs[3,0]  # control July 
        ax5  = axs[4,0]  # control August
        ax6  = axs[5,0]  # control September
        ax7  = axs[0,1]  # low     April 
        ax8  = axs[1,1]  # low     May 
        ax9  = axs[2,1]  # low     June 
        ax10 = axs[3,1]  # low     July 
        ax11 = axs[4,1]  # low     August
        ax12 = axs[5,1]  # low     September
        ax13 = axs[0,2]  # high    April 
        ax14 = axs[1,2]  # high    May 
        ax15 = axs[2,2]  # high    June 
        ax16 = axs[3,2]  # high    July 
        ax17 = axs[4,2]  # high    August
        ax18 = axs[5,2]  # high    September
        """
        orient_add = '_vert'

    # Plot the first row: control (65 - 87)
    # -------------------------------------
    ax_list = [ax1, ax2, ax3, ax4, ax5, ax6]
    first_min, first_max = \
        plot_type_forcing_v3_all_months_arctic_avg_manyrefice(\
            all_month_vals, \
            OMI_monthly_data, axs = ax_list, \
            minlat = 65, maxlat = 87., trend_type = 'standard', stype = stype, \
            ptype = ptype, vtype = 'v4', slope_type = slope_type, version4 = True)


    # Plot the second row: low (65 - 75)
    # ----------------------------------
    ax_list = [ax7, ax8, ax9, ax10, ax11, ax12]
    secnd_min, secnd_max = \
        plot_type_forcing_v3_all_months_arctic_avg_manyrefice(\
            all_month_vals, \
            OMI_monthly_data, axs = ax_list, \
            minlat = 65, maxlat = 75., trend_type = 'standard', stype = stype, \
            ptype = ptype, vtype = 'v4', slope_type = slope_type, version4 = True)

    # Plot the third row: high (75 - 87)
    # ----------------------------------
    ax_list = [ax13, ax14, ax15, ax16, ax17, ax18]
    third_min, third_max = \
        plot_type_forcing_v3_all_months_arctic_avg_manyrefice(\
            all_month_vals, \
            OMI_monthly_data, axs = ax_list, \
            minlat = 75, maxlat = 87., trend_type = 'standard', stype = stype, \
            ptype = ptype, vtype = 'v4', slope_type = slope_type, version4 = True)

    all_ranges = [first_min, first_max, secnd_min, secnd_max, \
        third_min, third_max]

    for pax in axs.flatten():
        pax.set_ylim(np.min(all_ranges), np.max(all_ranges))


    if(horiz_orient):
        ax1.set_ylabel('Forcing [Wm$^{-2}$]')
        ax7.set_ylabel('Forcing [Wm$^{-2}$]')
        ax13.set_ylabel('Forcing [Wm$^{-2}$]')

        ax1.set_title('April', weight = 'bold')
        ax2.set_title('May', weight = 'bold')
        ax3.set_title('June', weight = 'bold')
        ax4.set_title('July', weight = 'bold')
        ax5.set_title('August', weight = 'bold')
        ax6.set_title('September', weight = 'bold')

        plot_figure_text(ax1, '65$^{o}$ - 87$^{o}$ N', \
            location = 'lower_left', xval = None, yval = None, \
            fontsize = 12, color = 'red', weight = 'bold', backgroundcolor = 'white')
        plot_figure_text(ax7, '65$^{o}$ - 75$^{o}$ N', \
            location = 'lower_left', xval = None, yval = None, \
            fontsize = 12, color = 'red', weight = 'bold', \
            backgroundcolor = 'white')
        plot_figure_text(ax13, '75$^{o}$ - 87$^{o}$ N', \
            location = 'lower_left', xval = None, yval = None, \
            fontsize = 12, color = 'red', weight = 'bold', \
            backgroundcolor = 'white')

    else:
        ax1.set_ylabel('Forcing [Wm$^{-2}$]')
        ax2.set_ylabel('Forcing [Wm$^{-2}$]')
        ax3.set_ylabel('Forcing [Wm$^{-2}$]')
        ax4.set_ylabel('Forcing [Wm$^{-2}$]')
        ax5.set_ylabel('Forcing [Wm$^{-2}$]')
        ax6.set_ylabel('Forcing [Wm$^{-2}$]')

        ax1.set_title('65$^{o}$ - 87$^{o}$ N', weight = 'bold')
        ax7.set_title('65$^{o}$ - 75$^{o}$ N', weight = 'bold')
        ax13.set_title('75$^{o}$ - 87$^{o}$ N', weight = 'bold')

        fig.subplots_adjust(left = 0.2)

        row_label_size = 12 
        fig.text(0.07, 0.82, 'April', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.69, 'May', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.555, 'June', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.42, 'July', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.295, 'August', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.07, 0.16, 'September', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)

        norm = mc.Normalize(vmin=2005, vmax = 2020)
        bounds = np.arange(2005 - 0.5, 2020 + 1.5) 
        ticks = np.arange(2005, 2020 + 1) 
        cmap = 'turbo'

        cbar_ax = fig.add_axes([0.20, 0.07, 0.7, 0.01])
        if(stype == 'ice'):
            cbar_adder = 'sea ice'
        elif(stype == 'cld'):
            cbar_adder = 'cloud optical depth'
        cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),
                     cax=cbar_ax, orientation='horizontal', ticks = ticks[::3], \
                     boundaries = bounds, label='Year of ' + cbar_adder + ' used in forcing calculation')


    plt.suptitle('Arctic Regional Average Aerosol Forcing\n' + 
        'Sensitivity to ' + title_add,  \
        weight = 'bold') 
    #plt.suptitle('Arctic Regional Average Aerosol Forcing\n' + 
    #    'Sensitivity to ' + title_add + \
    #    '\nTop row: 65 - 87\nMiddle row: 65 - 75\nBottom row: 75 - 87',
    #    weight = 'bold') 

    ##row_label_size = 10 
    ##fig.text(0.05, 0.70, 'Climatology', ha='center', va='center', \
    ##    rotation='vertical',weight='bold',fontsize=row_label_size)
    ##fig.text(0.05, 0.30, 'Trend', ha='center', va='center', \
    ##    rotation='vertical',weight='bold',fontsize=row_label_size)
    
    #fig.tight_layout()
    if(save):
        file_add = ''
        trend_add = ''
        outname = 'calc_arctic_forcing_v4_monthly_arcticavg_aimin0_useintcpt_ref' + stype + \
            '_' + trend_type + file_add + trend_add + '_' + slope_type + \
        orient_add + '_combined.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


# This function gathers the refice and refcld simulation results,
# as well as the uncertainty analyses,
# calculates the trends from each refice and refcld simulation,
# and calculates the mean and standard error of those
# ref ice/cld simulations. The results are printed in a table
def calc_print_forcing_slope_error_v4(all_month_vals, \
        OMI_monthly_data, all_month_files, \
        min_year = 2005, max_year = 2020, \
        minlat = 65., maxlat = 87., \
        trend_type = 'standard', ptype = 'forcing', \
        vtype = '', version4 = True, slope_type = 'lin', 
        save = False):
    
    ice_trend_dict = \
        calc_forcing_slopes_v4_all_months_arctic_avg_manyrefice(all_month_vals, \
        OMI_monthly_data, minlat = 65., maxlat = 87., stype = 'ice', slope_type = slope_type)

    cld_trend_dict = \
        calc_forcing_slopes_v4_all_months_arctic_avg_manyrefice(all_month_vals, \
        OMI_monthly_data, minlat = 65., maxlat = 87., stype = 'cld', slope_type = slope_type)

    uncert_trend_dict = calc_forcing_slopes_v4_all_months_arctic_avg_uncert(all_month_files, \
            OMI_monthly_data, slope_type = slope_type)

    years = np.arange(min_year, max_year + 1)


    # = = = = = = = = = = = =
    #
    # refice results
    #
    # = = = = = = = = = = = =
    keys = ['65-87','65-75','75-87']
    months = ['Apr','May','Jun','Jul','Aug','Sep']
    print('\nsimul   latrng  month  orig_ΔF avg_ΔF_err std_ΔF_err ' + \
          'min_ΔF_err max_ΔF_err avg_%_err std_%_err min_%_err max_%_err')
    for key in keys:
        for ii, month in enumerate(months):
            trend_errs = (ice_trend_dict['calc_trends'][key][:,ii] - \
                           ice_trend_dict['orig_trends'][key][ii])
            pcnt_errs =  (trend_errs / \
                           ice_trend_dict['orig_trends'][key][ii]) * 100.
        

            avg_trend_err = np.mean(trend_errs)
            std_trend_err = np.std(trend_errs)
            min_trend_err = np.min(trend_errs)
            max_trend_err = np.max(trend_errs)
            avg_pcnt_err = np.mean(pcnt_errs)
            std_pcnt_err = np.std(pcnt_errs)
            min_pcnt_err = np.min(pcnt_errs)
            max_pcnt_err = np.max(pcnt_errs)

            #ser_trend = std_trend / np.sqrt(len(uncert_trend_dict['calc_trends'][key][:,ii]))
            # Format:
            # 0: orignal trend
            # 1: mean of the errors between the original trend and calculated trends
            # 2: standard error of the trend errors calculated from the refice simulation
            # 3: max error of the trend errors calculated from the refice simulation
            # 4: min error of the trend errors calculated from the refice simulation
            # 5: mean of the percent trend errors between the original trend and calculated trends
            # 6: standard error of the percent trend errors calculated from the refice simulation
            print('refice   ' + key + '   ' + month + ' {0:>9.1e} {1:>9.1e} {2:>10.1e} {3:>10.1e} {4:>10.1e} {5:>8.1f}% {6:>8.1f}% {7:>8.1f}% {8:>8.1f}%'.format(\
                ice_trend_dict['orig_trends'][key][ii], \
                avg_trend_err, std_trend_err, min_trend_err, max_trend_err, \
                avg_pcnt_err, std_pcnt_err, min_pcnt_err, max_pcnt_err, \
                ))

    # = = = = = = = = = = = =
    #
    # refcld results
    #
    # = = = = = = = = = = = =
    print('\nsimul   latrng  month  orig_ΔF avg_ΔF_err std_ΔF_err ' + \
          'min_ΔF_err max_ΔF_err avg_%_err std_%_err min_%_err max_%_err')
    for key in keys:
        for ii, month in enumerate(months):
            trend_errs = (cld_trend_dict['calc_trends'][key][:,ii] - \
                           cld_trend_dict['orig_trends'][key][ii])
            pcnt_errs =  (trend_errs / \
                           cld_trend_dict['orig_trends'][key][ii]) * 100.
        
            avg_trend_err = np.mean(trend_errs)
            std_trend_err = np.std(trend_errs)
            min_trend_err = np.min(trend_errs)
            max_trend_err = np.max(trend_errs)
            avg_pcnt_err = np.mean(pcnt_errs)
            std_pcnt_err = np.std(pcnt_errs)
            min_pcnt_err = np.min(pcnt_errs)
            max_pcnt_err = np.max(pcnt_errs)

            # Format:
            # 0: orignal trend
            # 1: mean of the errors between the original trend and calculated trends
            # 2: standard error of the trend errors calculated from the refice simulation
            # 3: max error of the trend errors calculated from the refice simulation
            # 4: min error of the trend errors calculated from the refice simulation
            # 5: mean of the percent trend errors between the original trend and calculated trends
            # 6: standard error of the percent trend errors calculated from the refice simulation
            #print('refcld   ' + key + '   ' + month + ' {0:>11.1e} {1:>11.1e} {2:>12.1e} {3:>12.1e} {4:>12.1e} {5:>8.1f}% {6:>8.1f}% {7:>8.1f}% {8:>8.1f}%'.format(\
            print('refcld   ' + key + '   ' + month + ' {0:>9.1e} {1:>9.1e} {2:>10.1e} {3:>10.1e} {4:>10.1e} {5:>8.1f}% {6:>8.1f}% {7:>8.1f}% {8:>8.1f}%'.format(\
                cld_trend_dict['orig_trends'][key][ii], \
                avg_trend_err, std_trend_err, min_trend_err, max_trend_err, \
                avg_pcnt_err, std_pcnt_err, min_pcnt_err, max_pcnt_err, \
                ))

    # = = = = = = = = = = = =
    #
    # uncertainty results
    #
    # = = = = = = = = = = = =
    print('\n\nCalculating the mean and stdev of the uncertainty runs. ' + \
          'No abs. errors here. ser = standard error (stdev / sqrt(n))\n' + \
           '\nsimul   latrng  month  orig_ΔF avg_unc_ΔF std_unc_ΔF ser_unc_ΔF ' + \
          'min_unc_ΔF max_unc_ΔF avg_%_err std_%_err min_%_err max_%_err')
    for key in keys:
        for ii, month in enumerate(months):
            #trend_errs = (uncert_trend_dict['calc_trends'][key][:,ii] - \
            #               uncert_trend_dict['orig_trends'][key][ii])
       
            avg_trend = np.mean(uncert_trend_dict['calc_trends'][key][:,ii])
            std_trend = np.std(uncert_trend_dict['calc_trends'][key][:,ii])
            min_trend = np.min(uncert_trend_dict['calc_trends'][key][:,ii])
            max_trend = np.max(uncert_trend_dict['calc_trends'][key][:,ii])

            ser_trend = std_trend / np.sqrt(len(uncert_trend_dict['calc_trends'][key][:,ii]))

            #print(avg_trend, uncert_trend_dict['calc_trends'][key][:,ii])
            pcnt_errs =  ((uncert_trend_dict['calc_trends'][key][:,ii] - avg_trend) / \
                           avg_trend) * 100.
 
            avg_pcnt_err = np.mean(pcnt_errs)
            std_pcnt_err = np.std(pcnt_errs)
            min_pcnt_err = np.min(pcnt_errs)
            max_pcnt_err = np.max(pcnt_errs)
    
            # Format:
            # 0: original trend
            # 1: mean of all the uncertainty runs
            # 2: stdev of all the uncertainty runs
            # 3: min of all the uncertainty runs
            # 4: max of all the uncertainty runs
            print('uncert   ' + key + '   ' + month + ' {0:>9.1e} {1:>9.1e} {2:>10.1e} {3:>10.1e} {4:>10.1e} {5:>10.1e} {6:>8.1f}% {7:>8.1f}% {8:>8.1f}% {9:>8.1f}%'.format(\
                uncert_trend_dict['orig_trends'][key][ii], \
                avg_trend, std_trend, ser_trend, min_trend, max_trend, \
                avg_pcnt_err, std_pcnt_err, min_pcnt_err, max_pcnt_err, \
                ))


    print('\n\nCalculating the mean and stdev of the errors between the original and uncertainty runs.\n' + \
          '\nsimul   latrng  month  orig_ΔF avg_ΔF_err std_ΔF_err ser_ΔF_err ' + \
          'min_ΔF_err max_ΔF_err avg_%_err std_%_err min_%_err max_%_err')
    for key in keys:
        for ii, month in enumerate(months):
            trend_errs = (uncert_trend_dict['calc_trends'][key][:,ii] - \
                           uncert_trend_dict['orig_trends'][key][ii])
            pcnt_errs =  (trend_errs / \
                           cld_trend_dict['orig_trends'][key][ii]) * 100.

 
            avg_trend_err = np.mean(trend_errs)
            std_trend_err = np.std(trend_errs)
            min_trend_err = np.min(trend_errs)
            max_trend_err = np.max(trend_errs)
            avg_pcnt_err = np.mean(pcnt_errs)
            std_pcnt_err = np.std(pcnt_errs)
            min_pcnt_err = np.min(pcnt_errs)
            max_pcnt_err = np.max(pcnt_errs)
    
            ser_trend_err = std_trend_err / np.sqrt(len(uncert_trend_dict['calc_trends'][key][:,ii]))

            # Format:
            # 0: orignal trend
            # 1: mean of the errors between the original trend and calculated trends
            # 2: standard error of the trend errors calculated from the refice simulation
            # 3: max error of the trend errors calculated from the refice simulation
            # 4: min error of the trend errors calculated from the refice simulation
            # 5: mean of the percent trend errors between the original trend and calculated trends
            # 6: standard error of the percent trend errors calculated from the refice simulation
            #print('uncert   ' + key + '   ' + month + ' {0:>9.1e} {1:>9.1e} {2:>10.1e} {3:>10.1e} {4:>10.1e} {5:>8.1f}% {6:>8.1f}% {7:>8.1f}% {8:>8.1f}%'.format(\
            print('uncert   ' + key + '   ' + month + ' {0:>9.1e} {1:>9.1e} {2:>10.1e} {3:>10.1e} {4:>10.1e} {5:>10.1e} {6:>8.1f}% {7:>8.1f}% {8:>8.1f}% {9:>8.1f}%'.format(\
                uncert_trend_dict['orig_trends'][key][ii], \
                avg_trend_err, std_trend_err, ser_trend_err, min_trend_err, max_trend_err, \
                avg_pcnt_err, std_pcnt_err, min_pcnt_err, max_pcnt_err, \
                ))
        



# Assume that the slopes are only going to be calculated using "forcing"
#   and not "error" or "pcnt_error"
#
# Assume that the first entry in all_month_files is the "original"
# forcing values
# ----------------------------------------------------------------
def calc_forcing_slopes_v4_all_months_arctic_avg_uncert(all_month_files, \
        OMI_monthly_data, \
        min_year = 2005, max_year = 2020, \
        minlat = 70., maxlat = 87., \
        trend_type = 'standard', ptype = 'forcing', stype = 'ice', \
        vtype = '', slope_type = 'lin', \
        save = False):

    years = np.arange(min_year, max_year + 1)
    months = np.arange(6)
    file_types = [ptitle.strip().split('/')[-1].split('_')[-1][3:].split('.')[0] \
        for ptitle in all_month_files][1:]

    calc_trend_dict = {}
    calc_trend_dict['orig_trends'] = {}
    calc_trend_dict['orig_force_vals'] = {}
    calc_trend_dict['calc_trends'] = {}
    calc_trend_dict['calc_force_vals'] = {}
    calc_trend_dict['simulations'] = file_types

    # Read in the original values, which should be housed in 
    # index 0
    # -------------------------------------------------------
    all_month_dict = read_daily_month_force_HDF5(all_month_files[0])
    all_month_vals = all_month_dict['FORCE_EST']

    # ----------------------------------
    #
    # Calculate the original slopes
    #
    # ----------------------------------
    lat_mins = [65., 65., 75.]
    lat_maxs = [87., 75., 87.]

    for kk in range(len(lat_mins)):

        # Mask the month values that are outside the desired
        # lat bounds
        # --------------------------------------------------
        lats = OMI_monthly_data['LAT'][:,0]
        keep_idxs = np.where((lats[:] >= (lat_mins[kk] + 0.5)) & \
            (lats[:] < (lat_maxs[kk] + 0.5)))[0]
        print('HERE', lats[keep_idxs])

        # Convert the base values to Arctic and monthly averages
        # ------------------------------------------------------
        local_month_vals =  np.array([np.nanmean(all_month_vals[idx::6,keep_idxs,:], \
            axis = (1,2)) for idx in range(6)])

        # Combined all the single-year stuff into one array
        # ------------------------------------------------- 
        combined_vals = np.full( (len(all_month_files), \
            6, local_month_vals.shape[1]), np.nan)

        # orig_slopes: the slopes of the original data, with the control
        #               ice/cld values.
        # calc_slopes: the slopes of the modified data, with the altered
        #               ice/cld values from each ref year.
        orig_slopes = np.full((6), np.nan)
        calc_slopes = np.full( (len(file_types), 6), np.nan)

        for ii in range(6):
            # Calculate the slopes here
            # -------------------------
            if(trend_type == 'theil-sen'):
                res = stats.theilslopes(local_month_vals[ii,:], years, 0.90)
                delta_flux = res[0] * len(years)
                #print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))
            elif((trend_type == 'linregress') | (trend_type == 'standard')):
                zdata = stats.linregress(years, local_month_vals[ii,:])
                delta_flux = len(years)* zdata.slope
   
            orig_slopes[ii] = delta_flux

        if(vtype != ''):
            vtype = vtype + '_'

        # Begin the index at 1 here to skip the original values. These will be
        # housed in a different key in the dictionary
        # --------------------------------------------------------------------
        for ii, infile in enumerate(all_month_files[1:]):

            print(infile)
            data_dict = read_daily_month_force_HDF5(infile)
            month_values = data_dict['FORCE_EST']

            # Calculate the monthly averages of the estimated forcings
            # over the entire Arctic region
            # --------------------------------------------------------
            combined_vals[ii,:,:] = np.array([np.nanmean(month_values[idx::6,keep_idxs,:], \
                axis = (1,2)) for idx in range(6)])

            # Now, loop over each month and plot accordingly
            # ----------------------------------------------
            for jj in range(6):

                # Calculate the slopes here
                # -------------------------
                if(trend_type == 'theil-sen'):
                    res = stats.theilslopes(combined_vals[ii,jj,:], years, 0.90)
                    delta_flux = res[0] * len(years)
                    #print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))
                elif((trend_type == 'linregress') | (trend_type == 'standard')):
                    zdata = stats.linregress(years, combined_vals[ii,jj,:])
                    delta_flux = len(years)* zdata.slope
   
                calc_slopes[ii,jj] = delta_flux

        # Insert these trends into the output dictionary
        # ----------------------------------------------
        dict_key = str(int(lat_mins[kk])) + '-' + str(int(lat_maxs[kk]))
        calc_trend_dict['orig_trends'][dict_key] = orig_slopes
        calc_trend_dict['orig_force_vals'][dict_key] = local_month_vals
        calc_trend_dict['calc_trends'][dict_key] = calc_slopes
        calc_trend_dict['calc_force_vals'][dict_key] = combined_vals
 
    #return orig_slopes, calc_slopes
    return calc_trend_dict


# Assume that the slopes are only going to be calculated using "forcing"
#   and not "error" or "pcnt_error"
def calc_forcing_slopes_v4_all_months_arctic_avg_manyrefice(all_month_vals, \
        OMI_monthly_data, \
        min_year = 2005, max_year = 2020, \
        minlat = 70., maxlat = 87., \
        trend_type = 'standard', ptype = 'forcing', stype = 'ice', \
        vtype = '', slope_type = 'lin', \
        save = False):

    years = np.arange(min_year, max_year + 1)
    months = np.arange(6)

    calc_trend_dict = {}
    calc_trend_dict['orig_trends'] = {}
    calc_trend_dict['orig_force_vals'] = {}
    calc_trend_dict['calc_trends'] = {}
    calc_trend_dict['calc_force_vals'] = {}
    calc_trend_dict['years'] = years

    # ----------------------------------
    #
    # Calculate the original slopes
    #
    # ----------------------------------
    lat_mins = [65., 65., 75.]
    lat_maxs = [87., 75., 87.]

    for kk in range(len(lat_mins)):

        # Mask the month values that are outside the desired
        # lat bounds
        # --------------------------------------------------
        lats = OMI_monthly_data['LAT'][:,0]
        keep_idxs = np.where((lats[:] >= (lat_mins[kk] + 0.5)) & \
            (lats[:] < (lat_maxs[kk] + 0.5)))[0]
        print('HERE', lats[keep_idxs])

        # Convert the base values to Arctic and monthly averages
        # ------------------------------------------------------
        local_month_vals =  np.array([np.nanmean(all_month_vals[idx::6,keep_idxs,:], \
            axis = (1,2)) for idx in range(6)])

        # Combined all the single-year stuff into one array
        # ------------------------------------------------- 
        combined_vals = np.full( (years.shape[0], \
            6, local_month_vals.shape[1]), np.nan)

        # orig_slopes: the slopes of the original data, with the control
        #               ice/cld values.
        # calc_slopes: the slopes of the modified data, with the altered
        #               ice/cld values from each ref year.
        orig_slopes = np.full((6), np.nan)
        calc_slopes = np.full( (years.shape[0], 6), np.nan)

        for ii in range(6):
            # Calculate the slopes here
            # -------------------------
            if(trend_type == 'theil-sen'):
                res = stats.theilslopes(local_month_vals[ii,:], years, 0.90)
                delta_flux = res[0] * len(years)
                #print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))
            elif((trend_type == 'linregress') | (trend_type == 'standard')):
                zdata = stats.linregress(years, local_month_vals[ii,:])
                delta_flux = len(years)* zdata.slope
   
            orig_slopes[ii] = delta_flux

        if(vtype != ''):
            vtype = vtype + '_'

        # ----------------------------------
        #
        # Calculate the modified slopes
        #
        # ----------------------------------
        # Loop over each reference ice year, reading in the 
        # individual forcing values.
        for ii, year in enumerate(years):
            # Read in the force vals with this ref ice
            # ----------------------------------------
            if(slope_type == 'thl'):
                filename = home_dir + '/Research/Arctic_compares/' + \
                    'arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_ref' + \
                    stype + str(year) + '.hdf5'
            else:
                filename = home_dir + '/Research/Arctic_compares/' + \
                    'arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_lin_ref' + \
                    stype + str(year) + '.hdf5'

            # Read in the force vals with this ref ice
            # ----------------------------------------
            print(filename)
            all_month_dict = read_daily_month_force_HDF5(filename)
            month_values = all_month_dict['FORCE_EST']

            # Calculate the monthly averages of the estimated forcings
            # over the entire Arctic region
            # --------------------------------------------------------
            #arctic_avgs = np.array([np.nanmean(month_values[idx::6,:,:], \
            combined_vals[ii,:,:] = np.array([np.nanmean(month_values[idx::6,keep_idxs,:], \
                axis = (1,2)) for idx in range(6)])

            #if(ptype == 'error'):
            #    combined_vals[ii,:,:] = (combined_vals[ii,:,:] - \
            #        all_month_vals[:,:])
            #elif(ptype == 'pcnt_error'):
            #    combined_vals[ii,:,:] = ((combined_vals[ii,:,:] - all_month_vals[:,:]) / \
            #        all_month_vals[:,:]) * 100.
   
            # Now, loop over each month and plot accordingly
            # ----------------------------------------------
            for jj in range(6):

                # Calculate the slopes here
                # -------------------------
                if(trend_type == 'theil-sen'):
                    res = stats.theilslopes(combined_vals[ii,jj,:], years, 0.90)
                    delta_flux = res[0] * len(years)
                    #print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))
                elif((trend_type == 'linregress') | (trend_type == 'standard')):
                    zdata = stats.linregress(years, combined_vals[ii,jj,:])
                    delta_flux = len(years)* zdata.slope
   
                calc_slopes[ii,jj] = delta_flux

        # Insert these trends into the output dictionary
        # ----------------------------------------------
        dict_key = str(int(lat_mins[kk])) + '-' + str(int(lat_maxs[kk]))
        calc_trend_dict['orig_trends'][dict_key] = orig_slopes
        calc_trend_dict['orig_force_vals'][dict_key] = local_month_vals
        calc_trend_dict['calc_trends'][dict_key] = calc_slopes
        calc_trend_dict['calc_force_vals'][dict_key] = combined_vals
 
    #return orig_slopes, calc_slopes
    return calc_trend_dict




# Assume that the slopes are only going to be calculated using "forcing"
#   and not "error" or "pcnt_error"
def calc_forcing_slopes_v3_all_months_arctic_avg_manyrefice(all_month_vals, \
        OMI_monthly_data, \
        min_year = 2005, max_year = 2020, \
        minlat = 70., maxlat = 87., \
        trend_type = 'standard', ptype = 'forcing', stype = 'ice', \
        vtype = '', version4 = False, \
        save = False):

    years = np.arange(min_year, max_year + 1)
    months = np.arange(6)

    # Convert the base values to Arctic and monthly averages
    # ------------------------------------------------------
    all_month_vals =  np.array([np.nanmean(all_month_vals[idx::6,:,:], \
        axis = (1,2)) for idx in range(6)])

    # Combined all the single-year stuff into one array
    # ------------------------------------------------- 
    combined_vals = np.full( (years.shape[0], \
        6, all_month_vals.shape[1]), np.nan)

    # orig_slopes: the slopes of the original data, with the control
    #               ice/cld values.
    # calc_slopes: the slopes of the modified data, with the altered
    #               ice/cld values from each ref year.
    orig_slopes = np.full((6), np.nan)
    calc_slopes = np.full( (years.shape[0], 6), np.nan)

    # ----------------------------------
    #
    # Calculate the original slopes
    #
    # ----------------------------------
    #arctic_avgs = np.array([np.nanmean(all_month_vals[idx::6,:,:], \
    #    axis = (1,2)) for idx in range(6)])

    for ii in range(6):
        # Calculate the slopes here
        # -------------------------
        if(trend_type == 'theil-sen'):
            res = stats.theilslopes(all_month_vals[ii,:], years, 0.90)
            delta_flux = res[0] * len(years)
            #print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))
        elif((trend_type == 'linregress') | (trend_type == 'standard')):
            zdata = stats.linregress(years, all_month_vals[ii,:])
            delta_flux = len(years)* zdata.slope
   
        orig_slopes[ii] = delta_flux

    if(vtype != ''):
        vtype = vtype + '_'

    # ----------------------------------
    #
    # Calculate the modified slopes
    #
    # ----------------------------------
    # Loop over each reference ice year, reading in the 
    # individual forcing values.
    for ii, year in enumerate(years):
        # Read in the force vals with this ref ice
        # ----------------------------------------
        if(version4):
            filename = home_dir + '/Research/Arctic_compares/' + \
                'arctic_month_est_forcing_dayaithresh07_v4_aimin0_useintcpt_ref' + \
                stype + str(year) + '.hdf5'
        else:
            filename = home_dir + '/Research/Arctic_compares/' + \
                'arctic_month_est_forcing_dayaithresh07_' + vtype + 'ref' + \
                stype + str(year) + '.hdf5'
        print(filename)
        all_month_dict = read_daily_month_force_HDF5(filename)

        combined_vals[ii,:,:] = np.array([np.nanmean(\
            all_month_dict['FORCE_EST'][idx::6,:,:], \
            axis = (1,2)) for idx in range(6)])

        #if(ptype == 'error'):
        #    combined_vals[ii,:,:] = (combined_vals[ii,:,:] - \
        #        all_month_vals[:,:])
        #elif(ptype == 'pcnt_error'):
        #    combined_vals[ii,:,:] = ((combined_vals[ii,:,:] - all_month_vals[:,:]) / \
        #        all_month_vals[:,:]) * 100.
   
        # Now, loop over each month and plot accordingly
        # ----------------------------------------------
        for jj in range(6):

            # Calculate the slopes here
            # -------------------------
            if(trend_type == 'theil-sen'):
                res = stats.theilslopes(combined_vals[ii,jj,:], years, 0.90)
                delta_flux = res[0] * len(years)
                #print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))
            elif((trend_type == 'linregress') | (trend_type == 'standard')):
                zdata = stats.linregress(years, combined_vals[ii,jj,:])
                delta_flux = len(years)* zdata.slope
   
            calc_slopes[ii,jj] = delta_flux
 
    return orig_slopes, calc_slopes

# This function gathers the refice and refcld simulation results,
# calculates the trends from each refice and refcld simulation,
# and calculates the mean and standard error of those
# ref ice/cld simulations. The results are printed in a table
def calc_print_forcing_slope_error_v3(all_month_vals, \
        OMI_monthly_data, \
        min_year = 2005, max_year = 2020, \
        minlat = 65., maxlat = 87., \
        trend_type = 'standard', ptype = 'forcing', \
        vtype = '', version4 = False, 
        save = False):
     
    orig_slopes, ice_slopes = \
        calc_forcing_slopes_v3_all_months_arctic_avg_manyrefice(all_month_vals, \
                    OMI_monthly_data, minlat = minlat, trend_type = trend_type, stype = 'ice', \
                    vtype = vtype, ptype = 'forcing', version4 = version4)

    orig_slopes, cld_slopes = \
        calc_forcing_slopes_v3_all_months_arctic_avg_manyrefice(all_month_vals, \
                    OMI_monthly_data, minlat = minlat, trend_type = trend_type, stype = 'cld', \
                    vtype = vtype, ptype = 'forcing', version4 = version4)

    years = np.arange(min_year, max_year + 1)

    mean_ice_slopes = np.mean(ice_slopes, axis = 0)
    mean_cld_slopes = np.mean(cld_slopes, axis = 0)
    err_ice_slopes = np.std(ice_slopes, axis = 0) / np.sqrt(len(years))
    err_cld_slopes = np.std(cld_slopes, axis = 0) / np.sqrt(len(years))

    if(vtype != ''):
        print("V3 FORCING: PERT, MINLAT 65")

    # Format:
    # 0: orignal trend
    # 1: mean of the trends calculated from the refice simulations
    # 2: standard error of the trends calculated from the refice simulations
    # 3: mean of the trends calculated from the refcld simulations
    # 4: standard error of the trends calculated from the refcld simulations
    for ii in range(6):
        print('{0:1.3e} {1:1.3e} {2:1.3e} {3:1.3e} {4:1.3e}   '.format(orig_slopes[ii], \
            mean_ice_slopes[ii], err_ice_slopes[ii], mean_cld_slopes[ii], err_cld_slopes[ii] ))


# comp_type: 'aerosol' or 'cloud'
# pvar: 'modis_ch7' or 'modis_cod'
def plot_compare_CH7_dist(combined_data, sza_min, sza_max, save = False, \
        ai_thresh = 1.0, normalize = False, comp_type = 'aerosol', \
        pvar = 'modis_ch7'):

    plt.close('all')
    fig = plt.figure(figsize = (9, 9))
    ax1 = fig.add_subplot(4,2,1)
    ax2 = fig.add_subplot(4,2,3)
    ax3 = fig.add_subplot(4,2,5)
    ax4 = fig.add_subplot(4,2,7)
    ax5 = fig.add_subplot(4,2,2)
    ax6 = fig.add_subplot(4,2,4)
    ax7 = fig.add_subplot(4,2,6)
    ax8 = fig.add_subplot(4,2,8)

    # = = = = = = = = = = = = = = = = = = = = =
    #
    # "Clear-sky" (AI < 1.0)
    #
    # = = = = = = = = = = = = = = = = = = = = =

    def plot_single_comp_cloud(ax, combined_data, ice_min, ice_max, \
        sza_min, sza_max, ai_thresh, normalize = False, ptype = 'nosmoke', \
        pvar = 'modis_ch7'):

        # Plot the distributions for ocean
        # --------------------------------
        if(ptype == 'nosmoke'):
            clear_idx = np.where( ((combined_data['nsidc_ice'] >= ice_min) & \
                                   (combined_data['nsidc_ice'] <= ice_max) & \
                                   (combined_data['nsidc_ice'] != -999.)) & \
                                  (combined_data['modis_cld'] == 3) & \
                                  ((combined_data['omi_sza'] > sza_min) & \
                                   (combined_data['omi_sza'] < sza_max)) & \
                                  (combined_data['omi_uvai_pert'] < ai_thresh)\
                                )
            cloud_idx = np.where( ((combined_data['nsidc_ice'] >= ice_min) & \
                                   (combined_data['nsidc_ice'] <= ice_max) & \
                                   (combined_data['nsidc_ice'] != -999.)) & \
                                  (combined_data['modis_cld'] < 3) & \
                                  (combined_data['modis_cld'] >= 0) & \
                                  ((combined_data['omi_sza'] > sza_min) & \
                                   (combined_data['omi_sza'] < sza_max)) & \
                                  (combined_data['omi_uvai_pert'] < ai_thresh) \
                            )
        else:
            clear_idx = np.where( ((combined_data['nsidc_ice'] >= ice_min) & \
                                   (combined_data['nsidc_ice'] <= ice_max) & \
                                   (combined_data['nsidc_ice'] != -999.)) & \
                                  (combined_data['modis_cld'] == 3) & \
                                  ((combined_data['omi_sza'] > sza_min) & \
                                   (combined_data['omi_sza'] < sza_max)) & \
                                  (combined_data['omi_uvai_pert'] >= ai_thresh)\
                                )
            cloud_idx = np.where( ((combined_data['nsidc_ice'] >= ice_min) & \
                                   (combined_data['nsidc_ice'] <= ice_max) & \
                                   (combined_data['nsidc_ice'] != -999.)) & \
                                  (combined_data['modis_cld'] < 3) & \
                                  (combined_data['modis_cld'] >= 0) & \
                                  ((combined_data['omi_sza'] > sza_min) & \
                                   (combined_data['omi_sza'] < sza_max)) & \
                                  (combined_data['omi_uvai_pert'] >= ai_thresh) \
                            )
        
        print('Clear size:', clear_idx[0].shape, 'Cloud size:', cloud_idx[0].shape)

        ax.hist(combined_data[pvar][clear_idx], bins = 100, \
            alpha = 0.5, label = 'Clear')
        if(normalize):
            ax11 = ax.twinx()
            ax11.hist(combined_data[pvar][cloud_idx], \
                bins = 100, alpha = 0.5, label = 'Cloud', color = 'tab:orange')
            ax.set_ylabel('Clear counts')
            ax11.set_ylabel('Cloud counts')
        else:
            ax.hist(combined_data[pvar][cloud_idx], bins = 100, \
                alpha = 0.5, label = 'Cloud')


    def plot_single_comp_aerosol(ax, combined_data, ice_min, ice_max, \
        cld_val, sza_min, sza_max, ai_thresh, normalize = False):

        if(cld_val == 3):
            # Plot the distributions for ocean
            # --------------------------------
            clear_idx = np.where( ((combined_data['nsidc_ice'] >= ice_min) & \
                                   (combined_data['nsidc_ice'] <= ice_max) & \
                                   (combined_data['nsidc_ice'] != -999.)) & \
                                  (combined_data['modis_cld'] == cld_val) & \
                                  ((combined_data['omi_sza'] > sza_min) & \
                                   (combined_data['omi_sza'] < sza_max)) & \
                                  (combined_data['omi_uvai_pert'] < ai_thresh)\
                                )
            smoke_idx = np.where( ((combined_data['nsidc_ice'] >= ice_min) & \
                                   (combined_data['nsidc_ice'] <= ice_max) & \
                                   (combined_data['nsidc_ice'] != -999.)) & \
                                  (combined_data['modis_cld'] == cld_val) & \
                                  ((combined_data['omi_sza'] > sza_min) & \
                                   (combined_data['omi_sza'] < sza_max)) & \
                                  (combined_data['omi_uvai_pert'] >= ai_thresh) \
                            )
        else:
            clear_idx = np.where( ((combined_data['nsidc_ice'] >= ice_min) & \
                                   (combined_data['nsidc_ice'] <= ice_max) & \
                                   (combined_data['nsidc_ice'] != -999.)) & \
                                  (combined_data['modis_cld'] < 3) & \
                                  (combined_data['modis_cld'] >= 0) & \
                                  ((combined_data['omi_sza'] > sza_min) & \
                                   (combined_data['omi_sza'] < sza_max)) & \
                                  (combined_data['omi_uvai_pert'] < ai_thresh)\
                                )
            smoke_idx = np.where( ((combined_data['nsidc_ice'] >= ice_min) & \
                                   (combined_data['nsidc_ice'] <= ice_max) & \
                                   (combined_data['nsidc_ice'] != -999.)) & \
                                  (combined_data['modis_cld'] < 3) & \
                                  (combined_data['modis_cld'] >= 0) & \
                                  ((combined_data['omi_sza'] > sza_min) & \
                                   (combined_data['omi_sza'] < sza_max)) & \
                                  (combined_data['omi_uvai_pert'] >= ai_thresh) \
                            )
        
        print('Clear size:', clear_idx[0].shape, 'Smoke size:', smoke_idx[0].shape)

        ax.hist(combined_data[pvar][clear_idx], bins = 100, \
            alpha = 0.5, label = 'Clear')
        if(normalize):
            ax11 = ax.twinx()
            ax11.hist(combined_data[pvar][smoke_idx], \
                bins = 100, alpha = 0.5, label = 'Cloud', color = 'tab:orange')
            ax.set_ylabel('Clear counts')
            ax11.set_ylabel('Smoky counts')
        else:
            ax.hist(combined_data[pvar][smoke_idx], bins = 100, \
                alpha = 0.5, label = 'Smoke')


    if(comp_type == 'aerosol'):
        # Ocean, clear-sky (no smoke)
        plot_single_comp_cloud(ax1, combined_data, 0, 20, sza_min, sza_max, \
            ai_thresh, normalize = normalize, ptype = 'nosmoke', pvar = pvar)

        # Ice, clear-sky (no smoke)
        plot_single_comp_cloud(ax2, combined_data, 80, 100, sza_min, sza_max, \
            ai_thresh, normalize = normalize, ptype = 'nosmoke', pvar = pvar)

        # Land, clear-sky (no smoke)
        plot_single_comp_cloud(ax3, combined_data, 253, 253, sza_min, sza_max, \
            ai_thresh, normalize = normalize, ptype = 'nosmoke', pvar = pvar)

        # Mix, clear-sky (no smoke)
        plot_single_comp_cloud(ax4, combined_data, 20.001, 79.999, sza_min, sza_max, \
            ai_thresh, normalize = normalize, ptype = 'nosmoke', pvar = pvar)

        # Ocean, smoky
        plot_single_comp_cloud(ax5, combined_data, 0, 20, sza_min, sza_max, \
            ai_thresh, normalize = normalize, ptype = 'smoke', pvar = pvar)
        
        # Ice, smoky
        plot_single_comp_cloud(ax6, combined_data, 80, 100, sza_min, sza_max, \
            ai_thresh, normalize = normalize, ptype = 'smoke', pvar = pvar)
        
        # Land, smoky
        plot_single_comp_cloud(ax7, combined_data, 253, 253, sza_min, sza_max, \
            ai_thresh, normalize = normalize, ptype = 'smoke', pvar = pvar)
        
        # Mix, smoky
        plot_single_comp_cloud(ax8, combined_data, 20.001, 79.999, sza_min, sza_max, \
            ai_thresh, normalize = normalize, ptype = 'smoke', pvar = pvar)

        ax1.set_title('Ocean Clear-sky')
        ax2.set_title('Ice Clear-sky')
        ax3.set_title('Land Clear-sky')
        ax4.set_title('Mix Clear-sky')
        ax5.set_title('Ocean Aerosol-sky')
        ax6.set_title('Ice Aerosol-sky')
        ax7.set_title('Land Aerosol-sky')
        ax8.set_title('Mix Aerosol-sky')


    elif(comp_type == 'cloud'):
        # Ocean, clear (show the distributions between clear and cloudy, pvar = pvar)
        plot_single_comp_aerosol(ax1, combined_data, 0, 20, 3, sza_min, sza_max, \
            ai_thresh, normalize = normalize, pvar = pvar)

        # Ice, clear (show the distributions between clear and cloudy, pvar = pvar)
        plot_single_comp_aerosol(ax2, combined_data, 80, 100, 3, sza_min, sza_max, \
            ai_thresh, normalize = normalize, pvar = pvar)

        # Land, clear (show the distributions between clear and cloudy, pvar = pvar)
        plot_single_comp_aerosol(ax3, combined_data, 253, 253, 3, sza_min, sza_max, \
            ai_thresh, normalize = normalize, pvar = pvar)

        # Mix, clear (show the distributions between clear and cloudy, pvar = pvar)
        plot_single_comp_aerosol(ax4, combined_data, 20.001, 79.999, 3, sza_min, sza_max, \
            ai_thresh, normalize = normalize, pvar = pvar)

        # Ocean, cloud (show the distributions between clear and cloudy, pvar = pvar)
        plot_single_comp_aerosol(ax5, combined_data, 0, 20, 0, sza_min, sza_max, \
            ai_thresh, normalize = normalize, pvar = pvar)

        # Ice, cloud (show the distributions between clear and cloudy, pvar = pvar)
        plot_single_comp_aerosol(ax6, combined_data, 80, 100, 0, sza_min, sza_max, \
            ai_thresh, normalize = normalize, pvar = pvar)

        # Land, cloud (show the distributions between clear and cloudy, pvar = pvar)
        plot_single_comp_aerosol(ax7, combined_data, 253, 253, 0, sza_min, sza_max, \
            ai_thresh, normalize = normalize, pvar = pvar)

        # Mix, cloud (show the distributions between clear and cloudy, pvar = pvar)
        plot_single_comp_aerosol(ax8, combined_data, 20.001, 79.999, 0, sza_min, sza_max, \
            ai_thresh, normalize = normalize, pvar = pvar)

        ax1.set_title('Ocean No Cloud')
        ax2.set_title('Ice No Cloud')
        ax3.set_title('Land No Cloud')
        ax4.set_title('Mix No Cloud')
        ax5.set_title('Ocean Cloud')
        ax6.set_title('Ice Cloud')
        ax7.set_title('Land Cloud')
        ax8.set_title('Mix Cloud')

    else:
        print("BAD COMP TYPE")
        return


    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax5.legend()
    ax6.legend()
    ax7.legend()
    ax8.legend()


    fig.tight_layout()
    plt.show()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# ERROR FUNCTIONS
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def test_error_calc(OMI_daily_data, OMI_monthly_data, coloc_dict, \
        date_str, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw',\
        sfc_type_idx = 0, \
        filter_bad_vals = True, return_modis_nsidc = True, debug = False):

    print("error stuff")

    # Need to know:
    # - Error of the slope for each individual bin (?????)
    # - Error of each SZA-meaned slope
    #     - can this be derived from the "n" slopes from the "n" cloud/clear
    #       bins?
    # - Error of the MODIS MYD08 daily cloud fraction
    #     - again, can this be derived from the standard deviation and the
    #       ob counts provided in the data?
    # - Error of each single-day forcing estimate
    # - Error of each daily-gridded monthly estimate
    # - Error of the trend 

    # FOR NOW, IGNORING THE STANDARD ERROR OF THE SLOPE FROM EACH
    # BIN. WILL COME BACK TO THIS AND INCORPORATE THIS INTO THE 
    # FOLLOWING ERROR ANALYSIS WORKFLOW. - 2023/12/04

    # Retrieve the clear and cloudy slopes
    # ------------------------------------
    clear_slopes, cloud_slopes = identify_slope_clear_clean_sfctype(\
        coloc_dict, sfc_type_idx, cld_idx, min_cloud = min_cloud, \
        maxerr = maxerr, data_type = data_type, use_error = False)

    # Calculate the standard error of the SZA-meaned slopes using
    #
    # σ_error = std_dev / sqrt(N)
    #
    # ------------------------------------------------------------
    clear_error, cloud_error = calc_slope_error(clear_slopes, cloud_slopes)

    fig = plt.figure(figsize = (9, 4))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    ax1.plot(clear_mean, color = 'tab:blue')
    ax1.plot(clear_mean + clear_error, linestyle = '--',  color = 'tab:blue')
    ax1.plot(clear_mean - clear_error, linestyle = '--',  color = 'tab:blue')
    ax1.plot(clear_mean + clear_std, linestyle = ':',  color = 'tab:blue')
    ax1.plot(clear_mean - clear_std, linestyle = ':',  color = 'tab:blue')
    ax1.plot(cloud_mean, color = 'tab:orange')
    ax1.plot(cloud_mean + cloud_error, linestyle = '--',  color = 'tab:orange')
    ax1.plot(cloud_mean - cloud_error, linestyle = '--',  color = 'tab:orange')
    ax1.plot(cloud_mean + cloud_std, linestyle = ':',  color = 'tab:orange')
    ax1.plot(cloud_mean - cloud_std, linestyle = ':',  color = 'tab:orange')
    ax1.axhline(0, color = 'k')
    ax1.grid(alpha = 0.5, color = 'grey')

    # Plot fractional uncertainties
    # -----------------------------
    clear_frac_uncert = clear_error / np.abs(clear_mean)
    cloud_frac_uncert = cloud_error / np.abs(cloud_mean)

    ax2.plot(clear_frac_uncert)
    ax2.plot(cloud_frac_uncert)
    ax2.axhline(0, color = 'k')
    ax2.grid(alpha = 0.5, color = 'grey')

    fig.tight_layout()
    plt.show()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Neural Network functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


# This function selects all the neural network output files 
# for the given simulation , while also applying some QC criteria,
# and returns a dictionary containing the observed and NN-estimated
# SWF, in addition to the binning variables.
def combine_NN_data(sim_name):

    files = glob('neuralnet_output/test_calc_out_' + sim_name + '*.hdf5')

    if(len(files) == 0):
        print("ERROR: NO FILES FOUND FOR SIM " + sim_name)
        return

    # Figure out the total size to insert the data
    # ---------------------------------------------
    #minlat = 70.
    minlat = 65.    # 2024/01/10: changed to 65.
    total_size = 0
    for ff in files:
        data = h5py.File(ff,'r')
        print(ff)
        #local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
        local_data = np.ma.masked_invalid(data['omi_uvai_pert'])
        local_data = np.ma.masked_where(abs(local_data) > 12, local_data)
        local_data = np.ma.masked_where(data['omi_lat'][:,:] < minlat, \
            local_data) 
        local_data = np.ma.masked_where((data['ceres_swf'][:,:] < -200.) | \
            (data['ceres_swf'][:,:] > 3000), \
            local_data) 
        local_data = np.ma.masked_where(np.isnan(data['calc_swf'][:,:]), local_data)
        local_size = local_data.compressed().shape[0]
        #print(ff, local_size)
        total_size += local_size
    
        data.close()
    
    
    # Set up the data structure to hold all the data
    # ----------------------------------------------
    combined_data = {}
    combined_data['omi_uvai_pert'] = np.full(total_size, np.nan)
    #combined_data['omi_uvai_raw']  = np.full(total_size, np.nan)
    combined_data['modis_cld_top_pres']     = np.full(total_size, np.nan)
    combined_data['modis_cod']     = np.full(total_size, np.nan)
    combined_data['calc_swf']      = np.full(total_size, np.nan)
    combined_data['ceres_swf']     = np.full(total_size, np.nan)
    #combined_data['modis_ch7']     = np.full(total_size, np.nan)
    combined_data['omi_sza']       = np.full(total_size, np.nan)
    combined_data['omi_lat']       = np.full(total_size, np.nan)
    combined_data['nsidc_ice']     = np.full(total_size, np.nan)
 
    print("Loading data")
    
    # Loop back over the files and insert the data into the structure
    # ---------------------------------------------------------------
    total_size = 0
    beg_idx = 0
    end_idx = 0
    for ff in files:
    
        data = h5py.File(ff,'r')
        # NOTE: Changed the omi variable here from "pert" to "raw" on 20230623.
        #       This move should allow for coloc data to be read after 2020
        #local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
        local_data = np.ma.masked_invalid(data['omi_uvai_pert'])
        local_data = np.ma.masked_where(abs(local_data) > 12, local_data)
        local_data = np.ma.masked_where(data['omi_lat'][:,:] < minlat, \
            local_data) 
        local_data = np.ma.masked_where((data['ceres_swf'][:,:] < -200.) | \
            (data['ceres_swf'][:,:] > 3000), \
            local_data) 
        local_data = np.ma.masked_where(np.isnan(data['calc_swf'][:,:]), local_data)
        local_size = local_data.compressed().shape[0]
    
        beg_idx = end_idx
        end_idx = beg_idx + local_size
    
        for tkey in combined_data.keys():
            combined_data[tkey][beg_idx:end_idx] = \
                data[tkey][~local_data.mask]
    
        print(local_size)
        total_size += local_size
    
        data.close()

    combined_data['sim_name']      = sim_name   

    return combined_data

# This function identifies the indices that give the NN data
# that meet the AI, ICE, COD, and SZA requirements
def select_idxs(test_dict, ai_min, ice_min, ice_max, \
        cod_min, cod_max, sza_min, sza_max):

    good_idxs = np.where(\
        (test_dict['omi_uvai_pert'] >= ai_min) & \
        (test_dict['nsidc_ice'] >= ice_min) & \
        (test_dict['nsidc_ice'] < ice_max) & \
        (test_dict['modis_cod'] >= cod_min) & \
        (test_dict['modis_cod'] < cod_max) & \
        ((test_dict['omi_sza'] >= sza_min) & \
         (test_dict['omi_sza'] < sza_max))
    )

    return good_idxs

# This function calculates the slope and intercepts for the 
# forcing vs AI values. NOTE: here, forcing is calculated as NN - OBS.
# 
# NOTE: This is also where the error stuff can be done by 
#       switching the slope found by Theil-Sen (res[0])
#       to either the lower bound slope (res[2]) or 
#       upper bound slope (res[3])
# --------------------------------------------------------------------
def calc_NN_force_slope_intcpt(test_dict, ice_bin_edges, \
        sza_bin_edges, cod_bin_edges, ai_min = 1, min_ob = 50, \
        trend_type = 'theil-sen'):

    slope_dict = {}
    
    calc_slopes = np.full(((ice_bin_edges.shape[0] - 1), \
                           (sza_bin_edges.shape[0] - 1), \
                           (cod_bin_edges.shape[0] - 1)), np.nan)
    calc_upper_slopes  = np.full(((ice_bin_edges.shape[0] - 1), \
                            (sza_bin_edges.shape[0] - 1), \
                            (cod_bin_edges.shape[0] - 1)), np.nan)
    calc_lower_slopes  = np.full(((ice_bin_edges.shape[0] - 1), \
                            (sza_bin_edges.shape[0] - 1), \
                            (cod_bin_edges.shape[0] - 1)), np.nan)
    calc_intcpt = np.full(((ice_bin_edges.shape[0] - 1), \
                           (sza_bin_edges.shape[0] - 1), \
                           (cod_bin_edges.shape[0] - 1)), np.nan)
    calc_intcpt_error = np.full(((ice_bin_edges.shape[0] - 1), \
                           (sza_bin_edges.shape[0] - 1), \
                           (cod_bin_edges.shape[0] - 1)), np.nan)
    calc_counts = np.full(((ice_bin_edges.shape[0] - 1), \
                           (sza_bin_edges.shape[0] - 1), \
                           (cod_bin_edges.shape[0] - 1)), np.nan)
    
    # = = = = = = = = =
    # 
    # NOTE: From Scipy documentation, the theilslopes returns the following 4
    #       variables: slope, intercept, low_slope, and high_slope. The last
    #       two give the lower and upper bounds of the confidence interval
    #       for the slope. COULD THIS BE USED FOR THE UNCERTAINTY ANALYSIS
    #       INSTEAD OF THE SLOPE STANDARD ERROR VALUE FROM LINREGRESS???
    #
    #       Values 3 and 4 are slopes, but the lower bound slope of the 90%
    #       confidence interval and the upper bound slope. They do not
    #       contain any information about the confidence intervals of 
    #       the intercept, so the result is not a "true" regression 
    #       confidence interval.
    #
    # = = = = = = = = =
    
    # Loop over each of the bins and find the slope and intercept
    # (either from linear regression or from Theil-Sen)
    #trend_type = 'theil-sen'
    dict_dims = ['ice','sza','cod']
    for ii in range(ice_bin_edges.shape[0] - 1):
        print(ice_bin_edges[ii])
        for jj in range(sza_bin_edges.shape[0] - 1):
            for kk in range(cod_bin_edges.shape[0] - 1):
                good_idxs = select_idxs(test_dict, ai_min, \
                    ice_bin_edges[ii], ice_bin_edges[ii + 1], \
                    cod_bin_edges[kk], cod_bin_edges[kk+1], \
                    sza_bin_edges[jj], sza_bin_edges[jj+1])
    
                calc_counts[ii,jj,kk] = good_idxs[0].shape[0]
                
                if(calc_counts[ii,jj,kk] > 2):
                    raw_ydata = test_dict['calc_swf'][good_idxs] - test_dict['ceres_swf'][good_idxs]
                    raw_xdata = test_dict['omi_uvai_pert'][good_idxs] 
    
                    if(trend_type == 'theil-sen'):
                        #if((len(raw_xdata) >= 20)):
                        res = stats.theilslopes(raw_ydata, raw_xdata, 0.90)
                        calc_slopes[ii,jj,kk]   = res[0]
                        calc_intcpt[ii,jj,kk]   = res[1]
                        # NOTE: If wanting to do an uncertainty estimate with
                        #       Theil-Sen, this is where it would be done.
                        calc_upper_slopes[ii,jj,kk]   = res[2] # lower bound slope
                        calc_lower_slopes[ii,jj,kk]   = res[3] # upper bound slope
    
                        #if((len(grid_xdata) > 20)):
                        #    res = stats.theilslopes(grid_ydata, grid_xdata, 0.90)
                        #    grid_slopes[ii,jj,kk]   = res[0]
                    elif(trend_type == 'linregress'):
                        #if((len(raw_xdata) >= 20)):
                        result1 = stats.linregress(raw_xdata,raw_ydata)
                        calc_slopes[ii,jj,kk] = result1.slope 
                        calc_intcpt[ii,jj,kk] = result1.intercept
                        calc_intcpt_error[ii,jj,kk] = result1.intercept_stderr
                        calc_upper_slopes[ii,jj,kk] = result1.stderr
                        calc_lower_slopes[ii,jj,kk] = result1.pvalue
   
    calc_slopes = np.ma.masked_where(calc_counts < min_ob, calc_slopes)
    calc_upper_slopes = np.ma.masked_where(calc_counts < min_ob, calc_upper_slopes)
    calc_lower_slopes = np.ma.masked_where(calc_counts < min_ob, calc_lower_slopes)
    calc_intcpt = np.ma.masked_where(calc_counts < min_ob, calc_intcpt)
    calc_intcpt_error = np.ma.masked_where(calc_counts < min_ob, calc_intcpt_error)
    calc_counts = np.ma.masked_where(calc_counts < min_ob, calc_counts)

    slope_dict['sim_name'] = test_dict['sim_name'] 
    slope_dict['counts'] = calc_counts
    slope_dict['slopes'] = calc_slopes
    slope_dict['intercepts'] = calc_intcpt
    slope_dict['min_ob'] = min_ob
    slope_dict['ai_min'] = ai_min
    slope_dict['dict_dims'] = dict_dims
    if(trend_type == 'theil-sen'):
        slope_dict['upper_slopes'] = calc_upper_slopes
        slope_dict['lower_slopes'] = calc_lower_slopes
    else:
        slope_dict['intcpt_stderr'] = calc_intcpt_error
        slope_dict['slope_stderr'] = calc_upper_slopes
        slope_dict['slope_pvals']  = calc_lower_slopes
    slope_dict['trend_type'] = trend_type

    return slope_dict

# Plot forcing slopes for 6 sfc types. For now, only plot
# the 6 slopes
# pvar: 'slopes', 'counts', 'errors'
# -----------------------------------
def plot_NN_bin_slopes_6types(slope_dict, bin_dict, pvar, min_ob = 50, \
            plot_error = False, save = False):

    plt.close('all')

    num_sfc_types = bin_dict['ice_bin_means'].shape[0]

    fig = plt.figure(figsize = (9, 7))
    axs = fig.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True)
    #nrows_ncols = (2, 3)
    #grid = ImageGrid(fig, 111, \
    #    nrows_ncols = (3, 4), \
    #    axes_pad = 0.20, \
    #    share_all = True, \
    #    cbar_location = 'right', \
    #    cbar_mode = 'edge', \
    #    cbar_size = '7%', \
    #    cbar_pad = 0.15)

    #grid = ImageGrid(fig, 111, \
    #    nrows_ncols = nrows_ncols, \
    #    axes_pad = 0.20, \
    #    share_all = True, \
    #    cbar_location = 'bottom', \
    #    cbar_mode = 'edge', \
    #    cbar_size = '7%', \
    #    cbar_pad = 0.15)

    #axs = [ax for ax in grid] 
    axs = axs.flatten()
    ax1 = axs[0]  # ocean force
    ax2 = axs[1]  # mix force
    ax3 = axs[2]  # ice force
    ax4 = axs[3]  # land force
    ax5 = axs[4]  # ocean cloud
    ax6 = axs[5]  # mix cloud
 
    #calc_slopes = np.ma.masked_where(abs(slope_dict['slopes']) > 50, slope_dict['slopes'])
    calc_slopes = np.ma.masked_where(slope_dict['counts'] < min_ob, slope_dict['slopes'])
    calc_counts = np.ma.masked_where(slope_dict['counts'] < min_ob, slope_dict['counts'])
    if(pvar == 'errors'):
        calc_errors = np.ma.masked_where(slope_dict['counts'] < min_ob, slope_dict['slope_stderr'])
  
    print('slope mins', np.min(calc_slopes, axis = (1,2)))
    print('slope maxs', np.max(calc_slopes, axis = (1,2)))

    print('count mins', np.min(calc_counts, axis = (1,2)))
    print('count maxs', np.max(calc_counts, axis = (1,2)))

    if(pvar == 'errors'):
        shrk = 1.0
        cmap2 = plt.get_cmap('jet')
        colorvals = np.arange(0, 6, 1)
        norm = mpl.colors.BoundaryNorm(colorvals, cmap2.N, extend = 'max')

    for ii in range(6): 
        if(pvar == 'slopes'):
            mesh = axs[ii].pcolormesh(calc_slopes[ii,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
        elif(pvar == 'counts'):
            mesh = axs[ii].pcolormesh(calc_counts[ii,:,:].T, cmap = 'viridis', vmin = 0, vmax = 5000)
        elif(pvar == 'errors'):
            mesh = axs[ii].pcolormesh(calc_errors[ii,:,:].T, norm = norm, shading = 'auto', cmap = 'jet')

        if(bin_dict['ice_bin_edges'][ii] == 0.):
            title_str = 'Ocean\n(Ice conc. < ' + \
                str(int(bin_dict['ice_bin_edges'][ii + 1])) + '%)'
        elif( (bin_dict['ice_bin_edges'][ii] >= 20.) & \
                (bin_dict['ice_bin_edges'][ii + 1] <= 80) ):
            title_str = 'Mix\n(' + str(int(bin_dict['ice_bin_edges'][ii])) + \
                '% <= Ice conc. < ' + \
                str(int(bin_dict['ice_bin_edges'][ii + 1])) + '%)'
        elif( (bin_dict['ice_bin_edges'][ii] <= 80.) & \
                (bin_dict['ice_bin_edges'][ii + 1] <= 100.2) ):
            title_str = 'Ice\n(Ice conc. >= ' + \
                str(int(bin_dict['ice_bin_edges'][ii])) + '%)'
        elif( (bin_dict['ice_bin_edges'][ii + 1] > 150.) ):
            title_str = 'Land'
        axs[ii].set_title(title_str)
        #axs[ii].set_title(title_str + '\nn = ' + str(total_count))


    """ 
    ax1.pcolormesh(calc_slopes[0,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    ax2.pcolormesh(calc_slopes[1,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    ax3.pcolormesh(calc_slopes[2,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    ax4.pcolormesh(calc_slopes[3,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    ax5.pcolormesh(calc_slopes[4,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    mesh = ax6.pcolormesh(calc_slopes[5,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    """ 
    #cbar = plt.colorbar(mesh, ax = ax4, orientation = 'vertical', \
    #    shrink = 0.8, extend = 'both')

    if(pvar == 'slopes'):
        plabel = 'Forcing Efficiency [ΔForcing / ΔAI]'
    elif(pvar == 'counts'):
        plabel = 'Counts]'
    elif(pvar == 'errors'):
        plabel = 'Forcing Efficiency Errors [ΔForcing / ΔAI]'
    
    cbar_ax = fig.add_axes([0.17, 0.09, 0.70, 0.02])
    cbar = fig.colorbar(mesh, \
        cax = cbar_ax, shrink = 0.8, orientation = 'horizontal', \
        label = plabel)
    fig.tight_layout(rect = [0,0.1,1,1])

    sza_labels = [str(edge) for edge in bin_dict['sza_bin_edges']]
    cod_labels = [str(edge) for edge in bin_dict['cod_bin_edges']]
 
    #if(slope_dict['trend_type'] == 'theil-sen'): 
    if((bin_dict['sza_bin_edges'][1] - bin_dict['sza_bin_edges'][0]) <= 5):
        xtick_setter = np.arange(0, len(bin_dict['sza_bin_means']) + 1, 2)
        label_setter = sza_labels[::2]
    else:
        xtick_setter = np.arange(0, len(bin_dict['sza_bin_means']) + 1)
        label_setter = sza_labels[::1]
    if(not plot_error): 
        ax5.set_xticks(xtick_setter)
        ax5.set_xticklabels(label_setter)
    else:
        ax9.set_xticks(xtick_setter)
        ax9.set_xticklabels(label_setter)
    ax1.set_yticks([0,1,2,3,4,5,6,7,8])
    ax1.set_yticklabels(cod_labels[::1])
    ax4.set_yticks([0,1,2,3,4,5,6,7,8])
    ax4.set_yticklabels(cod_labels[::1])
     
    #if(slope_dict['trend_type'] == 'theil-sen'): 
    if(not plot_error): 
        ax4.set_xlabel('OMI SZA [$^{o}$]')
        ax5.set_xlabel('OMI SZA [$^{o}$]')
        ax6.set_xlabel('OMI SZA [$^{o}$]')

    ax1.set_ylabel('MODIS COD')
    ax4.set_ylabel('MODIS COD')

    #axs[0].cax.colorbar(mesh, label = 'δForcing/δAI')
    #axs[5].cax.colorbar(mesh2, label = 'Counts')
    ##if(slope_dict['trend_type'] == 'linregress'):
    #if(plot_error):
    #    axs[11].cax.colorbar(mesh3, label = 'Slope Stderr')
        

    if(slope_dict['trend_type'] == 'theil-sen'):
        type_text = 'Theil-Sen'
    else:
        type_text = 'Linear'
    
    if(pvar == 'slopes'):
        plt.suptitle('Slopes of UVAI vs. ADRF ' + type_text + ' Regression Equations\n' + \
            'Plotted for Grids with Counts $\geq$ ' + str(int(min_ob)))
    elif(pvar == 'counts'): 
        plt.suptitle('Counts in each bin\nGrids with counts >= ' + str(int(min_ob)))
    elif(pvar == 'errors'): 
        plt.suptitle('Slope Errors of AI vs. Forcing Regression (' + type_text + \
            ')\nGrids with counts >= ' + str(int(min_ob)))

    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax4, 'd)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax5, 'e)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax6, 'f)', fontsize = 11, backgroundcolor = None)

 
    #fig.tight_layout()
    fig.tight_layout(rect = [0,0.1,1,1])
 
    #fig.tight_layout()
    if(save): 
        if(slope_dict['trend_type'] == 'theil-sen'):
            type_adder = '_theilsen'
        else:
            type_adder = '_linregress'
        if(plot_error):
            error_add = '_error'
        else:
            error_add = ''
        outname = 'force_bin_' + pvar + '_NN_numsfcbins6_' + slope_dict['sim_name']  + type_adder + error_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image:", outname)
    else:
        plt.show()


# Plot the binned slopes for the 4 surface types
# ----------------------------------------------
def plot_NN_bin_slopes(slope_dict, bin_dict, min_ob = 50, \
            plot_error = False, save = False):

    plt.close('all')
    #ax1 = fig.add_subplot(2,4,1)
    #ax2 = fig.add_subplot(2,4,2)
    #ax3 = fig.add_subplot(2,4,3)
    #ax4 = fig.add_subplot(2,4,4)
    #ax5 = fig.add_subplot(2,4,5)
    #ax6 = fig.add_subplot(2,4,6)
    #ax7 = fig.add_subplot(2,4,7)
    #ax8 = fig.add_subplot(2,4,8)

    #axs = fig.subplots(nrows = 2, ncols = 4, sharex = True, sharey = True)
    #ax1 = axs[0,0]  # ocean force
    #ax2 = axs[0,1]  # mix force
    #ax3 = axs[0,2]  # ice force
    #ax4 = axs[0,3]  # land force
    #ax5 = axs[1,0]  # ocean cloud
    #ax6 = axs[1,1]  # mix cloud
    #ax7 = axs[1,2]  # ice cloud
    #ax8 = axs[1,3]  # land cloud

    num_sfc_types = bin_dict['ice_bin_means'].shape[0]

    if((slope_dict['trend_type'] == 'theil-sen') | \
            ((slope_dict['trend_type'] == 'linregress') & \
             (not plot_error))):  
        fig = plt.figure(figsize = (10, 5.5))
        nrows_ncols = (2, num_sfc_types)
    else:
        fig = plt.figure(figsize = (10, 7.5))
        nrows_ncols = (3, 4)
        #grid = ImageGrid(fig, 111, \
        #    nrows_ncols = (3, 4), \
        #    axes_pad = 0.20, \
        #    share_all = True, \
        #    cbar_location = 'right', \
        #    cbar_mode = 'edge', \
        #    cbar_size = '7%', \
        #    cbar_pad = 0.15)

    grid = ImageGrid(fig, 111, \
        nrows_ncols = nrows_ncols, \
        axes_pad = 0.20, \
        share_all = True, \
        cbar_location = 'right', \
        cbar_mode = 'edge', \
        cbar_size = '7%', \
        cbar_pad = 0.15)

    axs = [ax for ax in grid] 
    ax1 = axs[0]  # ocean force
    ax2 = axs[1]  # mix force
    ax3 = axs[2]  # ice force
    ax4 = axs[3]  # land force
    ax5 = axs[4]  # ocean cloud
    ax6 = axs[5]  # mix cloud
    ax7 = axs[6]  # ice cloud
    ax8 = axs[7]  # land cloud
    #if(slope_dict['trend_type'] == 'linregress'):
    if(plot_error):
        ax9 = axs[8]  # ocean stderr
        ax10 = axs[9]  # mix stderr
        ax11 = axs[10]  # ice stderr
        ax12 = axs[11]  # land stderr
 
    #calc_slopes = np.ma.masked_where(abs(slope_dict['slopes']) > 50, slope_dict['slopes'])
    calc_slopes = np.ma.masked_where(slope_dict['counts'] < min_ob, slope_dict['slopes'])
    calc_counts = np.ma.masked_where(slope_dict['counts'] < min_ob, slope_dict['counts'])
  
    print('count mins', np.min(calc_counts, axis = (1,2)))
    print('count maxs', np.max(calc_counts, axis = (1,2)))

 
    ax1.pcolormesh(calc_slopes[0,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    ax2.pcolormesh(calc_slopes[1,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    ax3.pcolormesh(calc_slopes[2,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    mesh = ax4.pcolormesh(calc_slopes[3,:,:].T, cmap = 'bwr', vmin = -15, vmax = 15)
    #cbar = plt.colorbar(mesh, ax = ax4, orientation = 'vertical', \
    #    shrink = 0.8, extend = 'both')
    ax5.pcolormesh(slope_dict['counts'][0,:,:].T, cmap = 'viridis', vmin = 50, vmax = 5500)
    ax6.pcolormesh(slope_dict['counts'][1,:,:].T, cmap = 'viridis', vmin = 50, vmax = 5500)
    ax7.pcolormesh(slope_dict['counts'][2,:,:].T, cmap = 'viridis', vmin = 50, vmax = 5500)
    mesh2 = ax8.pcolormesh(slope_dict['counts'][3,:,:].T, cmap = 'viridis', vmin = 50, vmax = 5500)

    #if(slope_dict['trend_type'] == 'linregress'):
    if(plot_error):
        cmap = 'bwr' 
        shrk = 1.0
        cmap2 = plt.get_cmap('jet')
        colorvals = np.arange(0, 6, 1)
        norm = mpl.colors.BoundaryNorm(colorvals, cmap2.N, extend = 'max')

        ax9.pcolormesh(slope_dict['slope_stderr'][0,:,:].T, norm = norm, shading = 'auto', cmap = 'jet')
        ax10.pcolormesh(slope_dict['slope_stderr'][1,:,:].T, norm =norm, shading = 'auto', cmap = 'jet')
        ax11.pcolormesh(slope_dict['slope_stderr'][2,:,:].T, norm =norm, shading = 'auto', cmap = 'jet')
        mesh3 = ax12.pcolormesh(slope_dict['slope_stderr'][3,:,:].T, cmap = 'jet', \
            norm = norm, shading = 'auto')
        #cbar = plt.colorbar(ScalarMappable(norm = norm, cmap = cmap2),\
        #    ax = ax12, orientation='vertical',shrink = shrk, extend = 'both')
        

    #fig.colorbar(mesh, ax = axs.ravel().tolist())

    sza_labels = [str(edge) for edge in bin_dict['sza_bin_edges']]
    cod_labels = [str(edge) for edge in bin_dict['cod_bin_edges']]
 
    #if(slope_dict['trend_type'] == 'theil-sen'): 
    if((bin_dict['sza_bin_edges'][1] - bin_dict['sza_bin_edges'][0]) <= 5):
        xtick_setter = np.arange(0, len(bin_dict['sza_bin_means']) + 1, 2)
        label_setter = sza_labels[::2]
    else:
        xtick_setter = np.arange(0, len(bin_dict['sza_bin_means']) + 1)
        label_setter = sza_labels[::1]
    if(not plot_error): 
        ax5.set_xticks(xtick_setter)
        ax5.set_xticklabels(label_setter)
    else:
        ax9.set_xticks(xtick_setter)
        ax9.set_xticklabels(label_setter)
    ax1.set_yticks([0,1,2,3,4,5,6,7,8])
    ax1.set_yticklabels(cod_labels[::1])
     
    ax1.set_title('Ocean')
    ax2.set_title('Mix')
    ax3.set_title('Ice')
    ax4.set_title('Land')
    
    #if(slope_dict['trend_type'] == 'theil-sen'): 
    if(not plot_error): 
        ax5.set_xlabel('OMI SZA [$^{o}$]')
        ax6.set_xlabel('OMI SZA [$^{o}$]')
        ax7.set_xlabel('OMI SZA [$^{o}$]')
        ax8.set_xlabel('OMI SZA [$^{o}$]')
    else:
        ax9.set_xlabel('OMI SZA [$^{o}$]')
        ax10.set_xlabel('OMI SZA [$^{o}$]')
        ax11.set_xlabel('OMI SZA [$^{o}$]')
        ax12.set_xlabel('OMI SZA [$^{o}$]')
        ax9.set_ylabel('MODIS COD')
    #ax2.set_xlabel('SZA dim')
    #ax3.set_xlabel('SZA dim')
    #ax4.set_xlabel('SZA dim')
    ax1.set_ylabel('MODIS COD')
    ax5.set_ylabel('MODIS COD')
    #ax2.set_ylabel('COD dim')
    #ax3.set_ylabel('COD dim')
    #ax4.set_ylabel('COD dim')

    ax1.set_title('Ocean\n(Ice % <= ' + str(int(bin_dict['ice_bin_edges'][1])) + '%)')
    ax2.set_title('Mix\n(' + str(int(bin_dict['ice_bin_edges'][1])) + \
        '% < Ice % < ' + str(int(bin_dict['ice_bin_edges'][2])) + '%)')
    ax3.set_title('Ice\n(Ice % >= ' + str(int(bin_dict['ice_bin_edges'][2])) + '%)')
    ax4.set_title('Land')
  
    axs[0].cax.colorbar(mesh, label = 'δForcing/δAI')
    axs[5].cax.colorbar(mesh2, label = 'Counts')
    #if(slope_dict['trend_type'] == 'linregress'):
    if(plot_error):
        axs[11].cax.colorbar(mesh3, label = 'Slope Stderr')
        
 
    if(slope_dict['trend_type'] == 'theil-sen'):
        plt.suptitle('Neural Net Forcing Slopes - Theil-Sen\nGrids with counts >= ' + str(int(min_ob)))
    else:
        plt.suptitle('Neural Net Forcing Slopes - Linear Regression\nGrids with counts >= ' + str(int(min_ob)))
 
    #fig.tight_layout()
    if(save): 
        if(slope_dict['trend_type'] == 'theil-sen'):
            type_adder = 'theilsen'
        else:
            type_adder = 'linregress'
        if(plot_error):
            error_add = '_error'
        else:
            error_add = ''
        outname = 'force_slopes_NN_' + type_adder + error_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image:", outname)
    else:
        plt.show()

# Plot the binned slopes for the 4 surface types
# ----------------------------------------------
def plot_NN_bin_slopes_codmean(slope_dict, bin_dict, min_ob = 50, save = False):

    plt.close('all')
    fig = plt.figure(figsize = (10, 3))
    ax1 = fig.add_subplot(1,4,1)
    ax2 = fig.add_subplot(1,4,2)
    ax3 = fig.add_subplot(1,4,3)
    ax4 = fig.add_subplot(1,4,4)

    #calc_slopes = np.ma.masked_where(abs(slope_dict['slopes']) > 50, slope_dict['slopes'])
    calc_slopes = np.ma.masked_where(slope_dict['counts'] < min_ob, slope_dict['slopes'])
   
    # Calculate the average and standard deviation along the y axis
    # -------------------------------------------------------------
    mean_slopes_as_func_cod = np.mean(calc_slopes, axis = 1)
    std_slopes_as_func_cod  = np.std(calc_slopes, axis = 1)

    ax1.plot( mean_slopes_as_func_cod[0])
    ax2.plot(mean_slopes_as_func_cod[1])
    ax3.plot(mean_slopes_as_func_cod[2])
    ax4.plot(mean_slopes_as_func_cod[3])
    ax1.plot( mean_slopes_as_func_cod[0] - std_slopes_as_func_cod[0], color = 'tab:blue', alpha = 0.5, linestyle = '--')
    ax2.plot(mean_slopes_as_func_cod[1] - std_slopes_as_func_cod[1], color = 'tab:blue', alpha = 0.5, linestyle = '--')
    ax3.plot(mean_slopes_as_func_cod[2] - std_slopes_as_func_cod[2], color = 'tab:blue', alpha = 0.5, linestyle = '--')
    ax4.plot(mean_slopes_as_func_cod[3] - std_slopes_as_func_cod[3], color = 'tab:blue', alpha = 0.5, linestyle = '--')
    ax1.plot( mean_slopes_as_func_cod[0] + std_slopes_as_func_cod[0], color = 'tab:blue', alpha = 0.5, linestyle = '--')
    ax2.plot(mean_slopes_as_func_cod[1] + std_slopes_as_func_cod[1], color = 'tab:blue', alpha = 0.5, linestyle = '--')
    ax3.plot(mean_slopes_as_func_cod[2] + std_slopes_as_func_cod[2], color = 'tab:blue', alpha = 0.5, linestyle = '--')
    ax4.plot(mean_slopes_as_func_cod[3] + std_slopes_as_func_cod[3], color = 'tab:blue', alpha = 0.5, linestyle = '--')
    ax1.set_ylim( -20, 20)
    ax2.set_ylim(-20, 20)
    ax3.set_ylim(-20, 20)
    ax4.set_ylim(-20, 20)
    ax1.axhline(0, color = 'k', linestyle = '--')
    ax2.axhline(0, color = 'k', linestyle = '--')
    ax3.axhline(0, color = 'k', linestyle = '--')
    ax4.axhline(0, color = 'k', linestyle = '--')

    ax1.set_title('Ocean')
    ax2.set_title('Mix')
    ax3.set_title('Ice')
    ax4.set_title('Land')

    ax1.set_xlabel('MODIS COD')
    ax2.set_xlabel('MODIS COD')
    ax3.set_xlabel('MODIS COD')
    ax4.set_xlabel('MODIS COD')
    ax1.set_ylabel('Forcing Eff. (δForcing/δAI)')   
 
    sza_labels = [str(edge) for edge in bin_dict['sza_bin_edges']]
    cod_labels = [str(edge) for edge in bin_dict['cod_bin_edges']]
  
    ax1.set_xticks([0,2,4,6,8])
    ax1.set_xticklabels(cod_labels[::2])
    ax2.set_xticks([0,2,4,6,8])
    ax2.set_xticklabels(cod_labels[::2])
    ax3.set_xticks([0,2,4,6,8])
    ax3.set_xticklabels(cod_labels[::2])
    ax4.set_xticks([0,2,4,6,8])
    ax4.set_xticklabels(cod_labels[::2])

    fig.tight_layout()
    if(save):
        outname = 'force_slopes_NN_szamean.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image:", outname)
    else:
        plt.show()

# Overlays high AI OMI areas on true color imagery and forcing stuff
def plot_compare_NN_output_overlay(calc_data, auto_zoom = True, save = False):
    
    # Load in the JSON file string relations
    # --------------------------------------
    with open(json_time_database, 'r') as fin:
        file_date_dict = json.load(fin)

    file_name = calc_data.strip().split('/')[-1] 
    date_str = file_name.split('.')[0].split('_')[-1]
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    
    sim_name = file_name.split('_')[3]
    print(sim_name, date_str)
    
    in_calc = h5py.File(calc_data)
    
    mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                   (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
    mask_lwf  = np.ma.masked_where((in_calc['ceres_lwf'][:,:] == -999.) | \
                                   (in_calc['ceres_lwf'][:,:] > 3000), in_calc['ceres_lwf'][:,:])
    mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
    mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)
    mask_AI   = np.ma.masked_invalid(in_calc['omi_uvai_pert'])
    
    ##!#mask_ctp = np.ma.masked_invalid(in_calc['modis_cld_top_pres'][:,:])
    ##!#mask_ctp = np.ma.masked_where(mask_ctp == -999., mask_ctp)
    ##!#
    ##!#mask_cod = np.ma.masked_invalid(in_calc['modis_cod'][:,:])
    ##!#mask_cod = np.ma.masked_where(mask_cod == -999., mask_cod)
    ##!#
    ##!#mask_ice = np.ma.masked_invalid(in_calc['nsidc_ice'][:,:])
    ##!#mask_ice = np.ma.masked_where(mask_ice == -999., mask_ice)
    
    diff_calc = mask_calc - mask_orig
    
    both_orig = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_orig)
    both_calc = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_calc)

    modis_date = file_date_dict[date_str]['MODIS'][0]
   
    if(auto_zoom):

        # Figure out where the pixels containing AI are located
        # -----------------------------------------------------
        high_AI_idxs = np.where(mask_AI > 1.5)

        # Switch the longitudes to 0 - 360 rather than -180 - 180
        # -------------------------------------------------------
        work_lons = np.copy(in_calc['omi_lon'][:,:])
        work_lons[work_lons < 0] = 360 + work_lons[work_lons < 0]
        

        # Figure out the minimum/maximum latitude/longitudes in the AI ranges
        # -------------------------------------------------------------------
        min_lat = np.min(in_calc['omi_lat'][:,:][high_AI_idxs]) - 2
        max_lat = np.max(in_calc['omi_lat'][:,:][high_AI_idxs]) + 2
        
        if(max_lat > 90):
            max_lat = 90

        # Calculate the 5th and 95th percentiles of longitude to prevent
        # outliers from ruining the mean or max/min values
        # --------------------------------------------------------------
        min_lon = np.percentile(work_lons[high_AI_idxs], 5) - 2
        max_lon = np.percentile(work_lons[high_AI_idxs], 95) + 2

        center_lon = (min_lon + max_lon) / 2.
        workcrs = ccrs.NorthPolarStereo(central_longitude = center_lon)

        #fig = plt.figure(figsize = (10.5, 6))
        fig = plt.figure(figsize = (11, 7))
        ax1 = fig.add_subplot(2,3,1, projection = workcrs) # OMI AI
        ax2 = fig.add_subplot(2,3,2, projection = workcrs) # True color
        ax3 = fig.add_subplot(2,3,3, projection = workcrs) # True color with OMI hashing
        ax6 = fig.add_subplot(2,3,4, projection = workcrs) # CERES lwf
        ax4 = fig.add_subplot(2,3,5, projection = workcrs) # Forcing
        ax5 = fig.add_subplot(2,3,6, projection = workcrs) # Forcing with OMI hashing
        
    else:
        workcrs = ccrs.NorthPolarStereo(central_longitude = 0)

        #fig = plt.figure(figsize = (10.5, 6))
        fig = plt.figure(figsize = (8, 7))
        ax5 = fig.add_subplot(2,2,1, projection = workcrs)
        ax1 = fig.add_subplot(2,2,2, projection = workcrs)
        ax2 = fig.add_subplot(2,2,3, projection = workcrs)
        ax3 = fig.add_subplot(2,2,4, projection = workcrs)
 
    #ax4 = fig.add_subplot(2,3,5)
    
    #ax4 = fig.add_subplot(2,2,4, projection = ccrs.NorthPolarStereo())
    
    mesh = ax1.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmax = 5, cmap = 'jet')
    cbar = fig.colorbar(mesh, ax = ax1, label = 'OMI UVAI')
    #ax1.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax1.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax1.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_title('OMI UVAI')
    ax1.coastlines()
    
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = True, ax = ax2)
    #mesh = ax2.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_calc, \
    #    transform = ccrs.PlateCarree(), shading = 'auto')
    #cbar = fig.colorbar(mesh, ax = ax2, label = 'SWF [Wm$^{-2}$]')
    #ax2.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax2.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax2.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.set_title('MODIS True Color')
    ax2.coastlines()
    
    #mesh = ax3.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], diff_calc, \
    #    transform = ccrs.PlateCarree(), shading = 'auto', vmin = -80, vmax = 80, cmap = 'bwr')
    #cbar = fig.colorbar(mesh, ax = ax3, label = 'SWF [Wm$^{-2}$]')
    #ax3.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    #ax3.set_boundary(circle, transform=ax3.transAxes)
    #ax3.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = True, ax = ax3)
    if(auto_zoom):
        ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax3.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.set_title('MODIS True Color (hashing)')
    ax3.coastlines()
    
    mesh = ax6.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_lwf, \
        transform = ccrs.PlateCarree(), vmin = None, vmax = None, shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax6, label = 'LWF [Wm$^{-2}$]')
    #ax5.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    #ax5.set_boundary(circle, transform=ax5.transAxes)
    #ax5.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax6.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax6.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax6.set_boundary(circle, transform=ax6.transAxes)
    ax6.set_title('Forcing')
    ax6.coastlines()
   
    mesh = ax4.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], diff_calc, \
        transform = ccrs.PlateCarree(), vmin = -80, vmax = 80, shading = 'auto', cmap = 'bwr')
    cbar = fig.colorbar(mesh, ax = ax4, label = 'Forcing [W$^{-2}$]')
    #ax5.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    #ax5.set_boundary(circle, transform=ax5.transAxes)
    #ax5.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax4.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax4.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax4.set_boundary(circle, transform=ax4.transAxes)
    ax4.set_title('Forcing')
    ax4.coastlines()
   
    mesh = ax5.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], diff_calc, \
        transform = ccrs.PlateCarree(), vmin = -80, vmax = 80, shading = 'auto', cmap = 'bwr')
    cbar = fig.colorbar(mesh, ax = ax5, label = 'Forcing [W$^{-2}$]')
    #ax5.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    #ax5.set_boundary(circle, transform=ax5.transAxes)
    #ax5.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax5.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax5.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax5.set_boundary(circle, transform=ax5.transAxes)
    ax5.set_title('Forcing (hashing)')
    ax5.coastlines()

    # Figure out the high AI indices
    # ------------------------------
    hatch_style = '///'
    if(auto_zoom):

        # Figure out where the pixels containing AI are located
        # -----------------------------------------------------
        hasher = np.ma.masked_invalid(mask_AI)
        hasher = np.ma.masked_where(hasher < 1.0, \
            hasher)

        print("MIN MASK AI:", np.min(hasher), 'MAX MASK AI:', np.nanmax(hasher))

        #hash0 = ax1.pcolor(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], \
        #    hasher, hatch = hatch_style, alpha=0., color = 'k', \
        #    shading = 'auto', transform = datacrs)
        #ax1.pcolor(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:],\
        #    hasher, hatch = 'xx', alpha=0., transform = datacrs, shading = 'auto')
        ax3.pcolor(work_lons, in_calc['omi_lat'][:,:],\
            hasher, hatch = 'xx', alpha=0., transform = datacrs, shading = 'auto')
        ax5.pcolor(work_lons, in_calc['omi_lat'][:,:],\
            hasher, hatch = 'xx', alpha=0., transform = datacrs, shading = 'auto')
        ax6.pcolor(work_lons, in_calc['omi_lat'][:,:],\
            hasher, hatch = 'xx', alpha=0., transform = datacrs, shading = 'auto')

    
  
    """ 
    if(auto_zoom): 
        work_force = np.abs(diff_calc)
        work_force = np.ma.masked_where(mask_AI < 1.0, work_force)
        mask_AI    = np.ma.masked_where(mask_AI < 1.0, mask_AI)  
        work_force = np.ma.masked_invalid(work_force)
        mask_AI    = np.ma.masked_invalid(mask_AI)
        work_force = np.ma.masked_where(mask_AI.mask == True, work_force)
        mask_AI    = np.ma.masked_where(work_force.mask == True, mask_AI)  


        #print("HERE:", work_force.compressed().shape, mask_AI.compressed().shape)

        if(len(mask_AI.compressed()) > 2):
            r2 = r2_score(mask_AI.compressed(), work_force.compressed())
            print('R2:', r2)
            ax4.set_title('r$^{2}$ = ' + str(np.round(r2, 3)))

        ax4.scatter(mask_AI, work_force, s = 6, color = 'k')
        ax4.set_xlabel('OMI UVAI Perturbation')
        ax4.set_ylabel('Forcing [W/m2]')
    """ 
 
    """ 
    xy = np.vstack([both_orig.compressed(), both_calc.compressed()])
    if(len(both_orig.compressed()) > 1):
        r2 = r2_score(both_orig.compressed(), both_calc.compressed())
        z = stats.gaussian_kde(xy)(xy)       
        ax4.scatter(both_orig.compressed(), both_calc.compressed(), c = z, s = 1)
        #ax4.set_title('r$^{2}$ = ' + str(np.round(r2, 3)))

        # Calculate RMSE
        # --------------
        rmse = mean_squared_error(both_orig.compressed(), both_calc.compressed(), squared = False)
        ptext = 'r$^{2}$ = ' + str(np.round(r2, 3)) + \
            '\nRMSE = ' + str(np.round(rmse, 1))
        #plot_figure_text(ax3, ptext, location = 'lower_right', \
        plot_figure_text(ax4, ptext, xval = 400, yval = 25, \
            fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')
    lims = [\
            np.min([ax4.get_xlim(), ax4.get_ylim()]),\
            np.max([ax4.get_xlim(), ax4.get_ylim()]),
    ]
    ax4.plot(lims, lims, 'r', linestyle = ':', alpha = 0.75)
    ax4.set_xlim(lims)
    ax4.set_ylim(lims)
    ax4.set_xlabel('Observed SWF [Wm$^{-2}$]')
    ax4.set_ylabel('Calculated SWF [Wm$^{-2}$]')
    """ 
   
    ##!#print(np.min(mask_cod), np.max(mask_cod))
     
    ##!#good_idxs = np.where( (mask_AI.mask == False) & \
    ##!#                      (mask_calc.mask == False) & \
    ##!#                      (mask_cod.mask == False) & \
    ##!#                      (mask_ice > 80.) & \
    ##!#                      (mask_ice < 101.) & \
    ##!#                      (mask_AI > 2))
    ##!#
    ##!#scat_ai   = mask_AI[good_idxs].flatten()
    ##!#scat_diff = diff_calc[good_idxs].flatten()
    ##!#ax6.scatter(scat_ai, scat_diff, c = 'k', s = 6)
    ##!#ax6.set_xlabel('OMI UVAI PERT')
    ##!#ax6.set_ylabel('Calc Aer Forcing [W/m2]')
    
    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax4, 'd)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax5, 'e)', fontsize = 11, backgroundcolor = None)
    #plot_subplot_label(ax4, 'e)', fontsize = 11, backgroundcolor = None)
    
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    
    #in_base.close()
    in_calc.close()
    
    fig.tight_layout()

    if(save):    
        if(auto_zoom):
            zoom_add = '_zoom'
        else:
            zoom_add = ''
        outname = 'ceres_nn_compare_hashing_' + date_str + '_' + sim_name + zoom_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


# Overlays high AI OMI areas on true color imagery and forcing stuff
def plot_compare_NN_output_overlay_v2(calc_data, auto_zoom = True, \
        label_xloc = None, label_yloc = None, save = False):
    
    # Load in the JSON file string relations
    # --------------------------------------
    with open(json_time_database, 'r') as fin:
        file_date_dict = json.load(fin)

    file_name = calc_data.strip().split('/')[-1] 
    date_str = file_name.split('.')[0].split('_')[-1]
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    
    sim_name = file_name.split('_')[3]
    print(sim_name, date_str)
    
    in_calc = h5py.File(calc_data)
    
    mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                   (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
    mask_lwf  = np.ma.masked_where((in_calc['ceres_lwf'][:,:] == -999.) | \
                                   (in_calc['ceres_lwf'][:,:] > 3000), in_calc['ceres_lwf'][:,:])
    mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
    mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)
    mask_AI   = np.ma.masked_invalid(in_calc['omi_uvai_pert'])
    
    ##!#mask_ctp = np.ma.masked_invalid(in_calc['modis_cld_top_pres'][:,:])
    ##!#mask_ctp = np.ma.masked_where(mask_ctp == -999., mask_ctp)
    ##!#
    ##!#mask_cod = np.ma.masked_invalid(in_calc['modis_cod'][:,:])
    ##!#mask_cod = np.ma.masked_where(mask_cod == -999., mask_cod)
    ##!#
    ##!#mask_ice = np.ma.masked_invalid(in_calc['nsidc_ice'][:,:])
    ##!#mask_ice = np.ma.masked_where(mask_ice == -999., mask_ice)
    

    modis_date = file_date_dict[date_str]['MODIS'][0]
   
    if(auto_zoom):

        # Figure out where the pixels containing AI are located
        # -----------------------------------------------------
        high_AI_idxs = np.where(mask_AI > 1.5)

        # Switch the longitudes to 0 - 360 rather than -180 - 180
        # -------------------------------------------------------
        work_lons = np.copy(in_calc['omi_lon'][:,:])
        work_lons[work_lons < 0] = 360 + work_lons[work_lons < 0]
        

        # Figure out the minimum/maximum latitude/longitudes in the AI ranges
        # -------------------------------------------------------------------
        min_lat = np.min(in_calc['omi_lat'][:,:][high_AI_idxs]) - 2
        max_lat = np.max(in_calc['omi_lat'][:,:][high_AI_idxs]) + 2
        
        if(max_lat > 90):
            max_lat = 90

        # Calculate the 5th and 95th percentiles of longitude to prevent
        # outliers from ruining the mean or max/min values
        # --------------------------------------------------------------
        min_lon = np.percentile(work_lons[high_AI_idxs], 5) - 2
        max_lon = np.percentile(work_lons[high_AI_idxs], 95) + 2

        center_lon = (min_lon + max_lon) / 2.
        workcrs = ccrs.NorthPolarStereo(central_longitude = center_lon)

        #fig = plt.figure(figsize = (10.5, 6))
        fig = plt.figure(figsize = (7.5, 4))
        ax1 = fig.add_subplot(1,3,1, projection = workcrs) # OMI AI
        ax2 = fig.add_subplot(1,3,2, projection = workcrs) # True color
        ax3 = fig.add_subplot(1,3,3, projection = workcrs) # True color
        
    else:
        workcrs = ccrs.NorthPolarStereo(central_longitude = 0)

        #fig = plt.figure(figsize = (10.5, 6))
        fig = plt.figure(figsize = (9, 4))
        ax1 = fig.add_subplot(1,2,1, projection = workcrs)
        ax2 = fig.add_subplot(1,2,2, projection = workcrs)
        ax3 = fig.add_subplot(1,2,3, projection = workcrs)
 
    #ax4 = fig.add_subplot(2,3,5)
    
    #ax4 = fig.add_subplot(2,2,4, projection = ccrs.NorthPolarStereo())
    
    mesh = ax1.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmax = 5, cmap = 'jet')
    cbar = fig.colorbar(mesh, ax = ax1, label = 'OMI UVAI')
    #ax1.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax1.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax1.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_title('OMI UVAI')
    ax1.coastlines()
    
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = True, ax = ax2)
    #mesh = ax2.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_calc, \
    #    transform = ccrs.PlateCarree(), shading = 'auto')
    #cbar = fig.colorbar(mesh, ax = ax2, label = 'SWF [Wm$^{-2}$]')
    #ax2.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax2.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax2.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.set_title('Aqua MODIS\nTrue Color')
    #ax2.set_title('Aqua MODIS\nHatch: OMI UVAI > 1.5')
    ax2.coastlines()
   
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = True, ax = ax3)
    #mesh = ax2.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_calc, \
    #    transform = ccrs.PlateCarree(), shading = 'auto')
    #cbar = fig.colorbar(mesh, ax = ax2, label = 'SWF [Wm$^{-2}$]')
    #ax2.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax3.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax3.set_boundary(circle, transform=ax2.transAxes)
    ax3.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmax = 5, cmap = 'jet', alpha = 0.1)
    cbar = fig.colorbar(mesh, ax = ax3, label = 'OMI UVAI')
    #ax3.set_title('Aqua MODIS\nTrue Color')
    ax3.set_title('Aqua MODIS\nOMI UVAI Overlay')
    ax3.coastlines()

 
    # Figure out the high AI indices
    # ------------------------------
    hatch_style = '///'
    #if(auto_zoom):

    #    # Figure out where the pixels containing AI are located
    #    # -----------------------------------------------------
    #    hasher = np.ma.masked_invalid(mask_AI)
    #    hasher = np.ma.masked_where(hasher < 1.5, \
    #        hasher)

    #    print("MIN MASK AI:", np.min(hasher), 'MAX MASK AI:', np.nanmax(hasher))

    #    ax2.pcolor(work_lons, in_calc['omi_lat'][:,:],\
    #        hasher, hatch = 'xx', alpha=0., transform = datacrs, shading = 'auto')

    font_size = 10 
    if((label_xloc is not None) & (label_yloc is not None)):
        plot_figure_text(ax1, 'a)', \
            xval = label_xloc, yval = label_yloc, transform = ccrs.PlateCarree(), \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax2, 'b)', \
            xval = label_xloc, yval = label_yloc, transform = ccrs.PlateCarree(), \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax3, 'c)', \
            xval = label_xloc, yval = label_yloc, transform = ccrs.PlateCarree(), \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
    else:
        plot_figure_text(ax1, 'a)', \
            xval = None, yval = None, transform = None, \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax2, 'b)', \
            xval = None, yval = None, transform = None, \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
        plot_figure_text(ax3, 'c)', \
            xval = None, yval = None, transform = None, \
            color = 'black', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'left', weight = 'bold')
    
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    
    #in_base.close()
    in_calc.close()
    
    fig.tight_layout()

    if(save):    
        if(auto_zoom):
            zoom_add = '_zoom'
        else:
            zoom_add = ''
        outname = 'ceres_nn_compare_hashing_' + date_str + '_' + sim_name + zoom_add + '_v2.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()





def plot_L2_validate_regress_all(sim_name, slope_dict, bin_dict, \
        ai_thresh = 0.7, mod_slopes = None, mod_intercepts = None, \
        mod_cod = None, mod_ice = None, use_intercept = False, \
        min_sza = None, max_sza = None, \
        min_ice = None, max_ice = None, \
        min_cod = None, max_cod = None, \
        return_values = False, save = False, save_values = False):

    # Retrieve all the aerosol NN output files
    # ----------------------------------------
    files = glob('neuralnet_output/test_calc_out_' + sim_name + '*.hdf5')

    # Set up arrays to hold all the calculated forcing and actual 
    # forcing values
    # ------------------------------------------------------------
    direct_forcings = np.full((110000), np.nan)
    calc_forcings   = np.full((110000), np.nan)

    beg_idx = 0
    end_idx = 0
    for calc_data in files:

        file_name = calc_data.strip().split('/')[-1] 
        date_str = file_name.split('.')[0].split('_')[-1]

        if(date_str == '201807082244'):
            continue

        dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
        
        sim_name = file_name.split('_')[3]
        print(sim_name, date_str)
        
        in_calc = h5py.File(calc_data)

        # Calculate forcing values estimated from the intercepts
        # ------------------------------------------------------
        testers = perform_forcing_calculation_v4(in_calc, slope_dict, bin_dict, \
                ai_thresh = ai_thresh, mod_slopes = mod_slopes, \
                mod_intercepts = mod_intercepts, mod_cod = mod_cod, \
                mod_ice = mod_ice, use_intercept = use_intercept)
        
        mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                       (in_calc['ceres_swf'][:,:] > 3000), \
                                        in_calc['ceres_swf'][:,:])
        mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
        mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)

        title_str = ''
        if(min_cod is not None):
            title_str = str(min_cod) + ' < COD'
            mask_calc = np.ma.masked_where(in_calc['modis_cod'][:,:] < min_cod, mask_calc)
        if(max_cod is not None):
            if(title_str != None):
                title_str = title_str + ' < ' + str(max_cod)
            else:
                title_str = 'COD < ' + str(max_cod) + '\n'
            mask_calc = np.ma.masked_where(in_calc['modis_cod'][:,:] > max_cod, mask_calc)
        if(min_sza is not None):
            mask_calc = np.ma.masked_where(in_calc['omi_sza'][:,:] < min_sza, mask_calc)
        if(max_sza is not None):
            mask_calc = np.ma.masked_where(in_calc['omi_sza'][:,:] > max_sza, mask_calc)
        if(min_ice is not None):
            mask_calc = np.ma.masked_where(in_calc['nsidc_ice'][:,:] < min_ice, mask_calc)
        if(max_ice is not None):
            mask_calc = np.ma.masked_where(in_calc['nsidc_ice'][:,:] > max_ice, mask_calc)
   
        testers = np.ma.masked_where(mask_calc.mask == True, testers)
 
        ##!#mask_ctp = np.ma.masked_invalid(in_calc['modis_cld_top_pres'][:,:])
        ##!#mask_ctp = np.ma.masked_where(mask_ctp == -999., mask_ctp)
        ##!#
        ##!#mask_cod = np.ma.masked_invalid(in_calc['modis_cod'][:,:])
        ##!#mask_cod = np.ma.masked_where(mask_cod == -999., mask_cod)
        ##!#
        ##!#mask_ice = np.ma.masked_invalid(in_calc['nsidc_ice'][:,:])
        ##!#mask_ice = np.ma.masked_where(mask_ice == -999., mask_ice)
        
        diff_calc = mask_calc - mask_orig
        
        mask_AI = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:]  == -999., \
            in_calc['omi_uvai_pert'][:,:])

        testers = np.ma.masked_where(testers == 0, testers)
        allmask_testers  = np.ma.masked_where((mask_orig.mask == True) | \
                                              (mask_calc.mask == True) | \
                                              (mask_AI.mask == True) | \
                                              (testers.mask == True), testers)
        allmask_diffcalc = np.ma.masked_where((mask_orig.mask == True) | \
                                              (mask_calc.mask == True) | \
                                              (mask_AI.mask == True) | \
                                              (testers.mask == True), diff_calc)

        in_calc.close()

        print(allmask_testers.compressed().shape, testers.compressed().shape)

        if( (allmask_testers.compressed().shape[0] == 0) | \
            (allmask_testers.compressed().shape[0] == 0)):
            print("ERROR: No good compressed values")
            print("       Not making scatter plot")
        else:
       
            end_idx = beg_idx + allmask_testers.compressed().shape[0]
    
            print("BEG IDX:", beg_idx, "END IDX:", end_idx)

            direct_forcings[beg_idx:end_idx] = allmask_diffcalc.compressed()
            calc_forcings[beg_idx:end_idx]   = allmask_testers.compressed()
           
            beg_idx = end_idx 

    direct_forcings = np.ma.masked_invalid(direct_forcings).compressed()
    calc_forcings   = np.ma.masked_invalid(calc_forcings).compressed()

    if(save_values):
        write_L2_L3_validation_values(direct_forcings, calc_forcings, \
            sim_name, ai_min, bin_dict)

    # Use the plotting function to generate a figure of the results
    # -------------------------------------------------------------
    #plot_scatter_hist_L2_L3_errors(direct_forcings, calc_forcings, \
    #    sim_name, save = save)

    if(return_values):
        return direct_forcings, calc_forcings

def plot_scatter_L2_L3_errors(direct_forcings, calc_forcings, sim_name, \
        ax = None, delta_calc = 20, save = False):
   
    # Calculate error statistics between the two
    # ------------------------------------------
    error_dict = calc_L2_L3_error_stats(direct_forcings, calc_forcings, \
        delta_calc = delta_calc)

    in_ax = True
    if(ax is None):
        plt.close('all')
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
 
    ax.scatter(direct_forcings, calc_forcings, s = 0.5, color = 'k')
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    ax.plot([-300, 300], [-300, 300], color = 'tab:blue')
    ax.axhline(0, linestyle = '--', color = 'grey', alpha = 0.5)
    ax.axvline(0, linestyle = '--', color = 'grey', alpha = 0.5)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_title('r$^{2}$ = ' + str(np.round(r2, 3)))
    ax.set_xlabel('L2 Forcing (NN - obs)')
    ax.set_ylabel('L3-style Forcing (AI-based)')
    
    # Use these error statistics to plot error bars
    # ---------------------------------------------
    ax.errorbar(error_dict['mean_direct_values'], error_dict['bin_centers'], \
        xerr = error_dict['stdev_errors'], fmt = '-o', markersize = 4.5, color = 'tab:red') 

    if(not in_ax):
        if(save):
            outname = 'force_L2L3_error_scatter_' + sim_name + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show()

def plot_hist_L2_L3_errors(direct_forcings, calc_forcings, sim_name, \
        xmin = None, xmax = None, ax = None, num_bins = 75, \
        use_correct_error_calc = False, \
        astrofit = False, save = False):

    # direct = L2
    # calc   = L3
    if(use_correct_error_calc):
        print("Calculating L2/L3 errors using: calc (L3) - direct (L2)")
        xlabel_text = 'Forcing error (L3 - L2) [Wm$^{-2}$]'
        errors = calc_forcings - direct_forcings
    else:
        print("Calculating L2/L3 errors using: direct (L2) - calc (L3)")
        xlabel_text = 'Forcing error (L2 - L3) [Wm$^{-2}$]'
        errors = direct_forcings - calc_forcings

    mean_error = np.mean(errors)
    stdev_error  = np.std(errors)
    median_error = np.median(errors)
    mode_error   = np.median(errors)

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    ax.hist(errors, bins = num_bins)
    #ax.axvline( mean_error + stdev_error, color = 'gray', linestyle = '--')
    #ax.axvline( mean_error - stdev_error, color = 'gray', linestyle = '--')
    ax.axvline(0, linestyle = '--', color = 'k', alpha = 0.5)
    if( (xmin is None) & (xmax is None) ):
        ax.set_xlim(-150, 150)
    else:
        ax.set_xlim(xmin, xmax)
    ax.set_xlabel(xlabel_text)
    ax.set_ylabel('Counts')
    ax.set_title('Errors between L2 forcing (NN-based) and\n' + \
        'L3-style forcing (AI & LUT) applied to L2 data')

    # Test using astropy to fit a Gaussian function
    # ---------------------------------------------
    if(astrofit):
        bin_heights,bin_borders = np.histogram(errors,bins=num_bins)
        bin_widths = np.diff(bin_borders)
        bin_centers = bin_borders[:-1] + bin_widths / 2
        t_init = models.Gaussian1D()
        fit_t = fitting.LevMarLSQFitter()
        t = fit_t(t_init, bin_centers, bin_heights)
        print('Amplitude: ',np.round(t.amplitude.value,3))
        print('Mean:      ',np.round(t.mean.value,3))
        print('StDev:     ',np.round(t.stddev.value,3))
        x_interval_for_fit = np.linspace(bin_centers[0],bin_centers[-1],200)
        ax.plot(x_interval_for_fit,t(x_interval_for_fit),label='fit',c='tab:red')
        xdiff = ax.get_xlim()[1] - ax.get_xlim()[0]
        plot_xval = ax.get_xlim()[0]  + xdiff * 0.65
        plot_yval = ax.get_ylim()[1] * 0.8
        ptext = 'μ = ' + str(np.round(t.mean.value, 1)) + ' Wm$^{-2}$\n' + \
            'σ = ' + str(np.round(t.stddev.value, 1)) + ' Wm$^{-2}$'
        plot_figure_text(ax, ptext, xval = plot_xval, yval = plot_yval, \
            fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')

    if(not in_ax):
        if(save):
            fig.tight_layout()
            if(use_correct_error_calc):
                err_calc_add = '_correct'
            else:
                err_calc_add = ''
            outname = 'force_L2L3_error_hist_' + sim_name + err_calc_add + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show()



# delta_calc = the width of the calculated forcing bins 
def plot_scatter_hist_L2_L3_errors(direct_forcings, calc_forcings, 
            sim_name, num_bins = 75, delta_calc = 20, astrofit = False, screen_outliers = False, \
            save = False):

    print("\nCalculating r2 score")
    r2 = r2_score(direct_forcings, calc_forcings)
    print('R2:', r2)

    # Calculate the individual errors between the forcings
    # ----------------------------------------------------
    errors = direct_forcings - calc_forcings

    # Check if the user wants to screen outliers (with forcings 
    # above or below 100?)
    # ---------------------------------------------------------
    if(screen_outliers):
        print("MASKING ERRORS FOR FORCINGS WITH MAGNITUDES ABOVE 100")
        errors = np.ma.masked_where(abs(calc_forcings) > 100, errors).compressed()

    # Calculate the mean & standard deviation of these errors
    # -------------------------------------------------------
    mean_error = np.mean(errors)
    stdev_error  = np.std(errors)
    median_error = np.median(errors)
    mode_error   = np.median(errors)
    
    print("Total mean error between L2 and L3 values:", mean_error)
    print("Total stdev error between L2 and L3 values:", stdev_error)
    print("Total median error between L2 and L3 values:", median_error)
    print("Total mode error between L2 and L3 values:", mode_error)
    print('\n')

    # Calculate RMSE
    # --------------
    rmse = mean_squared_error(direct_forcings, calc_forcings, squared = False)
    #ax3.set_title('r$^{2}$ = ' + str(np.round(r2, 3)) + \
    #    'RMSE = ' + str(np.round(rmse, 1)))
    ptext = 'r$^{2}$ = ' + str(np.round(r2, 3)) + \
        '\nRMSE = ' + str(np.round(rmse, 1))

    plt.close('all')
    fig = plt.figure(figsize = (9, 4))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Plot the scatter of the direct (L2) vs calc (L3) forcings
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    plot_scatter_L2_L3_errors(direct_forcings, calc_forcings, sim_name, \
        ax = ax1, delta_calc = delta_calc, save = False)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Plot a histogram of the errors
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    plot_hist_L2_L3_errors(direct_forcings, calc_forcings, sim_name, \
        ax = ax2, num_bins = num_bins, save = False)

    fig.tight_layout()
    if(save):
        if(astrofit):
            outname = 'validation_L2_allscatter_histogram_astrofit_' + sim_name + '_numsfcbins6.png'
        else:
            outname = 'validation_L2_allscatter_histogram_' + sim_name + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def write_L2_L3_validation_values(direct_forcings, calc_forcings, \
        sim_name, ai_min, bin_dict, name_add = ''):

    outfile = 'validate_values_' + sim_name + name_add + '.hdf5'
 
    dset = h5py.File(outfile,'w')

    cdt = dset.create_dataset('direct_forcings', data = direct_forcings)
    cdt = dset.create_dataset('calc_forcings', data = calc_forcings)
    #if(maxerr is not None):
    #    cdt.attrs['maxerr'] = str(maxerr)
    #if(ai_thresh is not None):
    #    cdt.attrs['ai_thresh'] = str(ai_thresh)
    #if(minlat is not None):
    #    cdt.attrs['minlat'] = str(minlat)
    #if(dtype is not None):
    #    cdt.attrs['omi_data_type'] = dtype

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  

# This function loops over bins along the calculated forcing
# range and determines the mean, standard deviation, and 
# count of the errors between the direct forcing (L2-style) and
# the calculated forcing (L3-style) in each of the bins.
# -------------------------------------------------------------
def calc_L2_L3_error_stats(direct_forcings, calc_forcings, delta_calc = 20):

    bin_edges = np.arange(-150, 112, delta_calc)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2

    mean_errors   = np.full(bin_centers.shape[0], np.nan)
    stdev_errors  = np.full(bin_centers.shape[0], np.nan)
    stderr_errors = np.full(bin_centers.shape[0], np.nan)
    count_errors  = np.full(bin_centers.shape[0], np.nan)

    # This array will hold the average direct forcing value
    # in each of the calc forcing bins
    mean_direct_values = np.full(bin_centers.shape[0], np.nan)

    for ii in range(bin_centers.shape[0]):
        # Find the indices of all the calculated forcing values
        # within this range. This looks at a section of the y axis
        # of the scatter plot and grabs all those values in the
        # current range
        # --------------------------------------------------------
        match_idxs = np.where( (calc_forcings > bin_edges[ii]) & \
            (calc_forcings <= bin_edges[ii + 1]))

        # Calculate the errors between the calculated forcing values
        # and the corresponding direct forcing values
        # ----------------------------------------------------------
        errors = direct_forcings[match_idxs] - calc_forcings[match_idxs]

        # Calculate the mean direct forcing value for this
        # calculated forcing bin
        # ------------------------------------------------
        mean_direct_values[ii] = np.mean(direct_forcings[match_idxs])

        mean_errors[ii]   = np.mean(errors)
        stdev_errors[ii]  = np.std(errors)
        stderr_errors[ii] = sem(errors)
        count_errors[ii]  = errors.shape[0]

        print(bin_edges[ii], bin_centers[ii], bin_edges[ii + 1], \
            np.round(mean_errors[ii], 1), np.round(stdev_errors[ii], 1), \
            np.round(stderr_errors[ii], 1))

    out_dict = {}
    out_dict['bin_edges'] = bin_edges
    out_dict['bin_centers'] = bin_centers
    out_dict['mean_direct_values'] = mean_direct_values
    out_dict['mean_errors'] = mean_errors
    out_dict['stdev_errors']  = stdev_errors
    out_dict['stderr_errors']  = stderr_errors
    out_dict['count_errors'] = count_errors

    return out_dict
    

# Overlays high AI OMI areas on true color imagery and forcing stuff
def plot_compare_NN_output_L2_validate(calc_data, slope_dict, bin_dict, \
        mod_slopes = None, mod_intercepts = None, mod_cod = None, \
        mod_ice = None, use_intercept = False, auto_zoom = True, \
        ai_thresh = 0.7, save = False):
    
    # Load in the JSON file string relations
    # --------------------------------------
    with open(json_time_database, 'r') as fin:
        file_date_dict = json.load(fin)

    file_name = calc_data.strip().split('/')[-1] 
    date_str = file_name.split('.')[0].split('_')[-1]
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    
    sim_name = file_name.split('_')[3]
    print(sim_name, date_str)
    
    in_calc = h5py.File(calc_data)

    # Calculate forcing values estimated from the intercepts
    # ------------------------------------------------------
    testers = perform_forcing_calculation_v4(in_calc, slope_dict, bin_dict, \
            ai_thresh = ai_thresh, mod_slopes = mod_slopes, \
            mod_intercepts = mod_intercepts, mod_cod = mod_cod, \
            mod_ice = mod_ice, use_intercept = use_intercept)
    
    mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                   (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
    mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
    mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)
   
    testers = np.ma.masked_where(mask_calc.mask == True, testers)
 
    ##!#mask_ctp = np.ma.masked_invalid(in_calc['modis_cld_top_pres'][:,:])
    ##!#mask_ctp = np.ma.masked_where(mask_ctp == -999., mask_ctp)
    ##!#
    ##!#mask_cod = np.ma.masked_invalid(in_calc['modis_cod'][:,:])
    ##!#mask_cod = np.ma.masked_where(mask_cod == -999., mask_cod)
    ##!#
    ##!#mask_ice = np.ma.masked_invalid(in_calc['nsidc_ice'][:,:])
    ##!#mask_ice = np.ma.masked_where(mask_ice == -999., mask_ice)
    
    diff_calc = mask_calc - mask_orig
    
    both_orig = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_orig)
    both_calc = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_calc)

    modis_date = file_date_dict[date_str]['MODIS'][0]
   
    mask_AI = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:]  == -999., \
        in_calc['omi_uvai_pert'][:,:])

    if(auto_zoom):

        # Figure out where the pixels containing AI are located
        # -----------------------------------------------------
        high_AI_idxs = np.where(mask_AI > 1.5)
        #mask_AI = np.ma.masked_where(mask_AI < 1.5, mask_AI)

        # Switch the longitudes to 0 - 360 rather than -180 - 180
        # -------------------------------------------------------
        work_lons = np.copy(in_calc['omi_lon'][:,:])
        work_lons[work_lons < 0] = 360 + work_lons[work_lons < 0]
        

        # Figure out the minimum/maximum latitude/longitudes in the AI ranges
        # -------------------------------------------------------------------
        min_lat = np.min(in_calc['omi_lat'][:,:][high_AI_idxs]) - 2
        max_lat = np.max(in_calc['omi_lat'][:,:][high_AI_idxs]) + 2
        
        if(max_lat > 90):
            max_lat = 90

        # Calculate the 5th and 95th percentiles of longitude to prevent
        # outliers from ruining the mean or max/min values
        # --------------------------------------------------------------
        min_lon = np.percentile(work_lons[high_AI_idxs], 5) - 2
        max_lon = np.percentile(work_lons[high_AI_idxs], 95) + 2

        center_lon = (min_lon + max_lon) / 2.
        workcrs = ccrs.NorthPolarStereo(central_longitude = center_lon)
        
    else:
        workcrs = ccrs.NorthPolarStereo(central_longitude = 0)

    plt.close('all')
    fig = plt.figure(figsize = (11, 7))
    ax1 = fig.add_subplot(2,3,1, projection = workcrs) # OMI AI
    ax2 = fig.add_subplot(2,3,2, projection = workcrs) # True color
    ax3 = fig.add_subplot(2,3,4, projection = workcrs) # Forcing direct calc
    ax4 = fig.add_subplot(2,3,5, projection = workcrs) # Forcing from bin slopes
    ax5 = fig.add_subplot(2,3,6) # Regression of 2 forcings?


    # Plot OMI UVAI
    # ------------- 
    mesh = ax1.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmax = 5, cmap = 'jet')
    cbar = fig.colorbar(mesh, ax = ax1, label = 'OMI UVAI')
    if(auto_zoom):
        ax1.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax1.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_title('OMI UVAI')
    ax1.coastlines()
    
    # Plot MODIS True Color
    # ---------------------
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = True, ax = ax2)
    if(auto_zoom):
        ax2.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax2.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.set_title('MODIS True Color')
    ax2.coastlines()
   
    # Plot directly-calculated forcing (NN - obs)
    # -------------------------------------------
    mesh = ax3.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], diff_calc, \
        transform = ccrs.PlateCarree(), vmin = -80, vmax = 80, shading = 'auto', cmap = 'bwr')
    cbar = fig.colorbar(mesh, ax = ax3, label = 'Forcing [W$^{-2}$]')
    if(auto_zoom):
        ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax3.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.set_title('Forcing')
    ax3.coastlines()
  
    # Plot forcing calculated from the forcing slopes and intercepts
    # -------------------------------------------------------------- 
    mesh = ax4.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], testers, \
        transform = ccrs.PlateCarree(), vmin = -80, vmax = 80, shading = 'auto', cmap = 'bwr')
    cbar = fig.colorbar(mesh, ax = ax4, label = 'Forcing [W$^{-2}$]')
    if(auto_zoom):
        ax4.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax4.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax4.set_boundary(circle, transform=ax3.transAxes)
    ax4.set_title('Forcing Validated')
    ax4.coastlines()


    testers = np.ma.masked_where(testers == 0, testers)
    allmask_testers  = np.ma.masked_where((mask_orig.mask == True) | \
                                         (mask_calc.mask == True) | \
                                         (testers.mask == True), testers)
    allmask_diffcalc = np.ma.masked_where((mask_orig.mask == True) | \
                                          (mask_calc.mask == True) | \
                                          (testers.mask == True), diff_calc)
    #both_calc = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_calc)

 
    #diff_calc = np.ma.masked_where(testers.mask == True, diff_calc) 
    print(allmask_testers.compressed().shape, testers.compressed().shape)
    if( (allmask_testers.compressed().shape[0] == 0) | \
        (allmask_testers.compressed().shape[0] == 0)):
        print("ERROR: No good compressed values")
        print("       Not making scatter plot")
    else:
        r2 = r2_score(allmask_diffcalc.compressed(), allmask_testers.compressed())
        print('R2:', r2)
        ax5.set_title('r$^{2}$ = ' + str(np.round(r2, 3)))
        ax5.scatter(allmask_diffcalc.compressed(), allmask_testers.compressed(), s = 6, color = 'k')

 
    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax4, 'd)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax5, 'e)', fontsize = 11, backgroundcolor = None)
    
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    
    in_calc.close()
    
    fig.tight_layout()

    if(save):    
        if(auto_zoom):
            zoom_add = '_zoom'
        else:
            zoom_add = ''
        outname = 'ceres_nn_L2_validate_' + date_str + '_' + sim_name + zoom_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


# Zooms in the images on just the plume area
def plot_compare_NN_output_v2(calc_data, auto_zoom = True, save = False):
    
    # Load in the JSON file string relations
    # --------------------------------------
    with open(json_time_database, 'r') as fin:
        file_date_dict = json.load(fin)

    file_name = calc_data.strip().split('/')[-1] 
    date_str = file_name.split('.')[0].split('_')[-1]
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    
    sim_name = file_name.split('_')[3]
    print(sim_name, date_str)
    
    in_calc = h5py.File(calc_data)
    
    mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                   (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
    mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
    mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)
    
    if('ceres_lwf' in in_calc.keys()):
        l_contains_lwf = True
        mask_lwf  = np.ma.masked_where((in_calc['ceres_lwf'][:,:] == -999.) | \
                                       (in_calc['ceres_lwf'][:,:] > 3000), in_calc['ceres_lwf'][:,:])
    else:
        l_contains_lwf = False
    
    ##!#mask_ctp = np.ma.masked_invalid(in_calc['modis_cld_top_pres'][:,:])
    ##!#mask_ctp = np.ma.masked_where(mask_ctp == -999., mask_ctp)
    ##!#
    ##!#mask_cod = np.ma.masked_invalid(in_calc['modis_cod'][:,:])
    ##!#mask_cod = np.ma.masked_where(mask_cod == -999., mask_cod)
    ##!#
    ##!#mask_ice = np.ma.masked_invalid(in_calc['nsidc_ice'][:,:])
    ##!#mask_ice = np.ma.masked_where(mask_ice == -999., mask_ice)
    
    diff_calc = mask_calc - mask_orig
    
    both_orig = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_orig)
    both_calc = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_calc)

    modis_date = file_date_dict[date_str]['MODIS'][0]

    mask_AI = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:]  == -999., \
        in_calc['omi_uvai_pert'][:,:])
   
    if(auto_zoom):

        # Figure out where the pixels containing AI are located
        # -----------------------------------------------------
        high_AI_idxs = np.where(mask_AI > 1.5)
        #mask_AI = np.ma.masked_where(mask_AI < 1.5, mask_AI)

        # Switch the longitudes to 0 - 360 rather than -180 - 180
        # -------------------------------------------------------
        work_lons = np.copy(in_calc['omi_lon'][:,:])
        work_lons[work_lons < 0] = 360 + work_lons[work_lons < 0]
        

        # Figure out the minimum/maximum latitude/longitudes in the AI ranges
        # -------------------------------------------------------------------
        min_lat = np.min(in_calc['omi_lat'][:,:][high_AI_idxs]) - 2
        max_lat = np.max(in_calc['omi_lat'][:,:][high_AI_idxs]) + 2
        
        if(max_lat > 90):
            max_lat = 90

        # Calculate the 5th and 95th percentiles of longitude to prevent
        # outliers from ruining the mean or max/min values
        # --------------------------------------------------------------
        min_lon = np.percentile(work_lons[high_AI_idxs], 5) - 2
        max_lon = np.percentile(work_lons[high_AI_idxs], 95) + 2

        center_lon = (min_lon + max_lon) / 2.
        workcrs = ccrs.NorthPolarStereo(central_longitude = center_lon)

    else:
        workcrs = ccrs.NorthPolarStereo(central_longitude = 0)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # NEW: Print out the mean and standard deviation forcing value within
    #       the plume area
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    keep_idxs = np.where(mask_AI > 1.5)
    good_diff_calc = diff_calc[keep_idxs]
    mean_force = np.mean(good_diff_calc)
    max_force  = np.max(good_diff_calc)
    min_force  = np.min(good_diff_calc)
    std_force  = np.std(good_diff_calc) 
    print("\nFORCING STATS")
    print("    Mean forcing val: ", np.round(mean_force, 1))
    print("    StDv forcing val: ", np.round(std_force, 1))
    print("    Max  forcing val: ", np.round(max_force, 1))
    print("    Min  forcing val: ", np.round(min_force, 1))
    

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Make the figure
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #fig = plt.figure(figsize = (10.5, 6))
    fig = plt.figure(figsize = (11, 7))
    ax1 = fig.add_subplot(2,3,1, projection = workcrs) # OMI AI
    ax2 = fig.add_subplot(2,3,2, projection = workcrs) # True color
    if(l_contains_lwf):
        ax6 = fig.add_subplot(2,3,3, projection = workcrs) # CERES LWF obs
    ax3 = fig.add_subplot(2,3,4, projection = workcrs) # CERES SWF obs
    ax4 = fig.add_subplot(2,3,5, projection = workcrs) # NN SWF 
    ax5 = fig.add_subplot(2,3,6, projection = workcrs) # Forcing with OMI hashing
        
 
   
    # Plot OMI UVAI
    # ------------- 
    mesh = ax1.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmax = 5, cmap = 'jet')
    cbar = fig.colorbar(mesh, ax = ax1, label = 'OMI UVAI')
    if(auto_zoom):
        ax1.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax1.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_title('OMI UVAI')
    ax1.coastlines()
   
    # Plot MODIS True color
    # --------------------- 
    plot_MODIS_channel(modis_date, 'true_color', swath = True, \
        zoom = True, ax = ax2)
    if(auto_zoom):
        ax2.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax2.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.set_title('MODIS True Color')
    ax2.coastlines()
  
    max_swf = np.max([np.max(mask_orig), np.max(mask_calc)])
    min_swf = np.min([np.min(mask_orig), np.min(mask_calc)])
 
    if(l_contains_lwf):
        # Plot CERES LWF observations
        # --------------------------- 
        mesh = ax6.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_lwf, \
            transform = ccrs.PlateCarree(), shading = 'auto', vmin = None, vmax = None)
        cbar = fig.colorbar(mesh, ax = ax6, label = 'LWF [Wm$^{-2}$]')
        #ax1.set_extent([117, 170,70,  85], ccrs.PlateCarree())
        if(auto_zoom):
            ax6.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
        else:
            ax6.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
            ax6.set_boundary(circle, transform=ax6.transAxes)
        ax6.set_title('Observed LWF')
        ax6.coastlines()
   
    # Plot CERES SWF observations
    # --------------------------- 
    mesh = ax3.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_orig, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = min_swf, vmax = max_swf)
    cbar = fig.colorbar(mesh, ax = ax3, label = 'SWF [Wm$^{-2}$]')
    #ax1.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax3.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.set_title('Observed SWF')
    ax3.coastlines()
   
    # Plot neural net output
    # ---------------------- 
    mesh = ax4.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_calc, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = min_swf, vmax = max_swf)
        #transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax4, label = 'SWF [Wm$^{-2}$]')
    #ax2.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax4.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax4.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax4.set_boundary(circle, transform=ax4.transAxes)
    ax4.set_title('Calc. Aerosol-free SWF')
    ax4.coastlines()
   
    # Plot the estimated forcing
    # -------------------------- 
    mesh = ax5.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], diff_calc, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = -80, vmax = 80, cmap = 'bwr')
    cbar = fig.colorbar(mesh, ax = ax5, label = 'SWF [Wm$^{-2}$]')
    #ax3.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    #ax3.set_boundary(circle, transform=ax3.transAxes)
    #ax3.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax5.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax5.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax5.set_boundary(circle, transform=ax5.transAxes)
    ax5.set_title('Calculated - Observed')
    ax5.coastlines()
    
    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax4, 'd)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax5, 'e)', fontsize = 11, backgroundcolor = None)
    #plot_subplot_label(ax4, 'e)', fontsize = 11, backgroundcolor = None)
    
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    
    #in_base.close()
    in_calc.close()
    
    fig.tight_layout()

    if(save):    
        if(auto_zoom):
            zoom_add = '_zoom'
        else:
            zoom_add = ''
        outname = 'ceres_nn_compare_v2_' + date_str + '_' + sim_name + zoom_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

# Zooms in the images on just the plume area
def plot_compare_NN_output_double(infile1, infile2, auto_zoom = False, \
        save = False, include_scatter = False):
   
    file_name1 = infile1.strip().split('/')[-1] 
    file_name2 = infile2.strip().split('/')[-1] 
    date_str1 = file_name1.split('.')[0].split('_')[-1]
    date_str2 = file_name2.split('.')[0].split('_')[-1]
    dt_date_str1 = datetime.strptime(date_str1, '%Y%m%d%H%M')
    dt_date_str2 = datetime.strptime(date_str2, '%Y%m%d%H%M')
    
    sim_name1 = file_name1.split('_')[3]
    sim_name2 = file_name2.split('_')[3]
    if(sim_name1 != sim_name2):
        print("\nWARNING: SIM MISMATCH")
        print("\n    sim1 = ", sim_name1)
        print("    sim2 = ", sim_name2)

    print(sim_name1, sim_name2, date_str1, date_str2)
 

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    #
    # Process and plot output for the first file
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
 
    in_calc = h5py.File(infile1, 'r')
    mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                   (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
    mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
    mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)

    # Calculate forcing
    # -----------------
    diff_calc = mask_calc - mask_orig
    
    mask_AI = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:]  == -999., \
        in_calc['omi_uvai_pert'][:,:])
  
    if(auto_zoom):

        # Figure out where the pixels containing AI are located
        # -----------------------------------------------------
        high_AI_idxs = np.where(mask_AI > 1.5)
        #mask_AI = np.ma.masked_where(mask_AI < 1.5, mask_AI)

        # Switch the longitudes to 0 - 360 rather than -180 - 180
        # -------------------------------------------------------
        work_lons = np.copy(in_calc['omi_lon'][:,:])
        work_lons[work_lons < 0] = 360 + work_lons[work_lons < 0]
        

        # Figure out the minimum/maximum latitude/longitudes in the AI ranges
        # -------------------------------------------------------------------
        min_lat = np.min(in_calc['omi_lat'][:,:][high_AI_idxs]) - 2
        max_lat = np.max(in_calc['omi_lat'][:,:][high_AI_idxs]) + 2
        
        if(max_lat > 90):
            max_lat = 90

        # Calculate the 5th and 95th percentiles of longitude to prevent
        # outliers from ruining the mean or max/min values
        # --------------------------------------------------------------
        min_lon = np.percentile(work_lons[high_AI_idxs], 5) - 2
        max_lon = np.percentile(work_lons[high_AI_idxs], 95) + 2

        center_lon = (min_lon + max_lon) / 2.
        workcrs = ccrs.NorthPolarStereo(central_longitude = center_lon)

    else:
        workcrs = ccrs.NorthPolarStereo(central_longitude = 0)

    max_swf = np.max([np.max(mask_orig), np.max(mask_calc)])
    min_swf = np.min([np.min(mask_orig), np.min(mask_calc)])

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Make the figure
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #fig = plt.figure(figsize = (10.5, 6))
    if(include_scatter):
        #fig = plt.figure(figsize = (17.0, 7.5))
        #ax1 = fig.add_subplot(2,5,1, projection = workcrs) # date1 - OMI AI - date1
        #ax2 = fig.add_subplot(2,5,2, projection = workcrs) # date1 - CERES SWF
        #ax3 = fig.add_subplot(2,5,3, projection = workcrs) # date1 - NN SWF
        #ax4 = fig.add_subplot(2,5,4, projection = workcrs) # date1 - NN - CERES
        #ax5 = fig.add_subplot(2,5,6, projection = workcrs) # date2 - OMI AI - date1
        #ax6 = fig.add_subplot(2,5,7, projection = workcrs) # date2 - CERES SWF
        #ax7 = fig.add_subplot(2,5,8, projection = workcrs) # date2 - NN SWF
        #ax8 = fig.add_subplot(2,5,9, projection = workcrs) # date2 - NN - CERES
        #ax10= fig.add_subplot(2,5,10)
        #ax9 = fig.add_subplot(2,5,5, sharex = ax10)


        #fig = plt.figure(figsize = (14, 11))
        #ax1 = fig.add_subplot(3,4,1, projection = workcrs) # date1 - OMI AI - date1
        #ax2 = fig.add_subplot(3,4,2, projection = workcrs) # date1 - CERES SWF
        #ax3 = fig.add_subplot(3,4,3, projection = workcrs) # date1 - NN SWF
        #ax4 = fig.add_subplot(3,4,4, projection = workcrs) # date1 - NN - CERES
        #ax5 = fig.add_subplot(3,4,5, projection = workcrs) # date2 - OMI AI - date1
        #ax6 = fig.add_subplot(3,4,6, projection = workcrs) # date2 - CERES SWF
        #ax7 = fig.add_subplot(3,4,7, projection = workcrs) # date2 - NN SWF
        #ax8 = fig.add_subplot(3,4,8, projection = workcrs) # date2 - NN - CERES
        #ax10= fig.add_subplot(3,4,10)
        #ax9 = fig.add_subplot(3,4,9, sharey = ax10)

        fig = plt.figure(figsize = (10, 10), constrained_layout = True)
        subfigs = fig.subfigures(nrows = 2, ncols = 1, wspace = 0.01)
      
        ax1 = subfigs[0].add_subplot(2,4,1, projection = ccrs.NorthPolarStereo()) 
        ax2 = subfigs[0].add_subplot(2,4,2, projection = ccrs.NorthPolarStereo()) 
        ax3 = subfigs[0].add_subplot(2,4,3, projection = ccrs.NorthPolarStereo()) 
        ax4 = subfigs[0].add_subplot(2,4,4, projection = ccrs.NorthPolarStereo()) 
        ax5 = subfigs[0].add_subplot(2,4,5, projection = ccrs.NorthPolarStereo()) 
        ax6 = subfigs[0].add_subplot(2,4,6, projection = ccrs.NorthPolarStereo()) 
        ax7 = subfigs[0].add_subplot(2,4,7, projection = ccrs.NorthPolarStereo()) 
        ax8 = subfigs[0].add_subplot(2,4,8, projection = ccrs.NorthPolarStereo()) 

        axs2  = subfigs[1].subplots(nrows = 1, ncols = 2, sharex = True, sharey = True)
        ax9 = axs2[0]
        ax10 = axs2[1]
        #ax9  = subfigs[1].add_subplot(2,1,1)
        #ax10 = subfigs[1].add_subplot(2,1,2)

        #axs1 = [subfigs[0].add_subplot(6,1,ii, projection = ccrs.NorthPolarStereo()) for ii in range(1,7)]
        #axs2 = [subfigs[1].add_subplot(6,1,ii, projection = ccrs.NorthPolarStereo()) for ii in range(1,7)]
        #axs3 = [subfigs[2].add_subplot(6,1,ii, projection = ccrs.NorthPolarStereo()) for ii in range(1,7)]


    else:
        fig = plt.figure(figsize = (12.0, 7.5))
        #fig = plt.figure(figsize = (13, 10))
        ax1 = fig.add_subplot(2,4,1, projection = workcrs) # date1 - OMI AI - date1
        ax2 = fig.add_subplot(2,4,2, projection = workcrs) # date1 - CERES SWF
        ax3 = fig.add_subplot(2,4,3, projection = workcrs) # date1 - NN SWF
        ax4 = fig.add_subplot(2,4,4, projection = workcrs) # date1 - NN - CERES
        ax5 = fig.add_subplot(2,4,5, projection = workcrs) # date2 - OMI AI - date1
        ax6 = fig.add_subplot(2,4,6, projection = workcrs) # date2 - CERES SWF
        ax7 = fig.add_subplot(2,4,7, projection = workcrs) # date2 - NN SWF
        ax8 = fig.add_subplot(2,4,8, projection = workcrs) # date2 - NN - CERES

        fig2 = plt.figure(figsize = (9, 4))
        ax9 = fig2.add_subplot(1,2,1)
        ax10 = fig2.add_subplot(1,2,2)
   
    # Plot OMI UVAI
    # ------------- 
    ai_mesh = ax1.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmax = 5, cmap = 'jet')
    #cbar = fig.colorbar(mesh, ax = ax1, label = 'OMI UVAI')
    if(auto_zoom):
        ax1.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax1.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_title('OMI UVAI')
    ax1.coastlines()
   
    # Plot CERES SWF observations
    # --------------------------- 
    ceres_mesh = ax2.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_orig, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = min_swf, vmax = max_swf)
    #cbar = fig.colorbar(mesh, ax = ax2, label = 'SWF [Wm$^{-2}$]')
    #ax1.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax2.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax2.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.set_title('Aqua CERES SWF')
    ax2.coastlines()
   
    # Plot neural net output
    # ---------------------- 
    mesh = ax3.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_calc, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = min_swf, vmax = max_swf)
        #transform = ccrs.PlateCarree(), shading = 'auto')
    #cbar = fig.colorbar(mesh, ax = ax3, label = 'SWF [Wm$^{-2}$]')
    #ax2.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax3.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.set_title('Neural Network Aerosol-free SWF')
    ax3.coastlines()
   
    # Plot the estimated forcing
    # -------------------------- 
    force_mesh = ax4.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], diff_calc, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = -80, vmax = 80, cmap = 'bwr')
    #cbar = fig.colorbar(mesh, ax = ax4, label = 'SWF [Wm$^{-2}$]')
    #ax3.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    #ax3.set_boundary(circle, transform=ax3.transAxes)
    #ax3.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax4.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax4.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax4.set_boundary(circle, transform=ax4.transAxes)
    ax4.set_title('NN (aerosol-free) - CERES (obs.)')
    ax4.coastlines()
   

    # Plot scatter figures
    # --------------------
    #if(include_scatter):
    #mask_orig = np.ma.masked_where(np.isnan(in_calc['omi_uvai_pert'][:,:]), mask_orig)
    #mask_calc = np.ma.masked_where(np.isnan(in_calc['omi_uvai_pert'][:,:]), mask_calc)
    mask_orig = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:] > 1.0, mask_orig)
    mask_calc = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:] > 1.0, mask_calc)
    both_orig = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_orig)
    both_calc = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_calc)
    both_orig = both_orig.compressed()
    both_calc = both_calc.compressed()
    xy = np.vstack([both_orig, both_calc])
    print("HERE:", both_orig.shape[0])
    if(len(both_orig) > 1):
        #r2 = r2_score(both_orig, both_calc)
        pearsonr, pearson_pval = stats.pearsonr(both_orig, both_calc)
        print("pearsonr, pval = ", pearsonr, pearson_pval)
        z = stats.gaussian_kde(xy)(xy)       
        ax9.scatter(both_orig, both_calc, c = z, s = 1)
        lims = [\
                np.min([ax9.get_xlim(), ax9.get_ylim()]),\
                np.max([ax9.get_xlim(), ax9.get_ylim()]),
        ]
        ax9.plot(lims, lims, 'r', linestyle = ':', alpha = 1.00)
        ax9.set_xlim(lims)
        ax9.set_ylim(lims)
        ax9.grid(alpha = 0.5)
        ax9.set_xlabel('CERES SWF [Wm$^{-2}$]')
        ax9.set_ylabel('NN SWF [Wm$^{-2}$]')
        #ax9.set_xticklabels([])
        ax9.set_title(dt_date_str1.strftime('%Y-%m-%d %H:%M UTC\nAerosol-free Swath'))
        #ptext = 'r$^{2}$ = ' + str(np.round(r2, 3))
        ptext = 'r = ' + str(np.round(pearsonr, 3))
        plot_figure_text(ax9, ptext, xval = 450, yval = 100, \
            fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')

    in_calc.close()

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    #
    # Process and plot output for the second file
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
 
    in_calc = h5py.File(infile2, 'r')
    mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                   (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
    mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
    mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)

    # Calculate forcing
    # -----------------
    diff_calc = mask_calc - mask_orig
    
    mask_AI = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:]  == -999., \
        in_calc['omi_uvai_pert'][:,:])
   
    # Plot OMI UVAI
    # ------------- 
    mesh = ax5.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmax = 5, cmap = 'jet')
    #cbar = fig.colorbar(mesh, ax = ax5, label = 'OMI UVAI')
    if(auto_zoom):
        ax5.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax5.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax5.set_boundary(circle, transform=ax5.transAxes)
    #ax5.set_title('OMI UVAI')
    ax5.coastlines()
   
    # Plot CERES SWF observations
    # --------------------------- 
    mesh = ax6.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_orig, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = min_swf, vmax = max_swf)
    #cbar = fig.colorbar(mesh, ax = ax6, label = 'SWF [Wm$^{-2}$]')
    #ax1.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax6.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax6.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax6.set_boundary(circle, transform=ax6.transAxes)
    #ax6.set_title('Observed SWF')
    ax6.coastlines()
   
    # Plot neural net output
    # ---------------------- 
    mesh = ax7.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_calc, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = min_swf, vmax = max_swf)
        #transform = ccrs.PlateCarree(), shading = 'auto')
    #cbar = fig.colorbar(mesh, ax = ax7, label = 'SWF [Wm$^{-2}$]')
    #ax2.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax7.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax7.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax7.set_boundary(circle, transform=ax7.transAxes)
    #ax7.set_title('Calc. Aerosol-free SWF')
    ax7.coastlines()
   
    # Plot the estimated forcing
    # -------------------------- 
    mesh = ax8.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], diff_calc, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = -80, vmax = 80, cmap = 'bwr')
    #cbar = fig.colorbar(mesh, ax = ax8, label = 'SWF [Wm$^{-2}$]')
    #ax3.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    #ax3.set_boundary(circle, transform=ax3.transAxes)
    #ax3.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax8.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax8.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax8.set_boundary(circle, transform=ax8.transAxes)
    #ax8.set_title('Calculated - Observed')
    ax8.coastlines()
  

    #if(include_scatter):
    # Plot scatter figures
    # --------------------
    #mask_orig = np.ma.masked_where(np.isnan(in_calc['omi_uvai_pert'][:,:]), mask_orig)
    #mask_calc = np.ma.masked_where(np.isnan(in_calc['omi_uvai_pert'][:,:]), mask_calc)
    mask_orig = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:] > 1.0, mask_orig)
    mask_calc = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:] > 1.0, mask_calc)
    both_orig = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_orig)
    both_calc = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_calc)
    both_orig = both_orig.compressed()
    both_calc = both_calc.compressed()
    xy = np.vstack([both_orig, both_calc])
    print("HERE:", both_orig.shape[0])
    if(len(both_orig) > 1):
        #r2 = r2_score(both_orig, both_calc)
        pearsonr, pearson_pval = stats.pearsonr(both_orig, both_calc)
        print("pearsonr, pval = ", pearsonr, pearson_pval)
        z = stats.gaussian_kde(xy)(xy)       
        ax10.scatter(both_orig, both_calc, c = z, s = 1)
        lims = [\
                np.min([ax10.get_xlim(), ax10.get_ylim()]),\
                np.max([ax10.get_xlim(), ax10.get_ylim()]),
        ]
        ax10.plot(lims, lims, 'r', linestyle = ':', alpha = 1.00)
        ax10.set_xlim(lims)
        ax10.set_ylim(lims)
        ax10.grid(alpha = 0.5)
        ax10.set_xlabel('CERES SWF [Wm$^{-2}$]')
        #ax10.set_ylabel('NN SWF')
        ax10.set_yticklabels([])
        ax10.set_title(dt_date_str2.strftime('%Y-%m-%d %H:%M UTC\nAerosol Swath (Plotted for UVAI < 1)'))
        #ptext = 'r$^{2}$ = ' + str(np.round(r2, 3))
        ptext = 'r = ' + str(np.round(pearsonr, 3))
        plot_figure_text(ax10, ptext, xval = 450, yval = 100, \
            fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')

 
    in_calc.close()

    if(include_scatter):
        cbar_ax1 = fig.add_axes([0.02, 0.09, 0.15, 0.015])
        cbar = fig.colorbar(ai_mesh, \
            cax = cbar_ax1, shrink = 0.8, orientation = 'horizontal', \
            label = 'AI', extend = 'max')
        cbar_ax2 = fig.add_axes([0.25, 0.09, 0.30, 0.015])
        cbar = fig.colorbar(ceres_mesh, \
            cax = cbar_ax2, shrink = 0.8, orientation = 'horizontal', \
            label = 'SWF [Wm$^{-2}$]')
        cbar_ax3 = fig.add_axes([0.62, 0.09, 0.15, 0.015])
        cbar = fig.colorbar(force_mesh, \
            cax = cbar_ax3, shrink = 0.8, orientation = 'horizontal', \
            label = 'Forcing [Wm$^{-2}$]')
    else:
        cbar_ax1 = fig.add_axes([0.03, 0.08, 0.20, 0.015])
        cbar = fig.colorbar(ai_mesh, \
            cax = cbar_ax1, shrink = 0.8, orientation = 'horizontal', \
            label = 'UVAI', extend = 'max')
        cbar_ax2 = fig.add_axes([0.30, 0.08, 0.40, 0.015])
        cbar = fig.colorbar(ceres_mesh, \
            cax = cbar_ax2, shrink = 0.8, orientation = 'horizontal', \
            label = 'SWF [Wm$^{-2}$]')
        cbar_ax3 = fig.add_axes([0.77, 0.08, 0.20, 0.015])
        cbar = fig.colorbar(force_mesh, \
            cax = cbar_ax3, shrink = 0.8, orientation = 'horizontal', \
            label = 'Forcing [Wm$^{-2}$]')




 
    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax4, 'd)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax5, 'e)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax6, 'f)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax7, 'g)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax8, 'h)', fontsize = 11, backgroundcolor = None)

    plot_subplot_label(ax9, 'i)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax10,'j)', fontsize = 11, backgroundcolor = None)
    
    fig.suptitle(dt_date_str1.strftime('Top:       %Y-%m-%d %H:%M UTC') + '\n' + \
                 dt_date_str2.strftime('Bottom: %Y-%m-%d %H:%M UTC'))
   
    fig2.suptitle('Comparison of NN and CERES SWF')
 
    fig.tight_layout(rect = [0,0.05,1,1])
    fig2.tight_layout()

    #fig.subplots_adjust(wspace = 0.1, hspace = 0.00)

    if(save):    
        if(auto_zoom):
            zoom_add = '_zoom'
        else:
            zoom_add = ''
        outname = 'ceres_nn_compare_dual_v2_' + date_str1 + '_' + date_str2 + '_' + sim_name1 + zoom_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)

        outname = 'ceres_nn_compare_dual_scatter_' + date_str1 + '_' + date_str2 + '_' + sim_name1 + '.png'
        fig2.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.subplot_tool()
        plt.show()



def plot_compare_NN_output(calc_data, auto_zoom = False, save = False):
    
    file_name = calc_data.strip().split('/')[-1] 
    date_str = file_name.split('.')[0].split('_')[-1]
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    
    sim_name = file_name.split('_')[3]
    print(sim_name, date_str)
    
    in_calc = h5py.File(calc_data)
    
    mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                   (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
    mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
    mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)
    
    ##!#mask_ctp = np.ma.masked_invalid(in_calc['modis_cld_top_pres'][:,:])
    ##!#mask_ctp = np.ma.masked_where(mask_ctp == -999., mask_ctp)
    ##!#
    ##!#mask_cod = np.ma.masked_invalid(in_calc['modis_cod'][:,:])
    ##!#mask_cod = np.ma.masked_where(mask_cod == -999., mask_cod)
    ##!#
    ##!#mask_ice = np.ma.masked_invalid(in_calc['nsidc_ice'][:,:])
    ##!#mask_ice = np.ma.masked_where(mask_ice == -999., mask_ice)
    
    diff_calc = mask_calc - mask_orig
    
    both_orig = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_orig)
    both_calc = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_calc)
   
    mask_AI = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:]  == -999., \
        in_calc['omi_uvai_pert'][:,:])

    if(auto_zoom):

        # Figure out where the pixels containing AI are located
        # -----------------------------------------------------
        high_AI_idxs = np.where(mask_AI > 1.5)
        #mask_AI = np.ma.masked_where(mask_AI < 1.5, mask_AI)

        # Switch the longitudes to 0 - 360 rather than -180 - 180
        # -------------------------------------------------------
        work_lons = np.copy(in_calc['omi_lon'][:,:])
        work_lons[work_lons < 0] = 360 + work_lons[work_lons < 0]
        

        # Figure out the minimum/maximum latitude/longitudes in the AI ranges
        # -------------------------------------------------------------------
        min_lat = np.min(in_calc['omi_lat'][:,:][high_AI_idxs]) - 2
        max_lat = np.max(in_calc['omi_lat'][:,:][high_AI_idxs]) + 2
        
        if(max_lat > 90):
            max_lat = 90

        # Calculate the 5th and 95th percentiles of longitude to prevent
        # outliers from ruining the mean or max/min values
        # --------------------------------------------------------------
        min_lon = np.percentile(work_lons[high_AI_idxs], 5) - 2
        max_lon = np.percentile(work_lons[high_AI_idxs], 95) + 2

        center_lon = (min_lon + max_lon) / 2.
        workcrs = ccrs.NorthPolarStereo(central_longitude = center_lon)

        #fig = plt.figure(figsize = (10.5, 6))
        fig = plt.figure(figsize = (11, 7))
        ax5 = fig.add_subplot(2,3,1, projection = workcrs)
        ax1 = fig.add_subplot(2,3,2, projection = workcrs)
        ax2 = fig.add_subplot(2,3,4, projection = workcrs)
        ax3 = fig.add_subplot(2,3,5, projection = workcrs)
        ax4 = fig.add_subplot(2,3,3)
        
    else:
        workcrs = ccrs.NorthPolarStereo(central_longitude = 0)

        #fig = plt.figure(figsize = (10.5, 6))
        fig = plt.figure(figsize = (8, 7))
        ax5 = fig.add_subplot(2,2,1, projection = workcrs)
        ax1 = fig.add_subplot(2,2,2, projection = workcrs)
        ax2 = fig.add_subplot(2,2,3, projection = workcrs)
        ax3 = fig.add_subplot(2,2,4, projection = workcrs)
 
    #ax4 = fig.add_subplot(2,3,5)
    
    #ax4 = fig.add_subplot(2,2,4, projection = ccrs.NorthPolarStereo())
    
    max_swf = np.max([np.max(mask_orig), np.max(mask_calc)])
    min_swf = np.min([np.min(mask_orig), np.min(mask_calc)])
 
    mesh = ax1.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_orig, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = min_swf, vmax = max_swf)
        #transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax1, label = 'SWF [Wm$^{-2}$]')
    #ax1.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax1.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax1.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_title('Observed SWF')
    ax1.coastlines()
    
    mesh = ax2.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_calc, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = min_swf, vmax = max_swf)
        #transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax2, label = 'SWF [Wm$^{-2}$]')
    #ax2.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax2.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax2.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.set_title('Calc. Aerosol-free SWF')
    ax2.coastlines()
    
    mesh = ax3.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], diff_calc, \
        transform = ccrs.PlateCarree(), shading = 'auto', vmin = -80, vmax = 80, cmap = 'bwr')
    cbar = fig.colorbar(mesh, ax = ax3, label = 'SWF [Wm$^{-2}$]')
    #ax3.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    #ax3.set_boundary(circle, transform=ax3.transAxes)
    #ax3.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax3.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax3.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.set_title('Calculated - Observed')
    ax3.coastlines()
    
    mesh = ax5.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), vmin = -2, vmax = 4, shading = 'auto', cmap = 'jet')
    cbar = fig.colorbar(mesh, ax = ax5, label = 'UVAI Perturbation')
    #ax5.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    #ax5.set_boundary(circle, transform=ax5.transAxes)
    #ax5.set_extent([117, 170,70,  85], ccrs.PlateCarree())
    if(auto_zoom):
        ax5.set_extent([min_lon, max_lon , min_lat, max_lat], ccrs.PlateCarree())
    else:
        ax5.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
        ax5.set_boundary(circle, transform=ax5.transAxes)
    ax5.set_title('OMI UVAI')
    ax5.coastlines()
   
    if(auto_zoom): 
        work_force = np.abs(diff_calc)
        work_force = np.ma.masked_where(mask_AI < 1.0, work_force)
        mask_AI    = np.ma.masked_where(mask_AI < 1.0, mask_AI)  
        work_force = np.ma.masked_invalid(work_force)
        mask_AI    = np.ma.masked_invalid(mask_AI)
        work_force = np.ma.masked_where(mask_AI.mask == True, work_force)
        mask_AI    = np.ma.masked_where(work_force.mask == True, mask_AI)  


        #print("HERE:", work_force.compressed().shape, mask_AI.compressed().shape)

        if(len(mask_AI.compressed()) > 2):
            r2 = r2_score(mask_AI.compressed(), work_force.compressed())
            print('R2:', r2)
            ax4.set_title('r$^{2}$ = ' + str(np.round(r2, 3)))

        ax4.scatter(mask_AI, work_force, s = 6, color = 'k')
        ax4.set_xlabel('OMI UVAI Perturbation')
        ax4.set_ylabel('Forcing [W/m2]')
 
    """ 
    xy = np.vstack([both_orig.compressed(), both_calc.compressed()])
    if(len(both_orig.compressed()) > 1):
        r2 = r2_score(both_orig.compressed(), both_calc.compressed())
        z = stats.gaussian_kde(xy)(xy)       
        ax4.scatter(both_orig.compressed(), both_calc.compressed(), c = z, s = 1)
        #ax4.set_title('r$^{2}$ = ' + str(np.round(r2, 3)))

        # Calculate RMSE
        # --------------
        rmse = mean_squared_error(both_orig.compressed(), both_calc.compressed(), squared = False)
        ptext = 'r$^{2}$ = ' + str(np.round(r2, 3)) + \
            '\nRMSE = ' + str(np.round(rmse, 1))
        #plot_figure_text(ax3, ptext, location = 'lower_right', \
        plot_figure_text(ax4, ptext, xval = 400, yval = 25, \
            fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')
    lims = [\
            np.min([ax4.get_xlim(), ax4.get_ylim()]),\
            np.max([ax4.get_xlim(), ax4.get_ylim()]),
    ]
    ax4.plot(lims, lims, 'r', linestyle = ':', alpha = 0.75)
    ax4.set_xlim(lims)
    ax4.set_ylim(lims)
    ax4.set_xlabel('Observed SWF [Wm$^{-2}$]')
    ax4.set_ylabel('Calculated SWF [Wm$^{-2}$]')
    """ 
   
    ##!#print(np.min(mask_cod), np.max(mask_cod))
     
    ##!#good_idxs = np.where( (mask_AI.mask == False) & \
    ##!#                      (mask_calc.mask == False) & \
    ##!#                      (mask_cod.mask == False) & \
    ##!#                      (mask_ice > 80.) & \
    ##!#                      (mask_ice < 101.) & \
    ##!#                      (mask_AI > 2))
    ##!#
    ##!#scat_ai   = mask_AI[good_idxs].flatten()
    ##!#scat_diff = diff_calc[good_idxs].flatten()
    ##!#ax6.scatter(scat_ai, scat_diff, c = 'k', s = 6)
    ##!#ax6.set_xlabel('OMI UVAI PERT')
    ##!#ax6.set_ylabel('Calc Aer Forcing [W/m2]')
    
    plot_subplot_label(ax3, 'd)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax1, 'b)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax2, 'c)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax5, 'a)', fontsize = 11, backgroundcolor = None)
    #plot_subplot_label(ax4, 'e)', fontsize = 11, backgroundcolor = None)
    
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    
    #in_base.close()
    in_calc.close()
    
    fig.tight_layout()

    if(save):    
        if(auto_zoom):
            zoom_add = '_zoom'
        else:
            zoom_add = ''
        outname = 'ceres_nn_compare_' + date_str + '_' + sim_name + zoom_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

# Makes a plot of Observed SWF, NN SWF, and Regression
def plot_compare_NN_output_noaer(calc_data, num_bins = 100, astrofit = False, \
        save = False):
   
    file_name = calc_data.strip().split('/')[-1] 
    date_str = file_name.split('.')[0].split('_')[-1]
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    
    sim_name = file_name.split('_')[3]
    print(sim_name, date_str)
    
    in_calc = h5py.File(calc_data)
    
    mask_orig = np.ma.masked_where((in_calc['ceres_swf'][:,:] == -999.) | \
                                   (in_calc['ceres_swf'][:,:] > 3000), in_calc['ceres_swf'][:,:])
    mask_calc = np.ma.masked_invalid(in_calc['calc_swf'])
    mask_calc = np.ma.masked_where(mask_calc == -999., mask_calc)
    
    ##!#mask_ctp = np.ma.masked_invalid(in_calc['modis_cld_top_pres'][:,:])
    ##!#mask_ctp = np.ma.masked_where(mask_ctp == -999., mask_ctp)
    ##!#
    ##!#mask_cod = np.ma.masked_invalid(in_calc['modis_cod'][:,:])
    ##!#mask_cod = np.ma.masked_where(mask_cod == -999., mask_cod)
    ##!#
    ##!#mask_ice = np.ma.masked_invalid(in_calc['nsidc_ice'][:,:])
    ##!#mask_ice = np.ma.masked_where(mask_ice == -999., mask_ice)
    
    ##!#diff_calc = mask_calc - mask_orig
    
    both_orig = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_orig)
    both_calc = np.ma.masked_where((mask_orig.mask == True) | (mask_calc.mask == True), mask_calc)
    
    fig = plt.figure(figsize = (11, 6))
    ax4 = fig.add_subplot(2,3,1, projection = ccrs.NorthPolarStereo())
    ax1 = fig.add_subplot(2,3,2, projection = ccrs.NorthPolarStereo())
    ax2 = fig.add_subplot(2,3,4, projection = ccrs.NorthPolarStereo())
    ax3 = fig.add_subplot(2,3,5)
    ax5 = fig.add_subplot(2,3,6)
    
    mask_AI = np.ma.masked_where(in_calc['omi_uvai_pert'][:,:]  == -999., in_calc['omi_uvai_pert'][:,:])
    mesh = ax4.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_AI, \
        transform = ccrs.PlateCarree(), vmin = -2, vmax = 4, shading = 'auto', cmap = 'jet')
    cbar = fig.colorbar(mesh, ax = ax4, label = 'UVAI Perturbation')
    ax4.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    ax4.set_title('OMI UVAI')
    ax4.set_boundary(circle, transform=ax4.transAxes)
    ax4.coastlines()
    
    mesh = ax1.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_orig, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax1, label = 'SWF [Wm$^{-2}$]')
    ax1.set_extent([-180, 180,65,  90], ccrs.PlateCarree())
    ax1.set_title('Observed SWF')
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.coastlines()
    
    
    mesh = ax2.pcolormesh(in_calc['omi_lon'][:,:], in_calc['omi_lat'][:,:], mask_calc, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax2, label = 'SWF [Wm$^{-2}$]')
    ax2.set_extent([-180, 180,65, 90], ccrs.PlateCarree())
    ax2.set_title('NN Aerosol-free SWF')
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.coastlines()
    
    
    both_orig = both_orig.compressed()
    both_calc = both_calc.compressed()
    xy = np.vstack([both_orig, both_calc])
    if(len(both_orig) > 1):
        r2 = r2_score(both_orig, both_calc)
        z = stats.gaussian_kde(xy)(xy)       
        ax3.scatter(both_orig, both_calc, c = z, s = 1)

        # Calculate RMSE
        # --------------
        rmse = mean_squared_error(both_orig, both_calc, squared = False)
        #ax3.set_title('r$^{2}$ = ' + str(np.round(r2, 3)) + \
        #    'RMSE = ' + str(np.round(rmse, 1)))
        ptext = 'r$^{2}$ = ' + str(np.round(r2, 3)) + \
            '\nRMSE = ' + str(np.round(rmse, 1))
        #plot_figure_text(ax3, ptext, location = 'lower_right', \
        plot_figure_text(ax3, ptext, xval = 400, yval = 50, \
            fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')
            
        # Plot a histogram of the errors between these two
        # ------------------------------------------------
        errors = both_calc - both_orig
        ax5.hist(errors, bins = num_bins)
        if(astrofit):
            bin_heights,bin_borders = np.histogram(errors,bins=num_bins)
            bin_widths = np.diff(bin_borders)
            bin_centers = bin_borders[:-1] + bin_widths / 2
            t_init = models.Gaussian1D()
            fit_t = fitting.LevMarLSQFitter()
            t = fit_t(t_init, bin_centers, bin_heights)
            print('Amplitude: ',np.round(t.amplitude.value,3))
            print('Mean:      ',np.round(t.mean.value,3))
            print('StDev:     ',np.round(t.stddev.value,3))
            x_interval_for_fit = np.linspace(bin_centers[0],bin_centers[-1],200)
            ax5.plot(x_interval_for_fit,t(x_interval_for_fit),label='fit',c='tab:red')
            ax5.axvline(0, color = 'k', linestyle = '--')
            ax5.axvline(t.mean.value, color = 'r', linestyle = '--')
            ptext = 'μ = ' + str(np.round(t.mean.value, 1)) + \
                '\nσ = ' + str(np.round(t.stddev.value, 1))
            #plot_figure_text(ax3, ptext, location = 'lower_right', \
            ax5.set_title(ptext)
            #plot_figure_text(ax5, ptext, xval = 100, yval = 800, \
            #    fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')
        ax5.set_xlabel('Error (NN - obs)')
        ax5.set_ylabel('Counts')

    lims = [\
            np.min([ax3.get_xlim(), ax3.get_ylim()]),\
            np.max([ax3.get_xlim(), ax3.get_ylim()]),
    ]
    ax3.plot(lims, lims, 'r', linestyle = ':', alpha = 0.75)
    ax3.set_xlim(lims)
    ax3.set_ylim(lims)
    ax3.set_xlabel('Observed SWF [Wm$^{-2}$]')
    ax3.set_ylabel('Calculated SWF [Wm$^{-2}$]')
  
    
 
    plot_subplot_label(ax4, 'a)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax1, 'b)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax2, 'c)', fontsize = 11, backgroundcolor = None)
    plot_subplot_label(ax3, 'd)', fontsize = 11, backgroundcolor = None)
    
    plt.suptitle(dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    
    in_calc.close()
    
    fig.tight_layout()

    if(save):    
        outname = 'ceres_nn_compare_aerfree_' + date_str + '_' + sim_name + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

# Example COD bin edges:
#       cod_bin_edges = np.array([0,0.0001,2,4,8,12,20,30,50])
#
# This function plots the NN output for given COD bin edges
# and for the specified AI, ICE, and SZA min/max values. For
# each COD bin, the selected, scattered forcing 
# (calculated as NN - OBS) data are plotted as a function of AI.
# There is one plot for each COD bin.
# --------------------------------------------------------------
def plot_NN_scatter_multiCOD(test_dict, cod_bin_edges, \
        ai_min, ice_min, ice_max, sza_min, sza_max, trend_type = 'linregress', \
        plot_bounds = False, save = False):

    plt.close('all')
    fig = plt.figure(figsize = (12, 5))
    axs = fig.subplots(nrows = 2, ncols = 4)
    
    flat_axs = axs.flatten()
    
    #ice_min = 101
    #ice_max = 255
    #sza_min = 50
    #sza_max = 55
    for ii in range(flat_axs.shape[0]):
        print("COD MIN", cod_bin_edges[ii], "COD MAX", cod_bin_edges[ii+1])
        good_idxs = select_idxs(test_dict, ai_min, ice_min, ice_max, cod_bin_edges[ii], \
            cod_bin_edges[ii+1], sza_min, sza_max) 
    
        diff_calc = test_dict['calc_swf'][good_idxs] - test_dict['ceres_swf'][good_idxs]
    
        local_xdata =test_dict['omi_uvai_pert'][good_idxs] 
    
        flat_axs[ii].scatter(local_xdata, diff_calc, s = 3, color = 'k')
       
        flat_axs[ii].set_title('COD: ' + str(cod_bin_edges[ii]) + ' - ' + str(cod_bin_edges[ii + 1]))
    
        flat_axs[ii].set_xlabel('OMI AI')
        flat_axs[ii].set_ylabel('Forcing [W/m2]')
 
        #trend_type = 'theil-sen'
        #trend_type = 'linregress'
        if(len(local_xdata) > 10000):
            print("ERROR: Plotting trend line with too many points")
            print("       Preventing this for the sake of computer safety")
        else:
            plot_trend_line(flat_axs[ii], local_xdata, diff_calc, color='tab:red', \
                linestyle = '-',  slope = trend_type, plot_bounds = plot_bounds)
 
    if(ice_min >= 101):
        title_str = 'Land\n'
    elif( (ice_min < 20 ) & (ice_max <= 20) ):
        title_str = 'Ocean (' + str(ice_min) + '% < Sea Ice Conc. < ' + str(ice_max) + '%)\n'
    elif( (ice_min >= 20 ) & (ice_max <= 80) ):
        title_str = 'Mixed ice/ocean (' + str(ice_min) + ' < Sea Ice Conc. < ' + str(ice_max) + '%)\n'
    elif( (ice_min >= 80 ) & (ice_max <= 100) ):
        title_str = 'Ice (' + str(ice_min) + ' < Sea Ice Conc. < ' + str(ice_max) + '%)\n'
    title_str = title_str + 'SZA: ' + str(sza_min) + ' - ' + str(sza_max)
    plt.suptitle(title_str) 
     
    fig.tight_layout()
    if(save):
        outname = 'nn_force_scatter_multiCOD' + \
            '_ice' + str(int(ice_min)) + 'to' + str(int(ice_max)) + \
            '_sza' + str(int(sza_min)) + 'to' + str(int(sza_max)) + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:   
        plt.show()

# This function plots the NN output for given COD bin edges
# and for the specified AI, ICE, and SZA min/max values. For
# each COD bin, the selected, scattered forcing 
# (calculated as NN - OBS) data are plotted as a function of AI.
# There is one plot for each COD bin.
#
# show_specific_cod: the COD index 
# --------------------------------------------------------------
def plot_NN_scatter_combined_alltypes(test_dict, bin_dict, \
        ai_min, sza_min, sza_max, trend_type = 'linregress', \
        plot_bounds = False, show_specific_cod = None, \
        min_ai_for_stats = None, return_line_vals = False, \
        ai_calc_bins = None, save = False):

    plt.close('all')
    if(len(bin_dict['ice_bin_means']) == 4):
        fig = plt.figure(figsize = (9, 8))
        axs = fig.subplots(nrows = 2, ncols = 2, sharex = True, sharey = True)
        num_bins = 4
    elif(len(bin_dict['ice_bin_means']) == 6):
        fig = plt.figure(figsize = (9, 7.5))
        axs = fig.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True)
        num_bins = 6
    
    flat_axs = axs.flatten()

    xvals = np.arange(len(bin_dict['cod_bin_means']))

    plot_c = cm.turbo((xvals-np.min(xvals))/\
        (np.max(xvals)-np.min(xvals)))
 
    if(return_line_vals):
        if(ai_calc_bins is None):
            ai_calc_bins = np.array([0,2,4,6,30])
        print("CALCULATING AI AVGS FOR BINS:", ai_calc_bins)

        slope_vals = np.full( (len(bin_dict['ice_bin_means']), \
            len(bin_dict['cod_bin_means'])), np.nan)
        intpt_vals = np.full( (len(bin_dict['ice_bin_means']), \
            len(bin_dict['cod_bin_means'])), np.nan)
        mean_force_vals = np.full( (len(bin_dict['ice_bin_means']), \
            len(bin_dict['cod_bin_means']), len(ai_calc_bins) - 1), np.nan)
        std_force_vals = np.full( (len(bin_dict['ice_bin_means']), \
            len(bin_dict['cod_bin_means']), len(ai_calc_bins) - 1), np.nan)
        
 
    # Loop over surface types
    # -----------------------
    for ii in range(len(bin_dict['ice_bin_means'])):

        print('ICE BIN:', bin_dict['ice_bin_edges'][ii], bin_dict['ice_bin_edges'][ii + 1])

        total_count = 0

        # Loop over COD ranges
        # --------------------
        for jj in range(len(bin_dict['cod_bin_means'])):

            if( (show_specific_cod is None) | \
                ( (show_specific_cod is not None ) & \
                  (show_specific_cod == jj))):

                print('COD BIN:', bin_dict['cod_bin_edges'][jj], bin_dict['cod_bin_edges'][jj + 1])
                # Select the forcing values for these conditions
                # ----------------------------------------------
                good_idxs = select_idxs(test_dict, ai_min, \
                    bin_dict['ice_bin_edges'][ii], bin_dict['ice_bin_edges'][ii + 1], \
                    bin_dict['cod_bin_edges'][jj], bin_dict['cod_bin_edges'][jj + 1], \
                    sza_min, sza_max) 

                diff_calc = test_dict['calc_swf'][good_idxs] - test_dict['ceres_swf'][good_idxs]

                local_xdata =test_dict['omi_uvai_pert'][good_idxs] 
    
                flat_axs[ii].scatter(local_xdata, diff_calc, s = 1, color = plot_c[jj])
 
                print("SIZE:", len(local_xdata), len(diff_calc))
                if(len(local_xdata) > 30000):
                    print("ERROR: Plotting trend line with too many points")
                    print("       Preventing this for the sake of computer safety")
                else:
                    back_data = plot_trend_line(flat_axs[ii], local_xdata, diff_calc, color= plot_c[jj], \
                        linestyle = '-',  slope = trend_type, plot_bounds = plot_bounds)
                    if(return_line_vals):
                        slope_vals[ii,jj] = back_data.slope
                        intpt_vals[ii,jj] = back_data.intercept

                total_count += len(local_xdata)

                if(return_line_vals):
                    for kk in range(len(ai_calc_bins) - 1):
                        range_calc = np.ma.masked_where(\
                            (local_xdata < ai_calc_bins[kk]) | \
                            (local_xdata > ai_calc_bins[kk + 1]), \
                            diff_calc).compressed()
                        mean_force_vals[ii,jj,kk] = np.nanmean(range_calc)
                        std_force_vals[ii,jj,kk]  = np.nanstd(range_calc)

                if(min_ai_for_stats != None):
                    diff_calc = np.ma.masked_where(local_xdata < min_ai_for_stats, diff_calc).compressed()
                    mean_force = np.mean(diff_calc)
                    std_force  = np.std(diff_calc)
                    print("Mean forcing for AI > ", min_ai_for_stats, " = ", np.round(mean_force, 1))
                    print("StDv forcing for AI > ", min_ai_for_stats, " = ", np.round(std_force, 1))

        flat_axs[ii].axhline(0, linestyle = ':', color = 'gray', alpha = 1.0) 
        flat_axs[ii].grid(alpha = 0.5)

        # set up the subplot title
        if(bin_dict['ice_bin_edges'][ii] == 0.):
            title_str = 'Ocean\n(Ice conc. < ' + \
                str(int(bin_dict['ice_bin_edges'][ii + 1])) + '%)'
        elif( (bin_dict['ice_bin_edges'][ii] >= 20.) & \
                (bin_dict['ice_bin_edges'][ii + 1] <= 80) ):
            title_str = 'Mix\n(' + str(int(bin_dict['ice_bin_edges'][ii])) + \
                '% <= Ice conc. < ' + \
                str(int(bin_dict['ice_bin_edges'][ii + 1])) + '%)'
        elif( (bin_dict['ice_bin_edges'][ii] <= 80.) & \
                (bin_dict['ice_bin_edges'][ii + 1] <= 100.2) ):
            title_str = 'Ice\n(Ice conc. >= ' + \
                str(int(bin_dict['ice_bin_edges'][ii])) + '%)'
        elif( (bin_dict['ice_bin_edges'][ii + 1] > 150.) ):
            title_str = 'Land'
        flat_axs[ii].set_ylim(-225, 225)
        flat_axs[ii].set_title(title_str + '\nn = ' + str(total_count))
        #flat_axs[ii].set_title(str(bin_dict['ice_bin_edges'][ii]) + ' - ' + \
        #    str(bin_dict['ice_bin_edges'][ii + 1]) + '\nn = ' + str(total_count))
   
    #if(len(bin_dict['ice_bin_means']) == 4):
    plot_subplot_label(flat_axs[0], 'a)', fontsize = 12, backgroundcolor = None, location = 'lower_right')
    plot_subplot_label(flat_axs[1], 'b)', fontsize = 12, backgroundcolor = None, location = 'lower_right')
    plot_subplot_label(flat_axs[2], 'c)', fontsize = 12, backgroundcolor = None, location = 'lower_right')
    plot_subplot_label(flat_axs[3], 'd)', fontsize = 12, backgroundcolor = None, location = 'lower_right')
    if(num_bins == 4):
        flat_axs[2].set_xlabel('OMI UVAI Pert.')
        flat_axs[3].set_xlabel('OMI UVAI Pert.')
        flat_axs[0].set_ylabel('Direct Forcing [Wm$^{-2}$]')
        flat_axs[2].set_ylabel('Direct Forcing [Wm$^{-2}$]')
    #elif(len(bin_dict['ice_bin_means']) == 6):
    elif(num_bins == 6):
        flat_axs[3].set_xlabel('OMI UVAI Pert.')
        flat_axs[4].set_xlabel('OMI UVAI Pert.')
        flat_axs[5].set_xlabel('OMI UVAI Pert.')
        flat_axs[0].set_ylabel('Direct Forcing [Wm$^{-2}$]')
        flat_axs[3].set_ylabel('Direct Forcing [Wm$^{-2}$]')
        plot_subplot_label(flat_axs[4], 'e)', fontsize = 12, backgroundcolor = None, location = 'lower_right')
        plot_subplot_label(flat_axs[5], 'f)', fontsize = 12, backgroundcolor = None, location = 'lower_right')

    #cmap = cm.turbo
    #norm = mc.Normalize(vmin=np.min(xvals), \
    #    vmax = np.max(xvals))

    cmap = 'turbo' 
    shrk = 1.0
    cmap2 = plt.get_cmap('turbo')
    #colorvals = np.arange(0, len(xvals), 1)
    colorvals = bin_dict['cod_bin_edges']
    norm = mpl.colors.BoundaryNorm(colorvals, cmap2.N, extend = 'max')

    #cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
    #             ax=flat_axs[-1], orientation='vertical', label='COD')

    #mesh = ax.pcolormesh(lons, lats, data[tidx,:,:], transform = ccrs.PlateCarree(), \
    #    shading = 'auto', cmap = 'jet', vmax = value_max, vmin = -0.2)
    #if(plot_cbar):
    #fig.subplots_adjust(bottom = 0.8)
    #cbar_ax = fig.add_axes([0.17, 0.10, 0.70, 0.01])
    #fig.colorbar(cm.ScalarMappable(norm = norm, cmap = cmap), cax = cbar_ax, label = 'test', orientation = 'horizontal')


    cbar_ax = fig.add_axes([0.17, 0.09, 0.70, 0.01])
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), \
        cax = cbar_ax, shrink = 0.8, orientation = 'horizontal', \
        label = 'COD')
    fig.tight_layout(rect = [0,0.1,1,1])
    #cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax = axs, location = 'bottom', shrink = 0.8)


    if(save):
        if(show_specific_cod is not None):
            outname = 'nn_force_scatter_combined_' + test_dict['sim_name'] + '_numsfcbins' + str(num_bins) + \
                '_sza' + str(int(sza_min)) + 'to' + str(int(sza_max)) + '_codbin' + \
                str(show_specific_cod) + '.png'
        else:
            #outname = 'nn_force_scatter_combined_numsfcbins' + str(num_bins) + \
            outname = 'nn_force_scatter_combined_' + test_dict['sim_name'] + '_numsfcbins' + str(num_bins) + \
                '_sza' + str(int(sza_min)) + 'to' + str(int(sza_max)) + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:   
        plt.show()

    if(return_line_vals):
        return slope_vals, intpt_vals, mean_force_vals, std_force_vals, ai_calc_bins


def test_NN_forcing_daily_L2L3_errs(date_str, OMI_daily_data, \
        OMI_monthly_data, slope_dict, bin_dict, err_mean, err_std, error_type, \
        num_calcs, minlat = 65., maxlat = 87., \
        use_intercept = False, filter_bad_vals = True, \
        mod_L2_L3_error = None, \
        ai_thresh = -0.15, maxerr = 2., save = False):

    # Set up an array to hold all the results
    # ---------------------------------------
    full_array = np.full( (num_calcs, OMI_daily_data['grid_AI'].shape[1], \
        OMI_daily_data['grid_AI'].shape[2]), np.nan)
   
    for ii in range(num_calcs): 

        if(error_type == 'L2L3'):
            estimate_forcings, MYD08_data, NSIDC_data = \
                test_calculate_type_forcing_v4(OMI_daily_data, OMI_monthly_data, \
                slope_dict, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
                ai_thresh = ai_thresh, maxerr = maxerr,\
                filter_bad_vals = filter_bad_vals, \
                reference_ice = None, \
                reference_cld = None, \
                mod_slopes = None, \
                mod_L2_L3_error = mod_L2_L3_error, 
                L2L3_err_mean = err_mean, \
                L2L3_err_std = err_std, \
                return_modis_nsidc = True, \
                use_intercept = use_intercept)
        elif(error_type == 'ice'):
            estimate_forcings, MYD08_data, NSIDC_data = \
                test_calculate_type_forcing_v4(OMI_daily_data, OMI_monthly_data, \
                slope_dict, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
                ai_thresh = ai_thresh, maxerr = maxerr,\
                filter_bad_vals = filter_bad_vals, \
                reference_ice = None, \
                reference_cld = None, \
                mod_slopes = None, \
                ice_err_mean = err_mean,  \
                ice_err_std = err_std, \
                mod_L2_L3_error = mod_L2_L3_error, 
                L2L3_err_mean = None, \
                L2L3_err_std = None, \
                return_modis_nsidc = True, \
                use_intercept = use_intercept)

        full_array[ii,:,:] = estimate_forcings

    estimate_forcings, MYD08_data, NSIDC_data = \
        test_calculate_type_forcing_v4(OMI_daily_data, OMI_monthly_data, \
        slope_dict, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
        ai_thresh = ai_thresh, maxerr = maxerr,\
        filter_bad_vals = filter_bad_vals, \
        reference_ice = None, \
        reference_cld = None, \
        mod_slopes = None, \
        mod_L2_L3_error = None, 
        L2L3_err_mean = None, \
        L2L3_err_std = None, \
        return_modis_nsidc = True, \
        use_intercept = use_intercept)

    return estimate_forcings, full_array

# This also works for ice errors
# error_type: 'L2L3', 'ice', 'cod',...
def plot_NN_forcing_daily_L2L3_errors(date_str, OMI_daily_data, \
        OMI_monthly_data, slope_dict, bin_dict, err_mean, err_std, error_type, \
        num_calcs, minlat = 65., maxlat = 87., \
        use_intercept = False, filter_bad_vals = True, \
        mod_L2_L3_error = None, \
        ai_thresh = -0.15, maxerr = 2., save = False):


    estimate_forcings, full_array = \
        test_NN_forcing_daily_L2L3_errs(date_str, OMI_daily_data, \
        OMI_monthly_data, slope_dict, bin_dict, err_mean, err_std, error_type, \
        num_calcs, minlat = minlat, maxlat = maxlat, \
        use_intercept = use_intercept, filter_bad_vals = filter_bad_vals, \
        mod_L2_L3_error = mod_L2_L3_error, \
        ai_thresh = ai_thresh, maxerr = maxerr, save = False)
    
    full_masked = np.ma.masked_where(full_array == 0, full_array)
    estimate_mask = np.ma.masked_where(estimate_forcings == 0, estimate_forcings)
  
    #return estimate_forcings, full_masked

    plt.close('all') 
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax2 = ax.twinx()
    for ii in range(full_masked.shape[0]):
        local_array = full_masked[ii,:,:].compressed().flatten()
        hist = ax.hist(local_array, bins = 100, alpha = 0.25)
   
    hist = ax2.hist(estimate_mask.compressed().flatten(), 100, color = 'tab:blue')
    ax.set_title(date_str + ' - ' + error_type + ' Errors\nOriginal mean: ' + \
        str(np.round(np.mean(estimate_mask), 2)) + \
        '\nNew mean: ' + str(np.round(np.nanmean(full_masked), 2)) + \
        '\n# Runs = ' + str(num_calcs))
    ax.set_xlabel('Daily Forcing [W/m2]')
    ax.set_ylabel('Counts of Error Runs')
    ax2.set_ylabel('Counts of Daily Forcings (Blue)')
    fig.tight_layout()

    if(save):
        outname = 'daily_forcing_with_L2L3_err_numerr' + str(num_calcs) + '_test_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()
    

def plot_NN_forcing_daily(date_str, OMI_daily_data, OMI_monthly_data, \
        slope_dict, bin_dict, minlat = 65., maxlat = 87., \
        use_intercept = False, filter_bad_vals = True, \
        mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
        calc_from_bckgd = True, \
        sim_name = '', \
        ai_thresh = -0.15, maxerr = 2., save = False):

    # Load in the daily MODIS and NSIDC data
    # --------------------------------------
    file_strs = np.array([str(tval) for tval in OMI_daily_data['day_values']])
    if(not (date_str in file_strs)):
        # Return a nan array
        # ------------------
        print("WARNING: Date", date_str, "not in daily data. Returning nans")
        return np.full(OMI_daily_data['grid_AI'][10,:,:].shape, np.nan)
         
    match_idx = np.where(date_str == file_strs)[0][0]
    local_OMI_daily = np.ma.masked_where(\
        OMI_daily_data['count_AI'][match_idx,:,:] == 0, \
        OMI_daily_data['grid_AI'][match_idx,:,:])

    estimate_forcings, MYD08_data, NSIDC_data = \
        test_calculate_type_forcing_v4(OMI_daily_data, OMI_monthly_data, \
        slope_dict, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
        ai_thresh = ai_thresh, maxerr = maxerr,\
        filter_bad_vals = filter_bad_vals, \
        reference_ice = None, \
        reference_cld = None, \
        mod_slopes = None, \
        mod_L2_L3_error = mod_L2_L3_error, 
        L2L3_err_mean = L2L3_err_mean, \
        L2L3_err_std = L2L3_err_std, \
        return_modis_nsidc = True, \
        calc_from_bckgd = calc_from_bckgd, \
        use_intercept = use_intercept)

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')

    plt.close('all')
    fig = plt.figure(figsize = (7, 6))
    ax1 = fig.add_subplot(2,2,1, projection = ccrs.NorthPolarStereo(central_longitude = 0))
    ax2 = fig.add_subplot(2,2,2, projection = ccrs.NorthPolarStereo(central_longitude = 0))
    ax3 = fig.add_subplot(2,2,3, projection = ccrs.NorthPolarStereo(central_longitude = 0))
    ax4 = fig.add_subplot(2,2,4, projection = ccrs.NorthPolarStereo(central_longitude = 0))

    # Panel 1: OMI AI
    # ---------------
    mesh = ax1.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
        local_OMI_daily, shading = 'auto', transform = datacrs, \
        vmin = 0, vmax = 4, cmap = 'jet')
    cbar = fig.colorbar(mesh, ax = ax1, label = 'AI Perturbation')
    ax1.coastlines(color = 'lightgrey')
    ax1.set_extent([-180,180,65,90], datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_title('OMI UVAI Perturbation')

    # Panel 2: NSIDC ICE
    # ------------------
    #work_ice = np.ma.masked_where(NSIDC_data['grid_ice_conc'] > 80, NSIDC_data['grid_ice_conc'])
    work_ice = np.ma.masked_where(NSIDC_data['grid_ice_conc'] < 1, NSIDC_data['grid_ice_conc'])
    mesh = ax2.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
        work_ice, shading = 'auto', transform = datacrs, \
        vmin = 1, vmax = 100, cmap = 'ocean')
    cbar = fig.colorbar(mesh, ax = ax2, label = 'Sea Ice Conc [%]')
    ax2.coastlines()
    ax2.add_feature(cfeature.LAND, zorder = 100, edgecolor = 'k', facecolor = 'tab:grey', alpha = 1.0)
    ax2.set_extent([-180,180,65,90], datacrs)
    ax2.set_extent([-180,180,65,90], datacrs)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.set_title('SSMIS Sea Ice Conc.')

    ##!#keep_force = np.ma.masked_where(estimate_forcings == 0, estimate_forcings)
    ##!#ax3.hist(keep_force.compressed(), bins = 50)
    ##!#ax3.set_xlabel('Forcing [W/m2]')

    # Panel 3: MODIS COD
    # ------------------
    mesh = ax3.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
        MYD08_data['cod_mean'], shading = 'auto', transform = datacrs, \
        vmin = 0, vmax = 50, cmap = 'viridis')
    cbar = fig.colorbar(mesh, ax = ax3, label = 'Cloud Optical Depth')
    ax3.coastlines(color = 'lightgrey')
    ax3.set_extent([-180,180,65,90], datacrs)
    ax3.set_extent([-180,180,65,90], datacrs)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.set_title('Aqua MODIS COD')

    # Panel 4: Forcing
    # ----------------
    mesh = ax4.pcolormesh(OMI_monthly_data['LON'], OMI_monthly_data['LAT'], \
        estimate_forcings, shading = 'auto', transform = datacrs, \
        vmin = -50, vmax = 50, cmap = 'bwr')
    cbar = fig.colorbar(mesh, ax = ax4, label = 'Forcing [Wm$^{-2}$]')
    ax4.coastlines()
    ax4.set_extent([-180,180,65,90], datacrs)
    ax4.set_extent([-180,180,65,90], datacrs)
    ax4.set_boundary(circle, transform=ax4.transAxes)
    ax4.set_title('Aerosol Direct Forcing')

    plot_subplot_label(ax1, 'a)', fontsize = 12, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 12, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 12, backgroundcolor = None)
    plot_subplot_label(ax4, 'd)', fontsize = 12, backgroundcolor = None)
    
    plt.suptitle(dt_date_str.strftime('Daily Estimated Aerosol Direct ' + \
        'Radiative Forcing\n%Y-%m-%d'))

    fig.tight_layout()

    if(save):
        if(sim_name != ''):
            sim_name = '_' + sim_name
        if(calc_from_bckgd):
            bckgd_add = ''
        else:
            bckgd_add = '_noback'
        outname = 'test_forcing_v4_' + date_str + sim_name + bckgd_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()
    

#def test_calculate_type_forcing_v4(OMI_daily_data, OMI_monthly_data, slope_dict, \
#        bin_dict, date_str, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
#        maxerr = 2, 
#        reference_ice = None, reference_cld = None, mod_slopes = None, \
#        mod_intercepts = None, mod_ice = None, mod_cod = None, \
#        filter_bad_vals = True, return_modis_nsidc = True, use_intercept = False, \
#        debug = False):
def perform_forcing_calculation_v4(l2_data_dict, slope_dict, bin_dict, \
        ai_thresh = 0.7, mod_slopes = None, mod_intercepts = None, \
        mod_cod = None, mod_ice = None, use_intercept = False):

    local_OMI_data = l2_data_dict['omi_uvai_pert'][:,:]

    estimate_forcings = np.full(local_OMI_data.shape, np.nan)

    if(debug):
        print("OMI SHAPE = ", local_OMI_data.shape)
        #print("MOD SHAPE = ", MYD08_data['day_cld_frac_mean'].shape)
        #print("NSI SHAPE = ", NSIDC_data['grid_ice_conc'].shape)

    """
    final_pole_hole = np.where((grid_pole_hole_cc != -999.) & \
        ((grid_pole_hole_cc > grid_unused_cc) & \
         (grid_pole_hole_cc > grid_coastline_cc) & \
         (grid_pole_hole_cc > grid_land_cc)), \
        251, np.nan).squeeze()
    final_unused = np.where((grid_unused_cc != -999.) & \
        ((grid_unused_cc > grid_pole_hole_cc) & \
         (grid_unused_cc > grid_coastline_cc) & \
         (grid_unused_cc > grid_land_cc)), \
        252, np.nan).squeeze()
    final_coastline = np.where((grid_coastline_cc != -999.) & \
        ((grid_coastline_cc > grid_pole_hole_cc) & \
         (grid_coastline_cc > grid_unused_cc) & \
         (grid_coastline_cc > grid_land_cc)), \
        253, np.nan).squeeze()
    final_land = np.where((grid_land_cc != -999.) & \
        ((grid_land_cc >= grid_pole_hole_cc) & \
         (grid_land_cc >= grid_unused_cc) & \
         (grid_land_cc >= grid_coastline_cc)), \
        254, np.nan).squeeze()
    """

    # 2.  Loop over the x coordinate
    # ------------------------------
    for ii in range(local_OMI_data.shape[0]):

        # Loop over the y coordinate
        # --------------------------
        for jj in range(local_OMI_data.shape[1]):
    
            # 3.  Determine if a single gridbox has AI above the threshold
            if(local_OMI_data[ii,jj] < ai_thresh): 
                # 3a. If not above the threshold, set forcing to 0 and continue
                estimate_forcings[ii,jj] = 0.
            else:
                # 3b. If above the threshold,
                # 4.  Determine the difference between this pixel's AI and the clear-sky
                #     climatology.
                #delta_ai = local_OMI_daily[ii,jj] - \
                #           clear_sky_AI[ii,jj]

                # Extract the daily NSIDC and MODIS COD values here
                # Also, grab the pre-calculated 

                # 6.  Extract the MODIS MYD08 value and weigh the "clear" and "cloud"
                #     forcing values according to that cloud fraction.
                local_sza = l2_data_dict['omi_sza'][ii,jj]
                modis_cod = l2_data_dict['modis_cod'][ii,jj]
                nsidc_ice = l2_data_dict['nsidc_ice'][ii,jj]

                if(mod_cod is not None):
                    new_cod = modis_cod + mod_cod
                    if(new_cod < 0):
                        new_cod = 0
                    elif(new_cod > bin_dict['cod_bin_edges'][-1]):
                        new_cod = bin_dict['cod_bin_edges'][-1]

                    modis_cod = new_cod

                # NOTE. If the cloud fraction is very small (< 0.1, 0.2?) but the 
                # COD value is not 0, use a value of 0 for the COD
                #print('{0:5.1f} {1:5.1f} {2:5.1f} {3:5.2f} {4:5.2f} {5:5.2f} {6:5.2f}'.format(\
                #    OMI_daily_data['lat_values'][ii],OMI_daily_data['lon_values'][jj],\
                #    local_sza, modis_cod, cld_frac, local_OMI_daily[ii,jj], nsidc_ice))

                sza_idx = np.argmin(abs(local_sza - bin_dict['sza_bin_means']))
                cod_idx = np.argmin(abs(modis_cod - bin_dict['cod_bin_means']))

                #print(bin_dict['sza_bin_means'][sza_idx], bin_dict['cod_bin_means'][cod_idx])

                # 5.  Extract the NSIDC surface type to figure out which forcing value
                #     to use.
                ice_idx = -9
                if(nsidc_ice == -999.):
                #if( (pole_mask[ii,jj] == True) & (unused_mask[ii,jj] == True) & \
                #    (land_mask[ii,jj] == True) & (coast_mask[ii,jj] == True) & \
                #    (ice_mask[ii,jj] == True)):
                    if(debug):
                        print(pole_mask[ii,jj], unused_mask[ii,jj], \
                            land_mask[ii,jj], coast_mask[ii,jj], \
                            ice_mask[ii,jj], 'ALL MASKED')
                    #calc_forcing = 0.
                    ice_idx = -9
                else:
                    if( (nsidc_ice == 251.) | (nsidc_ice == 252.)):
                    #if((pole_mask[ii,jj] == False) | \
                    #   (unused_mask[ii,jj] == False)):
                        # Bad grid points. NO forcing
                        #calc_forcing =  0.
                        ice_idx = -9
                        if(debug):
                            print(pole_mask[ii,jj], unused_mask[ii,jj], \
                                land_mask[ii,jj], coast_mask[ii,jj], \
                                ice_mask[ii,jj], 'POLE/UNUSED')

                    elif((nsidc_ice == 253.) | \
                         (nsidc_ice == 254.)):
                    #elif((land_mask[ii,jj] == False) | \
                    #     (coast_mask[ii,jj] == False)):
                        # Use land forcing value
                        if(debug):
                            print(pole_mask[ii,jj], unused_mask[ii,jj], \
                                land_mask[ii,jj], coast_mask[ii,jj], \
                                ice_mask[ii,jj], 'LAND/COAST')

                        ice_idx = -1
                        #ice_idx = 3
                        # 2024/08/07 - currently unable to implement
                        #   this with the L2 validation stuff because
                        #   the GPQF data are not included in the
                        #   neuralnet output files.

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
                        ##!#if(OMI_daily_data['grid_GPQF'][match_idx,ii,jj] == 3):
                        ##!#    # Use the ice slope index
                        ##!#    ice_idx = 0

                        ##!#    #calc_forcing = cld_frac * \
                        ##!#    #    cloud_dict['ice_forcing'][ii] + \
                        ##!#    #    (1 - cld_frac) * clear_dict['ice_forcing'][ii]
                        ##!#else:
                        ##!#    # Use the land slope index 
                        ##!#    # CHANGED FROM EARLIER!!! 
                        ##!#    ice_idx = 3
                        ##!#    #calc_forcing = cld_frac * \
                        ##!#    #    cloud_dict['land_forcing'][ii] + \
                        ##!#    #    (1 - cld_frac) * clear_dict['land_forcing'][ii]
                        ##!##estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 

                    elif( (nsidc_ice >= 0) & (nsidc_ice <= 100)):
                    #elif(ice_mask[ii,jj] == False):
                        # Use either ice, mix, or ocean forcing value
                        if(debug):
                            print(pole_mask[ii,jj], unused_mask[ii,jj], \
                                land_mask[ii,jj], coast_mask[ii,jj], \
                                ice_mask[ii,jj], 'ICE/MIX/OCEAN')

                        # Account for if the user wants to test the sensitivity to ice
                        # ------------------------------------------------------------
                        if(mod_ice is not None):
                            new_ice = nsidc_ice + mod_ice
                            if(new_ice < 0):
                                new_ice = 0
                            elif(new_ice > 100):
                                new_ice = 100
                            #if((dt_date_str.year == 2015) & (dt_date_str.month == 7)):
                            #    print("OLD ICE:", nsidc_ice, "NEW ICE:", new_ice)

                            """
                            if((nsidc_ice < 20) & (new_ice >= 20)):
                                print("OCN TO MIX", np.round(nsidc_ice,1), np.round(new_ice,1), 'COD IDX = ', cod_idx)
                            elif((nsidc_ice >= 20) & (new_ice < 20)):
                                print("MIX TO OCN", np.round(nsidc_ice,1), np.round(new_ice,1), 'COD IDX = ', cod_idx)
                            elif((nsidc_ice < 80) & (new_ice >= 80)):
                                print("MIX TO ICE", np.round(nsidc_ice,1), np.round(new_ice,1), 'COD IDX = ', cod_idx)
                            elif((nsidc_ice >= 80) & (new_ice < 80)):
                                print("ICE TO MIX", np.round(nsidc_ice,1), np.round(new_ice,1), 'COD IDX = ', cod_idx)
                            """

                            nsidc_ice = new_ice

                        # USED FOR NEW MIX BINS
                        ice_idx = np.argmin(abs(nsidc_ice - bin_dict['ice_bin_means']))

                        """
                        if( (nsidc_ice < 20) ):
                            # Use ocean forcing
                            ice_idx = 0
                            #calc_forcing = cld_frac * cloud_dict['ocean_forcing'][ii] + \
                            #               (1 - cld_frac) * clear_dict['ocean_forcing'][ii]
                            #estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 
                            
                        elif( (nsidc_ice >= 20) & (nsidc_ice < 80)):
                            # Use mix forcing
                            ice_idx = 1
                            #calc_forcing = cld_frac * cloud_dict['mix_forcing'][ii] + \
                            #               (1 - cld_frac) * clear_dict['mix_forcing'][ii]
                            #estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 

                        elif( (nsidc_ice >= 80) ):
                            # Use ice forcing
                            ice_idx = 2
                            #calc_forcing = cld_frac * cloud_dict['ice_forcing'][ii] + \
                            #               (1 - cld_frac) * clear_dict['ice_forcing'][ii]
                            #estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 

                        else:
                            if(debug):
                                print("FAILED ICE")
                            #calc_forcing = 0.
                            ice_idx = -9
                        """

                    else:
                        if(debug):
                            print("FAILED EVERYTHING")
                        calc_forcing = 0.

                if(ice_idx != -9): 
                    #print("ICE INDEX = ", ice_idx)
                    #calc_forcing = cld_frac * cloud_dict['ice_forcing'][ii] + \
                    #               (1 - cld_frac) * clear_dict['ice_forcing'][ii]
                    #estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai
        
                    if(mod_slopes is None):
                        # Select the normal slope
                        work_slope = slope_dict['slopes'][ice_idx,sza_idx,cod_idx]
                    elif(mod_slopes == 'upper'):
                        if(slope_dict['trend_type'] == 'theil-sen'):
                            work_slope = slope_dict['upper_slopes'][ice_idx,sza_idx,cod_idx]
                        else:
                            work_slope = slope_dict['slopes'][ice_idx,sza_idx,cod_idx] + \
                                slope_dict['slope_stderr'][ice_idx,sza_idx,cod_idx]
                    elif(mod_slopes == 'lower'):
                        if(slope_dict['trend_type'] == 'theil-sen'):
                            work_slope = slope_dict['lower_slopes'][ice_idx,sza_idx,cod_idx]
                        else:
                            work_slope = slope_dict['slopes'][ice_idx,sza_idx,cod_idx] - \
                                slope_dict['slope_stderr'][ice_idx,sza_idx,cod_idx]
                    else:
                        print("WARNING: Invalid mod slopes. Must be " + \
                            "'upper' or 'lower'. Using default")
                        work_slope = slope_dict['slopes'][ice_idx,sza_idx,cod_idx]


                    if(use_intercept):
                        #slope_dict['slopes'][ice_idx,sza_idx,cod_idx] * \
                        if(mod_intercepts is None):
                            work_intercept = slope_dict['intercepts'][ice_idx,sza_idx,cod_idx]
                        else:
                            if(slope_dict['trend_type'] == 'linregress'):
                                if(mod_intercepts == 'upper'): 
                                    work_intercept = slope_dict['intercepts'][ice_idx,sza_idx,cod_idx] + \
                                                     slope_dict['intcpt_stderr'][ice_idx,sza_idx,cod_idx]
                                elif(mod_intercepts == 'lower'): 
                                    work_intercept = slope_dict['intercepts'][ice_idx,sza_idx,cod_idx] - \
                                                     slope_dict['intcpt_stderr'][ice_idx,sza_idx,cod_idx]
                                else:
                                    print("WARNING: Invalid mod intcpt. Must be " + \
                                        "'upper' or 'lower'. Using default")
                                    work_intercept = slope_dict['intercepts'][ice_idx,sza_idx,cod_idx]
                            else:
                                if(debug):
                                    print("ERROR: NO INTERCEPT ERROR FOR THEIL-SEN REGRESSION")
                                
                                
                        estimate_forcings[ii,jj] = \
                            work_slope * local_OMI_data[ii,jj] + \
                            work_intercept
                    else:
                        #slope_dict['slopes'][ice_idx,sza_idx,cod_idx] * \
                        estimate_forcings[ii,jj] = \
                            work_slope * local_OMI_data[ii,jj]

    estimate_forcings = np.ma.masked_invalid(estimate_forcings) 

    return estimate_forcings


# use_intercept: in the forcing estimation calculation, setting this to True 
#   makes the calculation use both the regression slope and regression to 
#   calculate the estiated forcing. Doing this makes a calculation of:
#
#       forcing = slope * AI + intercept
#
# Setting this to False means only the slope is used, following:
#
#       forcing = slope * AI 
#
#
# mod_slopes: 'upper', 'lower'
# mod_ice: -5, 5, -15, 15
# calc_from_bckgd: determines the rules by which forcing is calculated for
#                  a grid box. If set to True (default), then the code determines if
#                  the daily UVAI is greater than the OMI background by
#                  more than the threshold value. If set to False, then the
#                  code determines if the daily UVAI value is greater than
#                  the threshold value.
# ----------------------------------------------------------------------------
def test_calculate_type_forcing_v4(OMI_daily_data, OMI_monthly_data, slope_dict, \
        bin_dict, date_str, minlat = 70., maxlat = 87., ai_thresh = 0.7, \
        maxerr = 2, 
        reference_ice = None, reference_cld = None, mod_slopes = None, \
        mod_intercepts = None, calc_from_bckgd = True, \
        mod_ice = None, ice_err_mean = None, ice_err_std = None, \
        mod_cod = None, cod_err_mean = None, cod_err_std = None, \
        mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
        filter_bad_vals = True, return_modis_nsidc = True, use_intercept = False, \
        debug = False):

    ####################################
    # BEGIN CODE TRANSPLANT FROM calculate_type_forcing_v3
    ####################################

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')   
 
    # Load in the daily MODIS and NSIDC data
    # --------------------------------------
    file_strs = np.array([str(tval) for tval in OMI_daily_data['day_values']])
    if(not (date_str in file_strs)):
        # Return a nan array
        # ------------------
        print("WARNING: Date", date_str, "not in daily data. Returning nans")
        return np.full(OMI_daily_data['grid_AI'][10,:,:].shape, np.nan)
         
    match_idx = np.where(date_str == file_strs)[0][0]
    local_OMI_daily = np.ma.masked_where(\
        OMI_daily_data['count_AI'][match_idx,:,:] == 0, \
        OMI_daily_data['grid_AI'][match_idx,:,:])
    
    # tidx is the "month_idx"
    # -----------------------    
    tidx = int(date_str[4:6]) - 4
    
    if(reference_cld is None):
        cld_str = date_str
    else:
        cld_str = reference_cld + date_str[4:]
    
    MYD08_data = read_MODIS_MYD08_single(cld_str, minlat = minlat, \
        maxlat = maxlat)

    # It appears that the MODIS COD values are missing in clear-sky
    # conditions, so setting the missing COD value to 0 for now.    
    # Not sure if this can be done, but am testing here
    # -------------------------------------------------------------
    MYD08_data['cod_mean'] = np.where(MYD08_data['cod_mean'].mask == True, 0, MYD08_data['cod_mean'])
   
    print('max MYD08 cod', np.min(MYD08_data['cod_mean']), np.max(MYD08_data['cod_mean']))
 
    # Load in the single-day NSIDC ice concentration
    #
    # If the user wants to use a different year's ice values as a reference,
    # modify the date string here accordingly
    # ------------------------------------------------------------------------
    if(reference_ice is None):
        ice_str = date_str
    else:
        ice_str = reference_ice + date_str[4:]
    
    NSIDC_data =  readNSIDC_daily(ice_str, maxlat = maxlat)
    NSIDC_data = grid_data_conc(NSIDC_data, minlat = minlat, maxlat = maxlat)
    NSIDC_data['grid_ice_conc'] = np.ma.masked_where((NSIDC_data['grid_ice_conc'] < 0) | \
        (NSIDC_data['grid_ice_conc'] > 100), NSIDC_data['grid_ice_conc']).squeeze()

    clear_sky_AI = np.array([np.nanmean(\
        np.ma.masked_where(OMI_monthly_data['AI'][midx::6,:,:] > ai_thresh,\
        OMI_monthly_data['AI'][midx::6,:,:]), axis = 0) for midx in range(6)])
    clear_sky_AI = clear_sky_AI[tidx,:,:]

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

    print("In calc, calc_from_bckgd: ", calc_from_bckgd)

    # 2.  Loop over each individual month, over the latitudes first
    # -----------------------------------
    num_good_points = 0
    for ii in range(local_OMI_daily.shape[0]):

        # Calculate the solar zenith angle here
        # -------------------------------------
        day_idx = int(dt_date_str.strftime('%j'))
        del_angle = -23.45 * np.cos( np.radians((360 / 365) * (day_idx + 10)))
        local_sza = NSIDC_data['grid_lat'][ii,0] - del_angle
        #print('DATE:',date_str, ' LAT:', NSIDC_data['grid_lat'][ii,0],' SZA:',local_sza)
        #latitude = NSIDC_data['grid_lat'][ii,0]
        #lat_szas = np.array([tlat - del_angles for tlat in latitudes]).T

        ## Extract the calculated SZA for this latitude band
        #local_sza = lat_szas[day_idx,ii]

        # Loop over the longitudes
        for jj in range(local_OMI_daily.shape[1]):
    
            # 3b. If above the threshold,
            # 4.  Determine the difference between this pixel's AI and the clear-sky
            #     climatology.
            delta_ai = local_OMI_daily[ii,jj] - \
                       clear_sky_AI[ii,jj]

            # 3.  Determine if a single gridbox has AI that is above the clear-sky
            #       background by a threshold
            #if(local_OMI_daily[ii,jj] < ai_thresh): 
            #if(delta_ai < ai_thresh): 
            # NOTE: New for noland105_noback (2024/11/06)
            if( (calc_from_bckgd and (delta_ai < ai_thresh)) or \
                ((not calc_from_bckgd) and (local_OMI_daily[ii,jj] < ai_thresh)) ): 
                # 3a. If not above the threshold, set forcing to 0 and continue
                estimate_forcings[ii,jj] = 0.
            else:

                # Extract the daily NSIDC and MODIS COD values here
                # Also, grab the pre-calculated 

                # 6.  Extract the MODIS MYD08 value and weigh the "clear" and "cloud"
                #     forcing values according to that cloud fraction.
                #cld_frac = MYD08_data['day_cld_frac_mean'][ii,jj]
                modis_cod = MYD08_data['cod_mean'][ii,jj]
                #cld_frac  = MYD08_data['day_cld_frac_mean'][ii,jj]
                nsidc_ice = NSIDC_data['grid_ice_conc'][ii,jj]

                if(mod_cod is not None):
                    new_cod = modis_cod + mod_cod
                    if(new_cod < 0):
                        new_cod = 0
                    elif(new_cod > bin_dict['cod_bin_edges'][-1]):
                        new_cod = bin_dict['cod_bin_edges'][-1]

                    modis_cod = new_cod

                # Instead, see if the user wants to modify the COD value 
                # by COD errors that follow an error distribution
                # ----------------------------------------------------------
                elif(cod_err_mean is not None):
                    error_val = np.random.normal(cod_err_mean, cod_err_std)
                    new_cod = modis_cod + error_val

                    if(new_cod < 0):
                        new_cod = 0
                    elif(new_cod > bin_dict['cod_bin_edges'][-1]):
                        new_cod = bin_dict['cod_bin_edges'][-1]

                    print("ADDING COD ERROR OF", error_val, 'NEW COD = ', new_cod)

                    modis_cod = new_cod


                # NOTE. If the cloud fraction is very small (< 0.1, 0.2?) but the 
                # COD value is not 0, use a value of 0 for the COD
                #print('{0:5.1f} {1:5.1f} {2:5.1f} {3:5.2f} {4:5.2f} {5:5.2f} {6:5.2f}'.format(\
                #    OMI_daily_data['lat_values'][ii],OMI_daily_data['lon_values'][jj],\
                #    local_sza, modis_cod, cld_frac, local_OMI_daily[ii,jj], nsidc_ice))

                sza_idx = np.argmin(abs(local_sza - bin_dict['sza_bin_means']))
                cod_idx = np.argmin(abs(modis_cod - bin_dict['cod_bin_means']))

                #print(bin_dict['sza_bin_means'][sza_idx], bin_dict['cod_bin_means'][cod_idx])

                # 5.  Extract the NSIDC surface type to figure out which forcing value
                #     to use.
                ice_idx = -9
                if( (pole_mask[ii,jj] == True) & (unused_mask[ii,jj] == True) & \
                    (land_mask[ii,jj] == True) & (coast_mask[ii,jj] == True) & \
                    (ice_mask[ii,jj] == True)):
                    if(debug):
                        print(pole_mask[ii,jj], unused_mask[ii,jj], \
                            land_mask[ii,jj], coast_mask[ii,jj], \
                            ice_mask[ii,jj], 'ALL MASKED')
                    #calc_forcing = 0.
                    ice_idx = -9
                else:
                    if((pole_mask[ii,jj] == False) | \
                       (unused_mask[ii,jj] == False)):
                        # Bad grid points. NO forcing
                        #calc_forcing =  0.
                        ice_idx = -9
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
                            # Use the ice slope index
                            ice_idx = 0

                            #calc_forcing = cld_frac * \
                            #    cloud_dict['ice_forcing'][ii] + \
                            #    (1 - cld_frac) * clear_dict['ice_forcing'][ii]
                        else:
                            # Use the land slope index 
                            # CHANGED FROM EARLIER!!! 
                            #ice_idx = 3
                            ice_idx = -1
                            #calc_forcing = cld_frac * \
                            #    cloud_dict['land_forcing'][ii] + \
                            #    (1 - cld_frac) * clear_dict['land_forcing'][ii]
                        #estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 

                    elif(ice_mask[ii,jj] == False):
                        # Use either ice, mix, or ocean forcing value
                        if(debug):
                            print(pole_mask[ii,jj], unused_mask[ii,jj], \
                                land_mask[ii,jj], coast_mask[ii,jj], \
                                ice_mask[ii,jj], 'ICE/MIX/OCEAN')

                        # Account for if the user wants to test the sensitivity to ice
                        # ------------------------------------------------------------
                        if(mod_ice is not None):
                            new_ice = nsidc_ice + mod_ice
                            if(new_ice < 0):
                                new_ice = 0
                            elif(new_ice > 100):
                                new_ice = 100
                            #if((dt_date_str.year == 2015) & (dt_date_str.month == 7)):
                            #    print("OLD ICE:", nsidc_ice, "NEW ICE:", new_ice)

                            """
                            if((nsidc_ice < 20) & (new_ice >= 20)):
                                print("OCN TO MIX", np.round(nsidc_ice,1), np.round(new_ice,1), 'COD IDX = ', cod_idx)
                            elif((nsidc_ice >= 20) & (new_ice < 20)):
                                print("MIX TO OCN", np.round(nsidc_ice,1), np.round(new_ice,1), 'COD IDX = ', cod_idx)
                            elif((nsidc_ice < 80) & (new_ice >= 80)):
                                print("MIX TO ICE", np.round(nsidc_ice,1), np.round(new_ice,1), 'COD IDX = ', cod_idx)
                            elif((nsidc_ice >= 80) & (new_ice < 80)):
                                print("ICE TO MIX", np.round(nsidc_ice,1), np.round(new_ice,1), 'COD IDX = ', cod_idx)
                            """

                            nsidc_ice = new_ice

                        # Instead, see if the user wants to modify the ice value 
                        # by ice errors that follow an error distribution
                        # ----------------------------------------------------------
                        elif(ice_err_mean is not None):
                            error_val = np.random.normal(ice_err_mean, ice_err_std)
                            new_ice = nsidc_ice + error_val

                            if(new_ice < 0):
                                new_ice = 0
                            elif(new_ice > 100):
                                new_ice = 100

                            print("ADDING ICE ERROR OF", error_val, 'NEW ICE = ', new_ice)

                            nsidc_ice = new_ice

                        # USED FOR NEW MIX BINS

                        if(len(bin_dict['ice_bin_means']) == 6):
                            ice_idx = np.argmin(abs(nsidc_ice - bin_dict['ice_bin_means']))
                            #print("USING NEW STYLE ICE VALUES") 
                        else: 
                            #print("USING OLD STYLE ICE VALUES") 
                            if( (nsidc_ice < 20) ):
                                # Use ocean forcing
                                ice_idx = 0
                                #calc_forcing = cld_frac * cloud_dict['ocean_forcing'][ii] + \
                                #               (1 - cld_frac) * clear_dict['ocean_forcing'][ii]
                                #estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 
                                
                            elif( (nsidc_ice >= 20) & (nsidc_ice < 80)):
                                # Use mix forcing
                                ice_idx = 1
                                #calc_forcing = cld_frac * cloud_dict['mix_forcing'][ii] + \
                                #               (1 - cld_frac) * clear_dict['mix_forcing'][ii]
                                #estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 

                            elif( (nsidc_ice >= 80) ):
                                # Use ice forcing
                                ice_idx = 2
                                #calc_forcing = cld_frac * cloud_dict['ice_forcing'][ii] + \
                                #               (1 - cld_frac) * clear_dict['ice_forcing'][ii]
                                #estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai 

                            else:
                                if(debug):
                                    print("FAILED ICE")
                                #calc_forcing = 0.
                                ice_idx = -9

                    else:
                        if(debug):
                            print("FAILED EVERYTHING")
                        calc_forcing = 0.

                if(ice_idx != -9): 
                    num_good_points += 1
                    #print("ICE VAL:", nsidc_ice, "ICE INDEX = ", ice_idx)
                    #calc_forcing = cld_frac * cloud_dict['ice_forcing'][ii] + \
                    #               (1 - cld_frac) * clear_dict['ice_forcing'][ii]
                    #estimate_forcings[ii,jj] = 0. - calc_forcing * delta_ai
        
                    if(mod_slopes is None):
                        # Select the normal slope
                        work_slope = slope_dict['slopes'][ice_idx,sza_idx,cod_idx]
                    elif(mod_slopes == 'upper'):
                        if(slope_dict['trend_type'] == 'theil-sen'):
                            work_slope = slope_dict['upper_slopes'][ice_idx,sza_idx,cod_idx]
                        else:
                            work_slope = slope_dict['slopes'][ice_idx,sza_idx,cod_idx] + \
                                slope_dict['slope_stderr'][ice_idx,sza_idx,cod_idx]
                    elif(mod_slopes == 'lower'):
                        if(slope_dict['trend_type'] == 'theil-sen'):
                            work_slope = slope_dict['lower_slopes'][ice_idx,sza_idx,cod_idx]
                        else:
                            work_slope = slope_dict['slopes'][ice_idx,sza_idx,cod_idx] - \
                                slope_dict['slope_stderr'][ice_idx,sza_idx,cod_idx]
                    else:
                        print("WARNING: Invalid mod slopes. Must be " + \
                            "'upper' or 'lower'. Using default")
                        work_slope = slope_dict['slopes'][ice_idx,sza_idx,cod_idx]


                    if(use_intercept):
                        #slope_dict['slopes'][ice_idx,sza_idx,cod_idx] * \
                        if(mod_intercepts is None):
                            work_intercept = slope_dict['intercepts'][ice_idx,sza_idx,cod_idx]
                        else:
                            if(slope_dict['trend_type'] == 'linregress'):
                                if(mod_intercepts == 'upper'): 
                                    work_intercept = slope_dict['intercepts'][ice_idx,sza_idx,cod_idx] + \
                                                     slope_dict['intcpt_stderr'][ice_idx,sza_idx,cod_idx]
                                elif(mod_intercepts == 'lower'): 
                                    work_intercept = slope_dict['intercepts'][ice_idx,sza_idx,cod_idx] - \
                                                     slope_dict['intcpt_stderr'][ice_idx,sza_idx,cod_idx]
                                else:
                                    print("WARNING: Invalid mod intcpt. Must be " + \
                                        "'upper' or 'lower'. Using default")
                                    work_intercept = slope_dict['intercepts'][ice_idx,sza_idx,cod_idx]
                            else:
                                if(debug):
                                    print("ERROR: NO INTERCEPT ERROR FOR THEIL-SEN REGRESSION")
                                
                                
                        estimate_forcings[ii,jj] = \
                            work_slope * local_OMI_daily[ii,jj] + \
                            work_intercept
                    else:
                        #slope_dict['slopes'][ice_idx,sza_idx,cod_idx] * \
                        estimate_forcings[ii,jj] = \
                            work_slope * local_OMI_daily[ii,jj]

                    # Now, see if the user wants to modify the resulting
                    # forcings by the L2/L3 errors
                    # --------------------------------------------------
                    if(mod_L2_L3_error is not None):
                        estimate_forcings[ii,jj] = estimate_forcings[ii,jj] + mod_L2_L3_error

                    # Instead, see if the user wants to modify the resulting
                    # forcings by L2/L3 errors that follow an error distribution
                    # ----------------------------------------------------------
                    elif(L2L3_err_mean is not None):
                        error_val = np.random.normal(L2L3_err_mean, L2L3_err_std)
                        #print("ADDING ERROR OF", error_val)
                        estimate_forcings[ii,jj] = estimate_forcings[ii,jj] + error_val
    


    estimate_forcings = np.ma.masked_invalid(estimate_forcings) 

    # Here, figure out if there are less than a threshold (20?) points
    # with calculated forcings. If so, this may be noise that made it
    # through the system. Mask the daily values here if that's the case?
    #num_good_points = estimate_forcings.compressed().shape[0]

    print("BERFR: MIN MAX FORCINGS", np.round(np.nanmin(estimate_forcings),1), \
        np.round(np.nanmax(estimate_forcings), 1)) 
    if(filter_bad_vals):
        mean_val = np.nanmean(estimate_forcings)
        std_val  = np.nanstd(estimate_forcings)

        estimate_forcings = np.ma.masked_where(estimate_forcings > \
            (mean_val + 8.0 * std_val), estimate_forcings)

    print("AFTER: MIN MAX FORCINGS", np.round(np.nanmin(estimate_forcings),1), \
        np.round(np.nanmax(estimate_forcings), 1), \
        'MAX AI', np.round(np.max(local_OMI_daily),1), \
        'NUM GOOD POINTS', num_good_points)


    if(return_modis_nsidc):
        return estimate_forcings, MYD08_data, NSIDC_data
    else:
        return estimate_forcings

# This verison is set up to use daily data, but still needs the monthly
# OMI data to determine the clear-sky climatology.
# ---------------------------------------------------------------------
def calculate_type_forcing_v4_monthly(OMI_daily_data, OMI_monthly_data, \
        slope_dict, bin_dict, month_idx, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        maxerr = 2, \
        reference_ice = None, reference_cld = None, mod_slopes = None, \
        mod_intercepts = None, mod_ice = None, mod_cod = None, \
        mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
        filter_bad_vals = True, return_modis_nsidc = True, \
        use_intercept = False, debug = False):


    # Set up arrays to hold the monthly-averaged forcing calculations
    # ---------------------------------------------------------------
    xdim = OMI_daily_data['grid_AI'].shape[1]
    ydim = OMI_daily_data['grid_AI'].shape[2]
    daily_force_vals = np.full( (31, xdim, ydim), np.nan)

    if(str(month_idx) == 'all'):
        l_all_months = True
        month_force_vals = np.full(\
            (OMI_monthly_data['DATES'].shape[0], xdim, ydim), \
            np.nan)

        begin_date_str = datetime.strptime(\
            OMI_monthly_data['DATES'][0], '%Y%m')
        end_date_str = datetime.strptime(\
            OMI_monthly_data['DATES'][-1], '%Y%m') + \
            relativedelta(months = 1) - timedelta(days = 1)
    else:
        l_all_months = False
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
    #local_date_str = datetime(2015,7,1)
    #end_date_str = datetime(2015,7,30)
    while(local_date_str <= end_date_str):

        date_str = local_date_str.strftime('%Y%m%d')

        print(date_str, day_count, month_count)

        # Calculate the forcing value for this current day
        # ------------------------------------------------
        #estimate_forcing = \
        #    calculate_type_forcing_v3(OMI_daily_data, OMI_monthly_data, \
        #        coloc_dict, date_str, minlat = minlat, maxlat = maxlat, \
        #        ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr,\
        #        min_cloud = min_cloud, data_type = data_type,\
        #        filter_bad_vals = filter_bad_vals, \
        #        reference_ice = reference_ice, \
        #        reference_cld = reference_cld, \
        #        mod_slopes = mod_slopes, \
        #        return_modis_nsidc = False)

        estimate_forcing = \
            test_calculate_type_forcing_v4(OMI_daily_data, OMI_monthly_data, \
            slope_dict, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
            ai_thresh = ai_thresh, maxerr = maxerr,\
            filter_bad_vals = filter_bad_vals, \
            reference_ice = reference_ice, \
            reference_cld = reference_cld, \
            mod_slopes = mod_slopes, \
            mod_intercepts = mod_intercepts, \
            mod_ice = mod_ice, \
            mod_cod = mod_cod, \
            mod_L2_L3_error = mod_L2_L3_error, \
            L2L3_err_mean = L2L3_err_mean, \
            L2L3_err_std  = L2L3_err_std, \
            return_modis_nsidc = False,\
            use_intercept = use_intercept)

        # Insert the values into the daily holding array
        # ----------------------------------------------
        daily_force_vals[day_count,:,:] = estimate_forcing[:,:] 

        # Increment working date
        # ----------------------
        new_work_date = local_date_str + timedelta(days = 1)

        # If the new working date has a different month than
        # the previous working date, then average the daily values
        # and move the working date to the next desired month
        # --------------------------------------------------------
        if(new_work_date.month != local_date_str.month):
            month_force_vals[month_count,:,:] = np.nanmean(\
                daily_force_vals[:,:,:], axis = 0)
            day_count = 0
            month_count += 1            

            if(l_all_months):
                if(new_work_date.month == 10):
                    new_work_date = new_work_date + relativedelta(months = 6)
            else:
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

# This verison is set up to use daily data, but still needs the monthly
# OMI data to determine the clear-sky climatology.
# 
# Should be the same as calculate_type_forcing_v4_monthly, but
# keeps the daily averages and doesn't calculate monthly averages
# ---------------------------------------------------------------------
def calculate_type_forcing_v4_alldaily(OMI_daily_data, OMI_monthly_data, \
        slope_dict, bin_dict, month_idx, minlat = 70., maxlat = 87., ai_thresh = -0.15, \
        maxerr = 2, \
        reference_ice = None, reference_cld = None, mod_slopes = None, \
        mod_intercepts = None, \
        mod_ice = None, ice_err_mean = None, ice_err_std = None, \
        mod_cod = None,  cod_err_mean = None, cod_err_std = None, \
        calc_from_bckgd = True, \
        mod_L2_L3_error = None, L2L3_err_mean = None, L2L3_err_std = None, \
        filter_bad_vals = True, return_modis_nsidc = True, \
        use_intercept = False, debug = False):


    # Set up arrays to hold the monthly-averaged forcing calculations
    # ---------------------------------------------------------------
    xdim = OMI_daily_data['grid_AI'].shape[1]
    ydim = OMI_daily_data['grid_AI'].shape[2]
    tdim = OMI_daily_data['day_values'].shape[0]
    daily_force_vals = np.full( (tdim, xdim, ydim), np.nan)

    print("calc_from_bckgd: ", calc_from_bckgd)

    day_count = 0
    month_count = 0
    #local_date_str = datetime(2015,7,1)
    #end_date_str = datetime(2015,7,30)
    #while(local_date_str <= end_date_str):
    for ii in range(tdim):

        date_str = str(OMI_daily_data['day_values'][ii])

        print(date_str, day_count, month_count, 'DAILY CALC')

        # Calculate the forcing value for this current day
        # ------------------------------------------------
        #estimate_forcing = \
        #    calculate_type_forcing_v3(OMI_daily_data, OMI_monthly_data, \
        #        coloc_dict, date_str, minlat = minlat, maxlat = maxlat, \
        #        ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr,\
        #        min_cloud = min_cloud, data_type = data_type,\
        #        filter_bad_vals = filter_bad_vals, \
        #        reference_ice = reference_ice, \
        #        reference_cld = reference_cld, \
        #        mod_slopes = mod_slopes, \
        #        return_modis_nsidc = False)

        estimate_forcing = \
            test_calculate_type_forcing_v4(OMI_daily_data, OMI_monthly_data, \
            slope_dict, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
            ai_thresh = ai_thresh, maxerr = maxerr,\
            filter_bad_vals = filter_bad_vals, \
            reference_ice = reference_ice, \
            reference_cld = reference_cld, \
            mod_slopes = mod_slopes, \
            mod_intercepts = mod_intercepts, \
            mod_ice = mod_ice, \
            ice_err_mean = ice_err_mean, \
            ice_err_std  = ice_err_std, \
            mod_cod = mod_cod, \
            cod_err_mean = cod_err_mean, \
            cod_err_std  = cod_err_std, \
            mod_L2_L3_error = mod_L2_L3_error, \
            L2L3_err_mean = L2L3_err_mean, \
            L2L3_err_std  = L2L3_err_std, \
            calc_from_bckgd = calc_from_bckgd, \
            return_modis_nsidc = False,\
            use_intercept = use_intercept)

        # Insert the values into the daily holding array
        # ----------------------------------------------
        daily_force_vals[ii,:,:] = estimate_forcing[:,:] 

        """
        # Increment working date
        # ----------------------
        new_work_date = local_date_str + timedelta(days = 1)

        # If the new working date has a different month than
        # the previous working date, then average the daily values
        # and move the working date to the next desired month
        # --------------------------------------------------------
        if(new_work_date.month != local_date_str.month):
            month_force_vals[month_count,:,:] = np.nanmean(\
                daily_force_vals[:,:,:], axis = 0)
            day_count = 0
            month_count += 1            

            if(l_all_months):
                if(new_work_date.month == 10):
                    new_work_date = new_work_date + relativedelta(months = 6)
            else:
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
        """

    return daily_force_vals


def plot_NN_architecture(plot_lines = False, save = False):

    # Set up the figure
    # -----------------
    fig = plt.figure(figsize = (11, 7))
    ax = fig.add_subplot(1,1,1)

    # Set up the first layer, which is the input layer
    # ------------------------------------------------
    layer1_x  = [0] * 7
    layer1_y  = list(np.arange(7) + 28.5)
    ax.scatter(layer1_x, layer1_y, color = 'tab:blue')
    ax.text(0,70,'input',color = 'tab:blue', weight = 'bold', horizontalalignment = 'center')
    ax.text(0,68,'n = 7',color = 'tab:blue', weight = 'bold', horizontalalignment = 'center')

    # Set up the hidden layers
    # ------------------------
    delta_x = 1.0
    num_hidden = 11
    beg_hidden = 0 + delta_x
    end_hidden = beg_hidden + num_hidden * delta_x
    xvals = np.arange(beg_hidden, end_hidden, delta_x)
    yvals = np.array([8,12,16,24,32,64,32,24,16,12,8])

    print(xvals, yvals)
    print(len(xvals), len(yvals))

    for ii in range(len(xvals)):
        offset = int((64 - yvals[ii]) / 2)

        layerh_x = [xvals[ii]] * yvals[ii]
        layerh_y = list(np.arange(yvals[ii]) + offset)

        if(plot_lines):
            if(ii == 0):
                for jj in range(len(layer1_x)):
                    for kk in range(len(layerh_x)):
                        print(layer1_x[jj], layer1_y[jj], layerh_x[kk], layerh_y[kk])
                        ax.plot([layer1_x[jj], layerh_x[kk]], [layer1_y[jj], layerh_y[kk]], linewidth = 1, color = 'black', alpha = 0.1)
            else:
                prev_offset = int((64 - yvals[ii - 1]) / 2)
                prev_layerh_x = [xvals[ii - 1]] * yvals[ii - 1]
                prev_layerh_y = list(np.arange(yvals[ii - 1]) + prev_offset)

                for jj in range(len(prev_layerh_x)):
                    for kk in range(len(layerh_x)):
                        #print(prev_layerh_x[jj], layer1_y[jj], prev_layerh_x[kk], layerh_y[kk])
                        ax.plot([prev_layerh_x[jj], layerh_x[kk]], \
                            [prev_layerh_y[jj], layerh_y[kk]], linewidth = 1, color = 'black', alpha = 0.1)

        ax.scatter(layerh_x, layerh_y, color = 'tab:gray')

        ax.text(xvals[ii],70,'hidden' + str(int(ii + 1)),\
            color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
        ax.text(xvals[ii],68,'n = ' + str(int(yvals[ii])),\
            color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    
    ##!#layer2_x  = [xvals[0]] * 8
    ##!#layer3_x  = [xvals[1]] * 12
    ##!#layer4_x  = [xvals[2]] * 16
    ##!#layer5_x  = [xvals[3]] * 24
    ##!#layer6_x  = [xvals[4]] * 32
    ##!#layer7_x  = [xvals[5]] * 64
    ##!#layer8_x  = [xvals[6]] * 32
    ##!#layer9_x  = [xvals[7]] * 24
    ##!#layer10_x = [xvals[8]] * 16
    ##!#layer11_x = [xvals[9]] * 12
    ##!#layer12_x = [xvals[10]] * 8

    ##!#layer2_y  = list(np.arange(8) + 28)
    ##!#layer3_y  = list(np.arange(12) + 26)
    ##!#layer4_y  = list(np.arange(16) + 24)
    ##!#layer5_y  = list(np.arange(24) + 20)
    ##!#layer6_y  = list(np.arange(32) + 16)
    ##!#layer7_y  = list(np.arange(64))
    ##!#layer8_y  = list(np.arange(32) + 16)
    ##!#layer9_y  = list(np.arange(24) + 20)
    ##!#layer10_y = list(np.arange(16) + 24)
    ##!#layer11_y = list(np.arange(12) + 26)
    ##!#layer12_y = list(np.arange(8) + 28)

    ##!#total_x = layer2_x + layer3_x + layer4_x + layer5_x + \
    ##!#          layer6_x + layer7_x + layer8_x + layer9_x + layer10_x + \
    ##!#          layer11_x + layer12_x

    ##!#total_y = layer2_y + layer3_y + layer4_y + layer5_y + \
    ##!#          layer6_y + layer7_y + layer8_y + layer9_y + layer10_y + \
    ##!#          layer11_y + layer12_y

    #ax.scatter(total_x, total_y, color = 'tab:gray')


    ###ax.text(xvals[0],70,'hidden1',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[0],68,'n = 8',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[1],70,'hidden2',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[1],68,'n = 12',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[2],70,'hidden3',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[2],68,'n = 16',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[3],70,'hidden4',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[3],68,'n = 24',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[4],70,'hidden5',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[4],68,'n = 32',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[5],70,'hidden6',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[5],68,'n = 64',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[6],70,'hidden7',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[6],68,'n = 32',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[7],70,'hidden8',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[7],68,'n = 24',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[8],70,'hidden9',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[8],68,'n = 16',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[9],70,'hidden10',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[9],68,'n = 12',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[10],70,'hidden11',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')
    ###ax.text(xvals[10],68,'n = 8',color = 'tab:grey', weight = 'bold', horizontalalignment = 'center')

    # Set up the final layer, which is the output layer
    # -------------------------------------------------
    xval = end_hidden
    layer13_x = [xval] * 1
    layer13_y = [32]

    if(plot_lines):
        for jj in range(len(layerh_x)):
            for kk in range(len(layer13_x)):
                #print(prev_layerh_x[jj], layer1_y[jj], prev_layerh_x[kk], layerh_y[kk])
                ax.plot([layerh_x[jj], layer13_x[kk]], \
                    [layerh_y[jj], layer13_y[kk]], linewidth = 1, color = 'black', alpha = 0.1)

    ax.scatter(layer13_x, layer13_y, color = 'tab:red')
    ax.text(xval,70,'output',color = 'tab:red', weight = 'bold', horizontalalignment = 'center')
    ax.text(xval,68,'n = 1',color = 'tab:red', weight = 'bold', horizontalalignment = 'center')

    ax.axis('off')
    fig.tight_layout()
    if(save):
        if(plot_lines):
            line_add = '_lines'
        else:
            line_add = ''
        outname = 'nn_architecture' + line_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_test_forcing_v4(OMI_daily_data, OMI_month_data, date_str, \
        coloc_dict, minlat = 65., maxlat = 87., ai_thresh = -0.15, \
        cld_idx = 0, maxerr = 2, min_cloud = 0.95, data_type = 'raw', \
        save = False, filter_bad_vals = False):

    print("HERE1", date_str, minlat, maxlat, ai_thresh, cld_idx, maxerr, min_cloud, data_type)
    estimate_forcing, MYD08_data, NSIDC_data, clear_sky_AI = \
        calculate_type_forcing_v3(OMI_daily_data, OMI_month_data, \
            coloc_dict, date_str, minlat = minlat, maxlat = maxlat, \
            ai_thresh = ai_thresh, cld_idx = cld_idx, maxerr = maxerr,\
            min_cloud = min_cloud, data_type = data_type,\
            filter_bad_vals = filter_bad_vals, return_modis_nsidc = True)

    estimate_forcing, MYD08_data, NSIDC_data,  = \
        test_calculate_type_forcing_v4(OMI_daily_data, OMI_monthly_data, \
        slope_dict, bin_dict, date_str, minlat = minlat, maxlat = maxlat, \
        ai_thresh = ai_thresh, maxerr = maxerr,\
        filter_bad_vals = filter_bad_vals, \
        reference_ice = reference_ice, \
        reference_cld = reference_cld, \
        mod_slopes = mod_slopes, \
        return_modis_nsidc = False,\
        use_intercept = use_intercept)





    dt_date_str = datetime.strptime(date_str, '%Y%m%d')
    
    file_strs = np.array([str(tval) for tval in OMI_daily_data['day_values']])
    match_idx = np.where(date_str == file_strs)[0][0]
    
    plt.close('all')
    fig = plt.figure(figsize = (9, 5))
    ax1 = fig.add_subplot(2,3,1, projection = mapcrs)  # Map of original AI
    ax2 = fig.add_subplot(2,3,2, projection = mapcrs)  # Map of background clear AI
    ax3 = fig.add_subplot(2,3,3, projection = mapcrs)  # Map of AI above thresh
    ax4 = fig.add_subplot(2,3,4, projection = mapcrs)  # Map of forcing values
    ax5 = fig.add_subplot(2,3,5, projection = mapcrs)  # Map of cloud fraction
    ax6 = fig.add_subplot(2,3,6, projection = mapcrs)  # Map of ice concentration
    #ax6 = fig.add_subplot(2,3,2, projection = mapcrs)  # Map of GPQF values
    #ax3 = fig.add_subplot(2,3,6)  # histogram

    above_threshold = np.ma.masked_where(\
        OMI_daily_data['grid_AI'][match_idx,:,:] < ai_thresh, \
        OMI_daily_data['grid_AI'][match_idx,:,:])

    mesh = ax1.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], \
        OMI_daily_data['grid_AI'][match_idx,:,:], shading = 'auto', \
        transform = datacrs, cmap = 'jet', vmin = 0, vmax = 4.0)
    cbar = fig.colorbar(mesh, ax = ax1, label =  'UVAI')
    ax1.set_extent([-180,180,minlat,90], datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.coastlines()
    ax1.set_title('Daily OMI UVAI')

    mesh = ax2.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], \
        #OMI_daily_data['grid_GPQF'][match_idx,:,:], shading = 'auto', \
        clear_sky_AI, shading = 'auto', \
        transform = datacrs, cmap = 'jet', vmin = None, vmax = None)
    cbar = fig.colorbar(mesh, ax = ax2, label =  'UVAI')
    ax2.set_extent([-180,180,minlat,90], datacrs)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.coastlines()
    ax2.set_title('Daily OMI UVAI\nClear-sky Background')

    mesh = ax3.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], \
        #OMI_daily_data['grid_GPQF'][match_idx,:,:], shading = 'auto', \
        above_threshold, shading = 'auto', \
        transform = datacrs, cmap = 'jet', vmin = 0, vmax = 4.0)
    cbar = fig.colorbar(mesh, ax = ax3, label =  'UVAI')
    ax3.set_extent([-180,180,minlat,90], datacrs)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.coastlines()
    ax3.set_title('Daily OMI UVAI\nUVAI > ' + str(ai_thresh))

    mesh = ax4.pcolormesh(MYD08_data['lon'], \
        MYD08_data['lat'], MYD08_data['cld_frac_mean'][:,:], shading = 'auto',\
        transform = datacrs, \
        cmap = 'viridis')
    cbar = fig.colorbar(mesh, ax = ax4, label = 'Cloud Fraction')
    ax4.set_extent([-180,180,minlat,90], datacrs)
    ax4.set_boundary(circle, transform=ax4.transAxes)
    ax4.coastlines()
    ax4.set_title('Aqua MODIS\nCloud Fraction')

    mask_ice = np.ma.masked_where(NSIDC_data['grid_ice_conc'][:,:] == 0, \
        NSIDC_data['grid_ice_conc'][:,:])
    mesh = ax5.pcolormesh(NSIDC_data['grid_lon'], \
        NSIDC_data['grid_lat'], mask_ice, shading = 'auto',\
        transform = datacrs, \
        cmap = 'ocean', vmin = 1, vmax = 100)
    cbar = fig.colorbar(mesh, ax = ax5, label = 'Ice Concentration [%]')
    ax5.set_extent([-180,180,minlat,90], datacrs)
    ax5.set_boundary(circle, transform=ax5.transAxes)
    ax5.coastlines()
    ax5.set_title('SSMI/S Sea\nIce Concentration')

    min_force = np.nanmin(estimate_forcing[:,:])
    max_force = np.nanmax(estimate_forcing[:,:])
    if(abs(max_force) > abs(min_force)):
        lims = [-abs(max_force), abs(max_force)]
    else:
        lims = [-abs(min_force), abs(min_force)]
    mesh = ax6.pcolormesh(OMI_daily_data['lon_values'], \
        OMI_daily_data['lat_values'], estimate_forcing[:,:], shading = 'auto',\
        transform = datacrs, \
        cmap = 'bwr', vmin = lims[0], vmax = lims[1])
    cbar = fig.colorbar(mesh, ax = ax6, label = 'Aerosol Forcing [W/m2]')
    ax6.set_extent([-180,180,minlat,90], datacrs)
    ax6.set_boundary(circle, transform=ax6.transAxes)
    ax6.coastlines()
    ax6.set_title('Estimated\nAerosol Forcing')

    ###ax3.hist(np.ma.masked_where(estimate_forcing[:,:] == 0, \
    ###    estimate_forcing[:,:]).compressed(), bins = 'auto')
    ####ax6.set_yscale('log')
    ###ax3.set_xlabel('Estimated Forcing [W/m2]')
    ###ax3.set_ylabel('Counts')

    plt.suptitle(dt_date_str.strftime('%Y-%m-%d'))

    plot_subplot_label(ax1, 'a)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax2, 'b)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax3, 'c)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax4, 'd)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax5, 'e)', fontsize = 10, backgroundcolor = None)
    plot_subplot_label(ax6, 'f)', fontsize = 10, backgroundcolor = None)

    fig.tight_layout()
    if(save):
        #outname = 'test_calc_forcing_v3_' + date_str + '.png'
        dtype_add = ''
        if('data_type' in coloc_dict.keys()):
            if(coloc_dict['data_type'] == 'omi_uvai_pert'):
                dtype_add = '_pert'
        minlat_add = ''
        if('minlat' in coloc_dict.keys()):
            minlat_add = '_minlat' + str(int(coloc_dict['minlat']))
                
        outname = 'test_calc_forcing_v3_' + date_str + dtype_add + \
            minlat_add + '_v2.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

# This function grabs the OMI and CERES HDF5 files for given DTGs,
# grabs 10 (or other specified number) of OMI points, finds the closest
# CERES grid point to that OMI point, and figures out how many minutes
# are between the OMI and CERES observations at that point.
#
# file_list is a list of date-time groups from comp_data to analyze
def check_omi_ceres_time_offset(dir_list, num_points_per_file = 10):

    base_date = datetime(year=1970,month=1,day=1)
 
    # Loop over each directory
    # ------------------------
    for tdir in dir_list:
 
        print('\n' + tdir)

        day_str = tdir[:8]
 
        omi_dstr = tdir.strip().split('/')[-1]
      
        # Grab the OMI data file
        # ----------------------
        omi_data = h5py.File('comp_data/' + day_str + '/' + tdir + \
            '/omi_shawn_' + omi_dstr + '.hdf5')
    
        # Grab the CERES data file
        # ------------------------
        ceres_file = glob('comp_data/' + day_str + '/' + tdir + \
            '/ceres_*.hdf5')[0]
        ceres_data = h5py.File(ceres_file)

        # Flatten the OMI lat, lon, and time
        # ----------------------------------
        omi_lats_flat = omi_data['latitude'][:,:].flatten()
        omi_lons_flat = omi_data['longitude'][:,:].flatten()
        omi_time_flat = omi_data['time'][:,:].flatten()

        # Make sure the swaths only contain data north of 65
        # --------------------------------------------------
        check_minlat = np.min(omi_lats_flat)
        if(check_minlat < 60):
            print("ERROR: FULL SWATH")
        else:

            # Select "num_points_per_file" random indices from the OMI
            # file
            # --------------------------------------------------------
            test_idxs = (np.random.random(num_points_per_file) * \
                len(omi_time_flat)).astype('int')

            print('  omi_lat   omi_lon   cer_lat   cer_lon      Δkm     minutes')
            for tidx in test_idxs:

                omi_lat  = omi_lats_flat[tidx]
                omi_lon  = omi_lons_flat[tidx]
                omi_time = omi_time_flat[tidx]
    
                # Find the matching CERES grid points
                # -----------------------------------
                match_idx = nearest_gridpoint(omi_lat, omi_lon, \
                    ceres_data['latitude'][:,:], ceres_data['longitude'][:,:])

                ceres_lat = ceres_data['latitude'][:,:][match_idx][0]
                ceres_lon = ceres_data['longitude'][:,:][match_idx][0]
    
                # Calculate the OMI and CERES times based on base_date
                # ----------------------------------------------------
                omi_comp_time   = base_date + relativedelta(days = omi_time)
                ceres_comp_time = base_date + relativedelta(days = ceres_data['time'][:,:][match_idx][0])

                horiz_dist = find_distance_between_points(omi_lat, omi_lon, ceres_lat, ceres_lon)

                minutes_offset = (omi_comp_time - ceres_comp_time).seconds / 60

                print('{0:9.3f} {1:9.3f} {2:9.3f} {3:9.3f} {4:9.3f} {5:9.3f}'.format(\
                    omi_lat, omi_lon, ceres_lat, ceres_lon, horiz_dist, minutes_offset))

        omi_data.close()
        ceres_data.close()

# Plots the daily AI for a given date string as well as the corresponding
# monthly clear-sky ...
def plot_daily_OMI(daily_VSJ4, OMI_monthly_data, date_str, ai_thresh = 0.7):

    tidx = int(str(date_str)[4:6]) - 4

    clear_sky_AI = np.array([np.nanmean(\
        np.ma.masked_where(OMI_monthly_data['AI'][midx::6,:,:] > ai_thresh,\
        OMI_monthly_data['AI'][midx::6,:,:]), axis = 0) for midx in range(6)])
    clear_sky_AI = clear_sky_AI[tidx,:,:]

    match_idx = np.where(daily_VSJ4['day_values'] == int(date_str))[0][0]
    max_daily_AI = np.max(daily_VSJ4['grid_AI'][match_idx,:,:])
    plot_data = daily_VSJ4['grid_AI'][match_idx,:,:]

    delta_ai = plot_data - clear_sky_AI
   
    delta_ai = np.ma.masked_where(delta_ai < ai_thresh, delta_ai)
 
    plt.close('all')
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1, projection = ccrs.NorthPolarStereo())
    ax2 = fig.add_subplot(2,2,2, projection = ccrs.NorthPolarStereo())
    ax3 = fig.add_subplot(2,2,3, projection = ccrs.NorthPolarStereo())
    ax4 = fig.add_subplot(2,2,4)
    mesh = ax1.pcolormesh(daily_VSJ4['lon_values'][:], \
        daily_VSJ4['lat_values'][:], plot_data, \
        transform = ccrs.PlateCarree(), \
        shading = 'auto', cmap = 'jet', vmin = -2.0, vmax = 3)
    cbar = fig.colorbar(mesh, ax = ax1)
    ax1.coastlines()
    ax1.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax1.set_title(str(date_str) + '\nMax AI = ' + str(np.round(max_daily_AI, 1)))

    mesh = ax2.pcolormesh(OMI_monthly_data['LON'][:,:], \
        OMI_monthly_data['LAT'][:,:], clear_sky_AI[:,:], \
        transform = ccrs.PlateCarree(), \
        shading = 'auto', cmap = 'jet', vmin = -0.1, vmax = 0.1)
    cbar = fig.colorbar(mesh, ax = ax2)
    ax2.coastlines()
    ax2.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax2.set_title('Monthly clear-sky climo')

    mesh = ax3.pcolormesh(OMI_monthly_data['LON'][:,:], \
        OMI_monthly_data['LAT'][:,:], delta_ai, \
        transform = ccrs.PlateCarree(), \
        shading = 'auto', cmap = 'jet', vmin = 0, vmax = None)
    cbar = fig.colorbar(mesh, ax = ax3)
    ax3.coastlines()
    ax3.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax3.set_title('ΔAI (daily - clear_sky_bckd)')

    plotter = delta_ai.compressed()
    ax4.hist(plotter, bins = 20)
    ax4.set_xlabel('Delta AI')
    ax4.set_title('# valid points: ' + str(int(plotter.shape[0])))
    
    fig.tight_layout()
    plt.show()

def compare_nn_version_output(date_str, skip_version = None, save = False):

    # Grab all the desired simulation files
    # -------------------------------------
    infiles = glob('neuralnet_output/*' + date_str + '.hdf5')

    if(skip_version is not None):
        switch_array = np.full(len(infiles), True)


        for ii, tfile in enumerate(infiles):
            this_version = tfile.split('/')[-1].split('_')[-2] 
            print(this_version, this_version in skip_version)
            if(this_version in skip_version):
                switch_array[ii] = False

        infiles = np.array(infiles)[switch_array]

    # Figure out the size of the arrays needed
    # ----------------------------------------
    data = h5py.File(infiles[0])
    sizer = data['ceres_swf'].shape
    data.close()
    
    nn_force = np.full((len(infiles), sizer[0], sizer[1]), np.nan)
    nn_swf   = np.full((len(infiles), sizer[0], sizer[1]), np.nan)
    
    for ii, infile in enumerate(infiles):
        data = h5py.File(infile)
    
        calc_swf = data['calc_swf'][:,:]
        ceres_swf = data['ceres_swf'][:,:]
        ceres_swf = np.ma.masked_where( (ceres_swf == -999.) | (ceres_swf > 1000), ceres_swf)
    
        nn_swf[ii,:,:] = calc_swf
        nn_force[ii,:,:] = calc_swf - ceres_swf
        obs_swf = ceres_swf
        lats    = data['omi_lat'][:,:]
        lons    = data['omi_lon'][:,:]
        data.close()
        #plot_compare_NN_output_v2(infile, auto_zoom = False, save = False)
    
    nn_swf = np.ma.masked_invalid(nn_swf)
    nn_force = np.ma.masked_invalid(nn_force)
    
    obs_swf = np.ma.masked_where((obs_swf == -999.) | (obs_swf > 1000), obs_swf)
    
    mean_nn_swf = np.mean(nn_swf, axis = 0)
    std_nn_swf  = np.std(nn_swf, axis = 0)
    mean_nn_force = np.mean(nn_force, axis = 0)
    std_nn_force  = np.std(nn_force, axis = 0)
    
    fig = plt.figure(figsize = (9, 7))
    ax1 = fig.add_subplot(2,2,1, projection = ccrs.NorthPolarStereo())
    ax2 = fig.add_subplot(2,2,2, projection = ccrs.NorthPolarStereo())
    ax3 = fig.add_subplot(2,2,3, projection = ccrs.NorthPolarStereo())
    ax4 = fig.add_subplot(2,2,4)
    mesh = ax1.pcolormesh(lons, lats, obs_swf, transform = ccrs.PlateCarree(), \
        shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax1, label = "SWF [W/m2]")
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.coastlines()
    ax1.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax1.set_title('CERES Obs.')
    
    mesh = ax2.pcolormesh(lons, lats, mean_nn_swf, transform = ccrs.PlateCarree(), \
        shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax2, label = "SWF [W/m2]")
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ax2.coastlines()
    ax2.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax2.set_title('Average of NN versions')
    
    mesh = ax3.pcolormesh(lons, lats, std_nn_swf, transform = ccrs.PlateCarree(), \
        shading = 'auto', vmax = 20)
    cbar = fig.colorbar(mesh, ax = ax3, label = "SWF [W/m2]")
    ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.coastlines()
    ax3.set_extent([-180,180,65,90], ccrs.PlateCarree())
    ax3.set_title('Standard Dev. of NN versions')
    
    ax4.hist(std_nn_swf.compressed(), bins = 50)
    ax4.set_xlabel('STD of NN SWF')
    ax4.set_ylabel('Counts')
   
    plt.suptitle(date_str)
 
    fig.tight_layout()
    if(save):
        outname = 'nn_version_compare_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def compare_sza_bin_impact_on_slopes(test_dict, bin_dict, sfc_idx, \
        ai_min = 0, min_ob = 50, maxerr = 2, trend_type = 'linregress', \
        combined_plot = False, save = False):

    # Set up different ranges of SZA bins to test
    # -------------------------------------------
    delta_sza = np.array([5, 10, 15, 20])

    fig = plt.figure(figsize = (9, 7))
    if(combined_plot):
        axs = fig.subplots(nrows = 1, ncols = 1)
        flat_ax = axs
    else:
        axs = fig.subplots(nrows = 2, ncols = 2, sharex = True, sharey = True)
        flat_ax = axs.flatten()

    xvals = np.arange(bin_dict['cod_bin_means'].shape[0])
    cod_labels = [str(edge) for edge in bin_dict['cod_bin_edges']]
  
    linestyles = ['-','--',':','-.'] 
    
    min_force = 999
    max_force = -999
     
    for ii in range(len(delta_sza)):
        sza_bin_edges = np.arange(bin_dict['sza_bin_edges'][0], \
            bin_dict['sza_bin_edges'][-1] + delta_sza[ii], delta_sza[ii])
        sza_bin_means = (sza_bin_edges[1:] + sza_bin_edges[:-1]) / 2
        print(sza_bin_edges)

        # Calculate the regression slopes and intercepts for these SZA bins
        # -----------------------------------------------------------------
        slope_dict_lin = calc_NN_force_slope_intcpt(test_dict, bin_dict['ice_bin_edges'], \
                sza_bin_edges, bin_dict['cod_bin_edges'], ai_min = ai_min, min_ob = min_ob, \
                trend_type = 'linregress')

        plot_slopes = np.ma.masked_where(\
            slope_dict_lin['slope_stderr'] > maxerr, \
            slope_dict_lin['slopes'])

        if(np.nanmax(plot_slopes) > max_force):
            max_force = np.nanmax(plot_slopes)
        if(np.nanmin(plot_slopes) < min_force):
            min_force = np.nanmin(plot_slopes)

        # Plot these on the corresponding plot
        # ------------------------------------
        for jj in range(slope_dict_lin['slopes'].shape[1]): 
            if(combined_plot):
                flat_ax.plot(xvals, plot_slopes[sfc_idx,jj,:], \
                    linestyle = linestyles[ii], \
                    label = str(sza_bin_means[jj]))
            else:
                flat_ax[ii].plot(xvals, plot_slopes[sfc_idx,jj,:],\
                    linestyle = linestyles[ii], \
                    label = str(sza_bin_means[jj]))
            #flat_ax.plot(xvals, slope_dict_lin['slopes'][sfc_idx,jj,:],linestyle = linestyles[ii])

        if(combined_plot):
            flat_ax.axhline(0, linestyle = ':', color = 'gray')
            flat_ax.set_title('ΔSZA')
            flat_ax.set_xticks([0,1,2,3,4,5,6,7,8])
            flat_ax.set_xticklabels(cod_labels[::1])
            flat_ax.set_xlabel('COD')
            flat_ax.set_ylabel('Force. Eff. Slope [W m-2 AI-1]')
            flat_ax.legend()
        else:
            flat_ax[ii].axhline(0, linestyle = ':', color = 'gray')
            flat_ax[ii].set_title('ΔSZA = ' + str(delta_sza[ii]))
            flat_ax[ii].set_xticks([0,1,2,3,4,5,6,7,8])
            flat_ax[ii].set_xticklabels(cod_labels[::1])
            flat_ax[ii].grid(alpha = 0.50, linestyle = ':')
            flat_ax[ii].legend()

    if(not combined_plot):
        axs[0,0].set_ylabel('Force. Eff. Slope [W m-2 AI-1]')
        axs[1,0].set_ylabel('Force. Eff. Slope [W m-2 AI-1]')
        axs[1,0].set_xlabel('COD')
        axs[1,1].set_xlabel('COD')

    title_options = ['Ocean','Mix','Ice','Land']
    plt.suptitle(title_options[sfc_idx])
   
    fig.tight_layout() 
    if(save):
        outname = 'slope_sensitivity_to_sza_' + title_options[sfc_idx] + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()        

def plot_NN_error_dist_bulk(sim_name, num_bins = 100, xmin = None, \
        xmax = None, ax = None, astrofit = False, \
        use_correct_error_calc = False, save = False):

    # Read in the desired files
    # -------------------------
    #if( (sim_name == 'noland103') | (sim_name == 'noland104') | \
    #    (sim_name == 'noland105') | (sim_name == 'noland106') | |
    #    (sim_name == 'noland107') ):
    if( int(sim_name[6:]) >= 103) :
        files = glob('neuralnet_output_clear_newfiles/test_calc_out_' + sim_name + '*.hdf5')
    else:
        files = glob('neuralnet_output_clear/test_calc_out_' + sim_name + '*.hdf5')

    if(len(files) == 0):
        print("ERROR: NO CLEAR FILES FOUND FOR SIM " + sim_name)
        return

    # Figure out the total size to insert the data
    # ---------------------------------------------
    #minlat = 70.
    minlat = 65.    # 2024/01/10: changed to 65.
    total_size = 0
    for ff in files:
        data = h5py.File(ff,'r')
        #local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
        local_data = np.ma.masked_invalid(data['omi_uvai_pert'])
        local_data = np.ma.masked_where((local_data < -12) | (local_data > 1.0), local_data)
        local_data = np.ma.masked_where(data['omi_lat'][:,:] < minlat, \
            local_data) 
        local_data = np.ma.masked_where((data['ceres_swf'][:,:] < -200.) | \
            (data['ceres_swf'][:,:] > 3000), \
            local_data) 
        local_data = np.ma.masked_where(np.isnan(data['calc_swf'][:,:]), local_data)
        local_size = local_data.compressed().shape[0]
        max_ai = np.max(local_data)
        print(ff, local_size, np.round(max_ai, 3))
        total_size += local_size
    
        data.close()
    
    
    # Set up the data structure to hold all the data
    # ----------------------------------------------
    combined_data = {}
    combined_data['calc_swf']      = np.full(total_size, np.nan)
    combined_data['ceres_swf']     = np.full(total_size, np.nan)
 
    print("Loading data")
    
    # Loop back over the files and insert the data into the structure
    # ---------------------------------------------------------------
    total_size = 0
    beg_idx = 0
    end_idx = 0
    for ff in files:
    
        data = h5py.File(ff,'r')
        # NOTE: Changed the omi variable here from "pert" to "raw" on 20230623.
        #       This move should allow for coloc data to be read after 2020
        #local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
        local_data = np.ma.masked_invalid(data['omi_uvai_pert'])
        print("MAX AI FOR SWATH", ff, np.max(local_data))
        local_data = np.ma.masked_where((local_data < -12) | (local_data > 1.0), local_data)
        local_data = np.ma.masked_where(data['omi_lat'][:,:] < minlat, \
            local_data) 
        local_data = np.ma.masked_where((data['ceres_swf'][:,:] < -200.) | \
            (data['ceres_swf'][:,:] > 3000), \
            local_data) 
        local_data = np.ma.masked_where(np.isnan(data['calc_swf'][:,:]), local_data)
        local_size = local_data.compressed().shape[0]
    
        beg_idx = end_idx
        end_idx = beg_idx + local_size
    
        for tkey in combined_data.keys():
            combined_data[tkey][beg_idx:end_idx] = \
                data[tkey][~local_data.mask]
    
        #print(local_size)
        total_size += local_size
    
        data.close()

    if(use_correct_error_calc):
        print("Calculating NN errors using NN - CERES")
        errors = combined_data['calc_swf'] - combined_data['ceres_swf']
    else:
        print("Calculating NN errors using CERES - NN")
        errors = combined_data['ceres_swf'] - combined_data['calc_swf']

    mean_err = np.mean(errors)
    std_err  = np.std(errors)
    median_err = np.median(errors)

    print("Normal mean error:       ", np.round(mean_err, 1))
    print("Normal std dev error:    ", np.round(std_err, 1))
    print("Normal median dev error: ", np.round(median_err, 1))
   
    in_ax = True
    if(ax is None):
        in_ax = False 
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    if(astrofit):
        ax.hist(errors, bins = num_bins)
        bin_heights,bin_borders = np.histogram(errors,bins=num_bins)
        bin_widths = np.diff(bin_borders)
        bin_centers = bin_borders[:-1] + bin_widths / 2
        t_init = models.Gaussian1D()
        fit_t = fitting.LevMarLSQFitter()
        t = fit_t(t_init, bin_centers, bin_heights)
        print('Amplitude: ',np.round(t.amplitude.value,3))
        print('Mean:      ',np.round(t.mean.value,3))
        print('StDev:     ',np.round(t.stddev.value,3))
        x_interval_for_fit = np.linspace(bin_centers[0],bin_centers[-1],400)
        ax.plot(x_interval_for_fit,t(x_interval_for_fit),label='fit',c='tab:red')
        ax.set_xlim(xmin, xmax)

        # Add statistics to plot
        # ----------------------
        xdiff = ax.get_xlim()[1] - ax.get_xlim()[0]
        plot_xval = ax.get_xlim()[0]  + xdiff * 0.65
        plot_yval = ax.get_ylim()[1] * 0.8
        ptext = 'μ = ' + str(np.round(t.mean.value, 1)) + ' Wm$^{-2}$\n' + \
            'σ = ' + str(np.round(t.stddev.value, 1)) + ' Wm$^{-2}$'


        #plot_figure_text(ax3, ptext, location = 'lower_right', \
        #ax5.set_title(ptext)
        plot_figure_text(ax, ptext, xval = plot_xval, yval = plot_yval, \
            fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')

    else:
        ax.hist(errors, bins = num_bins, density = True)
        mu, std = statnorm.fit(errors)
        plot_xmin, plot_xmax = ax.get_xlim()
        xvals = np.linspace(plot_xmin, plot_xmax, 100)
        p = statnorm.pdf(xvals, mu, std)
        ax.plot(xvals, p, color = 'red')
        ax.set_xlim(xmin, xmax)

    ax.axvline(0, linestyle = '--', color = 'k', alpha = 0.5)

    if(use_correct_error_calc):
        ax.set_xlabel('TOA SWF Error (NN - CERES) [Wm$^{-2}$]')
    else:
        ax.set_xlabel('TOA SWF Error (CERES - NN) [Wm$^{-2}$]')
    #ax.set_xlabel('Forcing Error (NN - CERES) [Wm$^{-2}$]')
    ax.set_ylabel('Counts')
    ax.set_title('Errors between L2 CERES & NN output\n' + \
        'under aerosol-free conditions')

    if(not in_ax):
        fig.tight_layout()
        if(save):
            if(use_correct_error_calc):
                err_calc_add = '_correct'
            else:
                err_calc_add = ''
            outname = 'errors_aerfree_dist_' + sim_name + err_calc_add + \
                '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show()

def calc_monthly_force_from_daily(daily_dict, minlat = 65.5, maxlat = 90.5, \
        return_std = False):

    # Figure out the number of months 
    # -------------------------------
    months = np.array([int(str(tdate)[:6]) for tdate in daily_dict['dates']])
    unique_months = np.unique(months)
   
    test_lats = daily_dict['latitude'] 
    lat_idxs = np.where( (test_lats >= minlat) & (test_lats <= maxlat) )[0]
    min_lat_idx = lat_idxs[0]
    max_lat_idx = lat_idxs[-1] + 1
    test_lats = test_lats[lat_idxs]
    print("Working with latitudes",test_lats)

    full_arr = np.full( (unique_months.shape[0], \
        test_lats.shape[0], \
        #daily_dict['force_estimate_orig'].shape[1], \
        daily_dict['force_estimate_orig'].shape[2]), np.nan)
    if(return_std):
        full_std = np.full( (unique_months.shape[0], \
            test_lats.shape[0], \
            #daily_dict['force_estimate_orig'].shape[1], \
            daily_dict['force_estimate_orig'].shape[2]), np.nan)
    
    for ii in range(unique_months.shape[0]):
        print(unique_months[ii])
        time_idxs = np.where(months == unique_months[ii])[0]
    
        full_arr[ii,:,:] = np.nanmean(\
            daily_dict['force_estimate_orig'][time_idxs, \
            min_lat_idx:max_lat_idx,:], axis = (0))
        if(return_std):
            full_std[ii,:,:] = np.nanstd(\
                daily_dict['force_estimate_orig'][time_idxs, \
                min_lat_idx:max_lat_idx,:], axis = (0))
   
    if(return_std): 
        return full_arr, full_std       
    else:
        return full_arr


def calc_force_vals_bulk(daily_dict, err_mean, err_std, minlat = 65.5, \
        maxlat = 90.5, num_sims = 100, return_std = False):

    # Figure out the number of months 
    # -------------------------------
    months = np.array([int(str(tdate)[:6]) for tdate in daily_dict['dates']])
    unique_months = np.unique(months)
   
    test_lats = daily_dict['latitude'] 
    lat_idxs = np.where( (test_lats >= minlat) & (test_lats <= maxlat) )[0]
    min_lat_idx = lat_idxs[0]
    max_lat_idx = lat_idxs[-1] + 1
    test_lats = test_lats[lat_idxs]
    print("Working with latitudes",test_lats)

    full_arr = np.full( (num_sims, unique_months.shape[0], \
        test_lats.shape[0], \
        #daily_dict['force_estimate_orig'].shape[1], \
        daily_dict['force_estimate_orig'].shape[2]), np.nan)
    if(return_std):
        full_std = np.full( (num_sims, unique_months.shape[0], \
            test_lats.shape[0], \
            #daily_dict['force_estimate_orig'].shape[1], \
            daily_dict['force_estimate_orig'].shape[2]), np.nan)
    
    for ii in range(unique_months.shape[0]):
        print(unique_months[ii])
        time_idxs = np.where(months == unique_months[ii])[0]
       
        work_arr = np.full( (num_sims, time_idxs.shape[0], \
            test_lats.shape[0], \
            #daily_dict['force_estimate_orig'].shape[1], \
            daily_dict['force_estimate_orig'].shape[2]), np.nan)
    
        for jj in range(num_sims):
            work_arr[jj,:,:,:] = \
                daily_dict['force_estimate_orig'][time_idxs,min_lat_idx:max_lat_idx,:]
             
            work_arr[jj,:,:,:][work_arr[jj,:,:,:] != 0] += \
                np.random.normal(err_mean, err_std, work_arr[jj,:,:,:][work_arr[jj,:,:,:] != 0].shape)
    
        full_arr[:,ii,:,:] = np.nanmean(work_arr, axis = (1))
        if(return_std):
            full_std[:,ii,:,:] = np.nanstd(work_arr, axis = (1))
   
    if(return_std): 
        return full_arr, full_std       
    else:
        return full_arr

# pvar: 'mean', 'stdev', 'min', 'max'
def plot_bulk_sim_trends(daily_dict, forcing_trends, pvar, vmax = 1.0, \
        conf_window = 90, save = False):

    if(pvar == 'mean'):
        plot_trends = np.mean(forcing_trends, axis = 0)
        vmin = -vmax
        pcolor = 'bwr'
    elif(pvar == 'stdev'):
        plot_trends = np.std(forcing_trends, axis = 0)
        vmin = 0
        pcolor = 'jet'
    elif(pvar == 'min'):
        plot_trends = np.min(forcing_trends, axis = 0)
        vmin = -vmax
        pcolor = 'bwr'
    elif(pvar == 'max'):
        plot_trends = np.max(forcing_trends, axis = 0)
        vmin = -vmax
        pcolor = 'bwr'

    mean_trends = np.nanmean(forcing_trends, axis = 0)
    stdv_trends = np.std(forcing_trends, axis = 0)
    stderr_trends = sem(forcing_trends, axis = 0), 
    #stderr_trend = stdv_trend / np.sqrt(in_arr.shape[0])
    conf_intvl = statnorm.interval(alpha = 0.99, \
        loc = mean_trends, scale = stderr_trends)

    ## Calculate the standard error of the mean of the "n" monthly trends
    ## ------------------------------------------------------------------
    #num_sims = forcing_trends.shape[0]
    #stderr_trends = stdv_trends / np.sqrt(num_sims)

    if(conf_window == 90):
        conf_val = 1.645
    elif(conf_window == 95):
        conf_val = 1.960
    else:
        print("CURRENTLY INVALID CONFIDENCE INTERVAL. USING 90")
        conf_val = 1.645

    #window_min = plot_trends[month_idx,:,:] - (stdv_trends[month_idx,:,:] * conf_val)
    #window_max = plot_trends[month_idx,:,:] + (stdv_trends[month_idx,:,:] * conf_val)

    plt.close('all')
    fig = plt.figure(figsize = (9, 8.5))
    ax1  = fig.add_subplot(2,3,1, projection = ccrs.NorthPolarStereo())
    ax2  = fig.add_subplot(2,3,2, projection = ccrs.NorthPolarStereo())
    ax3  = fig.add_subplot(2,3,3, projection = ccrs.NorthPolarStereo())
    ax4  = fig.add_subplot(2,3,4, projection = ccrs.NorthPolarStereo())
    ax5  = fig.add_subplot(2,3,5, projection = ccrs.NorthPolarStereo())
    ax6  = fig.add_subplot(2,3,6, projection = ccrs.NorthPolarStereo())
    flat_axs = [ax1, ax2, ax3, ax4, ax5, ax6]

    titlers = ['Apr','May','June','July','Aug','Sep']
    for ii in range(6):

        #hash_idxs = np.where( np.abs(mean_trends[ii,:,:]) > stdv_trends[ii,:,:])

        mesh = flat_axs[ii].pcolormesh(daily_dict['longitude'][:], daily_dict['latitude'][:], plot_trends[ii,:,:], \
            cmap = pcolor, vmin = vmin, vmax = vmax, transform = ccrs.PlateCarree(), shading = 'auto')
        #cbar = fig.colorbar(mesh, ax = flat_axs[ii], pad = 0.03, fraction = 0.045)
        flat_axs[ii].set_boundary(circle, transform=flat_axs[ii].transAxes)
        flat_axs[ii].coastlines()
        flat_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        flat_axs[ii].set_title(str(ii))

        if(pvar == 'mean'):
            #window_min = plot_trends[ii,:,:] - (stdv_trends[ii,:,:] * 1.0)
            #window_max = plot_trends[ii,:,:] + (stdv_trends[ii,:,:] * 1.0)

            #hasher = np.ma.masked_where(  \
            #    (plot_trends[ii,:,:] == 0) |
            #    ((plot_trends[ii,:,:] < 0) & (window_max[:,:] > 0)) | \
            #    ((plot_trends[ii,:,:] > 0) & (window_min[:,:] < 0)), plot_trends[ii,:,:])

            hasher = np.ma.masked_where(  \
                 (plot_trends[ii,:,:] == 0) |
                 (np.abs(plot_trends[ii,:,:]) < \
                  (stdv_trends[ii,:,:] * conf_val)), plot_trends[ii, :,:])


            #hasher = np.ma.masked_where( (np.abs(plot_trends[ii,:,:]) <= 1.0 * stdv_trends[ii,:,:]), plot_trends[ii,:,:])
            flat_axs[ii].pcolor(daily_dict['longitude'][:], daily_dict['latitude'][:], \
                hasher, hatch = '....', alpha=0., shading = 'auto', transform = ccrs.PlateCarree())

            #hasher = np.ma.masked_where( (np.abs(plot_trends[ii,:,:]) <= 1.0 * stdv_trends[ii,:,:]), plot_trends[ii,:,:])
            #flat_axs[ii].pcolor(daily_dict['longitude'][:], daily_dict['latitude'][:], \
            #    hasher, hatch = '....', alpha=0., shading = 'auto', transform = ccrs.PlateCarree())


        #hasher = np.ma.masked_where(out_dict['raw_stderr'][0,:,:] > maxerr, \
        #    out_dict['raw_stderr'][0,:,:])
        #ax1.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        #    hasher, hatch = hatch_style, alpha=0., shading = 'auto')
        flat_axs[ii].set_title(titlers[ii])

    if(pvar == 'mean'):
        plt.suptitle('Mean of ' + str(int(forcing_trends.shape[0])) + \
            ' Observation-based Aerosol Direct Forcing Trends\nHashing = 90% Confidence in the trend being nonzero')
    elif(pvar == 'stdev'):
        plt.suptitle('Standard Deviation of ' + str(int(forcing_trends.shape[0])) + \
            ' Forcing Trends\nHashing = Trend mean is greater than standard deviation')


    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax3, 'c)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax4, 'd)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax5, 'e)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax6, 'f)', fontsize = 11, backgroundcolor = 'white')

    cbar_ax = fig.add_axes([0.17, 0.14, 0.70, 0.01])
    if( (pvar == 'mean') | (pvar == 'min') | (pvar == 'max') ):
        cbar = fig.colorbar(mesh, \
            cax = cbar_ax, shrink = 0.8, orientation = 'horizontal', \
            extend = 'both', \
            label = 'Forcing Trend [Wm$^{-2}$(study period)$^{-1}$]')
    elif(pvar == 'stdev'):
        cbar = fig.colorbar(mesh, \
            cax = cbar_ax, shrink = 0.8, orientation = 'horizontal', \
            extend = 'max', \
            label = 'Forcing Trend Std. Dev. [Wm$^{-2}$(study period)$^{-1}$]')
    
    fig.tight_layout(rect = [0,0.1,1,1])
    if(save):
        outname = 'bulk_trend_' + pvar + '_numsims' + \
            str(int(forcing_trends.shape[0])) + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_grid_OMI_trends_spatial(OMI_daily_data, min_AI = None, \
        max_AI = None, minlat = 65.5, maxlat = 90.5,  save = False):
    
    if(isinstance(OMI_daily_data, str)):
        if(min_AI is None):
            print("ERROR: If providing an OMI daily file, must provide min_AI")
            return

        print("Reading daily OMI data from", OMI_daily_data)
        print("Converting daily data to monthly omi_uvai_pert averages")
        OMI_daily_data = calcOMI_MonthAvg_FromDaily(OMI_daily_data, \
            min_AI = min_AI, max_AI = max_AI, minlat = minlat, maxlat = maxlat)

    # Calculate gridded trends
    #grid_trends = np.full( (6, OMI_daily_data['AI'].shape[1],  \
    #    OMI_daily_data['AI'].shape[2]), np.nan)
   
    plt.close('all') 
    fig = plt.figure(figsize = (9, 8.5))
    ax1 = fig.add_subplot(2,3,1, projection = ccrs.NorthPolarStereo())
    ax2 = fig.add_subplot(2,3,2, projection = ccrs.NorthPolarStereo())
    ax3 = fig.add_subplot(2,3,3, projection = ccrs.NorthPolarStereo())
    ax4 = fig.add_subplot(2,3,4, projection = ccrs.NorthPolarStereo())
    ax5 = fig.add_subplot(2,3,5, projection = ccrs.NorthPolarStereo())
    ax6 = fig.add_subplot(2,3,6, projection = ccrs.NorthPolarStereo())
    axs = [ax1, ax2, ax3, ax4, ax5, ax6]
    titlers = ['Apr','May','June','July','Aug','Sep']
    
    for ii in range(6):
        print(ii)
        grid_trends = np.full( (OMI_daily_data['AI'].shape[1],  \
            OMI_daily_data['AI'].shape[2]), np.nan)
        for jj in range(OMI_daily_data['AI'].shape[1]):
            for kk in range(OMI_daily_data['AI'].shape[2]):
                # Plot the individual runs
                # ------------------------
                local_sim_vals = OMI_daily_data['AI'][ii::6,jj,kk]
                x_vals = np.arange(local_sim_vals.shape[0])
    
                #if(trend_type=='standard'): 
                result = stats.linregress(x_vals, local_sim_vals[:])
                grid_trends[jj,kk] = result.slope * len(x_vals)
                #grid_trends[ii,jj,kk] = result.slope * len(x_vals)
    
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,work_mask.compressed())
                    #forcing_trends[i,j] = result.slope * len(x_vals)
                    #forcing_pvals[i,j]  = result.pvalue
                    #forcing_uncert[i,j] = result.stderr * len(x_vals)
                #else:
                #    res = stats.theilslopes(local_sim_vals[jj,:], x_vals, 0.90)
                #    trend_results[jj] = res[0] * len(x_vals)
                #    #forcing_trends[i,j] = res[0]*len(x_vals)
    
    
    #for ii in range(6):
        mesh = axs[ii].pcolormesh(OMI_daily_data['LON'], OMI_daily_data['LAT'], \
            grid_trends[:,:], transform = ccrs.PlateCarree(), \
            #grid_trends[ii,:,:], transform = ccrs.PlateCarree(), \
            shading = 'auto', cmap = 'bwr', vmin = -0.5, vmax = 0.5)
        axs[ii].set_boundary(circle, transform=axs[ii].transAxes)
        axs[ii].coastlines()
        axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        axs[ii].set_title(titlers[ii])

    plot_subplot_label(ax1, 'a)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax2, 'b)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax3, 'c)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax4, 'd)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax5, 'e)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax6, 'f)', fontsize = 12, backgroundcolor = 'white')

    plt.suptitle('OMI UVAI Trends\n2005 - 2020')
    cbar_ax = fig.add_axes([0.17, 0.13, 0.70, 0.01])
    cbar = fig.colorbar(mesh, \
        cax = cbar_ax, shrink = 0.8, orientation = 'horizontal', \
        extend = 'both', \
        label = 'Trend of Perturbed OMI UVAI [AI(study period)$^{-1}$]')
    
    fig.tight_layout(rect = [0,0.1,1,1])
   
    if(save):
        if(min_AI < 0):
            outname = 'omi_grid_trend_minAI_N' + str(int(abs(min_AI) * 100)).zfill(3) + '.png'
        else:
            outname = 'omi_grid_trend_minAI' + str(int(abs(min_AI) * 100)).zfill(3) + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show() 

def plot_bulk_force_AI_trend(daily_dict, forcing_trends, OMI_daily_data, \
        vmax = 1.0, min_AI = None, max_AI = None, minlat = 65.5, \
        maxlat = 90.5,  save = False):

    if(isinstance(OMI_daily_data, str)):
        if(min_AI is None):
            print("ERROR: If providing an OMI daily file, must provide min_AI")
            return

        print("Reading daily OMI data from", OMI_daily_data)
        print("Converting daily data to monthly omi_uvai_pert averages")
        OMI_daily_data = calcOMI_MonthAvg_FromDaily(OMI_daily_data, \
            min_AI = min_AI, max_AI = max_AI, minlat = minlat, maxlat = maxlat)


    mean_trends = np.nanmean(forcing_trends, axis = 0)
    stdv_trends = np.std(forcing_trends, axis = 0)
    stderr_trends = sem(forcing_trends, axis = 0), 
    #stderr_trend = stdv_trend / np.sqrt(in_arr.shape[0])
    conf_intvl = statnorm.interval(alpha = 0.99, \
        loc = mean_trends, scale = stderr_trends)

    plt.close('all')
    fig = plt.figure(figsize = (9, 12))
    #ax1  = fig.add_subplot(6,3,1,  projection = ccrs.NorthPolarStereo())
    #ax2  = fig.add_subplot(6,3,4,  projection = ccrs.NorthPolarStereo())
    #ax3  = fig.add_subplot(6,3,7,  projection = ccrs.NorthPolarStereo())
    #ax4  = fig.add_subplot(6,3,10, projection = ccrs.NorthPolarStereo())
    #ax5  = fig.add_subplot(6,3,13, projection = ccrs.NorthPolarStereo())
    #ax6  = fig.add_subplot(6,3,16, projection = ccrs.NorthPolarStereo())
    #ax7  = fig.add_subplot(6,3,2,  projection = ccrs.NorthPolarStereo())
    #ax8  = fig.add_subplot(6,3,5,  projection = ccrs.NorthPolarStereo())
    #ax9  = fig.add_subplot(6,3,8, projection = ccrs.NorthPolarStereo())
    #ax10 = fig.add_subplot(6,3,11, projection = ccrs.NorthPolarStereo())
    #ax11 = fig.add_subplot(6,3,14, projection = ccrs.NorthPolarStereo())
    #ax12 = fig.add_subplot(6,3,17, projection = ccrs.NorthPolarStereo())
    #ax13 = fig.add_subplot(6,3,3,  projection = ccrs.NorthPolarStereo())
    #ax14 = fig.add_subplot(6,3,6,  projection = ccrs.NorthPolarStereo())
    #ax15 = fig.add_subplot(6,3,9, projection = ccrs.NorthPolarStereo())
    #ax16 = fig.add_subplot(6,3,12, projection = ccrs.NorthPolarStereo())
    #ax17 = fig.add_subplot(6,3,15, projection = ccrs.NorthPolarStereo())
    #ax18 = fig.add_subplot(6,3,18, projection = ccrs.NorthPolarStereo())
    #omi_axs  = [ax1, ax2, ax3, ax4, ax5, ax6]
    #mean_axs = [ax7, ax8, ax9, ax10, ax11, ax12]
    #stdv_axs = [ax13, ax14, ax15, ax16, ax17, ax18]

    ax1  = fig.add_subplot(6,4,1,  projection = ccrs.NorthPolarStereo())
    ax2  = fig.add_subplot(6,4,5,  projection = ccrs.NorthPolarStereo())
    ax3  = fig.add_subplot(6,4,9,  projection = ccrs.NorthPolarStereo())
    ax4  = fig.add_subplot(6,4,13, projection = ccrs.NorthPolarStereo())
    ax5  = fig.add_subplot(6,4,17, projection = ccrs.NorthPolarStereo())
    ax6  = fig.add_subplot(6,4,21, projection = ccrs.NorthPolarStereo())
    ax7  = fig.add_subplot(6,4,2,  projection = ccrs.NorthPolarStereo())
    ax8  = fig.add_subplot(6,4,6,  projection = ccrs.NorthPolarStereo())
    ax9  = fig.add_subplot(6,4,10, projection = ccrs.NorthPolarStereo())
    ax10 = fig.add_subplot(6,4,14, projection = ccrs.NorthPolarStereo())
    ax11 = fig.add_subplot(6,4,18, projection = ccrs.NorthPolarStereo())
    ax12 = fig.add_subplot(6,4,22, projection = ccrs.NorthPolarStereo())
    ax13 = fig.add_subplot(6,4,3,  projection = ccrs.NorthPolarStereo())
    ax14 = fig.add_subplot(6,4,7,  projection = ccrs.NorthPolarStereo())
    ax15 = fig.add_subplot(6,4,11, projection = ccrs.NorthPolarStereo())
    ax16 = fig.add_subplot(6,4,15, projection = ccrs.NorthPolarStereo())
    ax17 = fig.add_subplot(6,4,19, projection = ccrs.NorthPolarStereo())
    ax18 = fig.add_subplot(6,4,23, projection = ccrs.NorthPolarStereo())
    ax19 = fig.add_subplot(6,4,4)
    ax20 = fig.add_subplot(6,4,8)
    ax21 = fig.add_subplot(6,4,12)
    ax22 = fig.add_subplot(6,4,16)
    ax23 = fig.add_subplot(6,4,20)
    ax24 = fig.add_subplot(6,4,24)
    omi_axs  = [ax1, ax2, ax3, ax4, ax5, ax6]
    mean_axs = [ax7, ax8, ax9, ax10, ax11, ax12]
    stdv_axs = [ax13, ax14, ax15, ax16, ax17, ax18]
    hist_axs = [ax19, ax20, ax21, ax22, ax23, ax24]

    titlers = ['Apr','May','June','July','Aug','Sep']
    for ii in range(6):

        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the OMI trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        grid_trends = np.full( (OMI_daily_data['AI'].shape[1],  \
            OMI_daily_data['AI'].shape[2]), np.nan)
        for jj in range(OMI_daily_data['AI'].shape[1]):
            for kk in range(OMI_daily_data['AI'].shape[2]):
                # Plot the individual runs
                # ------------------------
                local_sim_vals = OMI_daily_data['AI'][ii::6,jj,kk]
                x_vals = np.arange(local_sim_vals.shape[0])
    
                #if(trend_type=='standard'): 
                result = stats.linregress(x_vals, local_sim_vals[:])
                grid_trends[jj,kk] = result.slope * len(x_vals)
                #grid_trends[ii,jj,kk] = result.slope * len(x_vals)
    
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,work_mask.compressed())
                    #forcing_trends[i,j] = result.slope * len(x_vals)
                    #forcing_pvals[i,j]  = result.pvalue
                    #forcing_uncert[i,j] = result.stderr * len(x_vals)
                #else:
                #    res = stats.theilslopes(local_sim_vals[jj,:], x_vals, 0.90)
                #    trend_results[jj] = res[0] * len(x_vals)
                #    #forcing_trends[i,j] = res[0]*len(x_vals)
    
    
        mesh = omi_axs[ii].pcolormesh(OMI_daily_data['LON'], OMI_daily_data['LAT'], \
            grid_trends[:,:], transform = ccrs.PlateCarree(), \
            #grid_trends[ii,:,:], transform = ccrs.PlateCarree(), \
            shading = 'auto', cmap = 'bwr', vmin = -0.5, vmax = 0.5)
        omi_axs[ii].set_boundary(circle, transform=omi_axs[ii].transAxes)
        omi_axs[ii].coastlines()
        omi_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        #omi_axs[ii].set_title(titlers[ii])

 
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the mean forcing trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   

        mesh = mean_axs[ii].pcolormesh(daily_dict['longitude'][:], \
            daily_dict['latitude'][:], mean_trends[ii,:,:], \
            cmap = 'bwr', vmin = -vmax, vmax = vmax, \
            transform = ccrs.PlateCarree(), shading = 'auto')
        #cbar = fig.colorbar(mesh, ax = flat_axs[ii], pad = 0.03, fraction = 0.045)
        mean_axs[ii].set_boundary(circle, transform = mean_axs[ii].transAxes)
        mean_axs[ii].coastlines()
        mean_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        #omi_axs[ii].set_title(str(ii))

        pvar = 'mean'
        if(pvar == 'mean'):
            window_min = mean_trends[ii,:,:] - (stdv_trends[ii,:,:] * 1.0)
            window_max = mean_trends[ii,:,:] + (stdv_trends[ii,:,:] * 1.0)

            hasher = np.ma.masked_where(  \
                 (mean_trends[ii,:,:] == 0) |
                ((mean_trends[ii,:,:] < 0) & (window_max[:,:] > 0)) | \
                ((mean_trends[ii,:,:] > 0) & (window_min[:,:] < 0)), mean_trends[ii,:,:])

            #hasher = np.ma.masked_where( (np.abs(plot_trends[ii,:,:]) <= 1.0 * stdv_trends[ii,:,:]), plot_trends[ii,:,:])
            mean_axs[ii].pcolor(daily_dict['longitude'][:], \
                daily_dict['latitude'][:], hasher, hatch = '....', alpha=0., \
                shading = 'auto', transform = ccrs.PlateCarree())


        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the stdv forcing trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        mesh = stdv_axs[ii].pcolormesh(daily_dict['longitude'][:], \
            daily_dict['latitude'][:], stdv_trends[ii,:,:], \
            cmap = 'jet', vmin = 0, vmax = vmax, \
            transform = ccrs.PlateCarree(), shading = 'auto')
        #cbar = fig.colorbar(mesh, ax = flat_axs[ii], pad = 0.03, fraction = 0.045)
        stdv_axs[ii].set_boundary(circle, transform = stdv_axs[ii].transAxes)
        stdv_axs[ii].coastlines()
        stdv_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        #omi_axs[ii].set_title(str(ii))


    #row_label_size
    #fig.text(0.025, 0.65, 'Apr', ha='center', va='center', \
    #    rotation='vertical',weight='bold',fontsize=row_label_size + 0)
    #fig.text(0.025, 0.22, 'May', ha='center', va='center', \
    #    rotation='vertical',weight='bold',fontsize=row_label_size + 0)


    fig.tight_layout()
    plt.show()










def plot_grid_OMI_climo_spatial(OMI_daily_data, min_AI = None, \
        max_AI = None, minlat = 65.5, maxlat = 90.5,  save = False):
    
    if(isinstance(OMI_daily_data, str)):
        if(min_AI is None):
            print("ERROR: If providing an OMI daily file, must provide min_AI")
            return

        print("Reading daily OMI data from", OMI_daily_data)
        print("Converting daily data to monthly omi_uvai_pert averages")
        OMI_daily_data = calcOMI_MonthAvg_FromDaily(OMI_daily_data, \
            min_AI = min_AI, max_AI = max_AI, minlat = minlat, maxlat = maxlat)

    # Calculate gridded trends
    #grid_trends = np.full( (6, OMI_daily_data['AI'].shape[1],  \
    #    OMI_daily_data['AI'].shape[2]), np.nan)
   
    plt.close('all') 
    fig = plt.figure(figsize = (9, 8.5))
    ax1 = fig.add_subplot(2,3,1, projection = ccrs.NorthPolarStereo())
    ax2 = fig.add_subplot(2,3,2, projection = ccrs.NorthPolarStereo())
    ax3 = fig.add_subplot(2,3,3, projection = ccrs.NorthPolarStereo())
    ax4 = fig.add_subplot(2,3,4, projection = ccrs.NorthPolarStereo())
    ax5 = fig.add_subplot(2,3,5, projection = ccrs.NorthPolarStereo())
    ax6 = fig.add_subplot(2,3,6, projection = ccrs.NorthPolarStereo())
    axs = [ax1, ax2, ax3, ax4, ax5, ax6]
    titlers = ['Apr','May','June','July','Aug','Sep']
    
    for ii in range(6):
        print(ii)
        grid_trends = np.full( (OMI_daily_data['AI'].shape[1],  \
            OMI_daily_data['AI'].shape[2]), np.nan)

        grid_climo = np.nanmean( OMI_daily_data['AI'][ii::6,:,:], axis = 0)
        #for jj in range(OMI_daily_data['AI'].shape[1]):
        #    for kk in range(OMI_daily_data['AI'].shape[2]):
        #        # Plot the individual runs
        #        # ------------------------
        #        local_sim_vals = OMI_daily_data['AI'][ii::6,jj,kk]
        #        x_vals = np.arange(local_sim_vals.shape[0])
    
        #        #if(trend_type=='standard'): 
        #        result = stats.linregress(x_vals, local_sim_vals[:])
        #        grid_trends[jj,kk] = result.slope * len(x_vals)
        #        #grid_trends[ii,jj,kk] = result.slope * len(x_vals)
    
        #            #slope, intercept, r_value, p_value, std_err = \
        #            #    stats.linregress(x_vals,work_mask.compressed())
        #            #forcing_trends[i,j] = result.slope * len(x_vals)
        #            #forcing_pvals[i,j]  = result.pvalue
        #            #forcing_uncert[i,j] = result.stderr * len(x_vals)
        #        #else:
        #        #    res = stats.theilslopes(local_sim_vals[jj,:], x_vals, 0.90)
        #        #    trend_results[jj] = res[0] * len(x_vals)
        #        #    #forcing_trends[i,j] = res[0]*len(x_vals)
   
        mean_val = np.nanmean(grid_climo) 
    
    #for ii in range(6):
        mesh = axs[ii].pcolormesh(OMI_daily_data['LON'], OMI_daily_data['LAT'], \
            grid_climo[:,:], transform = ccrs.PlateCarree(), \
            #grid_trends[ii,:,:], transform = ccrs.PlateCarree(), \
            shading = 'auto', cmap = 'jet', vmin = -0.25, vmax = 0.25)
        axs[ii].set_boundary(circle, transform=axs[ii].transAxes)
        axs[ii].coastlines()
        axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        axs[ii].set_title(titlers[ii] + '\nμ = ' + str(np.round(mean_val, 1)))

    plot_subplot_label(ax1, 'a)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax2, 'b)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax3, 'c)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax4, 'd)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax5, 'e)', fontsize = 12, backgroundcolor = 'white')
    plot_subplot_label(ax6, 'f)', fontsize = 12, backgroundcolor = 'white')

    plt.suptitle('OMI UVAI Climo\n2005 - 2020')
    cbar_ax = fig.add_axes([0.17, 0.13, 0.70, 0.01])
    cbar = fig.colorbar(mesh, \
        cax = cbar_ax, shrink = 0.8, orientation = 'horizontal', \
        extend = 'both', \
        label = 'Perturbed OMI UVAI')
    
    fig.tight_layout(rect = [0,0.1,1,1])
   
    if(save):
        if(min_AI < 0):
            outname = 'omi_grid_climo_minAI_N' + str(int(abs(min_AI) * 100)).zfill(3) + '.png'
        else:
            outname = 'omi_grid_climo_minAI' + str(int(abs(min_AI) * 100)).zfill(3) + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show() 

def plot_grid_OMI_trends_arctic_avg(OMI_daily_data, min_AI = None, \
        max_AI = None, minlat = 65.5, maxlat = 90.5,  flat_axs = None, \
        save = False):
    
    if(isinstance(OMI_daily_data, str)):
        if(min_AI is None):
            print("ERROR: If providing an OMI daily file, must provide min_AI")
            return

        print("Reading daily OMI data from", OMI_daily_data)
        print("Converting daily data to monthly omi_uvai_pert averages")
        OMI_daily_data = calcOMI_MonthAvg_FromDaily(OMI_daily_data, \
            min_AI = min_AI, max_AI = max_AI, minlat = minlat, maxlat = maxlat)

        local_AI = OMI_daily_data['AI'][:,:,:]

    else:

        # See if the user wants different latitude ranges plotted
        # -------------------------------------------------------
        test_lats = OMI_daily_data['LAT'][:,0] 
        lat_idxs = np.where( (test_lats >= minlat) & (test_lats <= maxlat) )[0]
        min_lat_idx = lat_idxs[0]
        max_lat_idx = lat_idxs[-1] + 1
        test_lats = test_lats[lat_idxs]
        print("Working with latitudes",test_lats)
        #sim_values = sim_values[:,:,min_lat_idx:max_lat_idx,:]
        local_AI = OMI_daily_data['AI'][:,min_lat_idx:max_lat_idx,:]

    arctic_avgs = np.nanmean(local_AI, axis = (1, 2))

    in_ax = True
    if(flat_axs is None): 
        in_ax = False
        plt.close('all') 
        fig = plt.figure(figsize = (9, 6))
        axs = fig.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True)
        flat_axs = axs.flatten()

    for ii in range(flat_axs.shape[0]):

        # Plot the individual runs
        # ------------------------
        local_AI_vals = arctic_avgs[ii::6]
        x_vals = np.arange(local_AI_vals.shape[0])
        flat_axs[ii].plot(arctic_avgs[ii::6], color = 'k')

        plot_trend_line(flat_axs[ii], np.arange(local_AI_vals.shape[0]), \
            local_AI_vals, color='tab:red', linestyle = '-', \
            slope = 'linregress')
        #flat_axs[ii].axhline(0, linestyle = '--', color = 'k') 

       
    if(not in_ax): 
        fig.tight_layout() 
        plt.show()

def plot_grid_OMI_trends_arctic_avg_combined(OMI_daily_data, min_AI = None, \
        max_AI = None, save = False):

    plt.close('all')
    fig = plt.figure(figsize = (9, 12))
    axs = fig.subplots(nrows = 6, ncols = 3, sharex = True, sharey = True)

    # Plot regional-averaged trends over the entire Arctic
    # ----------------------------------------------------
    minlat = 65.5
    maxlat = 90.5
    plot_grid_OMI_trends_arctic_avg(OMI_daily_data, min_AI = min_AI, \
        max_AI = max_AI, minlat = minlat, maxlat = maxlat,  flat_axs = axs[:,0])

    # Plot regional-averaged trends over the lower Arctic
    # ----------------------------------------------------
    minlat = 65.5
    maxlat = 75.5
    plot_grid_OMI_trends_arctic_avg(OMI_daily_data, min_AI = min_AI, \
        max_AI = max_AI, minlat = minlat, maxlat = maxlat,  flat_axs = axs[:,1])

    # Plot regional-averaged trends over the lower Arctic
    # ----------------------------------------------------
    minlat = 75.5
    maxlat = 90.5
    plot_grid_OMI_trends_arctic_avg(OMI_daily_data, min_AI = min_AI, \
        max_AI = max_AI, minlat = minlat, maxlat = maxlat,  flat_axs = axs[:,2])

    axs[0,0].set_title('65 - 90')
    axs[0,1].set_title('65 - 75')
    axs[0,2].set_title('75 - 90')
    plt.suptitle('Perturbed OMI UVAI Regional Averages')

    fig.tight_layout()
    plt.show()


def plot_test_trends_stdevs(daily_dict, forcing_trends, meanmax = 1.0, \
        stdmax = 0.8, show_min_maxs = False):

    mean_trends = np.mean(forcing_trends, axis = 0)
    stdv_trends = np.std(forcing_trends, axis = 0)
    if(show_min_maxs):
        min_trends  = np.min(forcing_trends, axis = 0)
        max_trends  = np.max(forcing_trends, axis = 0)
   
    plt.close('all') 
    if(show_min_maxs):
        fig = plt.figure(figsize = (14, 11))
        ax1  = fig.add_subplot(4,6,1, projection = ccrs.NorthPolarStereo())
        ax2  = fig.add_subplot(4,6,2, projection = ccrs.NorthPolarStereo())
        ax3  = fig.add_subplot(4,6,3, projection = ccrs.NorthPolarStereo())
        ax4  = fig.add_subplot(4,6,4, projection = ccrs.NorthPolarStereo())
        ax5  = fig.add_subplot(4,6,5, projection = ccrs.NorthPolarStereo())
        ax6  = fig.add_subplot(4,6,6, projection = ccrs.NorthPolarStereo())
        ax7  = fig.add_subplot(4,6,7, projection = ccrs.NorthPolarStereo())
        ax8  = fig.add_subplot(4,6,8, projection = ccrs.NorthPolarStereo())
        ax9  = fig.add_subplot(4,6,9, projection = ccrs.NorthPolarStereo())
        ax10 = fig.add_subplot(4,6,10, projection = ccrs.NorthPolarStereo())
        ax11 = fig.add_subplot(4,6,11, projection = ccrs.NorthPolarStereo())
        ax12 = fig.add_subplot(4,6,12, projection = ccrs.NorthPolarStereo())
        ax13 = fig.add_subplot(4,6,13, projection = ccrs.NorthPolarStereo())
        ax14 = fig.add_subplot(4,6,14, projection = ccrs.NorthPolarStereo())
        ax15 = fig.add_subplot(4,6,15, projection = ccrs.NorthPolarStereo())
        ax16 = fig.add_subplot(4,6,16, projection = ccrs.NorthPolarStereo())
        ax17 = fig.add_subplot(4,6,17, projection = ccrs.NorthPolarStereo())
        ax18 = fig.add_subplot(4,6,18, projection = ccrs.NorthPolarStereo())
        ax19 = fig.add_subplot(4,6,19, projection = ccrs.NorthPolarStereo())
        ax20 = fig.add_subplot(4,6,20, projection = ccrs.NorthPolarStereo())
        ax21 = fig.add_subplot(4,6,21, projection = ccrs.NorthPolarStereo())
        ax22 = fig.add_subplot(4,6,22, projection = ccrs.NorthPolarStereo())
        ax23 = fig.add_subplot(4,6,23, projection = ccrs.NorthPolarStereo())
        ax24 = fig.add_subplot(4,6,24, projection = ccrs.NorthPolarStereo())
    else:
        fig = plt.figure(figsize = (14, 6))
        ax1  = fig.add_subplot(2,6,1, projection = ccrs.NorthPolarStereo())
        ax2  = fig.add_subplot(2,6,2, projection = ccrs.NorthPolarStereo())
        ax3  = fig.add_subplot(2,6,3, projection = ccrs.NorthPolarStereo())
        ax4  = fig.add_subplot(2,6,4, projection = ccrs.NorthPolarStereo())
        ax5  = fig.add_subplot(2,6,5, projection = ccrs.NorthPolarStereo())
        ax6  = fig.add_subplot(2,6,6, projection = ccrs.NorthPolarStereo())
        ax7  = fig.add_subplot(2,6,7, projection = ccrs.NorthPolarStereo())
        ax8  = fig.add_subplot(2,6,8, projection = ccrs.NorthPolarStereo())
        ax9  = fig.add_subplot(2,6,9, projection = ccrs.NorthPolarStereo())
        ax10 = fig.add_subplot(2,6,10, projection = ccrs.NorthPolarStereo())
        ax11 = fig.add_subplot(2,6,11, projection = ccrs.NorthPolarStereo())
        ax12 = fig.add_subplot(2,6,12, projection = ccrs.NorthPolarStereo())
    flat_axs = [ax1, ax2, ax3, ax4, ax5, ax6]
    flat_axs2 = [ax7, ax8, ax9, ax10, ax11, ax12]
    if(show_min_maxs):
        flat_axs3 = [ax13, ax14, ax15, ax16, ax17, ax18]
        flat_axs4 = [ax19, ax20, ax21, ax22, ax23, ax24]
    #axs = fig.subplots(nrows = 2, ncols = 3)
    #f8at_axs = axs.flatten()
    
    for ii in range(6):

        #hash_idxs = np.where( np.abs(mean_trends[ii,:,:]) > stdv_trends[ii,:,:])

        mesh = flat_axs[ii].pcolormesh(daily_dict['longitude'][:], daily_dict['latitude'][:], mean_trends[ii,:,:], \
            cmap = 'bwr', vmin = -meanmax, vmax = meanmax, transform = ccrs.PlateCarree(), shading = 'auto')
        cbar = fig.colorbar(mesh, ax = flat_axs[ii], pad = 0.03, fraction = 0.045)
        flat_axs[ii].set_boundary(circle, transform=flat_axs[ii].transAxes)
        flat_axs[ii].coastlines()
        flat_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        flat_axs[ii].set_title(str(ii))

        hasher = np.ma.masked_where( (np.abs(mean_trends[ii,:,:]) <= stdv_trends[ii,:,:]), mean_trends[ii,:,:])
        flat_axs[ii].pcolor(daily_dict['longitude'][:], daily_dict['latitude'][:], \
            hasher, hatch = '....', alpha=0., shading = 'auto', transform = ccrs.PlateCarree())
        #hasher = np.ma.masked_where(out_dict['raw_stderr'][0,:,:] > maxerr, \
        #    out_dict['raw_stderr'][0,:,:])
        #ax1.pcolor(out_dict['sza_mins'], out_dict[cloud_var + '_mins'], \
        #    hasher, hatch = hatch_style, alpha=0., shading = 'auto')

        mesh = flat_axs2[ii].pcolormesh(daily_dict['longitude'][:], daily_dict['latitude'][:], stdv_trends[ii,:,:], \
            cmap = 'jet', vmin = 0, vmax = stdmax, transform = ccrs.PlateCarree(), shading = 'auto')
        cbar = fig.colorbar(mesh, ax = flat_axs2[ii], pad = 0.03, fraction = 0.045)
        flat_axs2[ii].set_boundary(circle, transform=flat_axs2[ii].transAxes)
        flat_axs2[ii].coastlines()
        flat_axs2[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        flat_axs2[ii].set_title(str(ii))
   
        if(show_min_maxs): 
            mesh = flat_axs3[ii].pcolormesh(daily_dict['longitude'][:], daily_dict['latitude'][:], min_trends[ii,:,:], \
                cmap = 'bwr', vmin = -meanmax, vmax = meanmax, transform = ccrs.PlateCarree(), shading = 'auto')
            cbar = fig.colorbar(mesh, ax = flat_axs3[ii], pad = 0.03, fraction = 0.045)
            flat_axs3[ii].set_boundary(circle, transform=flat_axs3[ii].transAxes)
            flat_axs3[ii].coastlines()
            flat_axs3[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
            flat_axs3[ii].set_title(str(ii))
    
            mesh = flat_axs4[ii].pcolormesh(daily_dict['longitude'][:], daily_dict['latitude'][:], max_trends[ii,:,:], \
                cmap = 'bwr', vmin = -meanmax, vmax = meanmax, transform = ccrs.PlateCarree(), shading = 'auto')
            cbar = fig.colorbar(mesh, ax = flat_axs4[ii], pad = 0.03, fraction = 0.045)
            flat_axs4[ii].set_boundary(circle, transform=flat_axs4[ii].transAxes)
            flat_axs4[ii].coastlines()
            flat_axs4[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
            flat_axs4[ii].set_title(str(ii))
    
    fig.tight_layout()
    plt.show()
    

def plot_sim_errors_bulk_arctic_avg_combined(daily_filename, err_mean, err_std, \
        minlat = 65.5, maxlat = 90.5, num_sims = 100, sim_values = None, \
        plot_result_min_max_range = False, trend_type = 'linregress', \
        flat_axs = None):

    plt.close('all')
    fig = plt.figure(figsize = (9, 12))
    axs = fig.subplots(nrows = 6, ncols = 3, sharex = True, sharey = True)

    # Plot regional-averaged trends over the entire Arctic
    # ----------------------------------------------------
    minlat = 65.5
    maxlat = 90.5
    plot_sim_errors_bulk_arctic_avg(daily_filename, err_mean, \
        err_std, minlat = minlat, maxlat = maxlat, num_sims = num_sims, \
        sim_values = sim_values, plot_result_min_max_range = plot_result_min_max_range, \
        flat_axs = axs[:,0])

    # Plot regional-averaged trends over the lower Arctic
    # ----------------------------------------------------
    minlat = 65.5
    maxlat = 75.5
    plot_sim_errors_bulk_arctic_avg(daily_filename, err_mean, \
        err_std, minlat = minlat, maxlat = maxlat, num_sims = num_sims, \
        sim_values = sim_values, plot_result_min_max_range = plot_result_min_max_range, \
        flat_axs = axs[:,1])

    # Plot regional-averaged trends over the lower Arctic
    # ----------------------------------------------------
    minlat = 75.5
    maxlat = 90.5
    plot_sim_errors_bulk_arctic_avg(daily_filename, err_mean, \
        err_std, minlat = minlat, maxlat = maxlat, num_sims = num_sims, \
        sim_values = sim_values, plot_result_min_max_range = plot_result_min_max_range, \
        flat_axs = axs[:,2])

    fig.tight_layout()
    plt.show()


def plot_sim_errors_bulk_arctic_avg(daily_filename, err_mean, err_std, \
        minlat = 65.5, maxlat = 90.5, num_sims = 100, sim_values = None, \
        plot_result_min_max_range = False, trend_type = 'linregress', \
        flat_axs = None, return_std = False):

    daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)

    if(sim_values is None):
        if(return_std):
            sim_values, std_values = calc_force_vals_bulk(daily_dict, err_mean, err_std, \
                minlat = minlat, maxlat = maxlat, num_sims = num_sims, return_std = True)
        else:
            sim_values = calc_force_vals_bulk(daily_dict, err_mean, err_std, \
                minlat = minlat, maxlat = maxlat, num_sims = num_sims, return_std = False)
    else:
        # See if the user wants different latitude ranges plotted
        # -------------------------------------------------------
        test_lats = daily_dict['latitude'] 
        lat_idxs = np.where( (test_lats >= minlat) & (test_lats <= maxlat) )[0]
        min_lat_idx = lat_idxs[0]
        max_lat_idx = lat_idxs[-1] + 1
        test_lats = test_lats[lat_idxs]
        print("Working with latitudes",test_lats)
        sim_values = sim_values[:,:,min_lat_idx:max_lat_idx,:]

    # Calculate spatial averages
    # --------------------------
    arctic_avgs = np.nanmean(sim_values, axis = (2,3))
    #august_avgs = arctic_avgs[:,4::6]
  
    # Set up an array to hold the trend results for each run
    # ------------------------------------------------------
    trend_results = np.full( sim_values.shape[0], np.nan)

    in_ax = True
    if(flat_axs is None): 
        in_ax = False
        plt.close('all') 
        fig = plt.figure(figsize = (9, 6))
        axs = fig.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True)
        flat_axs = axs.flatten()

    for ii in range(flat_axs.shape[0]):

        # Plot the individual runs
        # ------------------------
        local_sim_vals = arctic_avgs[:,ii::6]
        x_vals = np.arange(local_sim_vals.shape[1])
        if(plot_result_min_max_range):
            min_line = np.min(local_sim_vals, axis = 0)
            max_line = np.max(local_sim_vals, axis = 0)
            flat_axs[ii].plot(min_line, alpha = 0.75, color = 'gray')
            flat_axs[ii].plot(max_line, alpha = 0.75, color = 'gray')
            flat_axs[ii].fill_between(x_vals, min_line, max_line, \
                alpha = 0.75, color = 'gray')
        else:
            flat_axs[ii].plot(arctic_avgs[:,ii::6].T, alpha = 0.5)

        # Calculate the trends for all the individual runs
        # ------------------------------------------------
        for jj in range(local_sim_vals.shape[0]):
            if(trend_type=='standard'): 
                result = stats.linregress(x_vals, local_sim_vals[jj,:])
                trend_results[jj] = result.slope * len(x_vals)

                #slope, intercept, r_value, p_value, std_err = \
                #    stats.linregress(x_vals,work_mask.compressed())
                #forcing_trends[i,j] = result.slope * len(x_vals)
                #forcing_pvals[i,j]  = result.pvalue
                #forcing_uncert[i,j] = result.stderr * len(x_vals)
            else:
                res = stats.theilslopes(local_sim_vals[jj,:], x_vals, 0.90)
                trend_results[jj] = res[0] * len(x_vals)
                #forcing_trends[i,j] = res[0]*len(x_vals)
            

        # Plot the mean of the runs
        # -------------------------
        sim_avg = np.nanmean(local_sim_vals, axis = 0)
        flat_axs[ii].plot(sim_avg, color = 'k')
        plot_trend_line(flat_axs[ii], np.arange(sim_avg.shape[0]), sim_avg, color='tab:red', linestyle = '-', \
            slope = 'linregress')
        flat_axs[ii].axhline(0, linestyle = '--', color = 'k') 

        # Calculate the trend of the mean time series
        # -------------------------------------------
        if((trend_type=='standard') | (trend_type == 'linregress')): 
            result = stats.linregress(x_vals, sim_avg)
            trend_of_mean = result.slope * len(x_vals)
        else:
            res = stats.theilslopes(sim_avg, x_vals, 0.90)
            trend_of_mean = res[0] * len(x_vals)

        # Calculate the mean of the (num_sims) trends
        # -------------------------------------------
        mean_of_trends = np.mean(trend_results)
        stdv_of_trends = np.std(trend_results)
        min_of_trends  = np.min(trend_results)
        max_of_trends  = np.max(trend_results)

        print("\nTrend of mean for month_idx",ii,":",trend_of_mean)
        print("\nMean of trend for month_idx",ii,":",mean_of_trends)
        print("StDv of trend for month_idx",ii,":",stdv_of_trends)
        print("Min  of trend for month_idx",ii,":",min_of_trends)
        print("Max  of trend for month_idx",ii,":",max_of_trends)

       
    if(not in_ax): 
        fig.tight_layout() 
        plt.show()


def test_plot_L2L3_err_sims(current_files, minlat = None, maxlat = None, save = False):

    plt.close('all')
    fig = plt.figure()
    axs = fig.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True)
    flat_axs = axs.flatten()
    titles = ['April','May','June','July','August','September']
    for tfile in current_files:
        data_dict = read_daily_month_force_L2L3_error_from_HDF5(tfile)
    
   
        if(minlat is not None):
            beg_idx = np.where(data_dict['latitude'] >= (minlat + 0.5))[0][0]
        else:
            beg_idx = 0
        if(maxlat is not None):
            end_idx = np.where(data_dict['latitude'] <= (maxlat + 0.5))[0][-1] + 1
        else:
            end_idx = None

        print('HERE3', data_dict['latitude'][beg_idx:end_idx])
 
        arctic_avgs = np.nanmean(data_dict['force_estimate'][:,:,beg_idx:end_idx,:], axis = (2, 3))
        if('force_estimate_orig' in data_dict.keys()):
            arctic_avgs_orig = np.nanmean(data_dict['force_estimate_orig'][:,beg_idx:end_idx,:], axis = (1, 2))
    
        # Plot the error results
        # ----------------------
        for ii in range(len(flat_axs)):
            flat_axs[ii].plot(arctic_avgs[:,ii::6].T)
    
            if('force_estimate_orig' in data_dict.keys()):
                flat_axs[ii].plot(arctic_avgs_orig[ii::6], color = 'k')
   
    for ii in range(6):     
        flat_axs[ii].set_title(titles[ii]) 
        flat_axs[ii].axhline(0, linestyle = ':', color = 'k')

    fig.tight_layout()
    plt.show()

def plot_ice_error_histogram(daily_filename, ice_filename, \
        log_scale = False, num_bins = 50, xmin = None, xmax = None, \
        use_correct_error_calc = False, \
        astrofit = False, ax = None, save = False):
    
    # Read the original daily values
    # ------------------------------
    daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)

    # Read the modified daily values
    # ------------------------------
    ice_dict = read_daily_month_force_L2L3_error_from_HDF5(ice_filename)

    # Make sure the versions are the same. If not, raise an error
    # -----------------------------------------------------------
    if( daily_dict['slope_dict']['sim_name'] != ice_dict['slope_dict']['sim_name']):
        print("# # # # # # # # # # # # # # # # # # # # # #")
        print("")
        print("    ERROR: MISMATCHED SIM TYPES")
        print("")
        print("    Daily dict sim name = ", daily_dict['slope_dict']['sim_name'])
        print("    Ice  dict sim name  = ", ice_dict['slope_dict']['sim_name'])
        print("")
        print("# # # # # # # # # # # # # # # # # # # # # #")
    

    # Flatten the arrays
    # ------------------
    flat_orig = daily_dict['force_estimate_orig'].flatten()
    flat_ice  = ice_dict['force_estimate_orig'].flatten()

    # Figure out the indices which contain non-zero forcing values and are not nans.
    # Only want to find differences when forcing values were actually calculated
    # ------------------------------------------------------------------------------
    good_idxs = np.where( (flat_orig != 0) & (np.isnan(flat_orig) == False) & \
        (flat_ice != 0) & (np.isnan(flat_ice) == False))

    # Calculate the differences between the original and modified L3 daily values
    # ---------------------------------------------------------------------------
    if(use_correct_error_calc):
        print("Calculating ice errors using: ice_err - original")
        xlabel_text = 'Forcing error (ADRF$_{ice    error}$ - ADRF$_{orig}$) [Wm$^{-2}$]'
        diffs = flat_ice[good_idxs] - flat_orig[good_idxs]
    else:
        print("Calculating ice errors using: original - ice_err")
        #xlabel_text = 'Forcing error (L2 - L3) [Wm$^{-2}$]'
        xlabel_text = 'Forcing error (ADRF$_{orig}$ - ADRF$_{ice    error}$) [Wm$^{-2}$]'
        diffs = flat_orig[good_idxs] - flat_ice[good_idxs]

    # Calculate statistics
    # --------------------
    mean_err = np.nanmean(diffs)
    std_err  = np.nanstd(diffs)

    # Plot the results
    # ----------------
    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    hist = ax.hist(diffs, bins = num_bins)
    if(log_scale):
        ax.set_yscale('log')
    ax.set_xlabel(xlabel_text)
    ax.set_ylabel('Counts')
    ax.set_title('L3 Daily Forcing Errors Due to Ice Errors' + \
        '\nIce errors: Mean = ' + str(ice_dict['attributes']['ice_err_mean']) + \
        '%, Std Dev = ' + str(ice_dict['attributes']['ice_err_std']) + \
        '%')

    # Test using astropy to fit a Gaussian function
    # ---------------------------------------------
    if(astrofit):
        bin_heights,bin_borders = np.histogram(diffs,bins=num_bins)
        bin_widths = np.diff(bin_borders)
        bin_centers = bin_borders[:-1] + bin_widths / 2
        t_init = models.Gaussian1D()
        fit_t = fitting.LevMarLSQFitter()
        t = fit_t(t_init, bin_centers, bin_heights)
        print('Amplitude: ',np.round(t.amplitude.value,3))
        print('Mean:      ',np.round(t.mean.value,3))
        print('StDev:     ',np.round(t.stddev.value,3))
        x_interval_for_fit = np.linspace(bin_centers[0],bin_centers[-1],200)
        ax.plot(x_interval_for_fit,t(x_interval_for_fit),label='fit',c='tab:red')
        #xdiff = ax.get_xlim()[1] - ax.get_xlim()[0]
        #plot_xval = ax.get_xlim()[0]  + xdiff * 0.8
        #plot_yval = ax.get_ylim()[1] * 0.5
        #ptext = 'μ = ' + str(np.round(t.mean.value, 1)) + \
        #    '\nσ = ' + str(np.round(t.stddev.value, 1))
        #plot_figure_text(ax, ptext, xval = plot_xval, yval = plot_yval, \
        #    fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')


        ax.set_xlim(xmin, xmax)

        xdiff = ax.get_xlim()[1] - ax.get_xlim()[0]
        plot_xval = ax.get_xlim()[0]  + xdiff * 0.65
        plot_yval = ax.get_ylim()[1] * 0.8
        ptext = 'μ = ' + str(np.round(t.mean.value, 1)) + ' Wm$^{-2}$\n' + \
            'σ = ' + str(np.round(t.stddev.value, 1)) + ' Wm$^{-2}$'


    else:

        ax.set_xlim(xmin, xmax)
        # Add statistics to plot
        # ----------------------
        ptext = 'μ = ' + str(np.round(mean_err, 1)) + ' Wm$^{-2}$' + \
            '\nσ = ' + str(np.round(std_err, 1)) + ' Wm$^{-2}$'

        xdiff = ax.get_xlim()[1] - ax.get_xlim()[0]
        plot_xval = ax.get_xlim()[0]  + xdiff * 0.65
        if(log_scale):
            plot_yval = ax.get_ylim()[1] * 0.05
        else:
            plot_yval = ax.get_ylim()[1] * 0.5
        #plot_figure_text(ax3, ptext, location = 'lower_right', \
        #ax5.set_title(ptext)
        plot_figure_text(ax, ptext, xval = plot_xval, yval = plot_yval, \
            fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')
   
        ax.axvline(0, linestyle = '--', color = 'k', alpha = 0.5)
 

    if(not in_ax):
        fig.tight_layout()
        if(save):
            num_ice_bins = len(daily_dict['bin_dict']['ice_bin_means'])
            if(log_scale):
                log_adder = '_logscale'
            else:
                log_adder = ''
            if(use_correct_NN_error_calc):
                err_calc_add = '_correctNNerr'
            else:
                err_calc_add = ''
            outname = 'ice_error_daily_forcing_impacts_numsfcbins' + \
                str(num_ice_bins) + '_' + daily_dict['slope_dict']['sim_name'] + \
                log_adder + err_calc_add + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show()

def plot_cod_error_histogram(daily_filename, cod_filename, num_bins = 50, \
        log_scale = False, xmin = None, xmax = None, ax = None, \
        use_correct_error_calc = False, \
        astrofit = False, save = False):
    
    # Read the original daily values
    # ------------------------------
    daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)

    # Read the modified daily values
    # ------------------------------
    cod_dict = read_daily_month_force_L2L3_error_from_HDF5(cod_filename)

    # Make sure the versions are the same. If not, raise an error
    # -----------------------------------------------------------
    if( daily_dict['slope_dict']['sim_name'] != cod_dict['slope_dict']['sim_name']):
        print("# # # # # # # # # # # # # # # # # # # # # #")
        print("")
        print("    ERROR: MISMATCHED SIM TYPES")
        print("")
        print("    Daily dict sim name = ", daily_dict['slope_dict']['sim_name'])
        print("    COD  dict sim name  = ", cod_dict['slope_dict']['sim_name'])
        print("")
        print("# # # # # # # # # # # # # # # # # # # # # #")

    # Flatten the arrays
    # ------------------
    flat_orig = daily_dict['force_estimate_orig'].flatten()
    flat_cod  = cod_dict['force_estimate_orig'].flatten()

    # Figure out the indices which contain non-zero forcing values and are not nans.
    # Only want to find differences when forcing values were actually calculated
    # ------------------------------------------------------------------------------
    good_idxs = np.where( (flat_orig != 0) & (np.isnan(flat_orig) == False) & \
        (flat_cod != 0) & (np.isnan(flat_cod) == False))

    # Calculate the differences between the original and modified L3 daily values
    # ---------------------------------------------------------------------------
    if(use_correct_error_calc):
        print("Calculating cod errors using: cod_err - original")
        xlabel_text = 'Forcing error (ADRF$_{cod error}$ - ADRF$_{orig}$) [Wm$^{-2}$]'
        diffs = flat_cod[good_idxs] - flat_orig[good_idxs]
    else:
        print("Calculating cod errors using: original - ice_err")
        #xlabel_text = 'Forcing error (L2 - L3) [Wm$^{-2}$]'
        xlabel_text = 'Forcing error (ADRF$_{orig}$ - ADRF$_{cod error}$) [Wm$^{-2}$]'
        diffs = flat_orig[good_idxs] - flat_cod[good_idxs]

    print("Num values = ", diffs.shape[0])

    # Calculate statistics
    # --------------------
    mean_err = np.nanmean(diffs)
    std_err  = np.nanstd(diffs)

    # Plot the results
    # ----------------
    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    hist = ax.hist(diffs, bins = num_bins)
    if(log_scale):
        ax.set_yscale('log')
    ax.set_xlabel(xlabel_text)
    ax.set_ylabel('Counts')
    ax.set_title('L3 Daily Forcing Errors Due to COD Errors' + \
        '\nCOD errors: Mean = ' + str(cod_dict['attributes']['cod_err_mean']) + \
        ', Std Dev = ' + str(cod_dict['attributes']['cod_err_std']))

    ax.set_xlim(xmin, xmax)
    ax.axvline(0, linestyle = '--', color = 'k', alpha = 0.5)

    # Test using astropy to fit a Gaussian function
    # ---------------------------------------------
    if(astrofit):
        bin_heights,bin_borders = np.histogram(diffs,bins=num_bins)
        bin_widths = np.diff(bin_borders)
        bin_centers = bin_borders[:-1] + bin_widths / 2
        t_init = models.Gaussian1D()
        fit_t = fitting.LevMarLSQFitter()
        t = fit_t(t_init, bin_centers, bin_heights)
        print('Amplitude: ',np.round(t.amplitude.value,3))
        print('Mean:      ',np.round(t.mean.value,3))
        print('StDev:     ',np.round(t.stddev.value,3))
        x_interval_for_fit = np.linspace(bin_centers[0],bin_centers[-1],200)
        ax.plot(x_interval_for_fit,t(x_interval_for_fit),label='fit',c='tab:red')
        xdiff = ax.get_xlim()[1] - ax.get_xlim()[0]
        plot_xval = ax.get_xlim()[0]  + xdiff * 0.8
        plot_yval = ax.get_ylim()[1] * 0.5
        #ptext = 'μ = ' + str(np.round(t.mean.value, 1)) + \
        #    '\nσ = ' + str(np.round(t.stddev.value, 1))
        #plot_figure_text(ax, ptext, xval = plot_xval, yval = plot_yval, \
        #    fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')



    else:

        # Add statistics to plot
        # ----------------------
        ptext = 'μ = ' + str(np.round(mean_err, 1)) + ' Wm$^{-2}$' + \
            '\nσ = ' + str(np.round(std_err, 1)) + ' Wm$^{-2}$'

        xdiff = ax.get_xlim()[1] - ax.get_xlim()[0]
        plot_xval = ax.get_xlim()[0]  + xdiff * 0.65
        if(log_scale):
            plot_yval = ax.get_ylim()[1] * 0.05
        else:
            plot_yval = ax.get_ylim()[1] * 0.5
        #plot_figure_text(ax3, ptext, location = 'lower_right', \
        #ax5.set_title(ptext)
        plot_figure_text(ax, ptext, xval = plot_xval, yval = plot_yval, \
            fontsize = 10, color = 'black', weight = None, backgroundcolor = 'white')
   
        ax.set_xlim(xmin, xmax)
        ax.axvline(0, linestyle = '--', color = 'k', alpha = 0.5)





    if(not in_ax):
        fig.tight_layout()
        if(save):
            num_ice_bins = len(daily_dict['bin_dict']['ice_bin_means'])
            if(log_scale):
                log_adder = '_logscale'
            else:
                log_adder = ''
            if(use_correct_NN_error_calc):
                err_calc_add = '_correctNNerr'
            else:
                err_calc_add = ''
            outname = 'cod_error_daily_forcing_impacts_numsfcbins' + \
                str(num_ice_bins) + '_' + daily_dict['slope_dict']['sim_name']+ \
                log_adder + err_calc_add + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show()

def plot_error_components_combined(direct_forcings, calc_forcings, sim_name, \
        daily_filename, ice_filename, cod_filename, num_bins, \
        use_correct_error_calc = False, \
        astrofit = False, log_scale = False, xmin = -100, xmax = 100, \
        run_type = None, save = False):

    plt.close('all')
    fig = plt.figure(figsize = (9, 8))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    plot_NN_error_dist_bulk(sim_name, ax = ax1, num_bins = num_bins * 4, \
        astrofit = astrofit, \
        xmin = xmin, xmax = xmax, \
        use_correct_error_calc = use_correct_error_calc, \
        save = False)
    plot_hist_L2_L3_errors(direct_forcings, calc_forcings, sim_name, \
        ax = ax2, num_bins = num_bins * 4, astrofit = astrofit, \
        use_correct_error_calc = use_correct_error_calc, \
        xmin = xmin, xmax = xmax, save = False)
    plot_ice_error_histogram(daily_filename, ice_filename, ax = ax3, \
        xmin = -100, xmax = 100, save = False, log_scale = log_scale, \
        use_correct_error_calc = use_correct_error_calc, \
        num_bins = num_bins, astrofit = False)
    plot_cod_error_histogram(daily_filename, cod_filename, ax = ax4, \
        xmin = None, xmax = None, save = False, log_scale = log_scale, \
        use_correct_error_calc = use_correct_error_calc, \
        num_bins = num_bins * 1, astrofit = False)

    plot_subplot_label(ax1, 'a)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax2, 'b)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax3, 'c)', fontsize = 11, backgroundcolor = 'white', yval = 20000)
    plot_subplot_label(ax4, 'd)', fontsize = 11, backgroundcolor = 'white', yval = 20000)

    fig.tight_layout()

    if(save):
        if(log_scale):
            log_adder = '_logscale'
        else:
            log_adder = ''
        if(run_type is not None):
            run_type = '_' + run_type
        else:
            run_type = ''
        if(use_correct_error_calc):
            err_calc_add = '_correcterr'
            #err_calc_add = '_correctNNerr'
        else:
            err_calc_add = ''
        outname = 'errors_L3_combined_' + sim_name + log_adder + run_type + err_calc_add + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_many_trend_test(daily_dict, OMI_data, forcing_trends, month_idx, lat_idx, lon_idx, \
        num_bins, conf_level = 90, save = False):

    plt.close('all')
    fig = plt.figure(figsize = (9, 4))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    x_vals = np.arange(2005, 2021)
    y_vals = OMI_data['AI'][month_idx::6,lat_idx, lon_idx]
    result = stats.linregress(x_vals, y_vals)
    ax1.plot(x_vals, y_vals)
    print(result.pvalue)

    test_error_dist(daily_dict, forcing_trends, month_idx, lat_idx, lon_idx, num_bins, \
        ax = ax2, conf_level = conf_level, save = False)

    fig.tight_layout()
    plt.show()

def test_error_dist(daily_dict, forcing_trends, month_idx, lat_idx, lon_idx, \
        num_bins, ax = None, conf_level = 90, sim_name = '', run_type = '', \
        save = False):

    in_arr = forcing_trends[:,month_idx,lat_idx,lon_idx]

    lat_val = daily_dict['latitude'][lat_idx]
    lon_val = daily_dict['longitude'][lon_idx]

    mean_trend = np.nanmean(in_arr[:])
    stdv_trend = np.nanstd(in_arr[:])
    #stderr_trend = stdv_trend / np.sqrt(in_arr.shape[0])
    stderr_trend = sem(in_arr[:])

    print('Mean trend:      ', mean_trend)
    print('StDv trend:      ', stdv_trend)
    print('StDv trend * 1.5:', stdv_trend * 1.5)
    print('StDv trend * 2.0:', stdv_trend * 2.0)
    print('StdErr trend:    ', stderr_trend)

    stdv_trend *= 1.0

    #conf_intvl = statnorm.interval(alpha = 0.99, \
    #    loc = mean_trend, scale = stderr_trend)

    # 90% conf = 1.645 * stdev
    # 95% conf = 1.965 * stdev
    # 99% conf = 2.576 * stdev

    str_conf_level = str(conf_level)
    if(conf_level > 1):
        conf_level = conf_level / 100.
    else:
        str_conf_level = str(conf_level * 100)
    conf_intvl = statnorm.interval(alpha = conf_level, \
        loc = mean_trend, scale = stdv_trend)
    print('Conf_intvl', conf_intvl)
    print('Othr_intvl', mean_trend - stdv_trend * 1.645, mean_trend + stdv_trend * 1.645)

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    month_names = ['April','May','June','July','August','September']

    hist = ax.hist(in_arr[:], bins = num_bins)
    ax.axvline(mean_trend, color = 'k', linestyle = '--', label = 'Mean')
    #ax.axvline( (mean_trend - stdv_trend), color = 'k', linestyle = ':', label = 'μ +/- ς')
    #ax.axvline( (mean_trend + stdv_trend), color = 'k', linestyle = ':')
    ax.axvline( conf_intvl[0], color = 'k', linestyle = ':', label = str_conf_level + '% Conf.')
    ax.axvline( conf_intvl[1], color = 'k', linestyle = ':')
    ax.legend()
    ax.set_xlabel('Forcing Trend [Wm$^{-2}$ (study period) $^{-1}$]')
    ax.set_ylabel('Counts')
    ax.set_title('Distribution of ' + str(int(forcing_trends.shape[0])) + \
        ' ' + month_names[month_idx] + ' ADRF Trends\n' + str(lat_val) + ' $^{o}$N, ' + str(lon_val) + ' $^{o}$E')

    if(not in_ax):
        if(save):
            if(lon_val < 0):
                loc_adder = '_' + str(int(lat_val)) + 'N' + str(int(lon_val)) + 'W'
            else:
                loc_adder = '_' + str(int(lat_val)) + 'N' + str(int(lon_val)) + 'E'
            if(sim_name != ''):
                sim_name = '_' + sim_name
            if(run_type != ''):
                run_type = '_' + run_type
            outname = 'trend_dist_monthidx' + str(month_idx) + '_numsims' + \
                str(int(in_arr.shape[0])) + loc_adder + sim_name + run_type + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show()    

def plot_bulk_trend_dist_multimonth(daily_dict, forcing_trends, bins = None, \
        save = False, ax = None, log_yscale = False, \
        plot_single_month = None, minlat = None, maxlat = None):

    in_ax = True 
    if(ax is None): 
        plt.close('all')
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    if(minlat is not None):
        beg_idx = np.where(daily_dict['latitude'] >= (minlat + 0.5))[0][0]
    else:
        beg_idx = 0
    if(maxlat is not None):
        end_idx = np.where(daily_dict['latitude'] <= (maxlat + 0.5))[0][-1] + 1
    else:
        end_idx = None

    print('HERE3', daily_dict['latitude'][beg_idx:end_idx])
 

    avg_trends = np.mean(forcing_trends[:,:,beg_idx:end_idx,:], axis = 0)  

    bins = np.linspace(np.min(avg_trends), np.max(avg_trends), 50)

    labels = ['Apr','May','Jun','Jul','Aug','Sep']

    if(plot_single_month is None):
        loop_vals = np.arange(6)
        histtype = 'step'
    else:
        loop_vals = [plot_single_month]
        histtype = 'bar'
     
    for ii in loop_vals:
        flat_trends = avg_trends[ii].flatten() 
        hist = ax.hist(flat_trends, bins = bins, histtype = histtype, label = labels[ii])
       
    ax.set_ylabel('Counts')
    ax.set_xlabel('Monthly Mean Trend [Wm$^{-2}$(study period)$^{-1}$]')
    if(log_yscale):
        ax.set_yscale('log')
    ax.set_title('Distributions of Monthly Gridded DRF Trends')
    ax.legend()

    if(not in_ax):
        fig.tight_layout()
        if(save):
            if(log_yscale):
                scale_adder = '_log'
            else:
                scale_adder = ''

            if(plot_single_month is not None):
                month_adder = '_monthidx' + str(int(plot_single_month))
            else:
                month_adder = ''

            outname = 'force_trend_monthly_dist' + scale_adder + month_adder + '.png'
            fig.savefig(outname, dpi = 200)
            print("Saved image", outname)
        else:
            plt.show()   

def plot_force_trend_mean_std_dist(daily_dict, forcing_trends, month_idx, \
        lat_idx, lon_idx, vmin = -1.5, vmax = 1.5, save = False, conf_window = 90):

    plt.close('all')
    fig = plt.figure(figsize = (10, 4))
    ax1 = fig.add_subplot(1,3,1, projection = ccrs.NorthPolarStereo())
    ax2 = fig.add_subplot(1,3,2, projection = ccrs.NorthPolarStereo())
    ax3 = fig.add_subplot(1,3,3)

    plot_trends = np.nanmean(forcing_trends, axis = 0)
    stdv_trends = np.nanstd(forcing_trends, axis = 0)

    lat_val = daily_dict['latitude'][lat_idx]
    lon_val = daily_dict['longitude'][lon_idx]
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Plot forcing trend
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    mesh = ax1.pcolormesh(daily_dict['longitude'][:], \
        daily_dict['latitude'][:], plot_trends[month_idx,:,:], cmap = 'bwr', \
        vmin = vmin, vmax = vmax, transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax1, pad = 0.03, fraction = 0.045, \
        label = 'Forcing Trend [Wm$^{-2}$]', orientation = 'horizontal')
    ax1.set_boundary(circle, transform = ax1.transAxes)
    ax1.coastlines()
    ax1.set_extent([-180,180,65,90], ccrs.PlateCarree())

    if(conf_window == 90):
        conf_val = 1.645
    elif(conf_window == 95):
        conf_val = 1.960
    else:
        print("CURRENTLY INVALID CONFIDENCE INTERVAL. USING 90")
        conf_val = 1.645

    window_min = plot_trends[month_idx,:,:] - (stdv_trends[month_idx,:,:] * conf_val)
    window_max = plot_trends[month_idx,:,:] + (stdv_trends[month_idx,:,:] * conf_val)

    hasher = np.ma.masked_where(  \
         (plot_trends[month_idx,:,:] == 0) |
         (np.abs(plot_trends[month_idx,:,:]) < \
          (stdv_trends[month_idx,:,:] * conf_val)), plot_trends[month_idx, :,:])

    #hasher = np.ma.masked_where(  \
    #     (plot_trends[month_idx,:,:] == 0) |
    #    ((plot_trends[month_idx,:,:] < 0) & (window_max[:,:] > 0)) | \
    #    ((plot_trends[month_idx,:,:] > 0) & (window_min[:,:] < 0)), plot_trends[month_idx, :,:])

    #hasher = np.ma.masked_where( (np.abs(plot_trends[ii,:,:]) <= 1.0 * stdv_trends[ii,:,:]), plot_trends[ii,:,:])
    ax1.pcolor(daily_dict['longitude'][:], daily_dict['latitude'][:], \
        hasher, hatch = '....', alpha=0., shading = 'auto', \
        transform = ccrs.PlateCarree())

    plot_point_on_map(ax1, lat_val, lon_val, markersize = 10, color = 'red', \
        alpha = 1.0, add_border = True)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Plot forcing trend standard deviation
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    mesh = ax2.pcolormesh(daily_dict['longitude'][:], \
        daily_dict['latitude'][:], stdv_trends[month_idx,:,:], cmap = 'jet', \
        vmin = 0, vmax = vmax, transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = fig.colorbar(mesh, ax = ax2, pad = 0.03, fraction = 0.045, \
        label = 'Forcing Trend Std. Dev. [Wm$^{-2}$]', orientation = 'horizontal')
    ax2.set_boundary(circle, transform = ax2.transAxes)
    ax2.coastlines()
    ax2.set_extent([-180,180,65,90], ccrs.PlateCarree())

    #hasher = np.ma.masked_where( (np.abs(plot_trends[ii,:,:]) <= 1.0 * stdv_trends[ii,:,:]), plot_trends[ii,:,:])
    ax2.pcolor(daily_dict['longitude'][:], daily_dict['latitude'][:], \
        hasher, hatch = '....', alpha=0., shading = 'auto', \
        transform = ccrs.PlateCarree())

    plot_point_on_map(ax2, lat_val, lon_val, markersize = 10, color = 'red', \
        alpha = 1.0, add_border = True)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Plot histogram of forcing trends at this point
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #test_error_dist(plot_trends, 4, 50)
    num_bins = 20
    test_error_dist(daily_dict, forcing_trends, month_idx, lat_idx, lon_idx, num_bins, \
        ax = ax3, conf_level = conf_window, save = False)
    
    fig.tight_layout()
    plt.show()



def plot_bulk_force_AI_trend_v2(daily_dict, forcing_trends, OMI_daily_data, \
        vmax = 1.0, min_AI = None, max_AI = None, minlat = 65.5, \
        maxlat = 90.5,  conf_level = 90., sim_name = '', save = False):

    if(isinstance(OMI_daily_data, str)):
        if(min_AI is None):
            print("ERROR: If providing an OMI daily file, must provide min_AI")
            return

        print("Reading daily OMI data from", OMI_daily_data)
        print("Converting daily data to monthly omi_uvai_pert averages")
        OMI_daily_data = calcOMI_MonthAvg_FromDaily(OMI_daily_data, \
            min_AI = min_AI, max_AI = max_AI, minlat = minlat, maxlat = maxlat)

    str_conf_level = str(conf_level)
    if(conf_level > 1):
        conf_level = conf_level / 100.
    else:
        str_conf_level = str(conf_level * 100)

    mean_trends = np.nanmean(forcing_trends, axis = 0)
    stdv_trends = np.std(forcing_trends, axis = 0)
    stderr_trends = sem(forcing_trends, axis = 0), 
    #stderr_trend = stdv_trend / np.sqrt(in_arr.shape[0])
    conf_intvl = statnorm.interval(alpha = conf_level, \
        loc = mean_trends, scale = stdv_trends)

    #conf_intvl = statnorm.interval(alpha = conf_level, \
    #    loc = mean_trend, scale = stdv_trend, )
   
    plt.close('all') 
    fig = plt.figure(figsize = (7, 14), constrained_layout = False)
    #subfigs = fig.subfigures(nrows = 1, ncols = 3, wspace = 0.04)
    #
    #axs1 = [subfigs[0].add_subplot(6,1,ii, projection = ccrs.NorthPolarStereo()) for ii in range(1,7)]
    #axs2 = [subfigs[1].add_subplot(6,1,ii, projection = ccrs.NorthPolarStereo()) for ii in range(1,7)]
    #axs3 = [subfigs[2].add_subplot(6,1,ii, projection = ccrs.NorthPolarStereo()) for ii in range(1,7)]
    
    ax01 = fig.add_subplot(6,3,1, projection = ccrs.NorthPolarStereo())
    ax02 = fig.add_subplot(6,3,4, projection = ccrs.NorthPolarStereo())
    ax03 = fig.add_subplot(6,3,7, projection = ccrs.NorthPolarStereo())
    ax04 = fig.add_subplot(6,3,10, projection = ccrs.NorthPolarStereo())
    ax05 = fig.add_subplot(6,3,13, projection = ccrs.NorthPolarStereo())
    ax06 = fig.add_subplot(6,3,16, projection = ccrs.NorthPolarStereo())

    ax11 = fig.add_subplot(6,3,2, projection = ccrs.NorthPolarStereo())
    ax12 = fig.add_subplot(6,3,5, projection = ccrs.NorthPolarStereo())
    ax13 = fig.add_subplot(6,3,8, projection = ccrs.NorthPolarStereo())
    ax14 = fig.add_subplot(6,3,11, projection = ccrs.NorthPolarStereo())
    ax15 = fig.add_subplot(6,3,14, projection = ccrs.NorthPolarStereo())
    ax16 = fig.add_subplot(6,3,17, projection = ccrs.NorthPolarStereo())

    ax21 = fig.add_subplot(6,3,3, projection = ccrs.NorthPolarStereo())
    ax22 = fig.add_subplot(6,3,6, projection = ccrs.NorthPolarStereo())
    ax23 = fig.add_subplot(6,3,9, projection = ccrs.NorthPolarStereo())
    ax24 = fig.add_subplot(6,3,12, projection = ccrs.NorthPolarStereo())
    ax25 = fig.add_subplot(6,3,15, projection = ccrs.NorthPolarStereo())
    ax26 = fig.add_subplot(6,3,18, projection = ccrs.NorthPolarStereo())

    omi_axs = [ax01,ax02, ax03, ax04, ax05, ax06]
    mean_axs = [ax11,ax12, ax13, ax14, ax15, ax16]
    stdv_axs = [ax21,ax22, ax23, ax24, ax25, ax26]
    #total_list = [axs1, axs2, axs3]
  
    hasher = np.ma.masked_where( \
        (mean_trends == 0 ) |
        ( ((mean_trends > 0) & (conf_intvl[0] < 0)) | \
          ((mean_trends < 0) & (conf_intvl[1] > 0)) ), mean_trends) 


    titlers = ['Apr','May','June','July','Aug','Sep']
    for ii in range(6):

        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the OMI trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        grid_trends = np.full( (OMI_daily_data['AI'].shape[1],  \
            OMI_daily_data['AI'].shape[2]), np.nan)
        for jj in range(OMI_daily_data['AI'].shape[1]):
            for kk in range(OMI_daily_data['AI'].shape[2]):
                # Plot the individual runs
                # ------------------------
                local_sim_vals = OMI_daily_data['AI'][ii::6,jj,kk]
                x_vals = np.arange(local_sim_vals.shape[0])
    
                #if(trend_type=='standard'): 
                result = stats.linregress(x_vals, local_sim_vals[:])
                grid_trends[jj,kk] = result.slope * len(x_vals)
                #grid_trends[ii,jj,kk] = result.slope * len(x_vals)
    
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,work_mask.compressed())
                    #forcing_trends[i,j] = result.slope * len(x_vals)
                    #forcing_pvals[i,j]  = result.pvalue
                    #forcing_uncert[i,j] = result.stderr * len(x_vals)
                #else:
                #    res = stats.theilslopes(local_sim_vals[jj,:], x_vals, 0.90)
                #    trend_results[jj] = res[0] * len(x_vals)
                #    #forcing_trends[i,j] = res[0]*len(x_vals)
    
    
        omi_mesh = omi_axs[ii].pcolormesh(OMI_daily_data['LON'], OMI_daily_data['LAT'], \
            grid_trends[:,:], transform = ccrs.PlateCarree(), \
            #grid_trends[ii,:,:], transform = ccrs.PlateCarree(), \
            shading = 'auto', cmap = 'bwr', vmin = -0.5, vmax = 0.5)
        omi_axs[ii].set_boundary(circle, transform=omi_axs[ii].transAxes)
        omi_axs[ii].coastlines()
        omi_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        #omi_axs[ii].set_title(titlers[ii])

 
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the mean forcing trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   

        mean_mesh = mean_axs[ii].pcolormesh(daily_dict['longitude'][:], \
            daily_dict['latitude'][:], mean_trends[ii,:,:], \
            cmap = 'bwr', vmin = -vmax, vmax = vmax, \
            transform = ccrs.PlateCarree(), shading = 'auto')
        #cbar = fig.colorbar(mesh, ax = flat_axs[ii], pad = 0.03, fraction = 0.045)
        mean_axs[ii].set_boundary(circle, transform = mean_axs[ii].transAxes)
        mean_axs[ii].coastlines()
        mean_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        #omi_axs[ii].set_title(str(ii))

        pvar = 'mean'
        if(pvar == 'mean'):
            #window_min = mean_trends[ii,:,:] - (stdv_trends[ii,:,:] * 1.0)
            #window_max = mean_trends[ii,:,:] + (stdv_trends[ii,:,:] * 1.0)

            #hasher = np.ma.masked_where(\
            #    (mean_trends[ii,:,:] == 0) |
            #    ((np.abs(mean_trends[ii,:,:]) < (stdv_trends[ii,:,:] + conf_inte)

            #hasher = np.ma.masked_where(  \
            #     (mean_trends[ii,:,:] == 0) |
            #    ((mean_trends[ii,:,:] < 0) & (window_max[:,:] > 0)) | \
            #    ((mean_trends[ii,:,:] > 0) & (window_min[:,:] < 0)), mean_trends[ii,:,:])

            #hasher = np.ma.masked_where( (np.abs(plot_trends[ii,:,:]) <= 1.0 * stdv_trends[ii,:,:]), plot_trends[ii,:,:])
            mean_axs[ii].pcolor(daily_dict['longitude'][:], \
                daily_dict['latitude'][:], hasher[ii,:,:], hatch = '....', alpha=0., \
                shading = 'auto', transform = ccrs.PlateCarree())


        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the stdv forcing trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        stdv_mesh = stdv_axs[ii].pcolormesh(daily_dict['longitude'][:], \
            daily_dict['latitude'][:], stdv_trends[ii,:,:], \
            cmap = 'jet', vmin = 0, vmax = vmax, \
            transform = ccrs.PlateCarree(), shading = 'auto')
        #cbar = fig.colorbar(mesh, ax = flat_axs[ii], pad = 0.03, fraction = 0.045)
        stdv_axs[ii].set_boundary(circle, transform = stdv_axs[ii].transAxes)
        stdv_axs[ii].coastlines()
        stdv_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        #omi_axs[ii].set_title(str(ii))
 
    omi_axs[0].set_title('OMI UVAI Trend')
    mean_axs[0].set_title('ADRF Trend Mean')
    stdv_axs[0].set_title('ADRF Trend St. Dev.')


    fig.tight_layout(rect = [0.04,0.05,1,1.00])

    row_label_size = 10
    plotloc = omi_axs[0].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'April', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    #fig.text( ((plotloc.xmax + plotloc.xmin) / 2), plotloc.ymax + 0.02, \
    #    'OMI UVAI Trend', ha = 'center', va = 'center', \
    #    rotation = 'horizontal', weight = 'bold', fontsize = row_label_size + 1)
    plotloc = omi_axs[1].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'May', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    plotloc = omi_axs[2].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'June', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    plotloc = omi_axs[3].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'July', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    plotloc = omi_axs[4].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'August', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    plotloc = omi_axs[5].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'September', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)



    #fig.text(axs1[5].get_position().xmin, axs1[5].get_position().ymin, '-')
    #fig.text(axs1[5].get_position().xmin, axs1[5].get_position().ymax, '+')
    #fig.text(axs1[5].get_position().xmax, axs1[5].get_position().ymin, '0')
    #fig.text(axs1[5].get_position().xmax, axs1[5].get_position().ymax, 'l')

    #fig.text(axs2[5].get_position().xmin, axs2[5].get_position().ymin, '-')
    #fig.text(axs2[5].get_position().xmin, axs2[5].get_position().ymax, '+')
    #fig.text(axs2[5].get_position().xmax, axs2[5].get_position().ymin, '0')
    #fig.text(axs2[5].get_position().xmax, axs2[5].get_position().ymax, 'l')

    #fig.text(axs3[5].get_position().xmin, axs3[5].get_position().ymin, '-')
    #fig.text(axs3[5].get_position().xmin, axs3[5].get_position().ymax, '+')
    #fig.text(axs3[5].get_position().xmax, axs3[5].get_position().ymin, '0')
    #fig.text(axs3[5].get_position().xmax, axs3[5].get_position().ymax, 'l')

    lowplot1 = omi_axs[5].get_position() 
    diff_val = 0.01
    cbar_ax1 = fig.add_axes([lowplot1.x0 + diff_val, lowplot1.y0 - 0.02, \
        ((lowplot1.x1 - lowplot1.x0) - (diff_val * 2)), 0.01])
    cbar = fig.colorbar(omi_mesh, \
        cax = cbar_ax1, shrink = 0.8, orientation = 'horizontal', \
        label = 'AI Trend', extend = 'max')
    lowplot1 = mean_axs[5].get_position() 
    #cbar_ax2 = fig.add_axes([0.40, 0.04, 0.23, 0.01])
    cbar_ax2 = fig.add_axes([lowplot1.x0 + diff_val, lowplot1.y0 - 0.02, \
        ((lowplot1.x1 - lowplot1.x0) - (diff_val * 2)), 0.01])
    cbar = fig.colorbar(mean_mesh, \
        cax = cbar_ax2, shrink = 0.8, orientation = 'horizontal', \
        label = 'ADRF Trend [Wm$^{-2}$]', extend = 'both')
    lowplot1 = stdv_axs[5].get_position() 
    #cbar_ax3 = fig.add_axes([0.70, 0.04, 0.23, 0.01])
    cbar_ax3 = fig.add_axes([lowplot1.x0 + diff_val, lowplot1.y0 - 0.02, \
        ((lowplot1.x1 - lowplot1.x0) - (diff_val * 2)), 0.01])
    cbar = fig.colorbar(stdv_mesh, \
        cax = cbar_ax3, shrink = 0.8, orientation = 'horizontal', \
        label = 'ADRF Trend Std Dev [Wm$^{-2}$]')

    if(save):
        if(sim_name != ''):
            sim_name = '_' + sim_name
        outname = 'ai_force_trend_combined' + sim_name + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


def plot_bulk_force_AI_trend_v3(daily_dict, forcing_trends, OMI_daily_data, \
        NSIDC_data, MODIS_data, modis_var = 'cld_frac_mean', vmax = 1.0, \
        min_AI = None, max_AI = None, minlat = 65.5, \
        maxlat = 90.5,  sim_name = '', run_type = '', conf_level = 90., save = False):

    if(isinstance(OMI_daily_data, str)):
        if(min_AI is None):
            print("ERROR: If providing an OMI daily file, must provide min_AI")
            return

        print("Reading daily OMI data from", OMI_daily_data)
        print("Converting daily data to monthly omi_uvai_pert averages")
        OMI_daily_data = calcOMI_MonthAvg_FromDaily(OMI_daily_data, \
            min_AI = min_AI, max_AI = max_AI, minlat = minlat, maxlat = maxlat)

    str_conf_level = str(conf_level)
    if(conf_level > 1):
        conf_level = conf_level / 100.
    else:
        str_conf_level = str(conf_level * 100)

    mean_trends = np.nanmean(forcing_trends, axis = 0)
    stdv_trends = np.std(forcing_trends, axis = 0)
    stderr_trends = sem(forcing_trends, axis = 0), 
    #stderr_trend = stdv_trend / np.sqrt(in_arr.shape[0])
    conf_intvl = statnorm.interval(alpha = conf_level, \
        loc = mean_trends, scale = stdv_trends)

    #conf_intvl = statnorm.interval(alpha = conf_level, \
    #    loc = mean_trend, scale = stdv_trend, )
   
    plt.close('all') 
    fig = plt.figure(figsize = (10.5, 14), constrained_layout = False)
    #subfigs = fig.subfigures(nrows = 1, ncols = 3, wspace = 0.04)
    #
    #axs1 = [subfigs[0].add_subplot(6,1,ii, projection = ccrs.NorthPolarStereo()) for ii in range(1,7)]
    #axs2 = [subfigs[1].add_subplot(6,1,ii, projection = ccrs.NorthPolarStereo()) for ii in range(1,7)]
    #axs3 = [subfigs[2].add_subplot(6,1,ii, projection = ccrs.NorthPolarStereo()) for ii in range(1,7)]
   
    nrows = 6
    ncols = 5
    ax01 = fig.add_subplot(nrows,ncols, (1 - 1) * ncols + 1, projection = ccrs.NorthPolarStereo())
    ax02 = fig.add_subplot(nrows,ncols, (2 - 1) * ncols + 1, projection = ccrs.NorthPolarStereo())
    ax03 = fig.add_subplot(nrows,ncols, (3 - 1) * ncols + 1, projection = ccrs.NorthPolarStereo())
    ax04 = fig.add_subplot(nrows,ncols, (4 - 1) * ncols + 1, projection = ccrs.NorthPolarStereo())
    ax05 = fig.add_subplot(nrows,ncols, (5 - 1) * ncols + 1, projection = ccrs.NorthPolarStereo())
    ax06 = fig.add_subplot(nrows,ncols, (6 - 1) * ncols + 1, projection = ccrs.NorthPolarStereo())

    ax11 = fig.add_subplot(nrows,ncols, (1 - 1) * ncols + 2, projection = ccrs.NorthPolarStereo())
    ax12 = fig.add_subplot(nrows,ncols, (2 - 1) * ncols + 2, projection = ccrs.NorthPolarStereo())
    ax13 = fig.add_subplot(nrows,ncols, (3 - 1) * ncols + 2, projection = ccrs.NorthPolarStereo())
    ax14 = fig.add_subplot(nrows,ncols, (4 - 1) * ncols + 2, projection = ccrs.NorthPolarStereo())
    ax15 = fig.add_subplot(nrows,ncols, (5 - 1) * ncols + 2, projection = ccrs.NorthPolarStereo())
    ax16 = fig.add_subplot(nrows,ncols, (6 - 1) * ncols + 2, projection = ccrs.NorthPolarStereo())

    ax21 = fig.add_subplot(nrows,ncols, (1 - 1) * ncols + 3, projection = ccrs.NorthPolarStereo())
    ax22 = fig.add_subplot(nrows,ncols, (2 - 1) * ncols + 3, projection = ccrs.NorthPolarStereo())
    ax23 = fig.add_subplot(nrows,ncols, (3 - 1) * ncols + 3, projection = ccrs.NorthPolarStereo())
    ax24 = fig.add_subplot(nrows,ncols, (4 - 1) * ncols + 3, projection = ccrs.NorthPolarStereo())
    ax25 = fig.add_subplot(nrows,ncols, (5 - 1) * ncols + 3, projection = ccrs.NorthPolarStereo())
    ax26 = fig.add_subplot(nrows,ncols, (6 - 1) * ncols + 3, projection = ccrs.NorthPolarStereo())

    ax31 = fig.add_subplot(nrows,ncols, (1 - 1) * ncols + 4, projection = ccrs.NorthPolarStereo())
    ax32 = fig.add_subplot(nrows,ncols, (2 - 1) * ncols + 4, projection = ccrs.NorthPolarStereo())
    ax33 = fig.add_subplot(nrows,ncols, (3 - 1) * ncols + 4, projection = ccrs.NorthPolarStereo())
    ax34 = fig.add_subplot(nrows,ncols, (4 - 1) * ncols + 4, projection = ccrs.NorthPolarStereo())
    ax35 = fig.add_subplot(nrows,ncols, (5 - 1) * ncols + 4, projection = ccrs.NorthPolarStereo())
    ax36 = fig.add_subplot(nrows,ncols, (6 - 1) * ncols + 4, projection = ccrs.NorthPolarStereo())

    ax41 = fig.add_subplot(nrows,ncols, (1 - 1) * ncols + 5, projection = ccrs.NorthPolarStereo())
    ax42 = fig.add_subplot(nrows,ncols, (2 - 1) * ncols + 5, projection = ccrs.NorthPolarStereo())
    ax43 = fig.add_subplot(nrows,ncols, (3 - 1) * ncols + 5, projection = ccrs.NorthPolarStereo())
    ax44 = fig.add_subplot(nrows,ncols, (4 - 1) * ncols + 5, projection = ccrs.NorthPolarStereo())
    ax45 = fig.add_subplot(nrows,ncols, (5 - 1) * ncols + 5, projection = ccrs.NorthPolarStereo())
    ax46 = fig.add_subplot(nrows,ncols, (6 - 1) * ncols + 5, projection = ccrs.NorthPolarStereo())

    #omi_axs = [ax01,ax02, ax03, ax04, ax05, ax06]
    #mean_axs = [ax11,ax12, ax13, ax14, ax15, ax16]
    #stdv_axs = [ax21,ax22, ax23, ax24, ax25, ax26]
    #ice_axs =  [ax31,ax32, ax33, ax34, ax35, ax36]
    #cld_axs =  [ax41,ax42, ax43, ax44, ax45, ax46]

    ice_axs = [ax01,ax02, ax03, ax04, ax05, ax06]
    cld_axs = [ax11,ax12, ax13, ax14, ax15, ax16]
    omi_axs = [ax21,ax22, ax23, ax24, ax25, ax26]
    mean_axs =  [ax31,ax32, ax33, ax34, ax35, ax36]
    stdv_axs =  [ax41,ax42, ax43, ax44, ax45, ax46]


    #total_list = [axs1, axs2, axs3]
  
    hasher = np.ma.masked_where( \
        (mean_trends == 0 ) |
        ( ((mean_trends > 0) & (conf_intvl[0] < 0)) | \
          ((mean_trends < 0) & (conf_intvl[1] > 0)) ), mean_trends) 

    if(modis_var == 'cld_frac_mean'):
        cld_lims = 40.
    elif(modis_var == 'cod_mean'):
        cld_lims = 10.

    titlers = ['Apr','May','June','July','Aug','Sep']
    for ii in range(6):

        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the OMI trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        grid_trends = np.full( (OMI_daily_data['AI'].shape[1],  \
            OMI_daily_data['AI'].shape[2]), np.nan)
        for jj in range(OMI_daily_data['AI'].shape[1]):
            for kk in range(OMI_daily_data['AI'].shape[2]):
                # Plot the individual runs
                # ------------------------
                local_sim_vals = OMI_daily_data['AI'][ii::6,jj,kk]
                x_vals = np.arange(local_sim_vals.shape[0])
    
                #if(trend_type=='standard'): 
                result = stats.linregress(x_vals, local_sim_vals[:])
                grid_trends[jj,kk] = result.slope * len(x_vals)
                #grid_trends[ii,jj,kk] = result.slope * len(x_vals)
    
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,work_mask.compressed())
                    #forcing_trends[i,j] = result.slope * len(x_vals)
                    #forcing_pvals[i,j]  = result.pvalue
                    #forcing_uncert[i,j] = result.stderr * len(x_vals)
                #else:
                #    res = stats.theilslopes(local_sim_vals[jj,:], x_vals, 0.90)
                #    trend_results[jj] = res[0] * len(x_vals)
                #    #forcing_trends[i,j] = res[0]*len(x_vals)
    
    
        omi_mesh = omi_axs[ii].pcolormesh(OMI_daily_data['LON'], OMI_daily_data['LAT'], \
            grid_trends[:,:], transform = ccrs.PlateCarree(), \
            #grid_trends[ii,:,:], transform = ccrs.PlateCarree(), \
            shading = 'auto', cmap = 'bwr', vmin = -0.3, vmax = 0.3)
        omi_axs[ii].set_boundary(circle, transform=omi_axs[ii].transAxes)
        omi_axs[ii].coastlines()
        omi_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        #omi_axs[ii].set_title(titlers[ii])

 
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the mean forcing trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   

        mean_mesh = mean_axs[ii].pcolormesh(daily_dict['longitude'][:], \
            daily_dict['latitude'][:], mean_trends[ii,:,:], \
            cmap = 'bwr', vmin = -vmax, vmax = vmax, \
            transform = ccrs.PlateCarree(), shading = 'auto')
        #cbar = fig.colorbar(mesh, ax = flat_axs[ii], pad = 0.03, fraction = 0.045)
        mean_axs[ii].set_boundary(circle, transform = mean_axs[ii].transAxes)
        mean_axs[ii].coastlines()
        mean_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        #omi_axs[ii].set_title(str(ii))

        pvar = 'mean'
        if(pvar == 'mean'):
            #window_min = mean_trends[ii,:,:] - (stdv_trends[ii,:,:] * 1.0)
            #window_max = mean_trends[ii,:,:] + (stdv_trends[ii,:,:] * 1.0)

            #hasher = np.ma.masked_where(\
            #    (mean_trends[ii,:,:] == 0) |
            #    ((np.abs(mean_trends[ii,:,:]) < (stdv_trends[ii,:,:] + conf_inte)

            #hasher = np.ma.masked_where(  \
            #     (mean_trends[ii,:,:] == 0) |
            #    ((mean_trends[ii,:,:] < 0) & (window_max[:,:] > 0)) | \
            #    ((mean_trends[ii,:,:] > 0) & (window_min[:,:] < 0)), mean_trends[ii,:,:])

            #hasher = np.ma.masked_where( (np.abs(plot_trends[ii,:,:]) <= 1.0 * stdv_trends[ii,:,:]), plot_trends[ii,:,:])
            mean_axs[ii].pcolor(daily_dict['longitude'][:], \
                daily_dict['latitude'][:], hasher[ii,:,:], hatch = '....', alpha=0., \
                shading = 'auto', transform = ccrs.PlateCarree())


        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the stdv forcing trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        stdv_mesh = stdv_axs[ii].pcolormesh(daily_dict['longitude'][:], \
            daily_dict['latitude'][:], stdv_trends[ii,:,:], \
            cmap = 'jet', vmin = 0, vmax = vmax, \
            transform = ccrs.PlateCarree(), shading = 'auto')
        #cbar = fig.colorbar(mesh, ax = flat_axs[ii], pad = 0.03, fraction = 0.045)
        stdv_axs[ii].set_boundary(circle, transform = stdv_axs[ii].transAxes)
        stdv_axs[ii].coastlines()
        stdv_axs[ii].set_extent([-180,180,65,90], ccrs.PlateCarree())
        #omi_axs[ii].set_title(str(ii))


        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the NSIDC ice trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        ice_mesh = plotNSIDC_MonthTrend(NSIDC_data,month_idx = ii,save=False,\
            trend_type='linregress',season='',minlat=65.,return_trend=False, \
            colorbar = False, colorbar_label_size = None,title = '', \
            pax = ice_axs[ii], show_pval = False, uncert_ax = None, return_mesh = True)

        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
        #
        # Plot the MODIS cloud frac trends
        #
        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        if(modis_var == 'cod_mean'):
            if(ii == 5):
                plotMODIS_MYD08_MonthTrend(MODIS_data,modis_var,\
                    month_idx = ii,\
                    trend_type='linregress',season= '', minlat=65.5,\
                    colorbar = False, title = '', \
                    ax = cld_axs[ii], show_pval = False, uncert_ax = None, \
                    norm_to_decade = False, vmin = -cld_lims, vmax = cld_lims, \
                    return_trend = False, return_mesh = False)
            else:
                cld_mesh = plotMODIS_MYD08_MonthTrend(MODIS_data,modis_var,\
                    month_idx = ii,\
                    trend_type='linregress',season= '', minlat=65.5,\
                    colorbar = False, title = '', \
                    ax = cld_axs[ii], show_pval = False, uncert_ax = None, \
                    norm_to_decade = False, vmin = -cld_lims, vmax = cld_lims, \
                    return_trend = False, return_mesh = True)
        else:
            cld_mesh = plotMODIS_MYD08_MonthTrend(MODIS_data,modis_var,\
                month_idx = ii,\
                trend_type='linregress',season= '', minlat=65.5,\
                colorbar = False, title = '', \
                ax = cld_axs[ii], show_pval = False, uncert_ax = None, \
                norm_to_decade = False, vmin = -cld_lims, vmax = cld_lims, \
                return_trend = False, return_mesh = True)

 
    omi_axs[0].set_title('OMI UVAI Trend')
    mean_axs[0].set_title('ADRF Trend Mean')
    stdv_axs[0].set_title('ADRF Trend St. Dev.')
    ice_axs[0].set_title('SSMIS\nIce Conc. Trend')
    if(modis_var == 'cld_frac_mean'):
        cld_axs[0].set_title('MODIS\nCld. Frac. Trend')
    else:
        cld_axs[0].set_title('MODIS COD Trend')


    fig.tight_layout(rect = [0.04,0.05,1,1.00])

    row_label_size = 10
    plotloc = ice_axs[0].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'April', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    #fig.text( ((plotloc.xmax + plotloc.xmin) / 2), plotloc.ymax + 0.02, \
    #    'OMI UVAI Trend', ha = 'center', va = 'center', \
    #    rotation = 'horizontal', weight = 'bold', fontsize = row_label_size + 1)
    plotloc = ice_axs[1].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'May', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    plotloc = ice_axs[2].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'June', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    plotloc = ice_axs[3].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'July', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    plotloc = ice_axs[4].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'August', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    plotloc = ice_axs[5].get_position() 
    fig.text(plotloc.xmin - 0.03, (plotloc.ymax + plotloc.ymin) / 2., \
        'September', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)



    #fig.text(axs1[5].get_position().xmin, axs1[5].get_position().ymin, '-')
    #fig.text(axs1[5].get_position().xmin, axs1[5].get_position().ymax, '+')
    #fig.text(axs1[5].get_position().xmax, axs1[5].get_position().ymin, '0')
    #fig.text(axs1[5].get_position().xmax, axs1[5].get_position().ymax, 'l')

    #fig.text(axs2[5].get_position().xmin, axs2[5].get_position().ymin, '-')
    #fig.text(axs2[5].get_position().xmin, axs2[5].get_position().ymax, '+')
    #fig.text(axs2[5].get_position().xmax, axs2[5].get_position().ymin, '0')
    #fig.text(axs2[5].get_position().xmax, axs2[5].get_position().ymax, 'l')

    #fig.text(axs3[5].get_position().xmin, axs3[5].get_position().ymin, '-')
    #fig.text(axs3[5].get_position().xmin, axs3[5].get_position().ymax, '+')
    #fig.text(axs3[5].get_position().xmax, axs3[5].get_position().ymin, '0')
    #fig.text(axs3[5].get_position().xmax, axs3[5].get_position().ymax, 'l')

    lowplot1 = omi_axs[5].get_position() 
    diff_val = 0.01
    cbar_ax1 = fig.add_axes([lowplot1.x0 + diff_val, lowplot1.y0 - 0.02, \
        ((lowplot1.x1 - lowplot1.x0) - (diff_val * 2)), 0.01])
    cbar = fig.colorbar(omi_mesh, \
        cax = cbar_ax1, shrink = 0.8, orientation = 'horizontal', \
        label = 'UVAI Trend\n[UVAI (study period) $^{-1}$]', extend = 'max')

    lowplot1 = mean_axs[5].get_position() 
    #cbar_ax2 = fig.add_axes([0.40, 0.04, 0.23, 0.01])
    cbar_ax2 = fig.add_axes([lowplot1.x0 + diff_val, lowplot1.y0 - 0.02, \
        ((lowplot1.x1 - lowplot1.x0) - (diff_val * 2)), 0.01])
    cbar = fig.colorbar(mean_mesh, \
        cax = cbar_ax2, shrink = 0.8, orientation = 'horizontal', \
        label = 'ADRF Trend\n[Wm$^{-2}$ (study period) $^{-1}$]', extend = 'both')

    lowplot1 = stdv_axs[5].get_position() 
    #cbar_ax3 = fig.add_axes([0.70, 0.04, 0.23, 0.01])
    cbar_ax3 = fig.add_axes([lowplot1.x0 + diff_val, lowplot1.y0 - 0.02, \
        ((lowplot1.x1 - lowplot1.x0) - (diff_val * 2)), 0.01])
    cbar = fig.colorbar(stdv_mesh, \
        cax = cbar_ax3, shrink = 0.8, orientation = 'horizontal', \
        label = 'ADRF Trend Std Dev\n[Wm$^{-2}$]')

    lowplot1 = ice_axs[5].get_position() 
    cbar_ax4 = fig.add_axes([lowplot1.x0 + diff_val, lowplot1.y0 - 0.02, \
        ((lowplot1.x1 - lowplot1.x0) - (diff_val * 2)), 0.01])
    cbar = fig.colorbar(ice_mesh, \
        cax = cbar_ax4, shrink = 0.8, orientation = 'horizontal', \
        label = 'Ice Conc. Trend\n[% (study period)$^{-1}$]')

    if(modis_var == 'cld_frac_mean'):
        label_text = ('Cld. Frac. Trend\n[% (study period$^{-1}$]')
    else:
        label_text = ('COD Trend\n[COD (study period)$^{-1}$]')
    lowplot1 = cld_axs[5].get_position() 
    cbar_ax5 = fig.add_axes([lowplot1.x0 + diff_val, lowplot1.y0 - 0.02, \
        ((lowplot1.x1 - lowplot1.x0) - (diff_val * 2)), 0.01])
    cbar = fig.colorbar(cld_mesh, \
        cax = cbar_ax5, shrink = 0.8, orientation = 'horizontal', \
        label = label_text)

    if(save):
        if(modis_var == 'cld_frac_mean'):
            modis_adder = '_cldfrac'
        elif(modis_var == 'cod_mean'):
            modis_adder = '_cod'
        if(sim_name != ''):
            sim_name = '_' + sim_name
        if(run_type != ''):
            run_type = '_' + run_type

        outname = 'ai_force_trend_ice' + modis_adder + '_combined' + sim_name + run_type + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()



def plot_arctic_avg_region_trends(sim_values, reg_idx):
    plt.close('all')
    fig = plt.figure(figsize = (10, 6))
    axs = fig.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True)
    flat_axs = axs.flatten()
  
    x_vals = np.arange(2005, 2021) 
     
    for ii in range(6):
        flat_axs[ii].plot(x_vals, sim_values[:,reg_idx,ii::6].T)
        flat_axs[ii].axhline(0, color = 'k', linestyle = ':')

        #result = stats.linregress(x_vals, work_mask.compressed())
        ##slope, intercept, r_value, p_value, std_err = \
        ##    stats.linregress(x_vals,work_mask.compressed())
        #forcing_trends[i,j] = result.slope * len(x_vals)
        
    fig.tight_layout()
    plt.show() 


def calc_arctic_avg_region_trends(sim_values):

    #trend_vals = np.full((sim_values.shape[0], 3, 6), np.nan)
    #trend_pval = np.full((sim_values.shape[0], 3, 6), np.nan)
   
    xvals = np.arange(2005, 2021)
 
    ## FOR USE WITH REGION-AVERAGED FORCING VALUES
    #for kk in range(3):
    #    print('region: ', kk)
    #    for jj in range(6):
    #
    #        # Trend analysis
    #        # --------------
    #        trend_vals = np.full( sim_values.shape[0], np.nan)
    #        for ii in range(trend_vals.shape[0]):
    #            result = stats.linregress(xvals, sim_values[ii,kk,jj::6])
    #            trend_vals[ii] = result.slope * len(xvals)
    #            #trend_pval[ii,kk,jj] = result.pvalue
    #        print(jj, np.round(np.mean(trend_vals[:]), 3), np.round(np.std(trend_vals[:]), 3))
    #
    #        # Mean analysis
    #        # -------------
    #        #mean_force = np.nanmean(sim_values[:,kk,jj::6])
    #        #std_force1 = np.nanstd(sim_values[:,kk,jj::6])
    #        #std_force2 = np.sqrt(np.sum(np.nanstd(sim_values[:,kk,jj::6], axis = 1)**2.) / \
    #        #    np.nanstd(sim_values[:,kk,jj::6], axis = 1).shape[0])
    #        #print(jj, np.round(mean_force, 3), np.round(std_force1, 3), np.round(std_force2, 3))

    trend_vals_all = np.full( (3, 6, sim_values.shape[0]), np.nan)
    trend_pval_all = np.full( (3, 6, sim_values.shape[0]), np.nan)
    # FOR USE WITH REGION-AVERAGED FORCING VALUES
    fmt_str = '    monthidx = {0:2d}, mean = {1:6.3f}, std = {2:6.3f}, ' + \
        'minpval = {3:6.3f}, meanpval = {4:6.3f}, maxpval = {5:6.3f}'
    for kk in range(3):
        print('region: ', kk)
        for jj in range(6):
    
            # Trend analysis
            # --------------
            #trend_vals = np.full( sim_values.shape[0], np.nan)
            for ii in range(trend_vals_all.shape[2]):
                result = stats.linregress(xvals, sim_values[ii,kk,jj::6])
                #trend_vals[ii] = result.slope * len(xvals)
                trend_vals_all[kk,jj,ii] = result.slope * len(xvals)
                trend_pval_all[kk,jj,ii] = result.pvalue
            print(fmt_str.format(jj, \
                np.round(np.mean(trend_vals_all[kk,jj,:]), 3), \
                np.round(np.std(trend_vals_all[kk,jj,:]), 3), \
                np.round(np.min(trend_pval_all[kk,jj,:]), 3), \
                np.round(np.mean(trend_pval_all[kk,jj,:]), 3), \
                np.round(np.max(trend_pval_all[kk,jj,:]), 3)))


    return trend_vals_all, trend_pval_all


def plot_reficecld_comps(daily_dict, refcld_dict, refice_dict, reg_idx, \
        trend_type = 'linregress', minlat = 65.5, maxlat = 90.5, save = False):

    daily_month_vals = calc_monthly_force_from_daily(daily_dict, minlat = 65.5, maxlat = 90.5, \
        return_std = False)
    
    refcld_month_vals = calc_monthly_force_from_daily(refcld_dict, minlat = 65.5, maxlat = 90.5, \
        return_std = False)
    
    refice_month_vals = calc_monthly_force_from_daily(refice_dict, minlat = 65.5, maxlat = 90.5, \
        return_std = False)
    
    daily_regions  = np.full( (3, daily_month_vals.shape[0]), np.nan)
    refcld_regions = np.full( (3, daily_month_vals.shape[0]), np.nan)
    refice_regions = np.full( (3, daily_month_vals.shape[0]), np.nan)
    
    lat_mins = [65.5, 65.5, 75.5]
    lat_maxs = [89.5, 75.5, 89.5]
    for jj in range(len(lat_mins)):
        lat_idxs = np.where( (daily_dict['latitude'][:] >= lat_mins[jj]) & (daily_dict['latitude'][:] < lat_maxs[jj])) 
        lat_beg_idx = lat_idxs[0][0]
        lat_end_idx = lat_idxs[0][-1] + 1
        daily_regions[jj,:] = \
            np.nanmean(daily_month_vals[:,lat_beg_idx:lat_end_idx,:], axis = (1,2))
        refcld_regions[jj,:] = \
            np.nanmean(refcld_month_vals[:,lat_beg_idx:lat_end_idx,:], axis = (1,2))
        refice_regions[jj,:] = \
            np.nanmean(refice_month_vals[:,lat_beg_idx:lat_end_idx,:], axis = (1,2))
    
    fig = plt.figure(figsize = (8, 5))
    axs = fig.subplots(2,3, sharex = True, sharey = True)
    flat_axs = axs.flatten()

    x_vals = np.arange(2005, 2021)   

    def plot_ref_time_series(x_data, y_data, pax, trend_type, label):
        
        flat_axs[ii].plot(x_data, y_data, label = label)
        if((trend_type == 'standard') | (trend_type == 'linregress')): 
            result = stats.linregress(x_vals, y_data)
            forcing_trend = result.slope * len(x_vals)
            forcing_pval  = result.pvalue
            forcing_uncert = result.stderr * len(x_vals)
        else:
            res = stats.theilslopes(y_data, x_vals, 0.90)
            forcing_trend = res[0]*len(x_vals)


        return forcing_trend    
 
    for ii in range(len(flat_axs)):

        # Do the daily trends first
        # -------------------------
        daily_trend = plot_ref_time_series(x_vals, \
            daily_regions[reg_idx,ii::6], flat_axs[ii], trend_type, 'Control')

        # Now, do ref cld trends
        # ----------------------
        refcld_trend = plot_ref_time_series(x_vals, \
            refcld_regions[reg_idx,ii::6], flat_axs[ii], trend_type, 'Cld2005')

        # Now, do ref ice trends
        # ----------------------
        refice_trend = plot_ref_time_series(x_vals, \
            refice_regions[reg_idx,ii::6], flat_axs[ii], trend_type, 'Ice2005')

        flat_axs[ii].axhline(0, color = 'k', linestyle = ':')

        print("DAILY:", np.round(daily_trend, 4), \
              "CLD:", np.round(refcld_trend, 4), \
              "ICE:", np.round(refice_trend, 4))

    flat_axs[5].legend()

    fig.tight_layout()       
    plt.show() 

def plot_reficecld_comps_many_allregions(daily_filename, comp_type, \
        file_start = 'arctic_daily_est_forcing_numsfcbins6', \
        trend_type = 'linregress', minlat = 65.5, maxlat = 90.5, \
        return_trends = False, save = False):

    daily_dict = read_daily_month_force_L2L3_error_from_HDF5(daily_filename)

    # Find all the ref*** simulation files
    # ------------------------------------
    ref_files = glob(file_start + '_ref' + comp_type + '*.hdf5')
   
    num_ref_sims = len(ref_files)
    
    # Calculate the values for the control run
    # ----------------------------------------
    daily_month_vals = calc_monthly_force_from_daily(daily_dict, minlat = 65.5, maxlat = 90.5, \
        return_std = False)
    
    daily_vals = np.full( (3, daily_month_vals.shape[0]), np.nan)
    ref_vals = np.full( (num_ref_sims, 3, daily_month_vals.shape[0]), np.nan)
    
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
    
        ref_dict = read_daily_month_force_L2L3_error_from_HDF5(ref_files[ii])
    
        ref_month_vals = calc_monthly_force_from_daily(ref_dict, minlat = 65.5, maxlat = 90.5, \
            return_std = False)
        
        for jj in range(len(lat_mins)):
            lat_idxs = np.where( (daily_dict['latitude'][:] >= lat_mins[jj]) & (daily_dict['latitude'][:] < lat_maxs[jj])) 
            lat_beg_idx = lat_idxs[0][0]
            lat_end_idx = lat_idxs[0][-1] + 1
            #daily_regions[jj,:] = \
            #    np.nanmean(daily_month_vals[:,lat_beg_idx:lat_end_idx,:], axis = (1,2))
            ref_vals[ii,jj,:] = \
                np.nanmean(ref_month_vals[:,lat_beg_idx:lat_end_idx,:], axis = (1,2))
   
    # Plot the results
    # ----------------
    fig = plt.figure(figsize = (7, 11))
    axs = fig.subplots(6,3, sharex = True, sharey = True)
    
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
    
    ref_trends = np.full( (6, 3, num_ref_sims), np.nan)
    daily_trends = np.full( (6, 3), np.nan)
  
    # Month loop  
    for ii in range(6):
 
        # Region loop 
        for kk in range(3):

            print("Region idx = ", kk)
 
            # Sim loop 
            for jj in range(ref_vals.shape[0]):
    
                # Now, do ref ice trends
                # ----------------------
                ref_trends[ii,kk,jj] = plot_ref_time_series(x_vals, \
                    ref_vals[jj,kk,ii::6], axs[ii,kk], trend_type, 'Ice2005')
    
            # Do the daily trends first
            # -------------------------
            daily_trends[ii,kk] = plot_ref_time_series(x_vals, \
                daily_vals[kk,ii::6], axs[ii,kk], trend_type, 'Control', color = 'k')
    
            axs[ii,kk].axhline(0, color = 'k', linestyle = ':')
    
            print("DAILY:", np.round(daily_trends[ii,kk], 4), \
                  "REF:",   np.round(np.mean(ref_trends[ii,kk,:]), 4))
    
    #flat_axs[5].legend()
    
    fig.tight_layout()       
    plt.show() 

    if(return_trends):
        return daily_trends, ref_trends

