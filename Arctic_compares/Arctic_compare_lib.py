"""
  NAME:

  PURPOSE:


"""
import sys
import json
sys.path.append('/home/bsorenson')
from python_lib import circle, plot_trend_line, nearest_gridpoint, \
    aerosol_event_dict, init_proj, plot_lat_circles, plot_figure_text, \
    plot_subplot_label
sys.path.append('/home/bsorenson/Research/OMI')
from OMILib import *
sys.path.append('/home/bsorenson/Research/CERES')
from gridCERESLib import *
sys.path.append('/home/bsorenson/Research/MODIS/obs_smoke_forcing')
from MODISLib import *
sys.path.append('/home/bsorenson/Research/NSIDC')
from NSIDCLib import *
sys.path.append('/home/bsorenson/Research/FuLiou')
from FuLiouLib import *

data_dir = '/home/bsorenson/Research/Arctic_compares/comp_data/'
datacrs = ccrs.PlateCarree()
mapcrs = ccrs.NorthPolarStereo()

json_time_database = '/home/bsorenson/Research/Arctic_compares/json_comp_times.txt'
json_file_database = '/home/bsorenson/Research/Arctic_compares/json_comp_files.txt'

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Automation
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# This function automates everything for a single OMI date string:
# the data downloading/matching, the first-look comparison-making, the 
# HDF5 file-making, and the data copying to Raindrop.
def automate_all_preprocess(date_str, download = True, images = True, \
        process = True):
    if(isinstance(date_str, str)):
        date_str = [date_str]

    if(download):
        # If desired, download all data and/or match up the files in the JSON
        # -------------------------------------------------------------------
        auto_all_download(date_str, download = download, rewrite_json = True)

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
            plot_compare_OMI_CERES_MODIS_NSIDC(dstr, 7, \
                omi_dtype = 'shawn', minlat = 65., zoom = True, save = True)

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
 
            # Run the data subsetter codes
            # ----------------------------
            write_shawn_to_HDF5(dstr, save_path = save_dir2, minlat = 65., \
                shawn_path = '/home/bsorenson/data/OMI/shawn_files/')

            write_CERES_hrly_grid_to_HDF5(file_date_dict[dstr]['CERES'], \
                save_path = save_dir2)

            MODIS_date = file_date_dict[dstr]['MODIS'][0]
            write_MODIS_to_HDF5(MODIS_date, channel = 2, swath = True, \
                save_path = save_dir2)
            write_MODIS_to_HDF5(MODIS_date, channel = 7, swath = True, \
                save_path = save_dir2)

            NSIDC_date = file_date_dict[dstr]['NSIDC'][:8]
            print(NSIDC_date)
            writeNSIDC_to_HDF5(NSIDC_date, save_path = save_dir2)

            # Finally, gzip the data
            # ---------------------
            os.chdir(data_dir)
            cmnd = dt_date_str.strftime('tar -cvzf ' + \
                'combined_subsets_%Y%m%d%H%M.tar.gz ' + short_save_dir2)
            print(cmnd)
            os.system(cmnd)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Downloading functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def auto_all_download(date_str, download = True, rewrite_json = False):
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

    for dstr in date_str:
        print(dstr)
        if(dstr not in out_time_dict.keys() or rewrite_json):
            OMI_base = readOMI_swath_shawn(dstr, latmin = 65., \
                shawn_path = '/home/bsorenson/data/OMI/shawn_files/')

            omi_file = subprocess.check_output('ls ' + \
                '/home/bsorenson/data/OMI/H5_files/' + \
                'OMI-*_' + OMI_base['date'][:4] + 'm' + OMI_base['date'][4:8] + 't' + \
                    OMI_base['date'][8:] + '*.he5', shell = True).decode(\
                'utf-8').strip().split('\n')[0]

            omi_shawn = '/home/bsorenson/data/OMI/shawn_files/' + OMI_base['date']

            CERES_date_str = np.min(OMI_base['TIME'][~OMI_base['UVAI_raw'].mask]).strftime('%Y%m%d%H')

            # Look for the CERES file
            # -----------------------
            try:
                ceres_file = subprocess.check_output('ls ' + \
                    '/home/bsorenson/data/CERES/SSF_Level2/Aqua/' + \
                    'CERES_SSF_Aqua-XTRK_Edition4A_Subset_' + \
                    CERES_date_str[:8] + '*', shell = True).decode(\
                    'utf-8').strip().split('\n')[0]
                print(ceres_file)
            except subprocess.CalledProcessError:
                print("CERES file not found. Must download")
                print(" Cycling the loop.")
                continue

            modis_date_list, file_list = download_MODIS_swath(CERES_date_str, \
                    dest_dir = '/home/bsorenson/data/MODIS/Aqua/', download = download)
   
            if(download): 
                download_NSIDC_daily(CERES_date_str[:8])

            nsidc_file = \
                '/home/bsorenson/data/NSIDC/NSIDC0051_SEAICE_PS_N25km_' + \
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

                out_file_dict[local_date][local_time]['omi']   = omi_file
                out_file_dict[local_date][local_time]['ceres'] = ceres_file
                out_file_dict[local_date][local_time]['modis'] = tfile
                out_file_dict[local_date][local_time]['ceres_time'] = CERES_date_str[8:10]
                out_file_dict[local_date][local_time]['swath']      = modis_date_list
                if('Lat' not in out_file_dict[local_date][local_time].keys()):
                    out_file_dict[local_date][local_time]['Lat'] = [60., 90.]
                if('Lon' not in out_file_dict[local_date][local_time].keys()):
                    out_file_dict[local_date][local_time]['Lon'] = [-180., 180.]
                if('modis_Lat' not in out_file_dict[local_date][local_time].keys()):
                    out_file_dict[local_date][local_time]['modis_Lat'] = [60., 90.]
                if('modis_Lon' not in out_file_dict[local_date][local_time].keys()):
                    out_file_dict[local_date][local_time]['modis_Lon'] = [-180., 180.]
    

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

def read_colocated(date_str, zoom = True):
   
    filename =  data_dir + 'colocated_subset_' + date_str + '.hdf5'
    print(filename)
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
    coloc_data['OMI'] = np.ma.masked_invalid(data['omi_uvai_pert'][:,:])
    coloc_data['MODIS_CH2'] = np.ma.masked_where((\
        data['modis_ch2'][:,:] == -999.) | (data['modis_ch2'][:,:] > 1.0), \
        data['modis_ch2'][:,:])
    coloc_data['MODIS_CH7'] = np.ma.masked_where((\
        data['modis_ch7'][:,:] == -999.) | (data['modis_ch7'][:,:] > 1.0), \
        data['modis_ch7'])
    coloc_data['NSIDC_ICE'] = np.ma.masked_where((\
        #data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.), \
        data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 100.) | \
        (data['nsidc_ice'][:,:] < 80.), \
        data['nsidc_ice'])
    coloc_data['NSIDC_OCEAN'] = np.ma.masked_where((\
        data['nsidc_ice'][:,:] == -999.) | (data['nsidc_ice'][:,:] > 0.), \
        data['nsidc_ice'])
    coloc_data['NSIDC_LAND'] = np.ma.masked_where((\
        data['nsidc_ice'][:,:] == -999.) |  (data['nsidc_ice'][:,:] != 254.), \
        data['nsidc_ice'])
    coloc_data['CERES_SWF'] = np.ma.masked_where((\
        data['ceres_swf'][:,:] == -999.) | (data['ceres_swf'][:,:] > 5000.), \
        data['ceres_swf'])
    coloc_data['CERES_LWF'] = np.ma.masked_where((\
        data['ceres_lwf'][:,:] == -999.) | (data['ceres_lwf'][:,:] > 5000.), \
        data['ceres_lwf'])
    coloc_data['LAT'] = data['omi_lat'][:,:]
    coloc_data['LON'] = data['omi_lon'][:,:]
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
            total_data = coloc_data
        else:
            for key in coloc_data.keys():
                if(key != 'date_str'):
                    total_data[key] = np.concatenate((total_data[key], coloc_data[key]))
        #test_data = np.concatenate((test_data, temp_data['omi_uvai_pert'][:,:]))
        #print(tfile, temp_data['omi_uvai_pert'].shape) 
        #temp_data.close()

    total_data['date_str'] = date_str

    total_data['OMI'] = np.ma.masked_invalid(total_data['OMI'][:,:])
    total_data['MODIS_CH2'] = np.ma.masked_where((\
        total_data['MODIS_CH2'][:,:] == -999.) | (total_data['MODIS_CH2'][:,:] > 1.0), \
        total_data['MODIS_CH2'][:,:])
    total_data['MODIS_CH7'] = np.ma.masked_where((\
        total_data['MODIS_CH7'][:,:] == -999.) | (total_data['MODIS_CH7'][:,:] > 1.0), \
        total_data['MODIS_CH7'])
    total_data['NSIDC_ICE'] = np.ma.masked_where((\
        #data['NSIDC_ICE'][:,:] == -999.) | (data['NSIDC_ICE'][:,:] > 100.), \
        total_data['NSIDC_ICE'][:,:] == -999.) | (total_data['NSIDC_ICE'][:,:] > 100.) | \
        (total_data['NSIDC_ICE'][:,:] < 80.), \
        total_data['NSIDC_ICE'])
    total_data['NSIDC_OCEAN'] = np.ma.masked_where((\
        total_data['NSIDC_ICE'][:,:] == -999.) | (total_data['NSIDC_ICE'][:,:] > 0.), \
        total_data['NSIDC_ICE'])
    total_data['NSIDC_LAND'] = np.ma.masked_where((\
        total_data['NSIDC_ICE'][:,:] == -999.) |  (total_data['NSIDC_ICE'][:,:] != 254.), \
        total_data['NSIDC_ICE'])
    total_data['CERES_SWF'] = np.ma.masked_where((\
        total_data['CERES_SWF'][:,:] == -999.) | (total_data['CERES_SWF'][:,:] > 5000.), \
        total_data['CERES_SWF'])
    total_data['CERES_LWF'] = np.ma.masked_where((\
        total_data['CERES_LWF'][:,:] == -999.) | (total_data['CERES_LWF'][:,:] > 5000.), \
        total_data['CERES_LWF'])

    return total_data
 
    #print(test_data.shape)    
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

    if('/home/bsorenson/Research/OMI/' not in sys.path):
        sys.path.append('/home/bsorenson/Research/OMI/')
    if('/home/bsorenson/Research/CERES/' not in sys.path):
        sys.path.append('/home/bsorenson/Research/CERES/')
    if('/home/bsorenson/Research/MODIS/obs_smoke_forcing/' not in sys.path):
        sys.path.append('/home/bsorenson/Research/MODIS/obs_smoke_forcing/')
    if('/home/bsorenson/Research/NSIDC/' not in sys.path):
        sys.path.append('/home/bsorenson/Research/NSIDC/')

    from OMILib import plotOMI_single_swath_figure
    from gridCERESLib import plotCERES_hrly_figure
    from MODISLib import plot_MODIS_channel
    from NSIDCLib import plotNSIDC_daily_figure

    

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
            zoom = zoom)
    
    # Plot the CERES data
    # -------------------
    plotCERES_hrly_figure(ceres_date, 'SWF',  \
        minlat = minlat, lat_circles = None, ax = ax5, title = 'SWF',\
        grid_data = True, zoom = zoom, vmax = 450, vmin = None, save = False)
    plotCERES_hrly_figure(ceres_date, 'ALB',  \
        minlat = minlat, lat_circles = None, ax = ax6, title = 'ALB',\
        grid_data = True, zoom = zoom, vmax = None, vmin = None, save = False)

    if(save):
        outname = 'omi_ceres_modis_nsidc_compare_' + omi_date + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_spatial(ax, lon, lat, data, date_str, cmap = 'Greys_r', \
        plabel = '', colorbar = True, zoom = False):
   
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


# cat is either "LAND", "ICE", "OCEAN", or "ALL"
def plot_compare_scatter_category(coloc_data, var1, var2, var3 = None, \
        cat = "ALL", minlat = 65., xmin = None, xmax = None, ymin = None, \
        ymax = None, ax = None, colorbar = True, trend = False, zoom = False, \
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
    
    if(var3 is None):
        ax.scatter(mask_xdata, mask_ydata, s = 6, color = pcolor)
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

    if(trend):
        print(cat)
        plot_trend_line(ax, mask_xdata.compressed(), mask_ydata.compressed(), color=trend_color, linestyle = '-', \
            slope = 'thiel-sen')

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

def plot_compare_combined_category(coloc_data, var1 = 'OMI', \
        var2 = 'CERES_SWF', var3 = None, cat = "ALL", minlat = 65., \
        xmin = 1.0, xmax = None, ymin = None, ymax = None, ax = None, \
        colorbar = True, trend = False, zoom = False, color = None, \
        save = False):

    if(isinstance(coloc_data, str)):
        dt_date_str = datetime.strptime(coloc_data, '%Y%m%d%H%M')
        coloc_data = read_colocated(coloc_data)
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
    plot_spatial(ax4, coloc_data['LON'], coloc_data['LAT'], coloc_data['OMI'], \
        pdate, cmap = 'jet', zoom = zoom)
    ax4.set_title('OMI AI Perturb')
    
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

    fig1.tight_layout()

    if(save):
        date_str = dt_date_str.strftime('%Y%m%d%H%M')
        outname = 'arctic_compare_combined_category_nomin_' + date_str + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()
