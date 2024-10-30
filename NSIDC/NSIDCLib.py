"""
  NAME:
    NSIDCLib.py

  PURPOSE:
    House all the functions used for reading and working with the ice data



"""
import numpy as np
import numpy.ma as ma
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import glob
import subprocess
import os
from scipy import stats
from datetime import datetime, timedelta
from netCDF4 import Dataset
import h5py
from dateutil.relativedelta import relativedelta

home_dir = os.environ['HOME']

sys.path.append(home_dir)
from python_lib import circle, plot_trend_line, nearest_gridpoint, \
    aerosol_event_dict, init_proj, plot_lat_circles, plot_figure_text, \
    plot_subplot_label, lat_lon_area, plot_point_on_map

data_dir = home_dir + '/data/NSIDC/'
datacrs = ccrs.PlateCarree()
mapcrs = ccrs.NorthPolarStereo()

rval_dict = {
    'Unchanged ocean':       0,   # 0
    'Unchanged mix':         1,   # 1
    'Unchanged ice':         2,   # 2 
    'Still primarily ocean': 3,   # 0
    'Mix to ocean':          4,   # 3
    'Ice to ocean':          5,   # 3
    'Still primarily mix':   6,   # 1
    'Ocean to mix':          7,   # 4
    'Ice to mix':            8,   # 4
    'Still primarily ice':   9,   # 2
    'Ocean to ice':         10,   # 5
    'Mix to ice':           11,   #
    'Pole Hole':            12,
    'Unused':               13,
    'Coastline':            14,
    'Land':                 15,
    'Other':                -1,
}
   



def plot_NSIDC_data(NSIDC_data,index):
    data = NSIDC_data['data'][index,:,:]
    plt.pcolormesh(data.T)
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.title(NSIDC_data['titles'][index])
    plt.show()

def trend_calc(NSIDC_data,x_ind,y_ind,thielsen=False):
    temp_data = NSIDC_data['data'][:,x_ind,y_ind]
    avg_ice = np.average(temp_data)
    interpx = np.arange(len(temp_data))
    #interpx = years
    ##interper = np.poly1d(np.polyfit(interpx,temp_data,1)) 
    ### Normalize trend by dividing by number of years
    ##trend = (interper(interpx[-1])-interper(interpx[0]))

    slope, intercept, r_value, p_value, std_err = stats.linregress(interpx,temp_data)
    ##print(slope/len(test_dict[dictkey].keys())
    regress_y = interpx*slope+intercept
    trend = regress_y[-1] - regress_y[0]

    if(thielsen==True):
        #The slope
        S=0
        sm=0
        nx = len(temp_data)
        num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
        Sn=np.zeros(num_d)
        for si in range(0,nx-1):
            for sj in range(si+1,nx):
                # Find the slope between the two points
                Sn[sm] = (temp_data[si]-temp_data[sj])/(si-sj) 
                sm=sm+1
            # Endfor
        # Endfor
        Snsorted=sorted(Sn)
        sm=int(num_d/2.)
        if(2*sm    == num_d):
            tsslope=0.5*(Snsorted[sm]+Snsorted[sm+1])
        if(2*sm+1 == num_d): 
            tsslope=Snsorted[sm+1]
        regress_ts = interpx*tsslope+intercept
        trend = regress_ts[-1]-regress_ts[0]


    pcnt_change = (trend/avg_ice)*100.
    # Find the percent change per decade
    return trend,pcnt_change

def read_ice(season,pre2001=False):
    spring=False
    summer=False
    autumn=False
    winter=False
    sunlight=False
    season_adder = ' '+season
    file_adder = '_'+season
    if(season=='spring'):
        spring=True
        ls_check = ['03_','04_','05_']
    elif(season=='summer'):
        summer=True
        ls_check = ['06_','07_','08_']
    elif(season=='autumn'):
        autumn=True
        ls_check = ['09_','10_','11_']
    elif(season=='winter'):
        winter=True
        ls_check = ['12_','01_','02_']
    elif(season=='sunlight'):
        sunlight=True
        ls_check = ['04_','05_','06_','07_','08_','09_']
    else:
        print("Analyzing all seasons. Option given was",season)
        season_adder = ''
        file_adder = ''
        ls_check=['01_','02_','03_','04_','05_','06_','07_','08_','09_','10_','11_','12_']
    
    
    # Grab all the ice files
    #file_names = glob.glob(data_loc+'*.bin')
    # Using this syntax, ignores 2000
    #cmnd = "ls /data/NSIDC/nt_*.bin"
    cmnd = "ls " + home_dir + "/data/NSIDC/nt_*.bin"
    #status,output = commands.getstatusoutput(cmnd)
    #file_initial = output.strip().split('\n')
    #
    #file_names = []
    #for fname in file_initial:
    #    if(fname[50:52] in ls_check):
    #        file_names.append(fname)
    
    
    file_initial = subprocess.check_output(cmnd,shell=True).decode('utf-8').strip().split('\n')
    file_names = []

    if(pre2001==True):
        cmnd = "ls " + home_dir + "/data/NSIDC/pre_2001/nt_*.bin"
        pre_file_initial = subprocess.check_output(cmnd,shell=True).decode('utf-8').strip().split('\n')
        for fname in pre_file_initial:
            if(fname[-17:-14] in ls_check):
                file_names.append(fname)
        
    for fname in file_initial:
        if(fname[-17:-14] in ls_check):
            file_names.append(fname)
    
    # Read in the latitude, longitude, and area data
    latfileo = open(home_dir + '/data/NSIDC/psn25lats_v3.dat','r')
    lats = np.reshape(np.fromfile(latfileo,dtype=np.uint32)/100000.,(448,304))
    lonfileo = open(home_dir + '/data/NSIDC/psn25lons_v3.dat','r')
    tlons = np.fromfile(lonfileo,dtype=np.uint32)/100000.
    tlons[tlons>180.] = tlons[tlons>180.]-42949.67296
    lons = np.reshape(tlons,(448,304))
    areafileo = open(home_dir + '/data/NSIDC/psn25area_v3.dat','r')
    areas = np.reshape(np.fromfile(areafileo,dtype=np.uint32)/1000.,(448,304))
    
    #lons = np.reshape(np.fromfile(lonfileo,dtype=np.uint32)/100000.,(448,304))
    
    num_files = len(file_names)
    NSIDC_data = {}
    NSIDC_data['data'] = np.full((num_files,448,304),-9.)
    NSIDC_data['lat']  = lats
    NSIDC_data['lon']  = lons
    NSIDC_data['area']  = areas
    NSIDC_data['trends'] = np.full((448,304),-9.)
    NSIDC_data['land_trends'] = np.full((448,304),-9.)
    NSIDC_data['titles'] = []
    NSIDC_data['season_adder'] = season_adder
    NSIDC_data['file_adder'] = file_adder
    
    count = 0
    for fname in file_names:
    
        total_name = fname
        #total_name = data_loc+fname
        print(total_name)
        fileo = open(total_name,'rb')
        test = bytearray(fileo.read())
        
        ## To print the header
        #print(test[0:299])
        
        # psg = polar stereographic grid
        
        missing_val             = int(test[0:5].decode('utf-8'))
        ncol_psg                = int(test[6:11].decode('utf-8'))
        nrow_psg                = int(test[12:17].decode('utf-8'))
        lat_enclosed_psg        = float(test[24:29].decode('utf-8'))
        greenwich_orient_psg    = float(test[30:35].decode('utf-8'))
        j_coord_pole            = float(test[42:47].decode('utf-8'))
        i_coord_pole            = float(test[48:53].decode('utf-8'))
        inst_desc               = test[54:59].decode('utf-8')
        data_desc               = test[60:65].decode('utf-8')
        julian_start_day        = int(test[66:71].decode('utf-8'))
        start_hour              = int(test[72:77].decode('utf-8'))
        start_min               = int(test[78:83].decode('utf-8'))
        julian_end_day          = int(test[84:89].decode('utf-8'))
        end_hour                = int(test[90:95].decode('utf-8'))
        end_min                 = int(test[96:101].decode('utf-8'))
        year                    = int(test[102:107].decode('utf-8'))
        julian_day              = int(test[108:113].decode('utf-8'))
        channel_desc            = int(test[114:119].decode('utf-8'))
        scaling_factor          = float(test[120:125].decode('utf-8'))
        file_name               = test[126:149].decode('utf-8')
        image_title             = test[150:229].decode('utf-8')
        data_info               = test[230:299].decode('utf-8')
        
        # Move the data into a np array
        data = np.reshape(np.array(test[300:]),(448,304))
        
        #for i in range(len(data)):
        #    templine = data[i,:].astype(np.float)
        #    templine[templine<251] = (templine[templine<251]/scaling_factor)*100.
        #    templine = templine.astype(int)
        #    print(templine)
        
        #data[np.where(data>=251)]=data[np.where(data>=251)]*0.
        data[np.where(data<251)]=  \
            (data[np.where(data<251)]/scaling_factor)*100.
    
        NSIDC_data['data'][count,:,:] = data[:,:]
        NSIDC_data['titles'].append(file_name)
    #    total_data[:,:,count] = data[:,:]
        count+=1

    return NSIDC_data

# This function downloads an NSIDC netCDF file of daily sea ice concentration
# NOTE: This function may not work if not already logged in to Earthdata.
#       May need to fix this to account for URS cookies.
# ---------------------------------------------------------------------------
def download_NSIDC_daily(date_str, save_dir = data_dir + 'daily/'):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')
    base_link = "https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0051.002/"
    file_str = 'NSIDC0051_SEAICE_PS_N25km_%Y%m%d_v2.0.nc'

    # First check if the file is already there
    # ----------------------------------------
    if(os.path.exists(dt_date_str.strftime(save_dir + file_str))):
        print("NSIDC file already exists.")
        return
    else:

        total_cmnd = dt_date_str.strftime('wget ' + base_link + '%Y.%m.%d/' + file_str)
        print(total_cmnd)
        os.system(total_cmnd)
        total_cmnd = dt_date_str.strftime('mv ' + file_str + ' ' + save_dir)
        print(total_cmnd)
        os.system(total_cmnd)

def download_NSIDC_monthly(date_str, save_dir = data_dir + 'monthly/'):

    dt_date_str = datetime.strptime(date_str, '%Y%m')
    base_link = "https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0051.002/"
    file_str = 'NSIDC0051_SEAICE_PS_N25km_%Y%m_v2.0.nc'

    # First check if the file is already there
    # ----------------------------------------
    if(os.path.exists(dt_date_str.strftime(save_dir + file_str))):
        print("NSIDC file already exists.")
        return
    else:

        total_cmnd = dt_date_str.strftime('wget ' + base_link + '%Y.%m.%d/' + file_str)
        print(total_cmnd)
        os.system(total_cmnd)
        total_cmnd = dt_date_str.strftime('mv ' + file_str + ' ' + save_dir)
        print(total_cmnd)
        os.system(total_cmnd)

def readNSIDC_monthly_grid_all(begin_date, end_date, season, minlat = 70., \
        maxlat = 90., calc_month = True):

    spring=False
    summer=False
    autumn=False
    winter=False
    sunlight=False
    if(season=='spring'):
        spring=True
    elif(season=='summer'):
        summer=True
    elif(season=='autumn'):
        autumn=True
    elif(season=='winter'):
        winter=True
    elif(season=='sunlight'):
        sunlight=True
   
    # Grab all the files
    base_path = home_dir + '/data/NSIDC/monthly/'
    total_list = sorted(glob.glob(base_path+'*.nc')) 

    # Loop over all files and find the ones that match with the desired times
    begin_date = datetime.strptime(str(begin_date),'%Y%m')
    end_date = datetime.strptime(str(end_date),'%Y%m')
    final_list = []
    for f in total_list:
        fdate = f.strip().split('/')[-1].split('_')[-2]
        fdate = datetime.strptime(str(fdate),'%Y%m')
        if((fdate >= begin_date) & (fdate <= end_date)):
            if(spring==True):
                if((fdate.month >= 3) & (fdate.month < 6)):
                    final_list.append(f)
            elif(summer==True):
                if((fdate.month >= 6) & (fdate.month < 9)):
                    final_list.append(f)
            elif(autumn==True):
                if((fdate.month >= 9) & (fdate.month < 12)):
                    final_list.append(f)
            elif(winter==True):
                if((fdate.month == 12) | (fdate.month < 3)):
                    final_list.append(f)
            elif(sunlight==True):
                if((fdate.month > 3) & (fdate.month < 10)):
                    final_list.append(f)
            else:
                final_list.append(f)
    time_dim = len(final_list)

    dates = [ffile.strip().split('/')[-1].split('_')[-2] for ffile in \
        final_list]

    NSIDC_data = {}
    NSIDC_data['begin_dt_date'] = begin_date
    NSIDC_data['end_dt_date']   = end_date
    NSIDC_data['season']        = season
    NSIDC_data['dates']         = dates
    NSIDC_data['data']          = np.full((time_dim, 448, 304), np.nan)
    #NSIDC_data['area']          = np.full((448, 304), np.nan)

    for ii, work_str in enumerate(dates):
        #work_str = ffile[0].strip().split('/')[-1].split('_')[-2]
        NSIDC_local = readNSIDC_monthly(work_str, maxlat = maxlat)

        NSIDC_data['data'][ii,:,:] = NSIDC_local['data']
        NSIDC_data['maxlat'] = maxlat
   
    NSIDC_data['lat'] = NSIDC_local['lat'] 
    NSIDC_data['lon'] = NSIDC_local['lon'] 
    NSIDC_data['area'] = NSIDC_local['area'] 

    NSIDC_data = grid_data_conc(NSIDC_data, minlat = minlat, maxlat = maxlat)

    if(calc_month == True):
        NSIDC_data = calcNSIDC_MonthClimo(NSIDC_data)
   
    # Process the original 25x25 data
    # ------------------------------- 
    NSIDC_data['pole_hole'] = np.ma.masked_where((NSIDC_data['data'] != 251), \
        NSIDC_data['data'])
    NSIDC_data['unused'] = np.ma.masked_where((NSIDC_data['data'] != 252), \
        NSIDC_data['data'])
    NSIDC_data['coastline'] = np.ma.masked_where((NSIDC_data['data'] != 253), \
        NSIDC_data['data'])
    NSIDC_data['land'] = np.ma.masked_where((NSIDC_data['data'] != 254), \
        NSIDC_data['data'])
    NSIDC_data['data'] = np.ma.masked_where((NSIDC_data['data'] < 0) | \
        (NSIDC_data['data'] > 100), NSIDC_data['data'])

    # Process the 1x1 gridded data
    # ---------------------------- 
    ##!#NSIDC_data['grid_pole_hole'] = np.ma.masked_where(\
    ##!#    (NSIDC_data['grid_ice_conc'] != 251), NSIDC_data['grid_ice_conc'])
    ##!#NSIDC_data['grid_unused'] = np.ma.masked_where(\
    ##!#    (NSIDC_data['grid_ice_conc'] != 252), NSIDC_data['grid_ice_conc'])
    ##!#NSIDC_data['grid_coastline'] = np.ma.masked_where(\
    ##!#    (NSIDC_data['grid_ice_conc'] != 253), NSIDC_data['grid_ice_conc'])
    ##!#NSIDC_data['grid_land'] = np.ma.masked_where(\
    ##!#    (NSIDC_data['grid_ice_conc'] != 254), NSIDC_data['grid_ice_conc'])
    NSIDC_data['grid_ice_conc'] = np.ma.masked_where((NSIDC_data['grid_ice_conc'] < 0) | \
        (NSIDC_data['grid_ice_conc'] > 100), NSIDC_data['grid_ice_conc'])
    NSIDC_data['MONTH_CLIMO'] = np.ma.masked_where((NSIDC_data['MONTH_CLIMO'] < 0) | \
        (NSIDC_data['MONTH_CLIMO'] > 100), NSIDC_data['MONTH_CLIMO'])

    NSIDC_data['minlat'] = minlat
    NSIDC_data['maxlat'] = maxlat

    return NSIDC_data


def readNSIDC_monthly(date_str, grid_data = False, maxlat = 90.):

    NSIDC_data = readNSIDC_daily(date_str, grid_data = grid_data, maxlat = maxlat)

    return NSIDC_data

def readNSIDC_daily(date_str, grid_data = False, maxlat = 90):

    if(len(date_str) == 8):
        str_fmt = '%Y%m%d'
        file_dir = data_dir + 'daily/'
    elif(len(date_str) == 6):
        str_fmt = '%Y%m'
        file_dir = data_dir + 'monthly/'
    else:
        print("ERROR: Invalid date format")
        print(date_str)
    
    dt_date_str = datetime.strptime(date_str, str_fmt) 
  
    # Read in the latitude, longitude, area and ice data
    # --------------------------------------------------
    latfileo = open(home_dir + '/data/NSIDC/psn25lats_v3.dat','r')
    lats = np.reshape(np.fromfile(latfileo,dtype=np.uint32)/100000.,(448,304))
    lonfileo = open(home_dir + '/data/NSIDC/psn25lons_v3.dat','r')
    tlons = np.fromfile(lonfileo,dtype=np.uint32)/100000.
    tlons[tlons>180.] = tlons[tlons>180.]-42949.67296
    lons = np.reshape(tlons,(448,304))
    areafileo = open(home_dir + '/data/NSIDC/psn25area_v3.dat','r')
    areas = np.reshape(np.fromfile(areafileo,dtype=np.uint32)/1000.,(448,304))

    # Look for both the binary file and the netCDF file.
    # --------------------------------------------------
    bin_filename = file_dir + 'nt_' + date_str + '_f17_v1.1_n.bin'
    net_filename = file_dir + 'NSIDC0051_SEAICE_PS_N25km_' + date_str + '_v2.0.nc'
    #print(bin_filename)

    bin_found = os.path.exists(bin_filename)
    net_found = os.path.exists(net_filename)

    if((bin_found) and not (net_found)):
 
        print('Reading ', bin_filename)
        fileo = open(filename,'rb')
        test = bytearray(fileo.read())
        
        # psg = polar stereographic grid
        
        missing_val             = int(test[0:5].decode('utf-8'))
        ncol_psg                = int(test[6:11].decode('utf-8'))
        nrow_psg                = int(test[12:17].decode('utf-8'))
        lat_enclosed_psg        = float(test[24:29].decode('utf-8'))
        greenwich_orient_psg    = float(test[30:35].decode('utf-8'))
        j_coord_pole            = float(test[42:47].decode('utf-8'))
        i_coord_pole            = float(test[48:53].decode('utf-8'))
        inst_desc               = test[54:59].decode('utf-8')
        data_desc               = test[60:65].decode('utf-8')
        julian_start_day        = int(test[66:71].decode('utf-8'))
        start_hour              = int(test[72:77].decode('utf-8'))
        start_min               = int(test[78:83].decode('utf-8'))
        julian_end_day          = int(test[84:89].decode('utf-8'))
        end_hour                = int(test[90:95].decode('utf-8'))
        end_min                 = int(test[96:101].decode('utf-8'))
        year                    = int(test[102:107].decode('utf-8'))
        julian_day              = int(test[108:113].decode('utf-8'))
        channel_desc            = int(test[114:119].decode('utf-8'))
        scaling_factor          = float(test[120:125].decode('utf-8'))
        file_name               = test[126:149].decode('utf-8')
        image_title             = test[150:229].decode('utf-8')
        data_info               = test[230:299].decode('utf-8')
        
        # Move the data into a np array
        data = np.reshape(np.array(test[300:]),(448,304))
        
        data[np.where(data<251)]=  \
            (data[np.where(data<251)]/scaling_factor)*100.
   
    else:
        print('Reading ', net_filename)
        # Read the netCDf file
        # --------------------
        in_data = Dataset(net_filename, 'r')
        if('F17_ICECON' in in_data.variables.keys()):
            data = in_data['F17_ICECON'][0,:,:].data
            data[data <= 1.0] = data[data <= 1.0] * 100.
            #data = in_data['F17_ICECON'][0,:,:] * 100.
        elif('F13_ICECON' in in_data.variables.keys()):
            data = in_data['F13_ICECON'][0,:,:].data
            data[data <= 1.0] = data[data <= 1.0] * 100.
        else:
            print("ERROR: incorrect NSIDC ice variable name")
            in_data.close()
            return

        in_data.close()

    data = np.ma.masked_where(lats >= maxlat, data)
    NSIDC_data = {}
    NSIDC_data['data'] = data
    NSIDC_data['lat']  = lats
    NSIDC_data['lon']  = lons
    NSIDC_data['area']  = areas
    NSIDC_data['date'] = date_str
    NSIDC_data['maxlat'] = maxlat

    if(grid_data):
        NSIDC_data = grid_NSIDC_data(NSIDC_data)
 
    return NSIDC_data

def write_NSIDC_daily_to_file(begin_date, end_date, minlat = 65.5):

    print("HIIII")

    dt_begin_str = datetime.strptime(begin_date, '%Y%m%d')
    dt_end_str   = datetime.strptime(end_date, '%Y%m%d')

    # Grab all the files
    base_path = home_dir + '/data/NSIDC/'
    total_list = sorted(glob.glob(base_path+'*.nc')) 

    # Loop over all files and find the ones that match with the desired times

    # Step 1: Find all the filenames to analyze
    # -----------------------------------------

    # Step 1b: Read in the daily OMI file to ensure that the file times
    # match those in the daily OMI product. There are some days missing
    # from the daily OMI product, so need to avoid including those 
    # missing times in the NSIDC daily products
    # ----------------------------------------------------------------

    # Step 2: From the number of daily filenames, allocate data arrays
    # to hold the final, combined data
    # ----------------------------------------------------------------

    # Step 3: Loop over each of the filenames
    # ---------------------------------------

    # Step 4: Load in the current day's data
    # --------------------------------------

    # Step 5: Grid the current data to the lat/lon grid
    # -------------------------------------------------

    # Step 6: Insert the current gridded data into the final array
    # ------------------------------------------------------------

def calcNSIDC_MonthClimo(NSIDC_data):

    if('grid_ice_conc' not in NSIDC_data.keys()):
        print("ERROR: no gridded data. Must grid data before month averaging")
        return   
 
    # Set up arrays to hold monthly climatologies
    month_grid_climo = np.full((6, NSIDC_data['grid_lat'].shape[0], \
        NSIDC_data['grid_lat'].shape[1]), np.nan)

    local_grid_data = np.copy(NSIDC_data['grid_ice_conc'][:,:,:])
    local_grid_mask = np.ma.masked_where((local_grid_data < 0) | \
        (local_grid_data > 100), local_grid_data)

    # Calculate monthly climatologies
    for m_i in range(6):
        month_grid_climo[m_i,:,:] = np.nanmean(\
            local_grid_mask[m_i::6,:,:],axis=0)
        print("Month: ",NSIDC_data['dates'][m_i][4:]) 
    ##!## Calculate monthly climatologies
    ##!## Assume all Aqua CERES data starting at 200207 is read in previously
    ##!#month_climo[0,:,:]  = np.nanmean(local_mask[6::12,:,:],axis=0)  # January 
    ##!#month_climo[1,:,:]  = np.nanmean(local_mask[7::12,:,:],axis=0)  # February
    ##!#month_climo[2,:,:]  = np.nanmean(local_mask[8::12,:,:],axis=0)  # March 
    ##!#month_climo[3,:,:]  = np.nanmean(local_mask[9::12,:,:],axis=0)  # April 
    ##!#month_climo[4,:,:]  = np.nanmean(local_mask[10::12,:,:],axis=0) # May 
    ##!#month_climo[5,:,:]  = np.nanmean(local_mask[11::12,:,:],axis=0) # June 
    ##!#month_climo[6,:,:]  = np.nanmean(local_mask[0::12,:,:],axis=0)  # July 
    ##!#month_climo[7,:,:]  = np.nanmean(local_mask[1::12,:,:],axis=0)  # August 
    ##!#month_climo[8,:,:]  = np.nanmean(local_mask[2::12,:,:],axis=0)  # September 
    ##!#month_climo[9,:,:]  = np.nanmean(local_mask[3::12,:,:],axis=0)  # October 
    ##!#month_climo[10,:,:] = np.nanmean(local_mask[4::12,:,:],axis=0)  # November
    ##!#month_climo[11,:,:] = np.nanmean(local_mask[5::12,:,:],axis=0)  # December

    # Insert data into dictionary
    NSIDC_data['MONTH_CLIMO'] = month_grid_climo

    return NSIDC_data

def calcNSIDC_grid_trend(NSIDC_data, month_idx, trend_type, mingrid_lat):

    if(month_idx == None):
        month_idx = 0
        index_jumper = 1
    else:
        month_adder = '_month'
        if(NSIDC_data['season'] == 'sunlight'):
            index_jumper = 6
        else:   
            index_jumper = 12

    lat_ranges = NSIDC_data['grid_lat'][:,0]
    #grid_lat_ranges = np.arange(mingrid_lat,90,1.0)
    lon_ranges = NSIDC_data['grid_lon'][0,:]
    #lon_ranges = np.arange(-180,180,1.0)

    # Make copy of NSIDC_data array
    print(NSIDC_data['dates'][month_idx::index_jumper])
    local_data   = np.copy(NSIDC_data['grid_ice_conc'][month_idx::index_jumper,:,:])
    local_data = np.ma.masked_invalid(local_data)
    #local_counts = np.copy(NSIDC_data['OB_COUNT'][month_idx::index_jumper,:,:])
    #local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where(((local_data < 0.) | \
        (local_data > 100) | 
        (NSIDC_data['grid_lat'] < mingrid_lat)), local_data)
    nsidc_trends = np.full(local_data.shape[1:], np.nan)
    nsidc_pvals  = np.full(local_data.shape[1:], np.nan)
    nsidc_uncert = np.full(local_data.shape[1:], np.nan)

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    for i in range(0,len(lat_ranges)):
        for j in range(0,len(lon_ranges)):
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
                    nsidc_trends[i,j] = result.slope * len(x_vals)
                    nsidc_pvals[i,j]  = result.pvalue
                    nsidc_uncert[i,j] = result.stderr * len(x_vals)
                else:
                    res = stats.theilslopes(work_mask.compressed(), x_vals, 0.90)
                    nsidc_trends[i,j] = res[0]*len(x_vals)
            #else:
            #    print('no data')

    #nsidc_trends = np.ma.masked_where(((NSIDC_data['grid_lat'] < mingrid_lat) | \
    #    (nsidc_trends == -999.)), nsidc_trends)
    nsidc_trends = np.ma.masked_where(NSIDC_data['grid_lat'] < mingrid_lat, nsidc_trends)
    nsidc_pvals  = np.ma.masked_where(NSIDC_data['grid_lat'] < mingrid_lat, nsidc_pvals)
    nsidc_uncert = np.ma.masked_where(NSIDC_data['grid_lat'] < mingrid_lat, nsidc_uncert)

    return nsidc_trends, nsidc_pvals, nsidc_uncert

# title is the plot title
# ptype is 'climo', 'trend'
def plotNSIDC_spatial(pax, plat, plon, pdata, ptype, ptitle = '', plabel = '', \
        vmin = None, vmax = None, colorbar_label_size = 16, minlat = 65., \
        colorbar = True, \
        pvals = None):

    if(vmin == None):
        vmin = np.nanmin(pdata)
    if(vmax == None):
        vmax = np.nanmax(pdata)

    if(ptype == 'trend'):
        colormap = plt.cm.bwr
    elif(ptype == 'uncert'):
        #colormap = plt.cm.plasma
        colormap = plt.cm.get_cmap('jet', 10)
        vmax = 30
    else:
        colormap = plt.cm.jet

    # Make copy of OMI_data array
    local_data  = np.copy(pdata)
    # Grid the data, fill in white space
    cyclic_data,cyclic_lons = add_cyclic_point(local_data,plon[0,:])
    plat,plon = np.meshgrid(plat[:,0],cyclic_lons)   
  
    # Mask any missing values
    mask_AI = np.ma.masked_where(cyclic_data < -998.9, cyclic_data)
    mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)

    # Plot lat/lon lines
    gl = pax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, \
        linewidth = 1, color = 'gray', alpha = 0.5, linestyle = '-',\
        y_inline = True, xlocs = range(-180, 180, 30), ylocs = range(70, 90, 10))
    #gl.top_labels = False
    #gl.bottom_labels = False
    #gl.left_labels = False
    #gl.right_labels = False
    #gl.xlines = False
    #gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 30))
    #gl.ylocator = mticker.FixedLocator(np.arange(70, 90, 10))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    #gl.xlabel_style = {'size': 5, 'color': 'gray'}
    gl.xlabel_style = {'color': 'gray', 'weight': 'bold'}
    gl.ylabel_style = {'size': 15, 'color': 'gray'}
    gl.ylabel_style = {'color': 'black', 'weight': 'bold'}
    #pax.gridlines()

    #pax.gridlines()
    pax.coastlines(resolution='50m')
    mesh = pax.pcolormesh(plon, plat,\
            mask_AI.T,transform = datacrs,shading = 'auto',\
            cmap = colormap,vmin=vmin,vmax=vmax)
    if(pvals is not None):
        cyclic_pvals, cyclic_lons = add_cyclic_point(pvals, plon[0,:])
        print(pvals.shape, plat.shape, plat2.shape)
        mask_pvals = np.ma.masked_where((plat < minlat) | \
            (pvals > 0.05), pvals)
        pax.pcolor(plon, plat, mask_pvals, hatch = '...', alpha = 0.0, \
            shading = 'auto', transform = datacrs)

    pax.set_extent([-180,180,minlat,90],datacrs)
    pax.set_boundary(circle, transform=pax.transAxes)
    #cbar = plt.colorbar(mesh,ticks = np.arange(-200.,400.,5.),\    
    if(colorbar):
        cbar = plt.colorbar(mesh,\
            ax = pax, orientation='vertical',shrink = 0.8, extend = 'both')
        cbar.set_label(plabel,fontsize=colorbar_label_size,weight='bold')
    print("USING ",ptitle)
    pax.set_title(ptitle)


# Plot a monthly climatology 
def plotNSIDC_MonthClimo(NSIDC_data,month_idx,minlat = 60, ax = None, \
        colorbar = True, title = None, save = False):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.ocean
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    month_obj = datetime(year = 1,month = int(NSIDC_data['dates'][month_idx][4:]),day = 1)
    str_month = month_obj.strftime('%B')

    # Make figure title
    if(title == None):
        title = 'NSIDC Sea Ice Conc. ' + str_month + \
            ' Climatology\n'+str_month+'. ' + NSIDC_data['dates'][0][:4] \
            + ' - '+str_month+'. ' + NSIDC_data['dates'][-1][:4]

    # Make figure
    in_ax = True 
    if(ax is None): 
        in_ax = False
        fig1 = plt.figure(figsize = (6,6))
        ax = plt.axes(projection = mapcrs)

    mask_data = np.ma.masked_where(\
        NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 10, 
        NSIDC_data['MONTH_CLIMO'][month_idx,:,:])
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'],\
            mask_data,transform = datacrs,\
            cmap = colormap)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
    #ax.set_xlim(-3430748.535086173,3430748.438879491)
    #ax.set_ylim(-3413488.8763307533,3443353.899053069)
    if(colorbar):
        #cbar = plt.colorbar(mesh,ax = pax, orientation='horizontal',pad=0,\
        cbar = plt.colorbar(mesh,ax = ax, orientation='vertical',\
            pad = 0.04, fraction = 0.040)
        cbar.set_label('Sea Ice Concentration', fontsize = 12, weight='bold')
    #cbar = plt.colorbar(mesh,ticks = np.arange(0.0,400.,25.),orientation='horizontal',pad=0,\
    #    aspect=50,shrink = 0.845,label=CERES_data['parm_name'])
    ax.set_title(title)

    if(not in_ax):
        fig1.tight_layout()
        if(save):
            outname = 'nsidc_month_climo_'+ month_obj.strftime('%b') + '.png'
            fig1.savefig(outname, dpi=300)
            print("Saved image",outname)
        else:
            plt.show()


# Designed to work with the netCDF data
def plotNSIDC_MonthTrend(NSIDC_data,month_idx=None,save=False,\
        trend_type='standard',season='',minlat=65.,return_trend=False, \
        colorbar = True, colorbar_label_size = None,title = None, \
        pax = None, show_pval = False, uncert_ax = None):

    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='_thielSen'

    if(month_idx == None):
        month_adder = ''
        month_idx = 0
        index_jumper = 1
        do_month = False
        v_max = 30.
        v_min = -30.
    else:
        month_adder = '_month'
        if(NSIDC_data['season'] == 'sunlight'):
            index_jumper = 6
        else:   
            index_jumper = 12
        do_month = True
        v_max = 30.
        v_min = -30.

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # --------------------------------------------------------------
    #
    # Use calcNSIDC_grid_trend to calculate the trends in the AI data
    #
    # --------------------------------------------------------------
    nsidc_trends, nsidc_pvals, nsidc_uncert = \
        calcNSIDC_grid_trend(NSIDC_data, month_idx, trend_type, \
        minlat)

    if(not show_pval):
        nsidc_pvals = None
    else:
        print('month_idx = ',month_idx,' PVAL nanmean = ', \
            np.nanmean(nsidc_pvals))

    if(uncert_ax is None):
        nsidc_uncert = None
    # --------------------------------------------------------------
    #
    # Plot the calculated trends on a figure
    #
    # --------------------------------------------------------------

    # Set up mapping variables 
    colormap = plt.cm.bwr

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(NSIDC_data['dates'][month_idx::index_jumper][0],'%Y%m')
    end_date   = datetime.strptime(NSIDC_data['dates'][month_idx::index_jumper][-1],'%Y%m')
    
    # Make figure title
    #date_month = datetime(year = 1,month = month_idx+1, day = 1).strftime('%B')
    month_string = ''
    if(do_month == True):
        month_string = start_date.strftime('%B') + ' '

    if(title is None):
        title = 'NSIDC ' +  month_string + 'Trends'\
            '\n'+start_date.strftime('%b. %Y') + ' - ' +\
            end_date.strftime('%b. %Y')

    # Call plotNSIDC_spatial to add the data to the figure

    #ii = 10
    #jj = 44
    #print(NSIDC_data['grid_lat'][ii,jj], NSIDC_data['grid_lon'][ii,jj])
    #print(NSIDC_data['grid_ice_conc'][month_idx::index_jumper,ii,jj])
    #fig2 = plt.figure()
    #ax4 = fig2.add_subplot(1,1,1)
    #ax4.plot(NSIDC_data['grid_ice_conc'][month_idx::index_jumper,ii,jj])
    

    if(pax is None):
        #plt.close('all')
        fig1 = plt.figure(figsize = (6,6))
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

        plotNSIDC_spatial(ax, NSIDC_data['grid_lat'], NSIDC_data['grid_lon'], \
            nsidc_trends, 'trend', ptitle = title, plabel = 'Sea Ice Conc. per study period', \
            vmin = v_min, vmax = v_max, colorbar_label_size = colorbar_label_size, \
            minlat = minlat, pvals = nsidc_pvals, colorbar = colorbar)

        fig1.tight_layout()

        if(save == True):
            month_adder = ''
            if(do_month == True):
                month_adder = '_' + start_date.strftime('%b') 
            out_name = 'nsidc_trend'+ month_adder + '_' + \
                start_date.strftime('%Y%m') + '_' + end_date.strftime('%Y%m') + \
                '_min' + str(int(minlat)) + '.png'
            plt.savefig(out_name,dpi=300)
            print("Saved image",out_name)
        else:
            plt.show()
    else:
        plotNSIDC_spatial(pax, NSIDC_data['grid_lat'], NSIDC_data['grid_lon'], \
            nsidc_trends, 'trend', ptitle = title, plabel = 'Sea Ice Conc. per study period', \
            vmin = v_min, vmax = v_max, colorbar_label_size = colorbar_label_size, \
            minlat = minlat, colorbar = colorbar)

    if(uncert_ax is not None):
        plotNSIDC_spatial(uncert_ax, NSIDC_data['grid_lat'], \
            NSIDC_data['grid_lon'], \
            nsidc_uncert, 'uncert', ptitle = title, plabel = 'Sea Ice Conc.', \
            colorbar = colorbar, colorbar_label_size = colorbar_label_size, \
            vmin = 0, vmax = 20.0, minlat = minlat)

    if(return_trend == True):
        return nsidc_trends

# Generate a 15-panel figure comparing the climatology and trend between 3
# versions of the CERES data for all months
def plotNSIDC_ClimoTrend_all(NSIDC_data,\
        trend_type = 'standard', minlat=65.,save=False):

    colormap = plt.cm.jet

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    index_jumper = 6 

    colorbar_label_size = 7
    axis_title_size = 8
    row_label_size = 10 
    #colorbar_label_size = 13
    #axis_title_size = 14.5
    #row_label_size = 14.5

    #fig = plt.figure()
    plt.close('all')
    fig = plt.figure(figsize=(9.2,13))
    #plt.suptitle('NAAPS Comparisons: '+start_date.strftime("%B"),y=0.95,\
    #    fontsize=18,fontweight=4,weight='bold')
    gs = gridspec.GridSpec(nrows=6, ncols=3, hspace = 0.001, wspace = 0.15)

    # - - - - - - - - - - - - - - - - - - - - -
    # Plot the climatologies along the top row
    # - - - - - - - - - - - - - - - - - - - - -
       
    # Plot DATA1 climos
    # -----------------
    ##!## Make copy of CERES_data array
    local_data1_Apr  = np.copy(NSIDC_data['MONTH_CLIMO'][0,:,:])
    local_data1_May  = np.copy(NSIDC_data['MONTH_CLIMO'][1,:,:])
    local_data1_Jun  = np.copy(NSIDC_data['MONTH_CLIMO'][2,:,:])
    local_data1_Jul  = np.copy(NSIDC_data['MONTH_CLIMO'][3,:,:])
    local_data1_Aug  = np.copy(NSIDC_data['MONTH_CLIMO'][4,:,:])
    local_data1_Sep  = np.copy(NSIDC_data['MONTH_CLIMO'][5,:,:])

    mask_AI1_Apr = np.ma.masked_where(local_data1_Apr == -999.9, local_data1_Apr)
    mask_AI1_May = np.ma.masked_where(local_data1_May == -999.9, local_data1_May)
    mask_AI1_Jun = np.ma.masked_where(local_data1_Jun == -999.9, local_data1_Jun)
    mask_AI1_Jul = np.ma.masked_where(local_data1_Jul == -999.9, local_data1_Jul)
    mask_AI1_Aug = np.ma.masked_where(local_data1_Aug == -999.9, local_data1_Aug)
    mask_AI1_Sep = np.ma.masked_where(local_data1_Sep == -999.9, local_data1_Sep)
    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,CERES_data1['LON'][0,:])
    ##!#plat,plon = np.meshgrid(CERES_data1['LAT'][:,0],cyclic_lons)   
  
    # Mask any missing values
    #mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)
    ax00 = plt.subplot(gs[0,0], projection=mapcrs)   # April climo original
    ax01 = plt.subplot(gs[0,1], projection=mapcrs)   # April climo screened
    ax02 = plt.subplot(gs[0,2], projection=mapcrs)   # April trend original
    ax10 = plt.subplot(gs[1,0], projection=mapcrs)   # May climo original
    ax11 = plt.subplot(gs[1,1], projection=mapcrs)   # May climo screened
    ax12 = plt.subplot(gs[1,2], projection=mapcrs)   # May trend original
    ax20 = plt.subplot(gs[2,0], projection=mapcrs)   # June climo original
    ax21 = plt.subplot(gs[2,1], projection=mapcrs)   # June climo screened
    ax22 = plt.subplot(gs[2,2], projection=mapcrs)   # June trend original
    ax30 = plt.subplot(gs[3,0], projection=mapcrs)   # July climo original
    ax31 = plt.subplot(gs[3,1], projection=mapcrs)   # July climo screened
    ax32 = plt.subplot(gs[3,2], projection=mapcrs)   # July trend original
    ax40 = plt.subplot(gs[4,0], projection=mapcrs)   # August climo original
    ax41 = plt.subplot(gs[4,1], projection=mapcrs)   # August climo screened
    ax42 = plt.subplot(gs[4,2], projection=mapcrs)   # August trend original
    ax50 = plt.subplot(gs[5,0], projection=mapcrs)   # September climo original
    ax51 = plt.subplot(gs[5,1], projection=mapcrs)   # September climo screened
    ax52 = plt.subplot(gs[5,2], projection=mapcrs)   # September trend original

    # Plot the figures in the first row: April
    # ---------------------------------------
    cbar_switch = False
    plotNSIDC_MonthClimo(NSIDC_data,0,minlat = minlat, ax = ax00, title = '', \
        colorbar = cbar_switch)
    plotNSIDC_MonthTrend(NSIDC_data,month_idx=0,save=False,\
        trend_type='standard',season='sunlight',minlat=minlat,\
        return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax01, show_pval = True, uncert_ax = ax02, title = '')
    #plotNSIDC_spatial(ax00, NSIDC_data['lats'], NSIDC_data['lons'], mask_AI1_Apr, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 1.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotNSIDC_MonthTrend(NSIDC_data,month_idx=0,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax01, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax02)
    #plotNSIDC_MonthTrend(NSIDC_data,month_idx=month_idx,save=False,\
    #    trend_type='standard',season='sunlight',minlat=65.,return_trend=False, \
    #    pax = ax1)

    # Plot the figures in the first row: May
    # ---------------------------------------
    #plotNSIDC_spatial(ax10, NSIDC_data['lats'], NSIDC_data['lons'], mask_AI1_May, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 2.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotNSIDC_MonthTrend(NSIDC_data,month_idx=1,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax11, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax12)
    plotNSIDC_MonthClimo(NSIDC_data,1,minlat = minlat, ax = ax10, title = '', \
        colorbar = cbar_switch)
    plotNSIDC_MonthTrend(NSIDC_data,month_idx=1,save=False,\
        trend_type='standard',season='sunlight',\
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax11, show_pval = True, uncert_ax = ax12, title = '')

    # Plot the figures in the first row: June
    # ---------------------------------------
    #plotNSIDC_spatial(ax20, NSIDC_data['lats'], NSIDC_data['lons'], mask_AI1_Jun, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 3.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotNSIDC_MonthTrend(NSIDC_data,month_idx=2,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax21, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax22)
    plotNSIDC_MonthClimo(NSIDC_data,2,minlat = minlat, ax = ax20, title = '', \
        colorbar = cbar_switch)
    plotNSIDC_MonthTrend(NSIDC_data,month_idx=2,save=False,\
        trend_type='standard',season='sunlight',
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax21, show_pval = True, uncert_ax = ax22, title = '')

    # Plot the figures in the second row: July
    # ----------------------------------------
    #plotNSIDC_spatial(ax30, NSIDC_data['lats'], NSIDC_data['lons'], mask_AI1_Jul, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 3.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotNSIDC_MonthTrend(NSIDC_data,month_idx=3,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax31, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax32)
    plotNSIDC_MonthClimo(NSIDC_data,3,minlat = minlat, ax = ax30, title = '', \
        colorbar = cbar_switch)
    plotNSIDC_MonthTrend(NSIDC_data,month_idx=3,save=False,\
        trend_type='standard',season='sunlight',\
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax31, show_pval = True, uncert_ax = ax32, title = '')

    # Plot the figures in the third row: August
    # -----------------------------------------
    #plotNSIDC_spatial(ax40, NSIDC_data['lats'], NSIDC_data['lons'], mask_AI1_Aug, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 3.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotNSIDC_MonthTrend(NSIDC_data,month_idx=4,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax41, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax42)
    plotNSIDC_MonthClimo(NSIDC_data,4,minlat = minlat, ax = ax40, title = '', \
        colorbar = cbar_switch)
    plotNSIDC_MonthTrend(NSIDC_data,month_idx=4,save=False,\
        trend_type='standard',season='sunlight',\
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax41, show_pval = True, uncert_ax = ax42, title = '')

    # Plot the figures in the third row: September
    # --------------------------------------------
    #plotNSIDC_spatial(ax50, NSIDC_data['lats'], NSIDC_data['lons'], mask_AI1_Sep, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 1.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotNSIDC_MonthTrend(NSIDC_data,month_idx=5,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax51, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax52)
    plotNSIDC_MonthClimo(NSIDC_data,5,minlat = minlat, ax = ax50, title = '', \
        colorbar = cbar_switch)
    plotNSIDC_MonthTrend(NSIDC_data,month_idx=5,save=False,\
        trend_type='standard',season='sunlight',\
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax51, show_pval = True, uncert_ax = ax52, title = '')

    fig.text(0.10, 0.82, 'April', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.692, 'May', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.565, 'June', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.435, 'July', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.305, 'August', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.18, 'September', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)

    fig.text(0.250, 0.90, 'Climatology', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.500, 0.90, 'Trend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.75, 0.90, 'Standard Error of\nScreened Trend (Slope)', \
        ha='center', va='center', rotation='horizontal',weight='bold',\
        fontsize=row_label_size)

    plt.suptitle('NSIDC')

    #cax = fig.add_axes([0.15, 0.09, 0.35, 0.01])
    #norm = mpl.colors.Normalize(vmin = -0.5, vmax = 0.5)
    #cb1 = mpl.colorbar.ColorbarBase(cax, cmap = plt.cm.bwr, norm = norm, \
    #    orientation = 'horizontal', extend = 'both')
    #cb1.set_label('Sfc Smoke Trend (Smoke / Study Period)', \
    #    weight = 'bold')

    #cax2 = fig.add_axes([0.530, 0.09, 0.35, 0.01])
    #norm2 = mpl.colors.Normalize(vmin = 0.0, vmax = 0.3)
    #cmap = plt.cm.get_cmap('jet', 6)
    #cb2 = mpl.colorbar.ColorbarBase(cax2, cmap = cmap, norm = norm2, \
    #    orientation = 'horizontal', extend = 'both')
    #cb2.set_label('Standard Error of Smoke Trend (Slope)', weight = 'bold')

    outname = '_'.join(['nsidc','ice','grid','comps','all6']) + '.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()


def plotNSIDC_daily(pax, NSIDC_data, minlat=65, \
        vmin = None, vmax = None, title = '', label = None, \
        labelsize = 11, labelticksize = 9, circle_bound = False, \
        gridlines = True, zoom = False):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.ocean
    ##!#if(pax is None):
    ##!#    if(minlat < 45):
    ##!#        mapcrs = ccrs.Miller()
    ##!#    else:
    ##!#        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    ##!#if(grid_data):
    ##!#    plot_lon  = NSIDC_data['lon']
    ##!#    plot_lat  = NSIDC_data['lat']
    ##!#    if(param == 'total'):
    ##!#        mask_flux = NSIDC_data['swf'] + NSIDC_data['lwf']
    ##!#    else:
    ##!#        mask_flux = NSIDC_data[param.lower()]
    ##!#else:
    mask_data = np.ma.masked_where((NSIDC_data['data'] > 250) | \
        (NSIDC_data['data'] == 0), \
        NSIDC_data['data'])
 
    if(title == ''):
        title =  'SSMI/S Sea Ice Concentration\n' + NSIDC_data['date']
    if(label == None):
        label = 'Concentration [%]'
    pax.set_title(title)

    #plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
    mesh = pax.pcolormesh(NSIDC_data['lon'], NSIDC_data['lat'], mask_data,\
        transform = datacrs, cmap = colormap, vmin = vmin, vmax = vmax, \
        shading = 'auto')
    cbar = plt.colorbar(mesh,ax = pax, orientation='vertical',\
        fraction = 0.046, pad = 0.04)

    cbar.set_label(label,fontsize = labelsize, weight='bold')
    cbar.ax.tick_params(labelsize=labelticksize)

    if(gridlines): 
        pax.gridlines()

    pax.coastlines(resolution = '50m')
    #pax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())

    if(circle_bound):
        pax.set_boundary(circle, transform=pax.transAxes)

# Plot just a single day 
# ----------------------
def plotNSIDC_daily_figure(date_str, minlat = 65., \
        lat_circles = None, grid_data = False, zoom = False, \
        vmin = None, vmax = None, circle_bound = True, 
        title = '', label = None, \
        ax = None, gridlines = False, save = False):

    dt_date_str = datetime.strptime(date_str, "%Y%m%d")
    
    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    ##!#if(grid_data):
    ##!#    NSIDC_data = readgridCERES_hrly_grid(date_str[:10], param, \
    ##!#        satellite = 'Aqua', minlat = minlat)

    ##!#    if(NSIDC_data is None):
    ##!#        print("ERROR: no data returned from readgridCERES_hrly_grid")
    ##!#        print("Quitting")
    ##!#        return
    ##!#else:
    NSIDC_data = readNSIDC_daily(date_str, grid_data = grid_data)

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    in_ax = True 
    if(ax is None): 
        in_ax = False
        plt.close('all')
        fig1 = plt.figure(figsize = (6,6))
        mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = init_proj(date_str)
        #mapcrs = ccrs.LambertConformal(central_longitude = -100.)
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot data 
    # ----------------------------------------------------
    if(zoom):
        circle_bound = False
    plotNSIDC_daily(ax, NSIDC_data, minlat = minlat, \
        vmin = vmin, vmax = vmax, title = title, label = label, \
        circle_bound = circle_bound, \
        gridlines = gridlines, zoom = False)

    if(zoom):
        ax.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lon'][0], \
                       aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lon'][1], \
                       aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lat'][0], \
                       aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lat'][1]],\
                       datacrs)
##!#        ax0.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
##!#                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
##!#                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
##!#                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
##!#                        datacrs)
    else:
        ax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())

    #ax.set_title(date_str)

    # ----------------------------------------------------
    # If the user wants circles along latitude lines,    
    # plot them here      
    # ----------------------------------------------------
    if(lat_circles is not None):
        plot_lat_circles(ax, lat_circles) 

    if(not in_ax):
        fig1.tight_layout()
        
        if(save == True):
            out_name = 'nsidc_daily_conc_' + NSIDC_data['date'] + '.png'
            plt.savefig(out_name,dpi=300)
            print('Saved image '+out_name)
        else:
            plt.show()

# Writes a MODIS channel dictionary to HDF5 for Fortran colocation
# ----------------------------------------------------------------
def writeNSIDC_to_HDF5(NSIDC_data, save_path = './', minlat = 65., \
        remove_empty_scans = False):

    if(isinstance(NSIDC_data, str)):
        NSIDC_data = readNSIDC_daily(NSIDC_data)

    if(remove_empty_scans):
        mask_data = np.ma.masked_where(NSIDC_data['lat'] < minlat, \
            NSIDC_data['data'])
        mask_dims = np.array([ (False in mask_data[ii,:].mask) for ii in \
            range(mask_data.shape[0])])
        keep_idxs = np.where(mask_dims == True)[0]
    else:
        keep_idxs = None 

    # Convert the filename object to datetime
    # ---------------------------------------
    file_date = NSIDC_data['date']
    dt_date_str = datetime.strptime(file_date, '%Y%m%d')

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'nsidc_comp_' + file_date + '.hdf5'
    dset = h5py.File(outfile,'w')
 
    dset.create_dataset('latitude',  data = \
        NSIDC_data['lat'][keep_idxs,:].squeeze())
    dset.create_dataset('longitude', data = \
        NSIDC_data['lon'][keep_idxs,:].squeeze())
    dset.create_dataset('SeaIceConcentration', data = \
        NSIDC_data['data'][keep_idxs,:].data.squeeze())

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  

def writeNSIDC_toNCDF(NSIDC_data, save_path = './'):

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'nsidc_'+ NSIDC_data['date'] + '.nc'
    nc = Dataset(outfile,'w',format='NETCDF4')
  
    # Dimensions = lon, lat
    # Create the sizes of each dimension in the file. In this case,
    # the dimensions are "# of latitude" and "# of longitude"
    # -------------------------------------------------------------
    num_x = NSIDC_data['data'].shape[0]
    num_y = NSIDC_data['data'].shape[1]
    
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_x  = nc.createDimension('nx',num_x)
    n_y  = nc.createDimension('ny',num_y)

    # Create variables for the three dimensions. Note that since these
    # are variables, they are still given 'dimensions' using 'dlat', 
    # and 'dlon'. Latitude and longitude are each 2-d grids, so they are 
    # given 2 dimensions (dlat, dlon).
    # ------------------------------------------------------------------
    #X_DIM=nc.createVariable('x_dim','i4',('nx','ny'))
    #X_DIM.description='Across track'
    #X_DIM.units='km?'
    #Y_DIM=nc.createVariable('Latitude','i4',('nx','ny'))
    #Y_DIM.description='Along track'
    #Y_DIM.units='km?'

    # Create a variable for the AI data, dimensioned using both 
    # dimensions.
    # ---------------------------------------------------------  
    LAT = nc.createVariable('Latitude','f4',('nx','ny'))
    LAT.description='Latitude'
    LAT.units      ='Degrees'
    LON = nc.createVariable('Longitude','f4',('nx','ny'))
    LON.description='Longitude'
    LON.units      ='Degrees'

    CONC = nc.createVariable('SeaIceConcentration','f4',('nx','ny'))
    CONC.description='25 x 25 km daily sea ice concentration.'
    CONC.units = 'Percent'

    AREA = nc.createVariable('GridBoxArea','f4',('nx','ny'))
    AREA.description='Area of each grid box.'
    AREA.units = 'km2'

    # Fill in dimension variables one-by-one.
    # NOTE: not sure if you can insert the entire data array into the
    # dimension (for example, doing :
    #      LAT, LON = np.meshgrid(lat_ranges,lon_ranges) ),
    # so this could be
    # something for you to try. Might make this faster if it works
    # ---------------------------------------------------------------
    LON[:,:]      = NSIDC_data['lat'][:,:]
    LAT[:,:]      = NSIDC_data['lon'][:,:]
    CONC[:,:]     = NSIDC_data['data'][:,:]
    AREA[:,:]     = NSIDC_data['area'][:,:]

    # Save, write, and close the netCDF file
    # --------------------------------------
    nc.close()

    print("Saved file ",outfile)  


def write_Ice(NSIDC_data):
    print("yay")

    # Print of a file for each month
    num_months = int(len(NSIDC_data['titles']))

    # Set up the format string for printing off the data
    # format_list includes the formats for all the yearly average concentrations
    # Separate each by a space
    formats = '\t'.join(['{:6.4}' for item in NSIDC_data['titles']])
    # Combine the formats string with the other format string
    total_ice_format = '\t'.join(['{:.3}\t{:6.4}',formats,'\n'])

    # Set up the header line format string
    # Extract the individual years from the titles parameter
    string_dates = [tstring.strip()[3:9] for tstring in NSIDC_data['titles']]
    header_formats = '\t'.join(['{:6}' for item in NSIDC_data['titles']])
    total_header_format = '\t'.join(['{:6}\t{:6}',header_formats,'\n'])   
   
    filename = "nsidc_grid_ice_values.txt"
    #for ti in range(num_months):
    #    str_date = NSIDC_data['titles'][ti].strip()[3:9]
    #    filename = "nsidc_grid_ice_"+str_date+'.txt'
    with open(filename,'w') as fout:
        fout.write(total_header_format.format('Lat','Lon',*string_dates))
        for xi in range(len(NSIDC_data['grid_lat'])):
            print("Printing latitude ",NSIDC_data['grid_lat'][xi])
            for yj in range(len(NSIDC_data['grid_lon'])):
                fout.write(total_ice_format.format(NSIDC_data['grid_lat'][xi],\
                        NSIDC_data['grid_lon'][yj],*NSIDC_data['grid_ice_conc'][:,xi,yj]))

    print("Saved file",filename) 
    


# ice_trendCalc calculates the trends over the time period at each
# grid point on the 25x25 km grid.
def ice_trendCalc(NSIDC_data,thielSen=False):
    # Loop over the data and calculate trends
    for i in range(448):
        print(i)
        max_trend = -99.
        min_trend = 99.
        for j in range(304):
            NSIDC_data['trends'][i,j],temp_pcnt = trend_calc(NSIDC_data,i,j,thielsen=thielSen)
            NSIDC_data['land_trends'][i,j],temp_pcnt = trend_calc(NSIDC_data,i,j,thielsen=thielSen)
            temp_trend = NSIDC_data['trends'][i,j]
            if(temp_trend>max_trend):
                max_trend = temp_trend
            if(temp_trend<min_trend):
                min_trend = temp_trend
        print("  max trend = ",max_trend)
        print("  min trend = ",min_trend)
        # Deal with land masks
        good_indices = np.where(NSIDC_data['data'][0,i,:]<251.)
        land_indices = np.where(NSIDC_data['data'][0,i,:]>=251.)
        NSIDC_data['trends'][i,land_indices] = np.nan
        NSIDC_data['land_trends'][i,good_indices] = np.nan

    return NSIDC_data

def ice_gridtrendCalc(NSIDC_data,area=True,thielSen=False):
    if(area==True):
        print("\nCalculating area trends\n")
    else:
        print("\nCalculating % concentration trends\n")
    NSIDC_data['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = NSIDC_data['grid_lat']
    lon_ranges = NSIDC_data['grid_lon']
    ## Create array to hold monthly averages over region
    #initial_avgs = np.full(len(CERES_dict['dates']),-9999.)
    #initial_years  = np.zeros(len(initial_avgs))
    #initial_months = np.zeros(len(initial_avgs))
   
    # Loop over the lat and lon dimensions
    for xi in range(len(lat_ranges)):
        max_trend = -999
        min_trend = 999
        for yj in range(len(lon_ranges)):
            # Calculate the trend at the current box
            if(area==True):
                interp_data = (NSIDC_data['grid_ice_conc'][:,xi,yj]/100.)*NSIDC_data['grid_total_area'][xi,yj]
                good_indices = np.where(np.isnan(interp_data)==False)
            else:
                interp_data = NSIDC_data['grid_ice_conc'][:,xi,yj]
                good_indices = np.where(interp_data!=-99.)
            if(len(good_indices[0])==0):
                total_trend = -9999.
            else:
                interpx = np.arange(len(interp_data))
                #print(CERES_dict['data'][:,xi,yj][good_indices])
                total_interper = np.poly1d(np.polyfit( \
                    interpx[good_indices],\
                    interp_data[good_indices],1)) 
                # Normalize trend by dividing by number of years
                total_trend = (total_interper(interpx[-1])-\
                               total_interper(interpx[0]))
                if(total_trend>max_trend):
                    max_trend = total_trend
                if(total_trend<min_trend):
                    min_trend = total_trend
            if(area==True):
                NSIDC_data['grid_ice_area_trend'][xi,yj] = total_trend
        print(xi)
        #print("max trend = ",max_trend)
        #print("min trend = ",min_trend)

    return NSIDC_data

# grid_data_conc grids the 25x25 km gridded ice concentration data into
# a 1x1 degree lat/lon grid
def grid_data_conc(NSIDC_data, minlat = 65., maxlat = 90):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(minlat,maxlat,1.0)

    #lon_ranges += 0.5
    #lat_ranges += 0.5

    if(NSIDC_data['data'].shape[0] == 448):
        NSIDC_data['data'] = np.expand_dims(NSIDC_data['data'], 0)

    grid_ice_conc    = np.full((len(NSIDC_data['data'][:,0,0]),\
        len(lat_ranges),len(lon_ranges)),-99.)
    grid_land_conc    = np.full((len(NSIDC_data['data'][:,0,0]),\
        len(lat_ranges),len(lon_ranges)),np.nan)
    grid_ice_conc_cc = np.full((len(NSIDC_data['data'][:,0,0]),\
        len(lat_ranges),len(lon_ranges)),-999.)
    grid_pole_hole_cc = np.full((len(NSIDC_data['data'][:,0,0]),\
        len(lat_ranges),len(lon_ranges)),-999.)
    grid_unused_cc = np.full((len(NSIDC_data['data'][:,0,0]),\
        len(lat_ranges),len(lon_ranges)),-999.)
    grid_coastline_cc = np.full((len(NSIDC_data['data'][:,0,0]),\
        len(lat_ranges),len(lon_ranges)),-999.)
    grid_land_cc = np.full((len(NSIDC_data['data'][:,0,0]),\
        len(lat_ranges),len(lon_ranges)),-999.)
    grid_ice_area    = np.zeros((len(lat_ranges),len(lon_ranges)))
    grid_ice_area_trend    = np.full((len(lat_ranges),len(lon_ranges)),\
        -999.)
    print("Size of grid array: ",grid_ice_conc.shape)
    for nt in range(grid_ice_conc.shape[0]):
        print(nt)
        for xi in range(448):
            # Don't grid the data if any portion of the lat/lon box is over land.
            # Don't include land data
            for yj in range(304):
                if((NSIDC_data['lat'][xi,yj] > minlat) & (NSIDC_data['lat'][xi,yj] < maxlat)):
                    lat_index = np.where(np.floor(NSIDC_data['lat'][xi,yj])>=\
                        lat_ranges)[-1][-1]
                    lon_index = np.where(np.floor(NSIDC_data['lon'][xi,yj])>=\
                        lon_ranges)[-1][-1]
#                    if(np.floor(NSIDC_data['lon'][xi,yj]) == -136
                    # Add the current pixel area into the correct grid box, no
                    # matter if the current box is missing or not.
                    #if((lat_index==20) & (lon_index==10)):
                    #    print("Current grid area = ",grid_ice_area[lat_index,lon_index])
                    #if((lat_index==20) & (lon_index==10)):
                    #    print("    New grid area = ",grid_ice_area[lat_index,lon_index])
                    if(NSIDC_data['data'][nt,xi,yj] == 251):
                        if(grid_pole_hole_cc[nt,lat_index,lon_index] == -999.):
                            grid_pole_hole_cc[nt,lat_index,lon_index] = 1
                        else:
                            grid_pole_hole_cc[nt,lat_index,lon_index] += 1
                    elif(NSIDC_data['data'][nt,xi,yj] == 252):
                        if(grid_unused_cc[nt,lat_index,lon_index] == -999.):
                            grid_unused_cc[nt,lat_index,lon_index] = 1
                        else:
                            grid_unused_cc[nt,lat_index,lon_index] += 1
                    elif(NSIDC_data['data'][nt,xi,yj] == 253):
                        if(grid_coastline_cc[nt,lat_index,lon_index] == -999.):
                            grid_coastline_cc[nt,lat_index,lon_index] = 1
                        else:
                            grid_coastline_cc[nt,lat_index,lon_index] += 1
                    elif(NSIDC_data['data'][nt,xi,yj] == 254):
                        if(grid_land_cc[nt,lat_index,lon_index] == -999.):
                            grid_land_cc[nt,lat_index,lon_index] = 1
                        else:
                            grid_land_cc[nt,lat_index,lon_index] += 1
                    elif(NSIDC_data['data'][nt,xi,yj] < 251):
                        if(nt==0): 
                            grid_ice_area[lat_index,lon_index] += \
                            NSIDC_data['area'][xi,yj]
                        if(grid_ice_conc_cc[nt,lat_index,lon_index]==-999.):
                            grid_ice_conc[nt,lat_index,lon_index] = \
                                NSIDC_data['data'][nt,xi,yj]
                            grid_ice_conc_cc[nt,lat_index,lon_index] = 1.
                        else:
                            grid_ice_conc[nt,lat_index,lon_index] += \
                                NSIDC_data['data'][nt,xi,yj]
                            grid_ice_conc_cc[nt,lat_index,lon_index] += 1
                            #grid_ice_conc[nt,lat_index,lon_index] = \
                            #((grid_ice_conc[nt,lat_index,lon_index]*\
                            #grid_ice_conc_cc[nt,lat_index,lon_index])+\
                            #NSIDC_data['data'][nt,xi,yj])/\
                            #(grid_ice_conc_cc[nt,lat_index,lon_index]+1.)
                    else:
                        print("WARNING: SHOULD NOT FIND NSIDC PIXELS HERE")
                        #grid_land_conc[nt,lat_index,lon_index] = \
                        #    NSIDC_data['data'][nt,xi,yj]
                        #if(nt==0): grid_ice_area[lat_index,lon_index] = np.nan

                    # end else
                # end if good ice check
            # end y grid loop
        # end x grid loop
    # end time loop 
       
    NSIDC_data['data'] = NSIDC_data['data'].squeeze()
 
    # Calc averages here
    final_grid_ice_conc = np.copy(grid_ice_conc)
    final_grid_ice_conc[grid_ice_conc_cc > 0] = \
        grid_ice_conc[grid_ice_conc_cc > 0] / \
        grid_ice_conc_cc[grid_ice_conc_cc > 0].squeeze()

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

    final_pole_hole = np.ma.masked_invalid(final_pole_hole)
    final_unused    = np.ma.masked_invalid(final_unused)
    final_coastline = np.ma.masked_invalid(final_coastline)
    final_land      = np.ma.masked_invalid(final_land)
       
    #print("HERE:", final_pole_hole.compressed().shape, final_unused.compressed().shape, \
    #    final_coastline.compressed().shape, final_land.compressed().shape)
 
    # Calculate the areas of each grid box
    # ------------------------------------
    lat_ranges += 0.5
    lon_ranges += 0.5
    grid_areas = np.array([[lat_lon_area(tlat+0.5, tlat-0.5, tlon+0.5, tlon-0.5) \
        for tlon in lon_ranges] for tlat in lat_ranges])

    xx, yy = np.meshgrid(lon_ranges, lat_ranges)
    NSIDC_data['grid_ice_conc']       = final_grid_ice_conc.squeeze()
    NSIDC_data['grid_pole_hole']      = final_pole_hole
    NSIDC_data['grid_unused']         = final_unused
    NSIDC_data['grid_coastline']      = final_coastline
    NSIDC_data['grid_land']           = final_land
    NSIDC_data['grid_pole_hole_cc']   = grid_pole_hole_cc.squeeze()
    NSIDC_data['grid_unused_cc']      = grid_unused_cc.squeeze()
    NSIDC_data['grid_coastline_cc']   = grid_coastline_cc.squeeze()
    NSIDC_data['grid_land_cc']        = grid_land_cc.squeeze()
    NSIDC_data['grid_total_area']     = grid_ice_area.squeeze()
    NSIDC_data['grid_ice_area_trend'] = grid_ice_area_trend.squeeze()
    NSIDC_data['grid_lat'] = yy
    NSIDC_data['grid_lon'] = xx
    NSIDC_data['grid_area'] = grid_areas

    #print("HERE AGAIN:", NSIDC_data['grid_pole_hole'].compressed().shape, NSIDC_data['grid_unused'].compressed().shape, \
    #    NSIDC_data['grid_coastline'].compressed().shape, NSIDC_data['grid_land'].compressed().shape)
    
    return NSIDC_data

# grid_data averages the 25x25 km trends into a 1x1 degree grid
# The 'grid_ice' paramter in the dictionary therefore contains
# the gridded trends.
def grid_data(NSIDC_data):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_ice    = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    grid_ice_cc = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    for xi in range(448):
        print(xi)
        for yj in range(304):
            if(not np.isnan(NSIDC_data['trends'][xi,yj])):
                lat_index = np.where(np.floor(NSIDC_data['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(NSIDC_data['lon'][xi,yj])>=lon_ranges)[-1][-1]
                if(grid_ice_cc[lat_index,lon_index]==-999.):
                    grid_ice[lat_index,lon_index] = NSIDC_data['trends'][xi,yj]
                    grid_ice_cc[lat_index,lon_index] = 1.
                else:
                    grid_ice[lat_index,lon_index] = ((grid_ice[lat_index,lon_index]*grid_ice_cc[lat_index,lon_index])+\
                                                     NSIDC_data['trends'][xi,yj])/(grid_ice_cc[lat_index,lon_index]+1.)
                    grid_ice_cc[lat_index,lon_index]+=1
    NSIDC_data['grid_ice'] = grid_ice
    NSIDC_data['grid_ice_cc'] = grid_ice_cc
    NSIDC_data['grid_lat'] = lat_ranges
    NSIDC_data['grid_lon'] = lon_ranges
    return NSIDC_data

# plot_grid_data generates a plot of the /
def plot_grid_data(NSIDC_data,t_ind,pvar,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    if(pvar=='ice'):
        plabel = "Percent Ice Concentration"

    local_grid_ice = np.copy(NSIDC_data['grid_ice_conc'][t_ind,:,:])
    local_grid_ice_bad = np.copy(NSIDC_data['grid_ice_conc'][t_ind,:,:])

    local_grid_ice[local_grid_ice==-999.] = np.nan
    local_grid_ice_bad[local_grid_ice!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_ice)
    plot_land_data = ma.masked_invalid(local_grid_ice_bad)

    colormap = plt.cm.ocean
    #coolwarm = plt.cm.coolwarm
    #ax = plt.axes(projection=ccrs.NorthPolarStereo())
    file_adder=''
    if(adjusted==True):
        file_adder='_adjusted'
        fig1 = plt.figure(figsize=(8,5))
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=45.))
    else:
        fig1 = plt.figure()
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
    ax.gridlines()
    mesh = plt.pcolormesh(lon_ranges,lat_ranges,plot_good_data,\
            transform=ccrs.PlateCarree(),vmin=0,vmax=100,cmap=colormap)
    if(adjusted==True):
        ax.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        axins = inset_axes(ax,width="5%",height="100%",loc='lower left',\
                    bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax.transAxes,\
                    borderpad=0)
        cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=plabel)
        ax.set_xlim(-5170748.535086173,5167222.438879491)
        ax.set_ylim(-3913488.8763307533,3943353.899053069)
        #ender = '_adjusted'+ender
    else:
        plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
        cbar = plt.colorbar(mesh,cmap=colormap,label=plabel)
    ax.coastlines()


    #fig1 = plt.figure(figsize=(9,5))
    #plt.pcolormesh(plot_good_data,cmap=plt.cm.bwr,vmin=-50,vmax=50)
    #plt.colorbar(label='Percent Ice Concentration Trend (%)')
    #plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
    #plt.gca().invert_xaxis()
    #plt.gca().invert_yaxis()
    ax.set_title('Gridded NSIDC Sea Ice Concentration\n'+NSIDC_data['titles'][t_ind])
    if(save==True):
        outname = 'ice_conc_gridded_200012_201812'+NSIDC_data['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
    

def plot_grid_trend(NSIDC_data,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    local_grid_ice = np.copy(NSIDC_data['grid_ice'])
    local_grid_ice_bad = np.copy(NSIDC_data['grid_ice'])

    local_grid_ice[local_grid_ice==-999.] = np.nan
    local_grid_ice_bad[local_grid_ice!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_ice)
    plot_land_data = ma.masked_invalid(local_grid_ice_bad)

    colormap = plt.cm.bwr
    #coolwarm = plt.cm.coolwarm
    #ax = plt.axes(projection=ccrs.NorthPolarStereo())
    file_adder=''
    if(adjusted==True):
        file_adder='_adjusted'
        fig1 = plt.figure(figsize=(8,5))
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=45.))
    else:
        fig1 = plt.figure()
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
    ax.gridlines()
    mesh = plt.pcolormesh(lon_ranges,lat_ranges,plot_good_data,\
            transform=ccrs.PlateCarree(),vmin=-50,vmax=50,cmap=colormap)
    if(adjusted==True):
        ax.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        axins = inset_axes(ax,width="5%",height="100%",loc='lower left',\
                    bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax.transAxes,\
                    borderpad=0)
        cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label='Ice Concentration Trend [%]')
        ax.set_xlim(-5170748.535086173,5167222.438879491)
        ax.set_ylim(-3913488.8763307533,3943353.899053069)
        #ender = '_adjusted'+ender
    else:
        plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
        cbar = plt.colorbar(mesh,cmap=colormap,label='Ice Concentration Trend [%]')
    ax.coastlines()


    #fig1 = plt.figure(figsize=(9,5))
    #plt.pcolormesh(plot_good_data,cmap=plt.cm.bwr,vmin=-50,vmax=50)
    #plt.colorbar(label='Percent Ice Concentration Trend (%)')
    #plt.pcolormesh(plot_land_data,cmap=plt.cm.Greys,vmin=-10,vmax=10)
    #plt.gca().invert_xaxis()
    #plt.gca().invert_yaxis()
    if(NSIDC_data['season_adder']!=''):
        ax.set_title('NSIDC Sea Ice Concentration'+NSIDC_data['season_adder'].title()+\
            ' Trends\nJan 2001 to Dec 2018')
    else:
        ax.set_title('NSIDC Sea Ice Concentration Trends\nJan 2001 to Dec 2018')
    if(save==True):
        outname = 'ice_trend_gridded_200101_201812'+NSIDC_data['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
    
def plot_grid_time_series(NSIDC_data,lat_ind,lon_ind,thielsen=False):
    inseason = NSIDC_data['season_adder'].strip()
    temp_data = NSIDC_data['grid_ice_conc'][:,lat_ind,lon_ind]
    interpx = np.arange(len(temp_data))
    #interpx = years
    ##interper = np.poly1d(np.polyfit(interpx,temp_data,1)) 
    ### Normalize trend by dividing by number of years
    ##trend = (interper(interpx[-1])-interper(interpx[0]))

    slope, intercept, r_value, p_value, std_err = stats.linregress(interpx,temp_data)
    ##print(slope/len(test_dict[dictkey].keys())
    regress_y = interpx*slope+intercept
    trend = regress_y[-1] - regress_y[0]

    if(thielsen==True):
        #The slope
        S=0
        sm=0
        nx = len(temp_data)
        num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
        Sn=np.zeros(num_d)
        for si in range(0,nx-1):
            for sj in range(si+1,nx):
                # Find the slope between the two points
                Sn[sm] = (temp_data[si]-temp_data[sj])/(si-sj) 
                sm=sm+1
            # Endfor
        # Endfor
        Snsorted=sorted(Sn)
        sm=int(num_d/2.)
        if(2*sm    == num_d):
            tsslope=0.5*(Snsorted[sm]+Snsorted[sm+1])
        if(2*sm+1 == num_d): 
            tsslope=Snsorted[sm+1]
        regress_ts = interpx*tsslope+intercept
        trend = regress_ts[-1]-regress_ts[0]

    fig1 = plt.figure() 
    plt.title('Gridded Ice Data: '+str(NSIDC_data['grid_lat'][lat_ind])+'x'+\
                str(NSIDC_data['grid_lon'][lon_ind])+'\n'+inseason.title()+\
                ' season of each year)')
    plt.plot(temp_data,label='observations') 
    plt.plot(regress_y,'--',label='trend')
    plt.ylabel('Ice Concentration [%]')
    plt.legend()
    outname = 'nsidc_grid_time_series_'+inseason+'_'+str(int(NSIDC_data['grid_lat'][lat_ind]))+'x'+\
                str(int(NSIDC_data['grid_lon'][lon_ind]))+'.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)   
    plt.show()

# This code replicates Jianglong's figure showing how the April - September
# ice concentration values have changed over the time period
# NOTE: To run, do something like
# >>> all_NSIDC_data = read_ice('all')
# >>> all_NSIDC_data = grid_data_conc(all_NSIDC_data)
# and then
# >>> plot_apr_sep_changes(all_NSIDC_data)
def plot_apr_sep_changes(NSIDC_data):
  
    inseason = NSIDC_data['season_adder'].strip()
 
    upper_vals = np.array([100,80,60,40])
    lower_vals = np.array([80,60,40,20])

    # Generate 3-year averages for 3 time periods:
    # 2001-2003
    # 2008-2011
    # 2016-2018
   
    # For now, look at April - September 
    if(inseason=='sunlight'):
        num_months = 6
        start_idx  = 0
    else:
        num_months = 12
        start_idx  = 4
    # Dimensions of ice_avgs are
    # 0 - year blocks (size = 3, see above)
    # 1 - months      (size = num_months)
    # 2 - ice values  (size = 4)
    ice_avgs = np.zeros((3,num_months,len(upper_vals)))
    indices = [0,8,15]
  
    # Base the locations for each ice concentration range on the average
    # April concentration between 2001 and 2003

    avg_Apr_old = \
        np.average(NSIDC_data['grid_ice_conc'][start_idx::num_months,:,:][:3,:,:],axis=0)
 
    # Use the locations during the first April average to base everything
    # on
    locations_80_100 = np.where((avg_Apr_old <= 100.) & (avg_Apr_old > 80.)) 
    locations_60_80  = np.where((avg_Apr_old <= 80.)  & (avg_Apr_old > 60.)) 
    locations_40_60  = np.where((avg_Apr_old <= 60.)  & (avg_Apr_old > 40.)) 
    locations_20_40  = np.where((avg_Apr_old <= 40.)  & (avg_Apr_old > 20.)) 

    # Fill the ice_avgs array with the data for each time period 
    for ri in range(len(indices)):
        for mi in range(num_months):
            print("Year block = ",ri,"  Month = ",mi)
            # Deal with 80 to 100
            ##temp_arr = np.zeros((3,len(locations_80_100)))
            ##temp_arr[0] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri]][locations_80_100]
            ##temp_arr[1] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 1][locations_80_100]
            ##temp_arr[2] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 2][locations_80_100]
            #ice_avgs[ri,mi,0] = np.average(temp_arr)
            ice_avgs[ri,mi,0] = \
                np.average(np.average(NSIDC_data['grid_ice_conc'][mi::num_months,:,:]
                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_80_100])

            ### Deal with 60 to 80
            ##temp_arr = np.zeros((3,len(locations_60_80)))
            ##temp_arr[0] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri]][locations_60_80]
            ##temp_arr[1] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 1][locations_60_80]
            ##temp_arr[2] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 2][locations_60_80]
            ##ice_avgs[ri,mi,1] = np.average(temp_arr)
            ice_avgs[ri,mi,1] = \
                np.average(np.average(NSIDC_data['grid_ice_conc'][mi::num_months,:,:]
                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_60_80])

            # Deal with 40 to 60
            ##temp_arr = np.zeros((3,len(locations_40_60)))
            ##temp_arr[0] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri]][locations_40_60]
            ##temp_arr[1] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 1][locations_40_60]
            ##temp_arr[2] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 2][locations_40_60]
            #ice_avgs[ri,mi,2] = np.average(temp_arr)
            ice_avgs[ri,mi,2] = \
                np.average(np.average(NSIDC_data['grid_ice_conc'][mi::num_months,:,:]
                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_40_60])

            # Deal with 40 to 60
            ##temp_arr = np.zeros((3,len(locations_20_40)))
            ##temp_arr[0] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri]][locations_20_40]
            ##temp_arr[1] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 1][locations_20_40]
            ##temp_arr[2] = NSIDC_data['grid_ice_conc'][mi::num_months][indices[ri] + 2][locations_20_40]
            ##ice_avgs[ri,mi,3] = np.average(temp_arr)
            ice_avgs[ri,mi,3] = \
                np.average(np.average(NSIDC_data['grid_ice_conc'][mi::num_months,:,:]
                    [indices[ri]:indices[ri] + 3,:,:],axis = 0)[locations_20_40])

            #ice_avgs[ri,mi] = \
                #np.average(NSIDC_data['grid_ice_conc'][0::num_months]
                #    [indices[ri]:indices[ri] + 3],axis = 0)

    # Generate the figure
    plt.close()
    fig1 = plt.figure(figsize=(8,6))
    ax = plt.subplot()
    if(inseason=='sunlight'):
        months = ['Apr','May','June','July','Aug','Sep']
    else:
        months = ['Dec','Jan','Feb','Mar','Apr','May',\
                  'June','July','Aug','Sep','Oct','Nov']
    # Plot the 2001 - 2003 data
    ax.plot(ice_avgs[0,:,0],color='black')
    ax.plot(ice_avgs[0,:,1],color='tab:blue')
    ax.plot(ice_avgs[0,:,2],color='tab:green')
    ax.plot(ice_avgs[0,:,3],color='tab:red')
    # Plot the 2009 - 2011 data
    ax.plot(ice_avgs[1,:,0],'--',color='black')
    ax.plot(ice_avgs[1,:,1],'--',color='tab:blue')
    ax.plot(ice_avgs[1,:,2],'--',color='tab:green')
    ax.plot(ice_avgs[1,:,3],'--',color='tab:red')
    # Plot the 2016 - 2018 data
    ax.plot(ice_avgs[2,:,0],linestyle='dotted',color='black')
    ax.plot(ice_avgs[2,:,1],linestyle='dotted',color='tab:blue')
    ax.plot(ice_avgs[2,:,2],linestyle='dotted',color='tab:green')
    ax.plot(ice_avgs[2,:,3],linestyle='dotted',color='tab:red')
    ax.set_xticks(np.arange(num_months),months)
    ax.set_ylabel('Ice Concentration [%]')

    # Shrink the current axis's height by 10% to make room for the legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,\
                    box.width, box.height * 0.9])

    # Make the legend
    custom_lines = [Line2D([0],[0],color='black'),\
                    Line2D([0],[0],color='black'),\
                    Line2D([0],[0],color='tab:blue'),\
                    Line2D([0],[0],linestyle='dashed',color='black'),\
                    Line2D([0],[0],color='tab:green'),\
                    Line2D([0],[0],linestyle='dotted',color='black'),\
                    Line2D([0],[0],color='tab:red')]
    ax.legend(custom_lines,['80 - 100%',\
                             '2001 - 2003',\
                             '60 - 80%',\
                             '2009 - 2011',\
                             '40 - 60%',\
                             '2016 - 2018',\
                             '20 - 40%'],\
              loc = 'upper center',bbox_to_anchor=(0.5,-0.05),\
              fancybox=True,shadow=True,ncol=4)
    plt.show()
    plt.close()

def calc_NSIDC_sfc_type_coverages(NSIDC_data, max_ice_for_ocean = 0, \
        min_ice_for_ice = 80, use_area = True, use_pcnt = True):

    total_array     = np.full(len(NSIDC_data['dates']), np.nan)
    other_array     = np.full(len(NSIDC_data['dates']), np.nan)
    coastline_array = np.full(len(NSIDC_data['dates']), np.nan)
    land_array      = np.full(len(NSIDC_data['dates']), np.nan)
    ocean_array     = np.full(len(NSIDC_data['dates']), np.nan)
    ice_array       = np.full(len(NSIDC_data['dates']), np.nan)
    mix_array       = np.full(len(NSIDC_data['dates']), np.nan)
    
    for tidx in range(len(NSIDC_data['dates'])):
    
        # Figure out the total number of grid boxes here. 
        # Land pixels also include pixels that are nonmissing in both land and ice data
        type_dict = calc_NSIDC_sfc_types(NSIDC_data, tidx, \
            max_ice_for_ocean = max_ice_for_ocean, \
            min_ice_for_ice   = min_ice_for_ice, \
            use_grid_data = use_area
            )
     
        if(use_area):
            total_array[tidx]    = sum(NSIDC_data['grid_area'][type_dict['total_idxs']])
            other_array[tidx]    = sum(NSIDC_data['grid_area'][type_dict['other_idxs']])
            land_array[tidx]     = sum(NSIDC_data['grid_area'][type_dict['land_idxs']])
            ocean_array[tidx]    = sum(NSIDC_data['grid_area'][type_dict['ocean_idxs']])
            ice_array[tidx]      = sum(NSIDC_data['grid_area'][type_dict['ice_idxs']])
            mix_array[tidx]      = sum(NSIDC_data['grid_area'][type_dict['mix_idxs']])
    
            calc_total = sum(NSIDC_data['grid_area'][type_dict['other_idxs']]) + \
                         sum(NSIDC_data['grid_area'][type_dict['land_idxs']]) + \
                         sum(NSIDC_data['grid_area'][type_dict['ocean_idxs']]) + \
                         sum(NSIDC_data['grid_area'][type_dict['ice_idxs']]) + \
                         sum(NSIDC_data['grid_area'][type_dict['mix_idxs']])
    
        else: 
            total_array[tidx]    = type_dict['total_pixels']
            other_array[tidx]    = type_dict['num_other']
            land_array[tidx]     = type_dict['num_land']
            ocean_array[tidx]    = type_dict['num_ocean_only']
            ice_array[tidx]      = type_dict['num_ice_only']
            mix_array[tidx]      = type_dict['num_mix_only']
    
            calc_total = type_dict['num_other'] + \
                         type_dict['num_land'] + \
                         type_dict['num_ocean_only'] + \
                         type_dict['num_ice_only'] + \
                         type_dict['num_mix_only']
    
        print(NSIDC_data['dates'][tidx], type_dict['total_pixels'], calc_total, type_dict['num_other'] + \
              type_dict['num_land'], type_dict['num_ocean_only'],type_dict['num_ice_only'],type_dict['num_mix_only'])
   
    if(use_pcnt): 
        pcnt_other    = ((other_array   / total_array) * 100.)
        pcnt_land     = ((land_array    / total_array) * 100.)
        pcnt_ocean    = ((ocean_array   / total_array) * 100.)
        pcnt_ice      = ((ice_array     / total_array) * 100.)
        pcnt_mix      = ((mix_array     / total_array) * 100.)
    else:
        pcnt_other = other_array
        pcnt_land  = land_array    
        pcnt_ocean = ocean_array   
        pcnt_ice   = ice_array     
        pcnt_mix   = mix_array     

    out_dict = {}
    out_dict['pcnt_other'] = pcnt_other    
    out_dict['pcnt_land']     = pcnt_land    
    out_dict['pcnt_ocean']    = pcnt_ocean
    out_dict['pcnt_ice']      = pcnt_ice
    out_dict['pcnt_mix']      = pcnt_mix
    out_dict['use_area']      = use_area
    out_dict['use_pcnt']      = use_pcnt
    
    return out_dict

def calc_NSIDC_sfc_types(NSIDC_data, tidx, max_ice_for_ocean = 0, \
        min_ice_for_ice = 80, use_grid_data = False):

    if(use_grid_data):
        total_idx = np.where(\
            (NSIDC_data['grid_pole_hole'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_unused'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_coastline'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_land'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))
        total_pixels = np.where(\
            (NSIDC_data['grid_pole_hole'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_unused'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_coastline'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_land'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))[0].shape[0]
        other_pixels = np.where(\
            (NSIDC_data['grid_pole_hole'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_unused'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_coastline'][tidx,:,:].mask == False))
            #(NSIDC_data['grid_land'][tidx,:,:].mask == False) & \
            #(NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))
        coastline_pixels = np.where(NSIDC_data['grid_coastline'][tidx,:,:].mask == False)
        num_other = other_pixels[0].shape[0]
        # The land only pixels are just "land"
        pole_unused = np.where(\
            ((NSIDC_data['grid_pole_hole'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_unused'][tidx,:,:].mask == False)) & \
            (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == True))
        num_pole_unused = pole_unused[0].shape[0]
        land_only = np.where(\
            ((NSIDC_data['grid_coastline'][tidx,:,:].mask == True) & \
             (NSIDC_data['grid_land'][tidx,:,:].mask == False)) & \
            (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == True))
        num_land = land_only[0].shape[0]
        oceanice_only = np.where(\
            ((NSIDC_data['grid_land'][tidx,:,:].mask == True) & \
            (NSIDC_data['grid_coastline'][tidx,:,:].mask == True)) & \
            (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))
        ocean_only = np.where( \
            ((((NSIDC_data['grid_land'][tidx,:,:].mask == True) & \
              (NSIDC_data['grid_coastline'][tidx,:,:].mask == True)) & \
            (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))) & 
                               (NSIDC_data['grid_ice_conc'][tidx,:,:] <= max_ice_for_ocean))
        num_ocean_only = ocean_only[0].shape[0]
        ice_only   = np.where( \
            ((((NSIDC_data['grid_land'][tidx,:,:].mask == True) & \
              (NSIDC_data['grid_coastline'][tidx,:,:].mask == True)) & \
            (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))) & 
                               (NSIDC_data['grid_ice_conc'][tidx,:,:] >= min_ice_for_ice))
        num_ice_only = ice_only[0].shape[0]
        mix_only   = np.where( \
            ((((NSIDC_data['grid_land'][tidx,:,:].mask == True) & \
              (NSIDC_data['grid_coastline'][tidx,:,:].mask == True)) & \
            (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))) & 
                               ((NSIDC_data['grid_ice_conc'][tidx,:,:] > max_ice_for_ocean) &
                               (NSIDC_data['grid_ice_conc'][tidx,:,:] < min_ice_for_ice))\
            )
        num_mix_only = mix_only[0].shape[0]

    else:
        total_pixels = np.where(\
            (NSIDC_data['pole_hole'][tidx,:,:].mask == False) | \
            (NSIDC_data['grid_unused'][tidx,:,:].mask == False) | \
            (NSIDC_data['coastline'][tidx,:,:].mask == False) | \
            (NSIDC_data['land'][tidx,:,:].mask == False) | \
            (NSIDC_data['data'][tidx,:,:].mask == False))[0].shape[0]
        other_pixels = np.where(\
            (NSIDC_data['pole_hole'][tidx,:,:].mask == False) & \
            (NSIDC_data['grid_unused'][tidx,:,:].mask == False) & \
            (NSIDC_data['coastline'][tidx,:,:].mask == False))
            #(NSIDC_data['land'][tidx,:,:].mask == False) & \
            #(NSIDC_data['data'][tidx,:,:].mask == False))
        num_other = other_pixels[0].shape[0]
        coastline_pixels = np.where(NSIDC_data['coastline'][tidx,:,:].mask == False)
        # The land only pixels are either "land" or "coastline"
        pole_unused = np.where(\
            ((NSIDC_data['pole_hole'][tidx,:,:].mask == False) | \
            (NSIDC_data['unused'][tidx,:,:].mask == False)) & \
            (NSIDC_data['data'][tidx,:,:].mask == True))
        num_pole_unused = pole_unused[0].shape[0]
        land_only = np.where(\
            #((NSIDC_data['coastline'][tidx,:,:].mask == False) | \
            ((NSIDC_data['land'][tidx,:,:].mask == False)) & \
            (NSIDC_data['data'][tidx,:,:].mask == True))
        num_land = land_only[0].shape[0]
        oceanice_only = np.where(\
            (NSIDC_data['land'][tidx,:,:].mask == True) & \
            (NSIDC_data['data'][tidx,:,:].mask == False))
        ocean_only = np.where( \
            (((NSIDC_data['land'][tidx,:,:].mask == True) & \
            (NSIDC_data['data'][tidx,:,:].mask == False))) & 
                               (NSIDC_data['data'][tidx,:,:] <= max_ice_for_ocean))
        num_ocean_only = ocean_only[0].shape[0]
        ice_only   = np.where( \
            (((NSIDC_data['land'][tidx,:,:].mask == True) & \
            (NSIDC_data['data'][tidx,:,:].mask == False))) & 
                               (NSIDC_data['data'][tidx,:,:] >= min_ice_for_ice))
        num_ice_only = ice_only[0].shape[0]
        mix_only   = np.where( \
            (((NSIDC_data['land'][tidx,:,:].mask == True) & \
            (NSIDC_data['data'][tidx,:,:].mask == False))) & 
                               ((NSIDC_data['data'][tidx,:,:] > max_ice_for_ocean) &
                                (NSIDC_data['data'][tidx,:,:] < min_ice_for_ice))\
            )
        num_mix_only = mix_only[0].shape[0]

    out_dict = {}
    out_dict['total_idxs']    = total_idx
    out_dict['total_pixels']  = total_pixels
    out_dict['other_idxs']    = other_pixels
    out_dict['num_other']     = num_other
    out_dict['coastline_idxs'] = coastline_pixels
    out_dict['pole_idxs']     = pole_unused
    out_dict['num_pole']      = num_pole_unused
    out_dict['land_idxs']     = land_only
    out_dict['num_land']      = num_land
    out_dict['ocean_idxs']    = ocean_only
    out_dict['num_ocean_only']     = num_ocean_only
    out_dict['ice_idxs']      = ice_only
    out_dict['num_ice_only']       = num_ice_only
    out_dict['mix_idxs']      = mix_only
    out_dict['num_mix_only']       = num_mix_only

    return out_dict

# Plots the surface types for a given month
# NSIDC_data: monthly data object
def plot_NSIDC_month_sfc_types(NSIDC_data, date_str, ax = None, save = False, \
        use_grid_data = False, max_ice_for_ocean = 0, min_ice_for_ice = 80):

    # Figure out the matching index
    # ------------------------------
    tidx = np.where(np.array(NSIDC_data['dates']) == date_str)[0][0]

    type_dict = calc_NSIDC_sfc_types(NSIDC_data, tidx, \
        use_grid_data = use_grid_data, max_ice_for_ocean = max_ice_for_ocean, \
        min_ice_for_ice = min_ice_for_ice)

    var_add = ''
    if(use_grid_data):
        var_add = 'grid_'
    ones   = np.full(NSIDC_data[var_add + 'lon'].shape, np.nan)
    twos   = np.full(NSIDC_data[var_add + 'lon'].shape, np.nan)
    threes = np.full(NSIDC_data[var_add + 'lon'].shape, np.nan)
    fours  = np.full(NSIDC_data[var_add + 'lon'].shape, np.nan)
    fives  = np.full(NSIDC_data[var_add + 'lon'].shape, np.nan)
    sixs   = np.full(NSIDC_data[var_add + 'lon'].shape, np.nan)
    sevens   = np.full(NSIDC_data[var_add + 'lon'].shape, np.nan)
  
    ones[type_dict['other_idxs']] = 1.
    twos[type_dict['land_idxs']]     = 2. 
    threes[type_dict['ocean_idxs']]  = 3. 
    fours[type_dict['ice_idxs']]     = 4. 
    fives[type_dict['mix_idxs']]     = 5.
    sixs[type_dict['pole_idxs']]     = 6.
    sevens[type_dict['coastline_idxs']]     = 7.
 
    #ones = np.ma.masked_where(type_dict['combined_idxs'].mask == False, \
    #    ones)
    #twos = np.ma.masked_where(type_dict['land_idxs'].mask == False, \
    #    twos)
    #three = np.ma.masked_where(type_dict['ocean_idxs'].mask == False, \
    #    threes)
    #fours = np.ma.masked_where(type_dict['ice_idxs'].mask == False, \
    #    fours)
    #fives = np.ma.masked_where(type_dict['mix_idxs'].mask == False, \
    #    fives)

    ones   = np.ma.masked_invalid(ones)
    twos   = np.ma.masked_invalid(twos)
    threes = np.ma.masked_invalid(threes)
    fours  = np.ma.masked_invalid(fours)
    fives  = np.ma.masked_invalid(fives)
    sixs   = np.ma.masked_invalid(sixs)
    #ones = np.ma.masked_where(mask_ice.mask == True, base_ones).T


    in_ax = True 
    if(ax is None): 
        in_ax = False
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1,projection = mapcrs)

    ax.pcolormesh(NSIDC_data[var_add + 'lon'], NSIDC_data[var_add + 'lat'], \
        ones, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10', label = 'Coastline')
    ax.pcolormesh(NSIDC_data[var_add + 'lon'], NSIDC_data[var_add + 'lat'], \
        twos, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10', label = 'Land')
    ax.pcolormesh(NSIDC_data[var_add + 'lon'], NSIDC_data[var_add + 'lat'], \
        threes, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10', label = 'Ocean')
    ax.pcolormesh(NSIDC_data[var_add + 'lon'], NSIDC_data[var_add + 'lat'], \
        fours, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10', label = 'Ice')
    ax.pcolormesh(NSIDC_data[var_add + 'lon'], NSIDC_data[var_add + 'lat'], \
        fives, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10', label = 'Mix')
    ax.pcolormesh(NSIDC_data[var_add + 'lon'], NSIDC_data[var_add + 'lat'], \
        sixs, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10', label = 'Pole')
    ax.pcolormesh(NSIDC_data[var_add + 'lon'], NSIDC_data[var_add + 'lat'], \
        sevens, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10', label = 'Pole')

    ax.coastlines()
    ax.set_extent([-180, 180, 65, 90], datacrs)
    
    if(not in_ax):
        fig.tight_layout()
        plt.show()


def plot_NSIDC_sfc_type_change_bar(NSIDC_data, max_ice_for_ocean = 0, \
        min_ice_for_ice = 80, use_area = True, use_pcnt = True, save = False):


    pcnt_coverages =  calc_NSIDC_sfc_type_coverages(NSIDC_data, \
            max_ice_for_ocean = max_ice_for_ocean, \
            min_ice_for_ice = min_ice_for_ice, \
            use_area = use_area, use_pcnt = True)
    
    pcnt_other    = pcnt_coverages['pcnt_other']
    pcnt_land     = pcnt_coverages['pcnt_land']
    pcnt_ocean    = pcnt_coverages['pcnt_ocean']
    pcnt_ice      = pcnt_coverages['pcnt_ice']
    pcnt_mix      = pcnt_coverages['pcnt_mix']
    
    int_years = np.array([int(tyear[:4]) for tyear in NSIDC_data['dates']])
    xvals = int_years[::6]

    if(use_area):
        ylabel = 'Percent Area Coverage [%]'
    else:
        ylabel = 'Percent Grid Box [%]'
  
    fig = plt.figure(figsize = (9, 6))
    axs = fig.subplots(2,3)
    flat_ax = axs.flatten()
    for ii, ax in enumerate(flat_ax):
        #ax = fig.add_subplot(1,1,1)
        ax.bar(xvals, pcnt_other[ii::6], label = 'Other', color = 'tab:purple')
        ax.bar(xvals, pcnt_land[ii::6],  bottom = pcnt_other[ii::6], label = 'Land', color = 'tab:red')
        ax.bar(xvals, pcnt_ocean[ii::6], bottom = pcnt_other[ii::6] + pcnt_land[ii::6], label = 'Ocean', color = 'tab:green')
        ax.bar(xvals, pcnt_mix[ii::6],   bottom = pcnt_other[ii::6] + pcnt_land[ii::6] + pcnt_ocean[ii::6], \
            label = 'Mix Ice/Ocn', color = 'tab:orange')
        ax.bar(xvals, pcnt_ice[ii::6],   bottom = pcnt_other[ii::6] + pcnt_land[ii::6] + pcnt_ocean[ii::6] +\
             pcnt_mix[ii::6], label = 'Ice', color = 'tab:blue')
        local_datetime = datetime(2000,4+ii,1)
        ax.set_title(local_datetime.strftime('%B'))
        ax.set_ylabel(ylabel)

    #axs[0,0].legend()
    custom_lines = [Line2D([0],[0],color='tab:blue'),\
                    Line2D([0],[0],color='tab:orange'),\
                    Line2D([0],[0],color='tab:green'),\
                    Line2D([0],[0],color='tab:red'),\
                    Line2D([0],[0],color='tab:purple')]
    fig.legend(custom_lines,['Ice',\
                             'Mix Ice/Ocn',\
                             'Ocean',\
                             'Land',\
                             'Other'],\
              loc = 'lower center',\
              fancybox=False,shadow=False,ncol=5)
    #fig.legend(\
    #          loc = 'lower center',\
    #          fancybox=True,shadow=True,ncol=5)


    if(use_area):
        plt.suptitle('Percent Area Coverage Per Surface Type\n' + \
            str(int(NSIDC_data['minlat'])) + '${^o}$ N - 87${^o}$ N\n2005 - 2020')
    else: 
        plt.suptitle('Percent 1$^{0}$ Grid Box Coverage Per Surface Type\n2005 - 2020')

    fig.tight_layout()

    plt.subplots_adjust(bottom = 0.11)

    if(save):
        outname = 'nsidc_sfc_type_change_bar.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_NSIDC_sfc_type_change_line(NSIDC_data, max_ice_for_ocean = 0, \
        min_ice_for_ice = 80, use_area = True, use_pcnt = True, save = False):

    pcnt_coverages =  calc_NSIDC_sfc_type_coverages(NSIDC_data, \
            max_ice_for_ocean = max_ice_for_ocean, \
            min_ice_for_ice = min_ice_for_ice, \
            use_area = use_area, use_pcnt = use_pcnt)
    
    pcnt_other    = pcnt_coverages['pcnt_other']
    pcnt_land     = pcnt_coverages['pcnt_land']
    pcnt_ocean    = pcnt_coverages['pcnt_ocean']
    pcnt_ice      = pcnt_coverages['pcnt_ice']
    pcnt_mix      = pcnt_coverages['pcnt_mix']
    
    int_years = np.array([int(tyear[:4]) for tyear in NSIDC_data['dates']])
    xvals = int_years[::6]
    
    fig2 = plt.figure(figsize = (9, 6))
    axs = fig2.subplots(2,3)
    flat_ax2 = axs.flatten()
    
    if(use_area):
        ylabel = 'Percent Area Coverage [%]'
    else:
        ylabel = 'Percent Grid Box [%]'
    
    for ii, ax in enumerate(flat_ax2):
        #ax.plot(xvals, pcnt_combined[ii::6], label = 'Combined') 
        #ax.plot(xvals, pcnt_land[ii::6], label = 'Land') 
        ax.plot(xvals, pcnt_ocean[ii::6], label = 'Ocean') 
        #ax.plot(xvals, pcnt_ice[ii::6], label = 'Ice') 
        #ax.plot(xvals, pcnt_mix[ii::6], label = 'Mix Ice/Ocn') 
        ax.plot(xvals, pcnt_mix[ii::6] + pcnt_ice[ii::6], label = 'Mix + Ice') 
        local_datetime = datetime(2000,4+ii,1)
        ax.set_title(local_datetime.strftime('%B'))
        
        first_out = plot_trend_line(ax, xvals, pcnt_ocean[ii::6], color='tab:green', linestyle = '--', \
            slope = 'thiel-sen')
        second_out = plot_trend_line(ax, xvals, pcnt_mix[ii::6] + pcnt_ice[ii::6], color='tab:brown', linestyle = '--', \
            slope = 'thiel-sen')
    
        first_slope = np.round(first_out[0], 3)
        first_change = np.round(first_slope * len(xvals), 1)
        second_slope  = np.round(second_out[0], 3)
        second_change = np.round(second_slope * len(xvals), 1)
    
        print("    Ocean:   {0}x, total = {1}".format(first_slope, first_change))
        print("    Mix+Ice: {0}x, total = {1}".format(second_slope, second_change))
    
        if(ii == 5):
            location = 'upper_right'
            align = 'right'
        else:
            location = 'mid_left'
            align = 'left'
        plot_figure_text(ax, 'ΔOcean = {0}%\nΔIce = {1}%'.format(first_change, second_change), \
            location = location, halign = align, fontsize = 10, weight = 'normal')
   
        ax.set_ylabel(ylabel)
 
    axs[0,0].legend()
    
    if(use_area):
        plt.suptitle('Change in Ocean vs. Ice Percent Area Coverage\n' + \
            str(int(NSIDC_data['minlat'])) + '${^o}$ N - 87${^o}$ N\n2005 - 2020')
    else: 
        plt.suptitle('Change in Ocean vs. Ice Percent 1$^{0}$ Grid Box Coverage\n2005 - 2020')
    
    fig2.tight_layout()
    
    if(save):
        outname = 'nsidc_sfc_type_change_line.png'
        fig2.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def check_type_change(NSIDC_data, month_idx, xidx, yidx, \
        max_ice_for_ocean = 20, min_ice_for_ice = 80, bin_size = 5, \
        plot_line = False, debug = False):

    data_before = NSIDC_data['grid_ice_conc'][month_idx::6,xidx,yidx][:bin_size]
    data_after  = NSIDC_data['grid_ice_conc'][month_idx::6,xidx,yidx][-bin_size:]

    # Check for other surface types
    # -----------------------------
    if(len(NSIDC_data[\
        'grid_pole_hole'][month_idx::6,xidx,yidx][:bin_size].compressed()\
         != 0)):
        type_str = 'Pole Hole' 
    elif(len(NSIDC_data[\
        'grid_unused'][month_idx::6,xidx,yidx][:bin_size].compressed()\
         != 0)):
        type_str = 'Unused' 
    elif(len(NSIDC_data[\
        'grid_coastline'][month_idx::6,xidx,yidx][:bin_size].compressed()\
         != 0)):
        type_str = 'Coastline' 
    elif(len(NSIDC_data[\
        'grid_land'][month_idx::6,xidx,yidx][:bin_size].compressed()\
         != 0)):
        type_str = 'Land' 
    else:    
        if(len(data_before.compressed()) == 0):
            type_str = 'Other'
        else:
            before_counts = np.full(3, np.nan)
            after_counts  = np.full(3, np.nan)

            before_count_ocean = len(data_before[data_before <= max_ice_for_ocean])
            before_count_mix   = len(data_before[(data_before > max_ice_for_ocean) & \
                                    (data_before < min_ice_for_ice)])
            before_count_ice   = len(data_before[data_before >= min_ice_for_ice])

            before_counts[0] = before_count_ocean
            before_counts[1] = before_count_mix
            before_counts[2] = before_count_ice

            after_count_ocean = len(data_after[data_after <= max_ice_for_ocean])
            after_count_mix   = len(data_after[(data_after > max_ice_for_ocean) & \
                                   (data_after < min_ice_for_ice)])
            after_count_ice   = len(data_after[data_after >= min_ice_for_ice])

            after_counts[0] = after_count_ocean
            after_counts[1] = after_count_mix
            after_counts[2] = after_count_ice

            delta_counts = after_counts - before_counts

            if(debug):        
                print(before_count_ocean, before_count_mix, before_count_ice) 
                print(after_count_ocean, after_count_mix, after_count_ice) 
                print(delta_counts[0], delta_counts[1], delta_counts[2]) 

            # NEED TO: add a check for the trend here. If the trend is very
            #           small, don't call for a change

            return_val = -9
            # First check: none of the counts per type changed
            if(delta_counts[0] == delta_counts[1] == delta_counts[2] == 0):
                if( np.argmax(before_counts) == 0):
                    type_str = "Unchanged ocean"
                elif( np.argmax(before_counts) == 1):
                    type_str = "Unchanged mix"
                else:
                    type_str = "Unchanged ice"
            else:
                # Check if this box has changed more to ocean
                if(np.argmax(delta_counts) == 0):

                    # First, see if it was also primarily ocean before
                    if(np.argmax(before_counts) == 0):
                        type_str = "Still primarily ocean"
                    # Now, see if the change is not significant
                    # -----------------------------------------
                    elif((np.argmax(before_counts) == 1) & \
                         (np.argmax(after_counts) == 1)):
                        type_str = 'Still primarily mix'
                    elif((np.argmax(before_counts) == 2) & \
                         (np.argmax(after_counts) == 2)):
                        type_str = "Still primarily ice"
                    else:
                        # Check if changing from mix or ice
            
                        # See if it was primarily mix before
                        if(before_counts[1] >= before_counts[2]):
                            type_str = "Mix to ocean"
                        else:               
                            type_str = "Ice to ocean"

                # Check if this box has changed more to mix
                elif(np.argmax(delta_counts) == 1):
                    # First, see if it was also primarily ocean before
                    if(np.argmax(before_counts) == 1):
                        type_str = "Still primarily mix"
                    # Now, see if the change is not significant
                    # -----------------------------------------
                    elif((np.argmax(before_counts) == 0) & \
                         (np.argmax(after_counts) == 0)):
                        type_str = "Still primarily ocean"
                    elif((np.argmax(before_counts) == 2) & \
                         (np.argmax(after_counts) == 2)):
                        type_str = "Still primarily ice"
                    else:
                        # See if it was primarily ocean before
                        if(before_counts[0] >= before_counts[2]):
                            type_str = "Ocean to mix"
                        else:               
                            type_str = "Ice to mix"

                # Check if this box has changed more to ice
                elif(np.argmax(delta_counts) == 2):
                    # First, see if it was also primarily ocean before
                    if(np.argmax(before_counts) == 2):
                        type_str = "Still primarily ice"
                    # Now, see if the change is not significant
                    # -----------------------------------------
                    elif((np.argmax(before_counts) == 0) & \
                         (np.argmax(after_counts) == 0)):
                        type_str = "Still primarily ocean"
                    elif((np.argmax(before_counts) == 1) & \
                         (np.argmax(after_counts) == 1)):
                        type_str = "Still primarily mix"
                    else:
                        # See if it was primarily ocean before
                        if(before_counts[0] > before_counts[1]):
                            type_str = "Ocean to ice"
                        else:               
                            type_str = "Mix to ice"
            if(plot_line): 
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.plot(NSIDC_data['grid_ice_conc'][month_idx::6,xidx,yidx])
                ax.axhline(max_ice_for_ocean, linestyle = ':', color = 'k')
                ax.axhline(min_ice_for_ice, linestyle = ':', color = 'k')
                ax.set_ylim(-5, 105)
                ax.set_title(type_str)
                fig.tight_layout()
                plt.show()

    if(debug): print(type_str)
    return_val = rval_dict[type_str]
    
    if(return_val == -9):
        print("SOMETHING WENT WRONG")
    return return_val
 
def plot_NSIDC_coverage_change(NSIDC_data, month_idx, xidx, yidx, \
        max_ice_for_ocean = 20, min_ice_for_ice = 80, bin_size = 5, \
        minlat = 70.5, save = False, cbar_switch = False):

    # Plotting types:
    # 0: unchanged ocean
    # 1: unchanged mix
    # 2: unchanged ice
    # 3: decrease in ice
    # 4: increase in ice

    #    'Unchanged ocean':       0,   # 0
    #    'Unchanged mix':         1,   # 1
    #    'Unchanged ice':         2,   # 2 
    #    'Still primarily ocean': 3,   # 0
    #    'Mix to ocean':          4,   # 3
    #    'Ice to ocean':          5,   # 3
    #    'Still primarily mix':   6,   # 1
    #    'Ocean to mix':          7,   # 4
    #    'Ice to mix':            8,   # 4
    #    'Still primarily ice':   9,   # 2
    #    'Ocean to ice':         10,   # 5
    #    'Mix to ice':           11,   #
    #    'Pole Hole':            12,
    #    'Unused':               13,
    #    'Coastline':            14,
    #    'Land':                 15,
    #    'Other':                -1,


    # ones:   all ocean all the time
    # twos:   all mix   all the time
    # threes: all ice   all the time
    # fours:  change from mix to ocean
    # fives:  change from ocean to mix
    # sixs:   change from ice to mix 
    base_ones   = np.full(NSIDC_data['grid_lat'].shape, 1)
    base_twos   = np.full(NSIDC_data['grid_lat'].shape, 2)
    base_threes = np.full(NSIDC_data['grid_lat'].shape, 3)
    base_fours  = np.full(NSIDC_data['grid_lat'].shape, 4)
    base_fives  = np.full(NSIDC_data['grid_lat'].shape, 5)
    base_sixs   = np.full(NSIDC_data['grid_lat'].shape, 6)

    # Go through and calculate the sfc type changes for each grid box
    # ---------------------------------------------------------------
    type_vals = np.array([[\
        check_type_change(NSIDC_data, month_idx, ii, jj, \
        min_ice_for_ice = min_ice_for_ice, \
        max_ice_for_ocean = max_ice_for_ocean, 
        bin_size = bin_size) \
        for jj in range(NSIDC_data['grid_lon'].shape[1])] \
        for ii in range(NSIDC_data['grid_lon'].shape[0])])

    ones    = np.ma.masked_where( ((type_vals != 0) & (type_vals != 3)), base_ones)
    twos    = np.ma.masked_where( ((type_vals != 1) & (type_vals != 6)), base_twos)
    threes  = np.ma.masked_where( ((type_vals != 2) & (type_vals != 9)), base_threes)
    fours   = np.ma.masked_where(  (type_vals != 4), base_fours)
    fives   = np.ma.masked_where(  (type_vals != 7), base_fives)
    sixs    = np.ma.masked_where(  (type_vals != 8), base_fives)
    #ax1.pcolormesh(day_data['longitude'][:], day_data['latitude'][:], \
    #    ones, transform = datacrs, shading = 'auto', \
    #    vmin = 1, vmax = 10, cmap = 'tab10')

    plt.close('all')
    fig = plt.figure(figsize = (9, 6))
    ax1 = fig.add_subplot(2,3,1, projection = mapcrs)
    ax2 = fig.add_subplot(2,3,2, projection = mapcrs)
    ax3 = fig.add_subplot(2,3,3, projection = mapcrs)
    ax4 = fig.add_subplot(2,3,4, projection = mapcrs)
    ax5 = fig.add_subplot(2,3,5, projection = mapcrs)
    ax6 = fig.add_subplot(2,3,6)
    
    plotNSIDC_MonthClimo(NSIDC_data, month_idx, minlat = minlat, ax = ax1, \
        title = '', colorbar = cbar_switch)
    plotNSIDC_MonthTrend(NSIDC_data,month_idx=month_idx,save=False,\
        trend_type='standard',season='sunlight',minlat=minlat,\
        return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = 10, \
        pax = ax2, show_pval = True, uncert_ax = ax3, title = '')

    avg1 = np.nanmean(NSIDC_data['grid_ice_conc'][month_idx::6,:,:][:5,:,:], \
        axis = 0)
    avg2 = np.nanmean(NSIDC_data['grid_ice_conc'][month_idx::6,:,:][-5:,:,:], \
        axis = 0)

    ax4.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
        avg1, transform = datacrs, shading = 'auto', cmap = 'ocean')
    ax4.coastlines()
    ax4.set_extent([-180,180,minlat,90], datacrs)

    #ax5.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    #    avg2, transform = datacrs, shading = 'auto', cmap = 'ocean')
    #ax5.coastlines()
    #ax5.set_extent([-180,180,minlat,90], datacrs)

    ax5.pcolormesh(NSIDC_data['grid_lon'][:,:], NSIDC_data['grid_lat'][:,:], \
        ones, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10')
    ax5.pcolormesh(NSIDC_data['grid_lon'][:,:], NSIDC_data['grid_lat'][:,:], \
        twos, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10')
    ax5.pcolormesh(NSIDC_data['grid_lon'][:,:], NSIDC_data['grid_lat'][:,:], \
        threes, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10')
    ax5.pcolormesh(NSIDC_data['grid_lon'][:,:], NSIDC_data['grid_lat'][:,:], \
        fours, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10')
    ax5.pcolormesh(NSIDC_data['grid_lon'][:,:], NSIDC_data['grid_lat'][:,:], \
        fives, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10')
    ax5.pcolormesh(NSIDC_data['grid_lon'][:,:], NSIDC_data['grid_lat'][:,:], \
        sixs, transform = datacrs, shading = 'auto', \
        vmin = 1, vmax = 10, cmap = 'tab10')
    ax5.coastlines()
    ax5.set_extent([-180,180,minlat,90], datacrs)

    plot_point_on_map(ax1, NSIDC_data['grid_lat'][xidx,yidx], \
        NSIDC_data['grid_lon'][xidx,yidx], markersize = 10, color = 'red')
    plot_point_on_map(ax2, NSIDC_data['grid_lat'][xidx,yidx], \
        NSIDC_data['grid_lon'][xidx,yidx], markersize = 10, color = 'red')
    plot_point_on_map(ax3, NSIDC_data['grid_lat'][xidx,yidx], \
        NSIDC_data['grid_lon'][xidx,yidx], markersize = 10, color = 'red')
    plot_point_on_map(ax4, NSIDC_data['grid_lat'][xidx,yidx], \
        NSIDC_data['grid_lon'][xidx,yidx], markersize = 10, color = 'red')
    plot_point_on_map(ax5, NSIDC_data['grid_lat'][xidx,yidx], \
        NSIDC_data['grid_lon'][xidx,yidx], markersize = 10, color = 'red')
 
    
    ax6.plot(NSIDC_data['grid_ice_conc'][month_idx::6,xidx,yidx])
    ax6.set_ylim(-5, 105)
    ax6.axhline(max_ice_for_ocean, linestyle = ':', color = 'k')
    ax6.axhline(min_ice_for_ice, linestyle = ':', color = 'k')
    ax6.set_title('Type idx = ' + str(type_vals[xidx, yidx]) )
 
    fig.tight_layout()
    plt.show()



