"""
  NAME:

  PURPOSE:
  
    NOTE: The MODIS channel and true color functions are designed to work with
    HDF MODIS files retriefed from 
    the data ordering website at this address:
    https://ladsweb.modaps.eosdis.nasa.gov/search/order/1/MODIS:Aqua


"""
import numpy as np
import numpy.ma as ma
import sys
from netCDF4 import Dataset
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
from scipy.stats import pearsonr,spearmanr
import subprocess
from scipy import stats
from pyhdf import SD
import pandas as pd
import h5py
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.tri import Triangulation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy.scene import Scene
from glob import glob

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Set up global variables
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
mapcrs = ccrs.NorthPolarStereo(central_longitude = 45.)
datacrs = ccrs.PlateCarree()

zoom_dict = {
    'Finland': [10,55,65,80]
}

proj_dict = {
    'Finland': ccrs.NorthPolarStereo(central_longitude = 35.) 
}

channel_dict = {
    '1': {
        'Name': 'EV_250_Aggr1km_RefSB',\
        'Index': 0,\
        'Bandwidth': [0.620, 0.670]
    },\
    '2': {
        'Name': 'EV_250_Aggr1km_RefSB',\
        'Index': 1,\
        'Bandwidth': [0.841, 0.876]
    },\
    '3': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 0,\
        'Bandwidth': [0.459, 0.479]
    },\
    '4': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 1,\
        'Bandwidth': [0.545, 0.565]
    },\
    '5': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 2,\
        'Bandwidth': [1.230, 1.250]
    },\
    '6': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 3,\
        'Bandwidth': [1.628, 1.652]
    },\
    '7': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 4,\
        'Bandwidth': [2.105, 2.155]
    },\
    '8': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 0,\
        'Bandwidth': [0.405, 0.420]
    },\
    '9': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 1,\
        'Bandwidth': [0.438, 0.448]
    },\
    '10': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 2,\
        'Bandwidth': [0.483, 0.493]
    },\
    '11': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 3,\
        'Bandwidth': [0.526, 0.536]
    },\
    '12': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 4,\
        'Bandwidth': [0.546, 0.556]
    },\
    '13lo': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 5,\
        'Bandwidth': [0.662, 0.672]
    },\
    '13hi': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 6,\
        'Bandwidth': [0.662, 0.672]
    },\
    '14lo': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 7,\
        'Bandwidth': [0.673, 0.683]
    },\
    '14hi': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 8,\
        'Bandwidth': [0.673, 0.683]
    },\
    '15': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 9,\
        'Bandwidth': [0.743, 0.753]
    },\
    '16': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 10,\
        'Bandwidth': [0.862, 0.877]
    },\
    '17': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 11,\
        'Bandwidth': [0.890, 0.920]
    },\
    '18': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 12,\
        'Bandwidth': [0.931, 0.941]
    },\
    '19': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 13,\
        'Bandwidth': [0.915, 0.965]
    },\
    '20': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 0,\
        'Bandwidth': [3.660, 3.840]
    },\
    '21': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 1,\
        'Bandwidth': [3.929, 3.989]
    },\
    '22': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 2,\
        'Bandwidth': [3.929, 3.989]
    },\
    '23': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 3,\
        'Bandwidth': [4.020, 4.080]
    },\
    '24': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 4,\
        'Bandwidth': [4.433, 4.498]
    },\
    '25': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 5,\
        'Bandwidth': [4.482, 4.549]
    },\
    '26': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 14,\
        'Bandwidth': [1.360, 1.390]
    },\
    '27': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 6,\
        'Bandwidth': [6.535, 6.895]
    },\
    '28': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 7,\
        'Bandwidth': [7.175, 7.475]
    },\
    '29': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 8,\
        'Bandwidth': [8.400, 8.700]
    },\
    '30': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 9,\
        'Bandwidth': [9.580, 9.880]
    },\
    '31': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 10,\
        'Bandwidth': [10.780, 11.280]
    },\
    '32': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 11,\
        'Bandwidth': [11.770, 12.270]
    },\
    '33': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 12,\
        'Bandwidth': [13.185, 13.485]
    },\
    '34': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 13,\
        'Bandwidth': [13.485, 13.785]
    },\
    '35': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 14,\
        'Bandwidth': [13.785, 14.085]
    },\
    '36': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 15,\
        'Bandwidth': [14.085, 14.385]
    }
}

plot_limits_dict = {
    "2021-07-20": {
        '2125': {
            'asos': 'asos_data_20210720.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021201.2125.061.2021202154814.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072010-2021072021.nc',
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5]
        }
    },
    "2021-07-22": {
        '2110': {
            'asos': 'asos_california_20210722.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5]
        }
    },
    "2021-08-05": {
        '2120': {
            'asos': 'asos_data_20210805.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021217.2120.061.2021218164201.hdf',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0]
        },
        '2125': {
            'asos': 'asos_california_20210805.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021217.2125.061.2021218161010.hdf',
            'Lat': [39.5, 42.5],
            'Lon': [-121.5, -119.5]
        }
    },
    "2021-08-06": {
        '2025': {
            'asos': 'asos_nevada_20210806.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0]
        }
    },
    "2021-08-17": {
        '2145': {
            'asos': 'asos_data_20210817.csv',
            'Lat': [38.0, 42.0],
            'Lon': [-122.0, -117.0]
        }
    },
    "2021-08-30": {
        '2115': {
            'asos': 'asos_data_20210830.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021242.2115.061.2021243183953.hdf',
            'Lat': [38.0, 40.0],
            'Lon': [-121.0, -118.5]
        }
    } 
}

def getCorners(centers):
    one = centers[:-1,:]
    two = centers[1:,:]
    d1 = (two - one) / 2.0
    one = one - d1
    two = two + d1
    stepOne = np.zeros((centers.shape[0] + 1, centers.shape[1]))
    stepOne[:-2,:] = one
    stepOne[-2:,:] = two[-2:,:]
    one = stepOne[:,:-1]
    two = stepOne[:,1:]
    d2 = (two - one) / 2.
    one = one - d2
    two = two + d2
    stepTwo = np.zeros((centers.shape[0] + 1, centers.shape[1] + 1))
    stepTwo[:,:-2] = one
    stepTwo[:,-2:] = two[:,-2:]
    return stepTwo

# Find the gridpoint in the gridded lat/lon data that 
# corresponds to the station at slat and slon
# ---------------------------------------------------- 
def nearest_gridpoint(slat, slon, grid_lat, grid_lon):
    fun_c = np.maximum(np.abs(grid_lat - slat), \
        np.abs(grid_lon - slon))
    m_idx = np.where(fun_c == np.min(fun_c))
    return m_idx
 
# Extract the MODIS information from a given channel at each ob point
# -------------------------------------------------------------------
def nearest_grid_values(MODIS_data):
    # Read in the correct ASOS file 
    asos_file = plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['asos']
    df = pd.read_csv(asos_file)
    df['valid'] = pd.to_datetime(df['valid'])
    df = df.set_index('valid')

    # Pull the event time from the plot_limits_dict
    event_date = datetime.strptime(MODIS_data['cross_date'], "%Y-%m-%d")
    first_time = MODIS_data['file_time']
    event_dtime = event_date + timedelta(hours = int(first_time[:2]))
    # Test what happens if looking at the data near the hour
    #event_dtime = event_date + timedelta(hours = int(first_time[:2]), \
    #    minutes = int(first_time[2:4]))
    print("Overpass time",event_dtime)

    begin_range = event_dtime - timedelta(minutes = 10)
    end_range   = event_dtime + timedelta(minutes = 10)

    compare_dict = {}
    station_names = sorted(set(df['station'].values))
    compare_dict['modis_time'] = event_dtime.strftime('%Y%m%d%H%M')
    compare_dict['stations'] = station_names
    compare_dict['stn_data'] = np.zeros((len(station_names)))
    compare_dict['mds_data'] = np.zeros((len(station_names)))
 
    for ii, station in enumerate(station_names):
        print(station)
        # Get the correct ob data
        stn_df = df[df['station'] == station]
        lat_stn = stn_df['lat'].values[0]
        lon_stn = stn_df['lon'].values[0]

        stn_df = stn_df[ begin_range : end_range]
        s_idx = np.argmin(np.abs((stn_df.index - event_dtime).total_seconds())) 
        if(stn_df['tmpc'][s_idx] == 'M'):
            stn_tmps = pd.to_numeric(stn_df['tmpc'], errors='coerce').values
            s_idx = np.where(~np.isnan(stn_tmps))[0][0]

        # Find the matching grid index
        m_idx = nearest_gridpoint(lat_stn, lon_stn, MODIS_data['lat'], \
            MODIS_data['lon'])
        m_data = MODIS_data['data'][m_idx][0]
      
        compare_dict['stn_data'][ii] = stn_df['tmpc'][s_idx]
        compare_dict['mds_data'][ii] = m_data
 
        print(station, lat_stn, MODIS_data['lat'][m_idx][0], \
            lon_stn, MODIS_data['lon'][m_idx][0], stn_df['tmpc'][s_idx], m_data )

    return compare_dict

 
# Plot the downloaded ASOS stations for each case
# ----------------------------------------------- 
def plot_ASOS_locs(pax,MODIS_data,color='red'):
    # Read in the correct ASOS file 
    asos_file = plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['asos']
    df = pd.read_csv(asos_file)

    station_names = set(df['station'].values)
    for ii, station in enumerate(station_names):
        lat_stn = df['lat'][df['station'] == station].values[0]
        lon_stn = df['lon'][df['station'] == station].values[0]
        pax.text(lon_stn, lat_stn, station, fontsize=10,weight='bold', \
            transform=datacrs, color=color)

# Determine areas of an image that are in smoke, defined by:
#    (ch1_refl - ch5_refl) < mean(ch1_refl - ch5_refl) * 0.25
def find_plume(filename):
    MODIS_ch1  = read_MODIS_channel(filename, 1,  zoom = True)
    MODIS_ch5  = read_MODIS_channel(filename, 5,  zoom = True)
    MODIS_ch31 = read_MODIS_channel(filename, 31, zoom = True)

    screen_limit = 0.40
    max_ch = 350.
    test_data = MODIS_ch1['data'] - MODIS_ch5['data']
    hash_data   = np.ma.masked_where(test_data <  \
        (np.nanmean(test_data) * screen_limit), MODIS_ch1['data'])
    nohash_data = np.ma.masked_where(test_data >= \
        (np.nanmean(test_data) * screen_limit), MODIS_ch1['data'])

    hash_data = np.ma.masked_where(MODIS_ch31['data'] > max_ch, hash_data)
    nohash_data = np.ma.masked_where(MODIS_ch31['data'] > max_ch, nohash_data)

    return hash_data, nohash_data

# Start_date and end_date must be formatted as 
# "YYYYMMDD"
# Variable must be either "CHL" or "PIC"
def read_MODIS(variable,start_date,latmin = 60):

    ##if(monthly == True):
    ##    dformat = '%Y%m' 
    ##    end_string_idx = -5
    ##else:
    ##    dformat = '%Y%m%d' 
    ##    end_string_idx = -3
    dformat = '%Y%m%d' 
    # Set up starting and ending datetime objects
    sdate = datetime.strptime(start_date,dformat)
    ##edate = datetime.strptime(end_date,dformat)
 
    
    # Grab all the modis files
    if(variable == 'CHL'):
        grabber = 'chlor-a'
        label = 'Chlorophyll Concentration, OCI Algorithm (mg m$^{-3}$)'
        ptitle = 'Chlorophyll-α'
        pticks = [0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0]
        ptick_labels = ['0.01','0.02','0.05','0.1','0.2','0.5','1','2','5',\
            '10','20']
    elif(variable == 'PIC'):
        grabber = 'pic'
        label = 'Calcite Concentration, Balch and Gordon (mol m$^{-3}$)'
        ptitle = 'Particulate Inorganic Carbon'
        pticks = [1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3,2e-3,5e-3,0.01,0.02,0.05]
        ptick_labels = ['1e-5','2e-5','5e-5','1e-4','2e-4','5e-4','1e-3',\
            '2e-3','5e-3','0.01','0.02','0.05']
    else:
        print("ERROR: variable must be CHL or PIC")
        return
    #cmnd = "ls /home/bsorenson/data/MODIS/CHL/A*.nc"
    cmnd = "ls /home/bsorenson/data/MODIS/"+variable+"/A*.nc"
    
    
    file_initial = subprocess.check_output(cmnd,shell=True).decode('utf-8').strip().split('\n')
    file_names = []

    for fname in file_initial:
        splitter = fname.split('/')[-1]
        fdate = datetime.strptime(splitter[1:5],'%Y') + relativedelta(days = int(splitter[5:8])-1)
        if(fdate == sdate):
            file_names.append(fname)
   
    lat_ranges = np.arange(latmin,90.25,0.25)
    lon_ranges = np.arange(-180.,180.25,0.25)
 
    # Read in the latitude, longitude, and area data
    num_files = len(file_names)
    modis_data = {}
    modis_data['data'] = np.full((num_files,len(lat_ranges),len(lon_ranges)),-9.)
    modis_data['lat']  = lat_ranges
    modis_data['lon']  = lon_ranges
    modis_data['variable'] = variable
    modis_data['ptitle'] = ptitle
    modis_data['label'] = label
    modis_data['grabber'] = grabber
    modis_data['pticks'] = pticks
    modis_data['ptick_labels'] = ptick_labels
    modis_data['titles']  = []
    
    count = 0
    for fname in file_names:
    
        total_name = fname
        #total_name = data_loc+fname
        print(total_name)
        data = Dataset(total_name,'r')
        
        for xi in range(len(lat_ranges)-1):       
            print(lat_ranges[xi])
            for yj in range(len(lon_ranges)-1):
                lat_indices = np.where((data.variables['lat'][:] >= lat_ranges[xi]) & (data.variables['lat'][:] < lat_ranges[xi+1]))
                lon_indices = np.where((data.variables['lon'][:] >= lon_ranges[yj]) & (data.variables['lon'][:] < lon_ranges[yj+1]))
                modis_data['data'][count,xi,yj] =\
                    np.average(data.variables[grabber][lat_indices[0],lon_indices[0]])
        data.close() 
        modis_data['titles'].append(fname)
        #modis_data['dates'].append(fname[-11:end_string_idx])
        count+=1

    # Convert the longitude values from 0 - 360 to -180 - 180
    return modis_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Trend calculating functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def trend_calc(modis_data,x_ind,y_ind,variable,thielsen=False):
    temp_data = modis_data[variable][:,x_ind,y_ind]
    temp_data[temp_data < -999.] = np.nan
    #temp_data = np.ma.masked_where(temp_data < -999., temp_data)
    # Don't calculate trends if there are less than 2 valid values
    if(np.count_nonzero(~np.isnan(temp_data)) < 2):
        trend = pcnt_change = np.nan
    else:
        avg_modis = np.nanmean(temp_data)
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


        pcnt_change = (trend/avg_modis)*100.
    # Find the percent change per decade
    return trend,pcnt_change


# modis_trendCalc calculates the trends over the time period at each
# grid point on the 25x25 km grid.
def modis_trendCalc(modis_data,thielSen=False):
    # Loop over the data and calculate trends
    for i in range(448):
        print(i)
        max_trend = -99.
        min_trend = 99.
        for j in range(304):
            modis_data['thick_trends'][i,j],temp_pcnt = \
                trend_calc(modis_data,i,j,'ice_thick',thielsen=thielSen)
            #modis_data['thick_land_trends'][i,j],temp_pcnt = \
            #    trend_calc(modis_data,i,j,'ice_thick',thielsen=thielSen)
            modis_data['con_trends'][i,j],temp_pcnt = \
                trend_calc(modis_data,i,j,'ice_con',thielsen=thielSen)
            #modis_data['con_land_trends'][i,j],temp_pcnt = \
            #    trend_calc(modis_data,i,j,'ice_con',thielsen=thielSen)
            ##temp_trend = modis_data['trends'][i,j]
            ##if(temp_trend>max_trend):
            ##    max_trend = temp_trend
            ##if(temp_trend<min_trend):
            ##    min_trend = temp_trend
        ##print("  max trend = ",max_trend)
        ##print("  min trend = ",min_trend)
        # Deal with land masks
        #good_indices = np.where(modis_data['data'][0,i,:]<251.)
        #land_indices = np.where(modis_data['data'][0,i,:]>=251.)
        #modis_data['trends'][i,land_indices] = np.nan
        #modis_data['land_trends'][i,good_indices] = np.nan

    return modis_data

# Calculate trends on the 1x1 degree lat/lon grid
def modis_gridtrendCalc(modis_data,area=True,thielSen=False):
    if(area==True):
        print("\nCalculating area trends\n")
    else:
        print("\nCalculating % concentration trends\n")
    modis_data['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = modis_data['grid_lat']
    lon_ranges = modis_data['grid_lon']
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
                interp_data = (modis_data['grid_modis_conc'][:,xi,yj]/100.)*modis_data['grid_total_area'][xi,yj]
                good_indices = np.where(np.isnan(interp_data)==False)
            else:
                interp_data = modis_data['grid_modis_conc'][:,xi,yj]
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
                modis_data['grid_modis_area_trend'][xi,yj] = total_trend
        print(xi)
        #print("max trend = ",max_trend)
        #print("min trend = ",min_trend)

    return modis_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Gridding functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# grid_data_conc grids the 25x25 km gridded modis concentration data into
# a 1x1 degree lat/lon grid
def grid_data_values(modis_data):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_con    = np.full((len(modis_data['ice_con'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_thick  = np.full((len(modis_data['ice_thick'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_con_cc = np.full((len(modis_data['ice_con'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    grid_thick_cc = np.full((len(modis_data['ice_thick'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    print("Size of grid array: ",grid_con.shape)
    for nt in range(len(modis_data['ice_con'][:,0,0])):
        print(nt)
        for xi in range(448):
            # Don't grid the data if any portion of the lat/lon box is over land.
            # Don't include land data
            for yj in range(304):
                lat_index = np.where(np.floor(modis_data['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(modis_data['lon'][xi,yj])>=lon_ranges)[-1][-1]
                # Add the current pixel area into the correct grid box, no
                # matter if the current box is missing or not.
                #if((lat_index==20) & (lon_index==10)):
                #    print("Current grid area = ",grid_modis_area[lat_index,lon_index])
                #if((lat_index==20) & (lon_index==10)):
                #    print("    New grid area = ",grid_modis_area[lat_index,lon_index])
                if(modis_data['ice_con'][nt,xi,yj] != -9999.0):
                    #if(nt==0): grid_modis_area[lat_index,lon_index] += modis_data['area'][xi,yj]
                    if(grid_con_cc[nt,lat_index,lon_index]==-999.):
                        grid_con[nt,lat_index,lon_index] = modis_data['ice_con'][nt,xi,yj]
                        grid_con_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_con[nt,lat_index,lon_index] = ((grid_con[nt,lat_index,\
                            lon_index]*grid_con_cc[nt,lat_index,lon_index])+\
                            modis_data['ice_con'][nt,xi,yj])/(grid_con_cc[nt,lat_index,\
                            lon_index]+1.)
                        grid_con_cc[nt,lat_index,lon_index]+=1
                if(modis_data['ice_thick'][nt,xi,yj] != -9999.0):
                    #if(nt==0): grid_modis_area[lat_index,lon_index] += modis_data['area'][xi,yj]
                    if(grid_thick_cc[nt,lat_index,lon_index]==-999.):
                        grid_thick[nt,lat_index,lon_index] = modis_data['ice_thick'][nt,xi,yj]
                        grid_thick_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_thick[nt,lat_index,lon_index] = ((grid_thick[nt,lat_index,\
                            lon_index]*grid_thick_cc[nt,lat_index,lon_index])+\
                            modis_data['ice_thick'][nt,xi,yj])/(grid_thick_cc[nt,lat_index,\
                            lon_index]+1.)
                        grid_thick_cc[nt,lat_index,lon_index]+=1
                #else:
                #    if(nt==0): grid_modis_area[lat_index,lon_index] = np.nan

                    # end else
                # end if good modis check
            # end y grid loop
        # end x grid loop
    # end time loop 
    modis_data['grid_con']   = grid_con
    modis_data['grid_thick'] = grid_thick
    modis_data['grid_con_cc'] = grid_con_cc
    modis_data['grid_thick_cc'] = grid_thick_cc
    modis_data['grid_lat'] = lat_ranges
    modis_data['grid_lon'] = lon_ranges
    
    return modis_data

# grid_data averages the 25x25 km trends into a 1x1 degree grid
# The 'grid_modis' paramter in the dictionary therefore contains
# the gridded trends.
def grid_data_trends(modis_dict):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_modis    = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    grid_modis_cc = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    for xi in range(448):
        print(xi)
        for yj in range(304):
            if(not np.isnan(modis_dict['trends'][xi,yj])):
                lat_index = np.where(np.floor(modis_dict['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(modis_dict['lon'][xi,yj])>=lon_ranges)[-1][-1]
                if(grid_modis_cc[lat_index,lon_index]==-999.):
                    grid_modis[lat_index,lon_index] = modis_dict['trends'][xi,yj]
                    grid_modis_cc[lat_index,lon_index] = 1.
                else:
                    grid_modis[lat_index,lon_index] = ((grid_modis[lat_index,lon_index]*grid_modis_cc[lat_index,lon_index])+\
                                                     modis_dict['trends'][xi,yj])/(grid_modis_cc[lat_index,lon_index]+1.)
                    grid_modis_cc[lat_index,lon_index]+=1
    modis_dict['grid_modis'] = grid_modis
    modis_dict['grid_modis_cc'] = grid_modis_cc
    modis_dict['grid_lat'] = lat_ranges
    modis_dict['grid_lon'] = lon_ranges
    return modis_dict

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def plot_true_color_satpy(date_str,zoom=True,save=False):

    # Determine the correct MODIS file associated with the date
    # ---------------------------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']
    print(filename)
    day_filenames = glob(filename[:50]+'*')

    # Use satpy (Scene) to open the file
    # ----------------------------------
    scn = Scene(reader = 'modis_l1b', filenames = day_filenames)

    scn.load(['1'])

    my_area = scn['1'].attrs['area'].compute_optimal_bb_area({\
        'proj': 'lcc', 'lon_0': -122., 'lat_0': 37., 'lat_1': 40., 'lat_2': 40.})
    new_scn = scn.resample(my_area)

    crs = new_scn['1'].attrs['area'].to_cartopy_crs()
    ax = plt.axes(projection=crs)
    print(new_scn)
    ax.coastlines()
    ax.gridlines()
    ax.set_global()
    plt.imshow(new_scn['1'], cmap='Greys_r',transform=crs, extent=crs.bounds, origin='upper')
    #cbar = plt.colorbar()
    #cbar.set_label("Kelvin")
    if(zoom):
        lat_lims = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']
        lon_lims = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']
        ax.set_extent([lat_lims[0],lat_lims[1],lon_lims[0],lon_lims[1]],\
                       crs)
    plt.show()

    """

    # Extract the true-color data
    # ---------------------------
    scn.load(['true_color'])

    scn.save_dataset('true_color','test_image_true'+date_str+'.png')
    scn.show('true_color')

    # Overlay the imagery on cartopy
    # ------------------------------
    ##my_area

    # Save the image
    # --------------
    """

def plot_true_color(filename,zoom=True):
    modis = SD.SD(filename)

    dat = modis.attributes().get('CoreMetadata.0').split()
    indx = dat.index('EQUATORCROSSINGDATE')+9
    cross_date = dat[indx][1:len(dat[indx])-1]
    cross_time = filename.strip().split('/')[-1].split('.')[2]

    lat5 = modis.select('Latitude').get()
    lon5 = modis.select('Longitude').get()

    red   = modis.select('EV_250_Aggr1km_RefSB').get()[0]
    green = modis.select('EV_500_Aggr1km_RefSB').get()[1]
    blue  = modis.select('EV_500_Aggr1km_RefSB').get()[0]

    red   = red[::5,::5]
    green = green[::5,::5]
    blue  = blue[::5,::5]

    red_scale    = modis.select('EV_250_Aggr1km_RefSB').attributes().get('reflectance_scales')[0]
    red_offset   = modis.select('EV_250_Aggr1km_RefSB').attributes().get('reflectance_offsets')[0]
    green_scale  = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_scales')[1]
    green_offset = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_offsets')[1]
    blue_scale   = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_scales')[0]
    blue_offset  = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_offsets')[0]

    modis.end()

    red   = (red - red_offset) * red_scale
    green = (green - green_offset) * green_scale
    blue  = (blue - blue_offset) * blue_scale

    red   = red*(255./1.1) 
    green = green*(255./1.1) 
    blue  = blue*(255./1.1) 

    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,110,160,210,240,255])
    red     = np.interp(red,old_val,ehn_val) / 255.
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,110,160,200,230,240])
    green   = np.interp(green,old_val,ehn_val) / 255.
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,100,150,210,240,255])
    blue    = np.interp(blue,old_val,ehn_val) / 255.
  
    image = np.zeros((red.shape[0],red.shape[1],3))
    image[:,:,0] = red
    image[:,:,1] = green
    image[:,:,2] = blue
   
    colortuple = tuple(np.array([image[:,:,0].flatten(), \
        image[:,:,1].flatten(), image[:,:,2].flatten()]).transpose().tolist())

    cornerLats = getCorners(lat5) ; cornerLons = getCorners(lon5)
    
    plt.close('all') 
    fig1 = plt.figure()
    ax = plt.axes(projection = ccrs.LambertConformal())

    image = np.ma.masked_where(np.isnan(image),image)

    print(lon5.shape, lat5.shape, image.shape)

    ax.pcolormesh(cornerLons,cornerLats,image[:,:,0],color= colortuple, shading='auto', \
        transform = ccrs.PlateCarree()) 
    
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES)
    ax.coastlines()
    if(zoom):
        ax.set_extent([plot_limits_dict[cross_date][cross_time]['Lon'][0], \
                       plot_limits_dict[cross_date][cross_time]['Lon'][1], \
                       plot_limits_dict[cross_date][cross_time]['Lat'][0], \
                       plot_limits_dict[cross_date][cross_time]['Lat'][1]], \
            ccrs.PlateCarree())

    plt.show()

def read_MODIS_channel(filename, channel, zoom = False):

    print("Reading MODIS channel",channel," from ",filename)

    MODIS_data = {}

    modis = SD.SD(filename)

    dat = modis.attributes().get('CoreMetadata.0').split()
    indx = dat.index('EQUATORCROSSINGDATE')+9
    cross_date = dat[indx][1:len(dat[indx])-1]

    print(cross_date)

    lat5 = modis.select('Latitude').get()
    lon5 = modis.select('Longitude').get()

    data  = modis.select(channel_dict[str(channel)]['Name']).get()[channel_dict[str(channel)]['Index']]

    data  = data[::5,::5]

    # Thermal emission data
    if((channel >= 20) & (channel != 26)):
        data_scale    = modis.select(channel_dict[str(channel)]['Name']\
            ).attributes().get('radiance_scales')[channel_dict[str(channel)]['Index']]
        data_offset   = modis.select(channel_dict[str(channel)]['Name']\
            ).attributes().get('radiance_offsets')[channel_dict[str(channel)]['Index']]

        data = (data - data_offset) * data_scale

        # Define constants for converting radiances to temperatures
        lmbda = (1e-6) * (np.average(channel_dict[str(channel)]['Bandwidth'])) # in m
        print("Average wavelength = ",np.average(channel_dict[str(channel)]['Bandwidth']))
        c_const = 3e8
        h_const = 6.626e-34 # J*s
        k_const = 1.381e-23 # J/K

        data = (h_const * c_const) / \
            (lmbda * k_const * np.log( ((2.0 * h_const * (c_const**2.0) ) / \
            ((lmbda**4.) * (lmbda / 1e-6) * data ) ) + 1.0 ) )
        #data = ((h_const * c_const)/(k_const * lmbda)) * (np.log((2.0 * h_const * (c_const ** 2.0) / \
        #    ((lmbda**5.0) * data)) + 1) ** -1.)

        colors = 'plasma'
        label = 'Blackbody Temperature [K]'

    # Reflectances
    else:
        data_scale    = modis.select(channel_dict[str(channel)]['Name']\
            ).attributes().get('reflectance_scales')[channel_dict[str(channel)]['Index']]
        data_offset   = modis.select(channel_dict[str(channel)]['Name']\
            ).attributes().get('reflectance_offsets')[channel_dict[str(channel)]['Index']]

        # Calculate reflectance using the scales and offsets
        # -------------------------------------------------
        data = ((data - data_offset) * data_scale)

        colors = 'Greys_r'
        label = 'Reflectance'
 
    modis.end()

    MODIS_data['data'] = data
    MODIS_data['lat']  = lat5
    MODIS_data['lon']  = lon5
    MODIS_data['variable']  = label
    MODIS_data['cross_date']  = cross_date
    MODIS_data['channel']  = channel
    MODIS_data['colors']  = colors
    MODIS_data['file_time'] = filename.strip().split('/')[-1].split('.')[2]

    if(zoom):
        # Mask MODIS_data['data'] that are outside the desired range
        # --------------------------------------------
        MODIS_data['data'][(((MODIS_data['lat'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
                             (MODIS_data['lat'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
                            ((MODIS_data['lon'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
                             (MODIS_data['lon'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.
        MODIS_data['lat'][ (((MODIS_data['lat'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
                             (MODIS_data['lat'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
                            ((MODIS_data['lon'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
                             (MODIS_data['lon'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.
        MODIS_data['lon'][ (((MODIS_data['lat'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
                             (MODIS_data['lat'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
                            ((MODIS_data['lon'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
                             (MODIS_data['lon'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.

        MODIS_data['data'] = np.ma.masked_where(MODIS_data['data'] == -999., MODIS_data['data'])
        MODIS_data['lat'] = np.ma.masked_where(MODIS_data['lat'] == -999., MODIS_data['lat'])
        MODIS_data['lon'] = np.ma.masked_where(MODIS_data['lon'] == -999., MODIS_data['lon'])


    return MODIS_data

def plot_MODIS_channel(date_str,channel,zoom=True,show_smoke=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel == 'red'):
        channel = 1
    elif(channel == 'green'):
        channel = 4
    elif(channel == 'blue'):
        channel = 3

    # Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # ---------------------------------------------------------------
    MODIS_data = read_MODIS_channel(filename, channel)

    print("Data max = ",np.max(MODIS_data['data']), "  Data min = ",np.min(MODIS_data['data']))

    plt.close('all')
    fig1 = plt.figure()
    datacrs = ccrs.PlateCarree()
    mapcrs = ccrs.LambertConformal()
    ax = plt.axes(projection = mapcrs)

    mesh = ax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
        MODIS_data['data'],cmap = MODIS_data['colors'], shading='auto', \
        transform = datacrs) 

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(filename) 
        hash0 = ax0.pcolor(MODIS_data['lon'],MODIS_data['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs)

    cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.850)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(MODIS_data['variable'],fontsize=16,weight='bold')
    
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES)
    ax.coastlines()
    if(zoom):
        ax.set_extent([plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0], \
                       plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1], \
                       plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0], \
                       plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1]],\
                       ccrs.PlateCarree())
    ax.set_title('Channel ' + str(channel) + '\n' + \
        str(channel_dict[str(channel)]['Bandwidth'][0]) + ' μm - ' + \
        str(channel_dict[str(channel)]['Bandwidth'][1]) + ' μm')

    plt.show()

# date_str of format "YYYYMMDDHHMM"
def compare_MODIS_3panel(date_str,channel1,channel2,channel3,zoom=True,save=False,\
        plot_ASOS_loc = False, show_smoke = True, compare_OMI = False, \
        compare_CERES = False, return_MODIS = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel1== 'red'):
        channel1= 1
    elif(channel1== 'green'):
        channel1= 4
    elif(channel1== 'blue'):
        channel1= 3

    # Step 1: Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    MODIS_data1 = read_MODIS_channel(filename, channel1, zoom = True)
    MODIS_data2 = read_MODIS_channel(filename, channel2, zoom = True)
    MODIS_data3 = read_MODIS_channel(filename, channel3, zoom = True)

    ##!#print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",np.min(MODIS_data1['data']))
    ##!#print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",np.min(MODIS_data2['data']))
    ##!#print("Data3 max = ",np.max(MODIS_data3['data']), "  Data3 min = ",np.min(MODIS_data3['data']))

    # Create copies to set the minimum and maximums for each plot
    # -----------------------------------------------------------
    cpy_1 = np.copy(MODIS_data1['data'])
    cpy_2 = np.copy(MODIS_data2['data'])
    cpy_3 = np.copy(MODIS_data3['data'])

    cpy_1 = np.ma.masked_where((((MODIS_data1['lat'] < plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][0]) | \
                         (MODIS_data1['lat'] > plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][1])) | \
                        ((MODIS_data1['lon'] < plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][0]) | \
                         (MODIS_data1['lon'] > plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][1]))), cpy_1)
    cpy_2 = np.ma.masked_where((((MODIS_data2['lat'] < plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][0]) | \
                         (MODIS_data2['lat'] > plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][1])) | \
                        ((MODIS_data2['lon'] < plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][0]) | \
                         (MODIS_data2['lon'] > plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][1]))), cpy_2)
    cpy_3 = np.ma.masked_where((((MODIS_data3['lat'] < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                         (MODIS_data3['lat'] > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                        ((MODIS_data3['lon'] < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                         (MODIS_data3['lon'] > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), cpy_3)


    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    datacrs = ccrs.PlateCarree() 
    mapcrs = ccrs.LambertConformal()

    plt.close('all')
    if(compare_OMI and compare_CERES):
        fig = plt.figure(figsize=(12,9))
        ax0 = fig.add_subplot(2,3,1,projection = mapcrs)
        ax1 = fig.add_subplot(2,3,2,projection = mapcrs)
        ax2 = fig.add_subplot(2,3,3,projection = mapcrs)
        axo = fig.add_subplot(2,3,4,projection = mapcrs)
        axcs = fig.add_subplot(2,3,5,projection = mapcrs)
        axcl = fig.add_subplot(2,3,6,projection = mapcrs)
    elif((compare_CERES and not compare_OMI)):
        fig = plt.figure(figsize=(12,9))
        ax0 = fig.add_subplot(2,3,1,projection = mapcrs)
        ax1 = fig.add_subplot(2,3,2,projection = mapcrs)
        ax2 = fig.add_subplot(2,3,3,projection = mapcrs)
        axcs = fig.add_subplot(2,3,4,projection = mapcrs)
        axcl = fig.add_subplot(2,3,5,projection = mapcrs)
    elif((compare_OMI and not compare_CERES)):
        fig = plt.figure(figsize=(9,9))
        ax0 = fig.add_subplot(2,2,1,projection = mapcrs)
        ax1 = fig.add_subplot(2,2,2,projection = mapcrs)
        ax2 = fig.add_subplot(2,2,3,projection = mapcrs)
        axo = fig.add_subplot(2,2,4,projection = mapcrs)
    else:
        fig = plt.figure(figsize=(14.5,5))
        ax0 = fig.add_subplot(1,3,1,projection = mapcrs)
        ax1 = fig.add_subplot(1,3,2,projection = mapcrs)
        ax2 = fig.add_subplot(1,3,3,projection = mapcrs)

    # Step 3: Plot the MODIS channel data in the first 2 panels
    # ---------------------------------------------------------
    mesh0 = ax0.pcolormesh(MODIS_data1['lon'],MODIS_data1['lat'],\
        MODIS_data1['data'],cmap = MODIS_data1['colors'], shading='auto', \
        vmin = np.nanmin(cpy_1), vmax = 330, transform = datacrs) 

    cbar0 = plt.colorbar(mesh0,ax=ax0,orientation='vertical',\
        pad=0.03,label=MODIS_data1['variable'])

    ax0.add_feature(cfeature.BORDERS)
    ax0.add_feature(cfeature.STATES)
    ax0.coastlines()
    if(zoom):
        ax0.set_extent([plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][0], \
                        plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][1], \
                        plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][0], \
                        plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][1]],\
                        datacrs)
    ax0.set_title('MODIS Ch. ' + str(channel1) + '\n' + \
        str(channel_dict[str(channel1)]['Bandwidth'][0]) + ' μm - ' + \
        str(channel_dict[str(channel1)]['Bandwidth'][1]) + ' μm')

    # Plot channel 2
    mesh1 = ax1.pcolormesh(MODIS_data2['lon'],MODIS_data2['lat'],\
        MODIS_data2['data'],cmap = MODIS_data2['colors'], shading='auto', \
        vmin = np.nanmin(cpy_2), vmax = np.nanmax(cpy_2), transform = datacrs) 


    cbar1 = plt.colorbar(mesh1,ax=ax1,orientation='vertical',\
        pad=0.03,label=MODIS_data2['variable'])
    
    ax1.add_feature(cfeature.BORDERS)
    ax1.add_feature(cfeature.STATES)
    ax1.coastlines()
    if(zoom):
        ax1.set_extent([plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][0], \
                        plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][1], \
                        plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][0], \
                        plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][1]],\
                        datacrs)
    ax1.set_title('MODIS Ch. ' + str(channel2) + '\n' + \
        str(channel_dict[str(channel2)]['Bandwidth'][0]) + ' μm - ' + \
        str(channel_dict[str(channel2)]['Bandwidth'][1]) + ' μm')

    # Plot channel 3
    mesh2 = ax2.pcolormesh(MODIS_data3['lon'],MODIS_data3['lat'],\
        MODIS_data3['data'],cmap = MODIS_data3['colors'], shading='auto', \
        vmin = np.nanmin(cpy_3), vmax = np.nanmax(cpy_3), transform = datacrs) 
    cbar2 = plt.colorbar(mesh2,ax=ax2,orientation='vertical',\
        pad=0.03,label=MODIS_data3['variable'])
    
    ax2.add_feature(cfeature.BORDERS)
    ax2.add_feature(cfeature.STATES)
    ax2.coastlines()
    if(zoom):
        ax2.set_extent([plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0], \
                        plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1], \
                        plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0], \
                        plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]],\
                        datacrs)
    ax2.set_title('MODIS Ch. ' + str(channel3) + '\n' + \
        str(channel_dict[str(channel3)]['Bandwidth'][0]) + ' μm - ' + \
        str(channel_dict[str(channel3)]['Bandwidth'][1]) + ' μm')


    if(compare_OMI):
        print("Reading OMI data")
        data = h5py.File(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['omi'],'r')
        LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
        LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
        UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
        mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), UVAI)
        mask_UVAI = np.ma.masked_where((((LAT < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), mask_UVAI)

        mesh3 = axo.pcolormesh(LON,LAT,mask_UVAI, cmap = 'plasma', shading='auto', \
            vmin = np.nanmin(mask_UVAI), vmax = np.nanmax(mask_UVAI), transform = datacrs) 
        axo.add_feature(cfeature.BORDERS)
        axo.add_feature(cfeature.STATES)
        axo.coastlines()
        if(zoom):
            axo.set_extent([plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0], \
                            plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1], \
                            plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0], \
                            plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]],\
                            datacrs)
        cbar3 = plt.colorbar(mesh3,ax=axo,orientation='vertical',\
            pad=0.03,label='OMI UVAI')
        axo.set_title('OMI UVAI')
        data.close()

    if(compare_CERES):
        print("Reading CERES data")

        # NOTE: need to use the time variable to screen out any data
        # that are not close to the event time.
        base_date = datetime(year=1970,month=1,day=1)
        start_date = dt_date_str - timedelta(hours = 1)
        end_date   = dt_date_str + timedelta(hours = 2)

        print(start_date, end_date)

        data = Dataset(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['ceres'],'r')
        LAT   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
        LON   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
        LON[LON>179.99] = -360.+LON[LON>179.99]
        swflux  = data.variables['CERES_SW_TOA_flux___upwards'][:]
        lwflux  = data.variables['CERES_LW_TOA_flux___upwards'][:]
        time  = data.variables['time'][:]
        local_time = np.array([base_date + relativedelta(days = ttime) for ttime in time])

        mask_LAT = LAT[ \
            (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_LON = LON[ \
            (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_swf = swflux[ \
            (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_lwf = lwflux[ \
            (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_time = local_time[ \
            (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]

        return mask_LAT, mask_LON, mask_swf, mask_lwf, mask_time
   
        # Removed masked data
        triMesh = Triangulation(mask_LON,mask_LAT)

        ##!## Reshape the data to make it work with pcolormesh
        ##!#first_size = int(np.sqrt(mask_swf.shape)) 
        ##!#max_size = int(first_size ** 2.)
        ##!#mask_swf = mask_swf[:max_size].reshape((first_size, first_size))
        ##!#mask_lwf = mask_lwf[:max_size].reshape((first_size, first_size))
        ##!#LAT = LAT[:max_size].reshape((first_size, first_size))
        ##!#LON = LON[:max_size].reshape((first_size, first_size))
      
        print(np.nanmax(mask_LAT.compressed()), np.min(mask_LAT.compressed()))
        print(np.nanmax(LAT), np.nanmin(LAT))
        print(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'])
 
        #scat3 = axcs.scatter(mask_LON, mask_LAT,mask_swf, transform = datacrs)
        mesh3 = axcs.scatter(mask_LON.compressed(), mask_LAT.compressed(),s = 120,marker='s',c = mask_swf.compressed(),cmap='plasma', transform = datacrs)
        #mesh3 = axcs.tricontourf(triMesh,mask_swf, cmap = 'plasma', shading='auto', \
        #    vmin = np.nanmin(mask_swf), vmax = np.nanmax(mask_swf), transform = datacrs) 
#        mesh3 = axcs.pcolormesh(LON,LAT,mask_swf, cmap = 'plasma', shading='auto', \
#            vmin = np.nanmin(mask_swf), vmax = np.nanmax(mask_swf), transform = datacrs) 
        axcs.add_feature(cfeature.BORDERS)
        axcs.add_feature(cfeature.STATES)
        axcs.coastlines()
        if(zoom):
            axcs.set_extent([plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0], \
                             plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1], \
                             plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0], \
                             plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]],\
                             datacrs)
        cbar3 = plt.colorbar(mesh3,ax=axcs,orientation='vertical',\
            pad=0.03,label='TOA SWF [W/m2]')
        axcs.set_title('CERES SWF')

        mesh4 = axcl.scatter(mask_LON.compressed(), mask_LAT.compressed(),s = 120,marker = 's',c = mask_lwf.compressed(),cmap='plasma', transform = datacrs)
        #mesh4 = axcl.tricontourf(triMesh,mask_lwf, cmap = 'plasma', shading='auto', \
        #    vmin = np.nanmin(mask_lwf), vmax = np.nanmax(mask_lwf), transform = datacrs) 
        ##!#mesh4 = axcl.tricontourf(LON,LAT,mask_lwf, cmap = 'plasma', shading='auto', \
        ##!#    vmin = np.nanmin(mask_lwf), vmax = np.nanmax(mask_lwf), transform = datacrs) 
        axcl.add_feature(cfeature.BORDERS)
        axcl.add_feature(cfeature.STATES)
        axcl.coastlines()
        if(zoom):
            axcl.set_extent([plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0], \
                            plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1], \
                            plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0], \
                            plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]],\
                            datacrs)
            cbar4 = plt.colorbar(mesh4,ax=axcl,orientation='vertical',\
                pad=0.03,label='TOA LWF [W/m2]')
        axcl.set_title('CERES LWF')

        data.close()

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(filename) 
        hash0 = ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs)
        hash1 = ax1.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs)
        hash2 = ax2.pcolor(MODIS_data3['lon'],MODIS_data3['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs)
        if(compare_OMI): 
            hash3 = axo.pcolor(MODIS_data3['lon'],MODIS_data3['lat'],\
                hash_data1, hatch = '///', alpha=0., transform = datacrs)
        if(compare_CERES): 
            hash4 = axcs.pcolor(MODIS_data3['lon'],MODIS_data3['lat'],\
                hash_data1, hatch = '///', alpha=0., transform = datacrs)
            hash5 = axcl.pcolor(MODIS_data3['lon'],MODIS_data3['lat'],\
                hash_data1, hatch = '///', alpha=0., transform = datacrs)
    

    if(plot_ASOS_loc):
        print("Plotting ASOS site")
        plot_ASOS_locs(ax0,MODIS_data1,color='black')
        plot_ASOS_locs(ax1,MODIS_data1)
        plot_ASOS_locs(ax2,MODIS_data1)
        if(compare_OMI): plot_ASOS_locs(axo,MODIS_data1,color='black')
        if(compare_CERES): 
            plot_ASOS_locs(axcs,MODIS_data1,color='black')
            plot_ASOS_locs(axcl,MODIS_data1,color='black')

    cross_date = MODIS_data1['cross_date']
    file_time  = MODIS_data1['file_time']
    plt.suptitle(cross_date + ' ' + file_time)
    if(save):
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel1) + '_vs_ch' + \
            str(channel2) + '_ch' + str(channel3) + '_' + pdate + '_3panel.png'
        if(compare_OMI):
            outname = 'modis_compare_ch' + str(channel1) + '_vs_ch' + \
                str(channel2) + '_ch' + str(channel3) + '_omi_' + pdate + '_3panel.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
        plt.show()

    if(return_MODIS):
        return MODIS_data1
    

def compare_MODIS_channels(date_str,channel1,channel2,zoom=True,save=False,\
        plot_ASOS_loc = False,show_smoke = True):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel1 == 'red'):
        channel1 = 1
    elif(channel1 == 'green'):
        channel1 = 4
    elif(channel1 == 'blue'):
        channel1 = 3

    # -----------------------------------------------------------------------
    # Step 1: Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    MODIS_data1 = read_MODIS_channel(filename, channel1, zoom = zoom)
    MODIS_data2 = read_MODIS_channel(filename, channel2, zoom = zoom)

    print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",np.min(MODIS_data1['data']))
    print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",np.min(MODIS_data2['data']))

    print(MODIS_data1['data'].shape, MODIS_data2['data'].shape)

    # --------------------------------------
    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    datacrs = ccrs.PlateCarree() 
    mapcrs = ccrs.LambertConformal()

    plt.close('all')
    fig = plt.figure(figsize=(14.5,5))
    ax0 = fig.add_subplot(1,3,2,projection = mapcrs)
    ax1 = fig.add_subplot(1,3,3,projection = mapcrs)
    ax2 = fig.add_subplot(1,3,1)

    # Step 3: Plot the MODIS channel data in the first 2 panels
    # ---------------------------------------------------------
    mesh0 = ax0.pcolormesh(MODIS_data1['lon'],MODIS_data1['lat'],\
        MODIS_data1['data'],cmap = MODIS_data1['colors'], shading='auto', \
        transform = datacrs)

    cbar0 = plt.colorbar(mesh0,ax=ax0,orientation='vertical',\
        pad=0.03,label=MODIS_data1['variable'])
    #cbar0 = plt.colorbar(mesh0,cax = ax0, orientation='horizontal',pad=0,\
    #    aspect=50,shrink = 0.850)
    #cbar0.ax.tick_params(labelsize=14)
    #cbar0.set_label(MODIS_data1['variable'],fontsize=16,weight='bold')
   
    ax0.add_feature(cfeature.BORDERS)
    ax0.add_feature(cfeature.STATES)
    ax0.coastlines()
    if(zoom):
        ax0.set_extent([plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][0], \
                        plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][1], \
                        plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][0], \
                        plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][1]],\
                        datacrs)
    ax0.set_title('MODIS Ch. ' + str(channel1) + '\n' + \
        str(channel_dict[str(channel1)]['Bandwidth'][0]) + ' μm - ' + \
        str(channel_dict[str(channel1)]['Bandwidth'][1]) + ' μm')

    # Plot channel 2
    mesh1 = ax1.pcolormesh(MODIS_data2['lon'],MODIS_data2['lat'],\
        MODIS_data2['data'],cmap = MODIS_data2['colors'], shading='auto', \
        transform = datacrs) 

    cbar1 = plt.colorbar(mesh1,ax=ax1,orientation='vertical',\
        pad=0.03,label=MODIS_data2['variable'])
    #cbar1 = plt.colorbar(mesh1,cax = ax1,orientation='horizontal',pad=1,\
    #    aspect=51,shrink = 1.851)
    #cbar1.ax.tick_params(labelsize=14)
    #cbar1.set_label(MODIS_data2['variable'],fontsize=16,weight='bold')
    
    ax1.add_feature(cfeature.BORDERS)
    ax1.add_feature(cfeature.STATES)
    ax1.coastlines()
    if(zoom):
        ax1.set_extent([plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][0], \
                        plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][1], \
                        plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][0], \
                        plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][1]],\
                        datacrs)
    ax1.set_title('MODIS Ch. ' + str(channel2) + '\n' + \
        str(channel_dict[str(channel2)]['Bandwidth'][0]) + ' μm - ' + \
        str(channel_dict[str(channel2)]['Bandwidth'][1]) + ' μm')

    ##!## Shade areas that are inside the plume, defined currently as areas
    ##!## with reflectance below the mean
    ##!## -----------------------------------------------------------------
    ##!#if(MODIS_data1['variable'] == 'Reflectance'):
    ##!#    hash_data1 = np.ma.masked_where(MODIS_data1['data'] < \
    ##!#        np.nanmean(MODIS_data1['data']),MODIS_data1['data'])
    ##!#    hash_data2 = np.ma.masked_where(MODIS_data1['data'] < \
    ##!#        np.nanmean(MODIS_data1['data']),MODIS_data2['data'])

    ##!#    nohash_data1 = np.ma.masked_where(MODIS_data1['data'] >= \
    ##!#        np.nanmean(MODIS_data1['data']),MODIS_data1['data'])
    ##!#    nohash_data2 = np.ma.masked_where(MODIS_data1['data'] >= \
    ##!#        np.nanmean(MODIS_data1['data']),MODIS_data2['data'])

    ##!#    hash1 = ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
    ##!#        hash_data1, hatch = '///', alpha=0., transform = datacrs),
    ##!#    hash2 = ax1.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
    ##!#        hash_data2, hatch = '///', alpha=0., transform = datacrs),
    ##!#if(MODIS_data2['variable'] == 'Reflectance'):
    ##!#    hash_data2 = np.ma.masked_where(MODIS_data2['data'] < \
    ##!#        np.nanmean(MODIS_data2['data']),MODIS_data2['data'])
    ##!#    hash_data1 = np.ma.masked_where(MODIS_data2['data'] < \
    ##!#        np.nanmean(MODIS_data2['data']),MODIS_data1['data'])

    ##!#    nohash_data2 = np.ma.masked_where(MODIS_data2['data'] >= \
    ##!#        np.nanmean(MODIS_data2['data']),MODIS_data2['data'])
    ##!#    nohash_data1 = np.ma.masked_where(MODIS_data2['data'] >= \
    ##!#        np.nanmean(MODIS_data2['data']),MODIS_data1['data'])
    ##!#

    ##!#    hash2 = ax1.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
    ##!#        hash_data2, hatch = '///', alpha=0., transform = datacrs),
    ##!#    hash1 = ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
    ##!#        hash_data1, hatch = '///', alpha=0., transform = datacrs),
    
    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(filename) 
        hash2 = ax1.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs),
        hash1 = ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs),
    
    # Step 4: Plot scatter MODIS channel comparison in third panel
    # ------------------------------------------------------------

    # Figure out which pixels are in the plume
    # ----------------------------------------
    max_ch = 350.

    tmp_data1 = np.copy(MODIS_data1['data'])
    tmp_data2 = np.copy(MODIS_data2['data'])

    tmp_data1 = np.ma.masked_where( (MODIS_data1['data'] > max_ch) | \
        (MODIS_data2['data'] > max_ch), tmp_data1)
    tmp_data2 = np.ma.masked_where( (MODIS_data1['data'] > max_ch) | \
        (MODIS_data2['data'] > max_ch), tmp_data2)

    ###hash_data1 = np.ma.masked_where((hash_data1 > max_ch) | \
    ###    (hash_data2 > max_ch), hash_data1)
    ###hash_data2 = np.ma.masked_where((hash_data1 > max_ch) | \
    ###    (hash_data2 > max_ch), hash_data2)
    ###nohash_data1 = np.ma.masked_where((nohash_data1 > max_ch) | \
    ###    (nohash_data2 > max_ch), nohash_data1)
    ###nohash_data2 = np.ma.masked_where((nohash_data1 > max_ch) | \
    ###    (nohash_data2 > max_ch), nohash_data2)

    plot_data1 = tmp_data1.compressed()
    plot_data2 = tmp_data2.compressed()

    ###print("Hash 1 size:", hash_data1.compressed().shape, "  Hash 2 size: ", hash_data2.compressed().shape)
    ###print("No Hash 1 size:", nohash_data1.compressed().shape, "  No Hash 2 size: ", nohash_data2.compressed().shape)

    rval_p = pearsonr(plot_data1,plot_data2)[0]
    #rval_s = spearmanr(tmp_data1,tmp_data2)[0]
    print("Pearson:  ",rval_p)
    #print("Spearman: ",rval_s)

    xy = np.vstack([plot_data1,plot_data2])
    z = stats.gaussian_kde(xy)(xy)

    #ax2.scatter(plot_data1,plot_data2,c=z,s=6)
    ax2.scatter(MODIS_data1['data'][np.where(~hash_data1.mask)], \
        MODIS_data2['data'][np.where(~hash_data1.mask)], s=6, \
        color='tab:blue', label='Inside Plume')
    ax2.scatter(MODIS_data1['data'][np.where(hash_data1.mask)], \
        MODIS_data2['data'][np.where(hash_data1.mask)], s=6, \
        color='tab:orange', label='Outside Plume')
    #ax2.scatter(nohash_data1.compressed(), nohash_data2.compressed(), s=6, \
    #    color='tab:orange', label='Outside Plume')
    ax2.set_xlabel('Ch. ' + str(MODIS_data1['channel']) +' ' + MODIS_data1['variable'])
    ax2.set_ylabel('Ch. ' + str(MODIS_data2['channel']) +' ' + MODIS_data2['variable'])
    ax2.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))
    ax2.legend()

    if(plot_ASOS_loc):
        print("Plotting ASOS site")
        plot_ASOS_locs(ax0,MODIS_data1)
        plot_ASOS_locs(ax1,MODIS_data1)

    if(save):
        cross_date = MODIS_data1['cross_date']
        file_time  = MODIS_data1['file_time']
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel1) + '_ch' + str(channel2) + '_' + pdate + '_2pannelscatter.png'
        plt.savefig(outname,dpi=300)
        print('Saved image',outname)
    else:
        plt.show()

    #return hash_data1, nohash_data1, MODIS_data1, MODIS_data2

def compare_MODIS_3scatter(date_str,channel0,channel1,channel2,channel3,\
        save=False, compare_OMI = False, avg_pixel = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel1 == 'red'):
        channel1 = 1
    elif(channel1 == 'green'):
        channel1 = 4
    elif(channel1 == 'blue'):
        channel1 = 3

    # Step 1: Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    MODIS_data0 = read_MODIS_channel(filename, channel0, zoom = True)
    MODIS_data1 = read_MODIS_channel(filename, channel1, zoom = True)
    MODIS_data2 = read_MODIS_channel(filename, channel2, zoom = True)
    MODIS_data3 = read_MODIS_channel(filename, channel3, zoom = True)

    print("Data0 max = ",np.max(MODIS_data0['data']), "  Data0 min = ",np.min(MODIS_data0['data']))
    print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",np.min(MODIS_data1['data']))
    print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",np.min(MODIS_data2['data']))
    print("Data3 max = ",np.max(MODIS_data3['data']), "  Data3 min = ",np.min(MODIS_data3['data']))

    max_ch = 350.

    # Determine where the smoke is located
    # ------------------------------------
    hash_data1, nohash_data1 = find_plume(filename) 

    tmp_data0 = np.copy(MODIS_data0['data'])
    tmp_data1 = np.copy(MODIS_data1['data'])
    tmp_data2 = np.copy(MODIS_data2['data'])
    tmp_data3 = np.copy(MODIS_data3['data'])

    tmp_data0 = np.ma.masked_where( (MODIS_data0['data'] > max_ch) | \
        (MODIS_data1['data'] > max_ch) | (MODIS_data2['data'] > max_ch) | \
        (MODIS_data3['data'] > max_ch), tmp_data0)
    tmp_data1 = np.ma.masked_where( (MODIS_data0['data'] > max_ch) | \
        (MODIS_data1['data'] > max_ch) | (MODIS_data2['data'] > max_ch) | \
        (MODIS_data3['data'] > max_ch), tmp_data1)
    tmp_data2 = np.ma.masked_where( (MODIS_data0['data'] > max_ch) | \
        (MODIS_data1['data'] > max_ch) | (MODIS_data2['data'] > max_ch) | \
        (MODIS_data3['data'] > max_ch), tmp_data2)
    tmp_data3 = np.ma.masked_where( (MODIS_data0['data'] > max_ch) | \
        (MODIS_data1['data'] > max_ch) | (MODIS_data2['data'] > max_ch) | \
        (MODIS_data3['data'] > max_ch), tmp_data3)
    tmp_lat0  = np.ma.masked_where( (MODIS_data0['data'] > max_ch) | \
        (MODIS_data1['data'] > max_ch) | (MODIS_data2['data'] > max_ch) | \
        (MODIS_data3['data'] > max_ch), MODIS_data0['lat'])
    tmp_lon0  = np.ma.masked_where( (MODIS_data0['data'] > max_ch) | \
        (MODIS_data1['data'] > max_ch) | (MODIS_data2['data'] > max_ch) | \
        (MODIS_data3['data'] > max_ch), MODIS_data0['lon'])

    ##!#plot_data0 = tmp_data0.compressed()
    ##!#plot_data1 = tmp_data1.compressed()
    ##!#plot_data2 = tmp_data2.compressed()
    ##!#plot_data3 = tmp_data3.compressed()
    ##!#plot_lat0  = tmp_lat0.compressed()
    ##!#plot_lon0  = tmp_lon0.compressed()

    cross_date = MODIS_data0['cross_date']
    file_time  = MODIS_data0['file_time']

    plt.close('all')
    if(compare_OMI):
        fig = plt.figure(figsize=(9,9))
        ax0 = fig.add_subplot(2,2,1)
        ax1 = fig.add_subplot(2,2,2)
        ax2 = fig.add_subplot(2,2,3)
        ax3 = fig.add_subplot(2,2,4)
    else:
        fig = plt.figure(figsize=(17,5))
        ax0 = fig.add_subplot(1,3,1)
        ax1 = fig.add_subplot(1,3,2)
        ax2 = fig.add_subplot(1,3,3)

    #xy = np.vstack([plot_data0,plot_data1])
    #z = stats.gaussian_kde(xy)(xy)
    #ax0.scatter(plot_data0,plot_data1,c=z,s=6)

    ax0.scatter(tmp_data0[np.where(~hash_data1.mask)], \
        tmp_data1[np.where(~hash_data1.mask)], s=6, \
        color='tab:blue')
    ax0.scatter(tmp_data0[np.where(hash_data1.mask)], \
        tmp_data1[np.where(hash_data1.mask)], s=6, \
        color='tab:orange')

    ax0.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
        str(np.average(channel_dict[str(channel0)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data0['variable'])
    ax0.set_ylabel('Ch. ' + str(MODIS_data1['channel']) +' [' + \
        str(np.average(channel_dict[str(channel1)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data1['variable'])
    plt.suptitle('Aqua MODIS ' + cross_date + ' ' + file_time)
    #ax0.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    #xy = np.vstack([plot_data0,plot_data2])
    #z = stats.gaussian_kde(xy)(xy)
    #ax1.scatter(plot_data0,plot_data2,c=z,s=6)
    ax1.scatter(tmp_data0[np.where(~hash_data1.mask)], \
        tmp_data2[np.where(~hash_data1.mask)], s=6, \
        color='tab:blue')
    ax1.scatter(tmp_data0[np.where(hash_data1.mask)], \
        tmp_data2[np.where(hash_data1.mask)], s=6, \
        color='tab:orange')

    ax1.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
        str(np.average(channel_dict[str(channel0)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data0['variable'])
    ax1.set_ylabel('Ch. ' + str(MODIS_data2['channel']) +' [' + \
        str(np.average(channel_dict[str(channel2)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data2['variable'])
    #ax0.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    #xy = np.vstack([plot_data0,plot_data3])
    #z = stats.gaussian_kde(xy)(xy)
    #ax2.scatter(plot_data0,plot_data3,c=z,s=6)
    ax2.scatter(tmp_data0[np.where(~hash_data1.mask)], \
        tmp_data3[np.where(~hash_data1.mask)], s=6, \
        color='tab:blue',label = 'Inside plume')
    ax2.scatter(tmp_data0[np.where(hash_data1.mask)], \
        tmp_data3[np.where(hash_data1.mask)], s=6, \
        color='tab:orange', label = 'Outside plume')
    ax2.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
        str(np.average(channel_dict[str(channel0)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data0['variable'])
    ax2.set_ylabel('Ch. ' + str(MODIS_data3['channel']) +' [' + \
        str(np.round(np.average(channel_dict[str(channel3)]['Bandwidth']),3)) \
        + ' μm] ' + MODIS_data3['variable'])
    #ax0.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    if(compare_OMI):
        print("Reading OMI data")
        data = h5py.File(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['omi'],'r')
        LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
        LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
        UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
        mask_LAT = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), LAT)
        mask_LON = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), LON)
        mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), UVAI)
        mask_LAT  = np.ma.masked_where((((LAT < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), mask_LAT)
        mask_LON  = np.ma.masked_where((((LAT < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), mask_LON)
        mask_UVAI = np.ma.masked_where((((LAT < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), mask_UVAI)


        #hash_data1, nohash_data1 = find_plume(filename) 
        # Colocate the masked OMI data with the MODIS channel0 data
        # ---------------------------------------------------------
        print("Colocating OMI data")

        # Remove the nans from the MODIS data, lats, and lons for both the
        # in-the-plume (hashed) and outside-the-plume (nohashed) data
        hash_plot_data0   = tmp_data0[np.where(~hash_data1.mask)].compressed()
        nohash_plot_data0 = tmp_data0[np.where(hash_data1.mask)].compressed()
        hash_plot_lat0   = tmp_lat0[np.where(~hash_data1.mask)].compressed()
        nohash_plot_lat0 = tmp_lat0[np.where(hash_data1.mask)].compressed()
        hash_plot_lon0   = tmp_lon0[np.where(~hash_data1.mask)].compressed()
        nohash_plot_lon0 = tmp_lon0[np.where(hash_data1.mask)].compressed()

        if(avg_pixel):
            print('averaging MODIS pixels')

            ##!## Determine the extreme bounds of the current 
            ##!#max_MODIS_lat = np.nanmax(tmp_lat0)
            ##!#min_MODIS_lat = np.nanmin(tmp_lat0)
            ##!#max_MODIS_lon = np.nanmax(tmp_lon0)
            ##!#min_MODIS_lon = np.nanmin(tmp_lon0)

            # Extract the OMI corner latitudes and longitudes
            # ----------------------------------------------- 
            crnr_LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/'+\
                'Geolocation Fields/FoV75CornerLatitude'][:,:,:]
            crnr_LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/'+\
                'Geolocation Fields/FoV75CornerLongitude'][:,:,:]

            # Pull out the OMI lat/lon corner pixels that are associated with
            # the masked OMI data
            # ---------------------------------------------------------------
            work_LAT   = mask_LAT.compressed()
            work_LON   = mask_LON.compressed()
            work_UVAI  = mask_UVAI.compressed()

            work_crnrLAT = crnr_LAT[np.where(~np.isnan(mask_UVAI)==True)[0],\
                np.where(~np.isnan(mask_UVAI)==True)[1],:] 
            work_crnrLON = crnr_LON[np.where(~np.isnan(mask_UVAI)==True)[0],\
                np.where(~np.isnan(mask_UVAI)==True)[1],:] 

            # Declare full arrays to hold the averaged MODIS data
            # Must be dimensioned to the shape of the OMI data
            # ---------------------------------------------------
            print("work UVAI shape = ",work_UVAI.shape)
            hash_avg_modis   = np.full(work_UVAI.shape,np.nan)
            hash_avg_mlat    = np.full(work_UVAI.shape,np.nan)
            hash_avg_mlon    = np.full(work_UVAI.shape,np.nan)
            nohash_avg_modis = np.full(work_UVAI.shape,np.nan)
            nohash_avg_mlat  = np.full(work_UVAI.shape,np.nan)
            nohash_avg_mlon  = np.full(work_UVAI.shape,np.nan)

            # Loop over screened out OMI pixels
            # ---------------------------------
            print("corners shape:",work_crnrLAT.shape)
            for ii in range(work_crnrLAT.shape[0]):

                # For each OMI pixel, find the in-plume and out-of-plume MODIS
                # pixels that are close to the OMI corner values
                # -------------------------------------------------------------
                hash_in = np.where((hash_plot_lat0 >= np.nanmin(work_crnrLAT[ii,:])) & \
                                   (hash_plot_lat0 <= np.nanmax(work_crnrLAT[ii,:])) & \
                                   (hash_plot_lon0 >= np.nanmin(work_crnrLON[ii,:])) & \
                                   (hash_plot_lon0 <= np.nanmax(work_crnrLON[ii,:])) )
                nohash_in = np.where((nohash_plot_lat0 >= np.nanmin(work_crnrLAT[ii,:])) & \
                                     (nohash_plot_lat0 <= np.nanmax(work_crnrLAT[ii,:])) & \
                                     (nohash_plot_lon0 >= np.nanmin(work_crnrLON[ii,:])) & \
                                     (nohash_plot_lon0 <= np.nanmax(work_crnrLON[ii,:])) )

                # Use the indices to pull out the MODIS pixels
                hash_modis_data   = hash_plot_data0[hash_in]
                hash_modis_lat    = hash_plot_lat0[hash_in]
                hash_modis_lon    = hash_plot_lon0[hash_in]
                nohash_modis_data = nohash_plot_data0[nohash_in]
                nohash_modis_lat  = nohash_plot_lat0[nohash_in]
                nohash_modis_lon  = nohash_plot_lon0[nohash_in]

                #print(work_crnrLAT[ii,:], work_crnrLON[ii,:], hash_modis_lat, \
                #    nohash_modis_lat)
            
                # Use the corner values to make a Polygon and determine if
                # each close MODIS pixel is within the OMI pixel
                # --------------------------------------------------------

                # Generate a Polygon using the OMI corner points
                points = []
                for x, y in zip(work_crnrLON[ii,:], work_crnrLAT[ii,:]):
                    points.append((x,y))
                omi_poly = Polygon(points)
              
                inside_hash = [omi_poly.contains(Point(hlon, hlat)) for \
                    hlon, hlat in zip(hash_modis_lon, hash_modis_lat)]
                inside_nohash = [omi_poly.contains(Point(nlon, nlat)) for \
                    nlon, nlat in zip(nohash_modis_lon, nohash_modis_lat)]
       
                if((len(inside_hash) > 0) & (True in inside_hash)):
                    hash_avg_modis[ii] = np.average(hash_modis_data[inside_hash])
                    hash_avg_mlat[ii]  = np.average(hash_modis_lat[inside_hash])
                    hash_avg_mlon[ii]  = np.average(hash_modis_lon[inside_hash])

                if((len(inside_nohash) > 0) & (True in inside_nohash)):
                    nohash_avg_modis[ii] = np.average(nohash_modis_data[inside_nohash])
                    nohash_avg_mlat[ii]  = np.average(nohash_modis_lat[inside_nohash])
                    nohash_avg_mlon[ii]  = np.average(nohash_modis_lon[inside_nohash])

                ##!## Loop over the close MODIS pixels and determine which ones
                ##!## are within the polygon
                ##!## ---------------------------------------------------------
                ##!#for jj in range(len(hash_modis_data)):

                ##!#    # Average each MODIS value that is within the polygon
                ##!#    # into a value, and insert the value into the avg
                ##!#    # modis array for hash and nohash.
                ##!#    # --------------------------------------------------- 
                ##!#    mpoint = Point(hash_modis_lon[jj], hash_modis_lat[jj])

                ##!#    inside_vals = [omi_poly.contains(mpoint)  ] 

            # end screened corner data loop

            # Plot the screened OMI data vs the averaged MODIS data
            # -----------------------------------------------------

            #print('hashed averaged MODIS data = ',hash_avg_modis)
            #print('nohashed averaged MODIS data = ',nohash_avg_modis)

            print(work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask])
            
            ax3.scatter(np.ma.masked_invalid(hash_avg_modis).compressed(), \
                work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask], s=6, \
                color='tab:blue')
            #ax3.scatter(np.ma.masked_invalid(nohash_avg_modis).compressed(), \
            #    work_UVAI[~np.ma.masked_invalid(nohash_avg_modis).mask], s=6, \
            #    color='tab:orange')
            #ax3.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
            #    color='tab:orange', label='Outside Plume')

            ax3.set_xlabel('Averaged Ch. ' + str(MODIS_data0['channel']) +' [' + \
                str(np.average(channel_dict[str(channel0)]['Bandwidth'])) \
                + ' μm] ' + MODIS_data0['variable'])
            ax3.set_ylabel('OMI UVAI')

            # Plot hash_avg_modis against mask_UVAI 

            # Plot nohash_avg_modis against mask_UVAI 
 
        else:
            hash_match_OMI   = np.full(hash_plot_data0.shape,-9.)
            hash_match_LAT   = np.full(hash_plot_lat0.shape,-9.)
            hash_match_LON   = np.full(hash_plot_lon0.shape,-9.)
            nohash_match_OMI = np.full(nohash_plot_data0.shape,-9.)
            nohash_match_LAT = np.full(nohash_plot_lat0.shape,-9.)
            nohash_match_LON = np.full(nohash_plot_lon0.shape,-9.)

            print(hash_plot_data0.shape)
            for ii in range(hash_match_OMI.shape[0]):
                # Find the gridpoint in the gridded lat/lon data that 
                # corresponds to the station at slat and slon
                # ---------------------------------------------------- 
                o_idx = nearest_gridpoint(hash_plot_lat0[ii], hash_plot_lon0[ii],\
                    mask_LAT, mask_LON)

                if(len(o_idx[0]) > 1):
                    o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
                hash_match_OMI[ii] = mask_UVAI[o_idx]
                hash_match_LAT[ii] = mask_LAT[o_idx] 
                hash_match_LON[ii] = mask_LON[o_idx] 

            print(nohash_plot_data0.shape)
            for ii in range(nohash_match_OMI.shape[0]):
                # Find the gridpoint in the gridded lat/lon data that 
                # corresponds to the station at slat and slon
                # ---------------------------------------------------- 
                o_idx = nearest_gridpoint(nohash_plot_lat0[ii], nohash_plot_lon0[ii],\
                    mask_LAT, mask_LON)

                if(len(o_idx[0]) > 1):
                    o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
                nohash_match_OMI[ii] = mask_UVAI[o_idx]
                nohash_match_LAT[ii] = mask_LAT[o_idx] 
                nohash_match_LON[ii] = mask_LON[o_idx] 
 
            #xy = np.vstack([plot_data0,match_OMI])
            #z = stats.gaussian_kde(xy)(xy)
            #ax3.scatter(plot_data0,match_OMI,c=z,s=6)

            ax3.scatter(hash_plot_data0, hash_match_OMI, s=6, \
                color='tab:blue')
            #ax3.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
            #    color='tab:orange')

            ax3.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
                str(np.average(channel_dict[str(channel0)]['Bandwidth'])) \
                + ' μm] ' + MODIS_data0['variable'])
            ax3.set_ylabel('OMI UVAI')
        
        data.close()

    # End compare_OMI
       
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    if(compare_OMI):
        plt.legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0, 0.01, 1, 1),\
            bbox_transform = plt.gcf().transFigure, ncol=2)
    else:
        plt.legend(lines, labels, loc = 'right', bbox_to_anchor = (0, 0., 1.0, 1),\
            bbox_transform = plt.gcf().transFigure, ncol=1)

    if(save):
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel0) + '_vs_ch' + \
            str(channel1) + '_ch' + str(channel2) + '_ch' + \
            str(channel3) + '_' + pdate + '_3scatter.png'
        if(compare_OMI):
            outname = 'modis_compare_ch' + str(channel0) + '_vs_ch' + \
                str(channel1) + '_ch' + str(channel2) + '_ch' + \
                str(channel3) + '_omi_' + pdate + '_3scatter.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
        plt.show()


# Compare colocated MODIS and ob data for two dates
def colocate_comparison(date1, date2, channel = 31):
    # Read in MODIS data for both cases
    # ---------------------------------
    dt_date_str1 = datetime.strptime(date1,"%Y%m%d%H%M")
    dt_date_str2 = datetime.strptime(date2,"%Y%m%d%H%M")
    filename1 = plot_limits_dict[dt_date_str1.strftime('%Y-%m-%d')][dt_date_str1.strftime('%H%M')]['modis']
    filename2 = plot_limits_dict[dt_date_str2.strftime('%Y-%m-%d')][dt_date_str2.strftime('%H%M')]['modis']

    MODIS_data1 = read_MODIS_channel(filename1, channel, zoom = True)
    MODIS_data2 = read_MODIS_channel(filename2, channel, zoom = True)

    # Use the colocation code to extract matching data
    # ------------------------------------------------
    compare_data1 = nearest_grid_values(MODIS_data1)
    compare_data2 = nearest_grid_values(MODIS_data2)

    # Plot the data
    # -------------

    plt.close('all')
    fig, ax = plt.subplots()
    print(compare_data1['stations'])
    ax.plot(np.arange(len(compare_data1['stations'])), compare_data1['stn_data'], color='tab:blue', label = '05 Stations')
    ax.plot(np.arange(len(compare_data1['stations'])), compare_data2['stn_data'], color = 'tab:blue', linestyle = '--', label = '06 Stations')

    ax2 = ax.twinx()
    ax2.plot(np.arange(len(compare_data1['stations'])),compare_data1['mds_data'], color='tab:orange', label = '05 Modis')
    ax2.plot(np.arange(len(compare_data1['stations'])),compare_data2['mds_data'], color='tab:orange', linestyle = '--', label = '06 Modis')
  
    ax.set_xticks(np.arange(len(compare_data1['stations']))) 
    ax.set_xticklabels(compare_data1['stations'])
 
    plt.show()
 
def plot_modis_data(modis_data,minlat=60,tind=0,zoom = None,save=False):

    data = modis_data['data'][tind,:,:]
    colormap = plt.cm.jet
    mask_data = np.ma.masked_where(data == -9., data)
    mask_data = np.ma.masked_invalid(mask_data)

    print(modis_data['titles'][0])
    splitter = modis_data['titles'][0].split('/')[-1]
    plot_date = datetime.strptime(splitter[1:5],'%Y') + \
        relativedelta(days = int(splitter[5:8])-1)

    if(minlat > 50):
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0)
    else:
        mapcrs = ccrs.Robinson()

    plot_lat, plot_lon = np.meshgrid(modis_data['lat'],modis_data['lon'])

    plt.close()
    fig1 = plt.figure(figsize=(8,8))
    if(zoom == None):
        ax = plt.axes(projection = mapcrs)
        ax.set_extent([-180,180,minlat,90],datacrs)
        saver = ''
    else:
        ax = plt.axes(projection = proj_dict[zoom])
        #ax.set_extent([17,49,68,75],datacrs)
        ax.set_extent(zoom_dict[zoom],datacrs)
        saver = '_'+zoom 
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(plot_lon,plot_lat,mask_data.T,\
            transform = datacrs, norm = cm.LogNorm(vmin = modis_data['pticks'][0],\
            vmax = modis_data['pticks'][-1]),cmap = colormap)
    #CS = ax.contour(longitude,latitude,smooth_thick,[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],transform = datacrs)
    
    # Adjust and make it look good
    #ax.add_feature(cfeature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
    ax.set_title('MODIS '+modis_data['ptitle']+'\n'+plot_date.strftime('%Y%m%d'))
    cbar = plt.colorbar(mesh,ticks = modis_data['pticks'],orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905, label=modis_data['label'])
    cbar.ax.set_xticklabels(modis_data['ptick_labels'])
    #ax.set_xlim(-4170748.535086173,4167222.438879491)
    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
    #ax.set_title(datetime.strftime(base_dtm,'%B %Y') + ' CryoSat-2 Data')

    if(save == True):
        outname = 'modis_'+modis_data['grabber'] + '_' + plot_date.strftime('%Y%m%d') + saver + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Plot a histogram of 25 x 25 km grid data
def plot_modis_hist(modis_data,tind,variable,bins = 100,save = False):
    data = modis_data[variable][tind,:,:]
    mask_data = np.ma.masked_where(data < -999., data)

    # Make a datetime object from the structure date
    # ----------------------------------------------
    if(len(modis_data['dates'][tind]) == 6):
        str_fmt = '%Y%m'
    elif(len(modis_data['dates'][tind]) == 8):
        str_fmt = '%Y%m%d'

    base_dtm = datetime.strptime(modis_data['dates'][tind],str_fmt)

    # Set up the x axis label
    # -----------------------
    if(variable == 'ice_thick'):
        xlabel = 'Derived Sea Ice Thickness [m]' 
    elif(variable == 'ice_con'):
        xlabel = 'Sea Ice Concentration [%]' 

    bin_heights,bin_borders = np.histogram(mask_data.compressed(),bins=bins)
    bin_widths = np.diff(bin_borders)
    bin_centers = bin_borders[:-1] + bin_widths / 2

    t_init = models.Gaussian1D()
    fit_t = fitting.LevMarLSQFitter()
    t = fit_t(t_init, bin_centers, bin_heights)

    print('Amplitude: ',np.round(t.amplitude.value,3))
    print('Mean:      ',np.round(t.mean.value,3))
    print('StDev:     ',np.round(t.stddev.value,3))

    x_interval_for_fit = np.linspace(bin_centers[0],bin_centers[-1],100)
    plt.close()
    plt.figure()
    plt.bar(bin_centers,bin_heights,width=bin_widths,label='histogram')
    plt.plot(x_interval_for_fit,t(x_interval_for_fit),label='fit',c='tab:red')
    plt.xlabel(xlabel)
    plt.ylabel('Counts')
    plt.title(datetime.strftime(base_dtm,'%B %Y') + ' CryoSat-2 Data')
    plt.legend()
    #plt.close()
    #plt.hist(mask_data.compressed(),bins=bins)
    #plt.title(modis_data['dates'][tind] + ' '+variable)
    #plt.xlabel(xlabel)
    #plt.ylabel('Counts')
    #plt.title(datetime.strftime(base_dtm,'%B %Y') + ' CryoSat-2 Data')

    if(save == True):
        outname = 'modissat2_' + variable + '_' + datetime.strftime(base_dtm,'%Y%m%d') + '_hist.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# plot_grid_data generates a plot of the /
def plot_grid_data(modis_data,t_ind,pvar,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    if(pvar=='grid_con'):
        plabel = "Percent Ice Concentration"
    elif(pvar=='grid_thick'):
        plabel = "Sea Ice Thickness"

    local_grid_modis = np.copy(modis_dict['grid_modis_conc'][t_ind,:,:])
    local_grid_modis_bad = np.copy(modis_dict['grid_modis_conc'][t_ind,:,:])

    local_grid_modis[local_grid_modis==-999.] = np.nan
    local_grid_modis_bad[local_grid_modis!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_modis)
    plot_land_data = ma.masked_invalid(local_grid_modis_bad)

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
    ax.set_title('Gridded NSIDC Sea Ice Concentration\n'+modis_dict['titles'][t_ind])
    if(save==True):
        outname = 'modis_conc_gridded_200012_201812'+modis_dict['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
   
# Plot the trends for the 25x25 km grid data 
def plot_trend(modis_data,variable):
    data = modis_data[variable][:,:]
    colormap = plt.cm.bwr
    mask_data = np.ma.masked_where(data == np.nan, data)

    plt.close()
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180,180,60,90])
    mesh = ax.pcolormesh(modis_data['lon'],modis_data['lat'],mask_data,\
            transform = datacrs, cmap = colormap,vmin=min_dict[variable],\
            vmax=max_dict[variable])
    #CS = ax.contour(modis_data['lon'],modis_data['lat'],np.ma.masked_where(modis_data['thick_trends'][:,:] == np.nan,modis_data['thick_trends'][:,:]),\
    #        np.linspace(-0.5,0.5,5),transform = datacrs)
    
    # Adjust and make it look good
    ax.add_feature(cfeature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
    cbar = plt.colorbar(mesh,ticks = tick_dict[variable],orientation='horizontal',pad=0,aspect=50,label=variable)
    cbar.ax.set_xticklabels(tick_label_dict[variable])
    ax.set_xlim(-4170748.535086173,4167222.438879491)
    ax.set_ylim(-2913488.8763307533,2943353.899053069)
    ax.set_title(variable)
    plt.show()

def plot_grid_trend(modis_dict,adjusted=False,save=False):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)

    local_grid_modis = np.copy(modis_dict['grid_modis'])
    local_grid_modis_bad = np.copy(modis_dict['grid_modis'])

    local_grid_modis[local_grid_modis==-999.] = np.nan
    local_grid_modis_bad[local_grid_modis!=-999.] = np.nan
    plot_good_data = ma.masked_invalid(local_grid_modis)
    plot_land_data = ma.masked_invalid(local_grid_modis_bad)

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
    if(modis_dict['season_adder']!=''):
        ax.set_title('NSIDC Sea Ice Concentration'+modis_dict['season_adder'].title()+\
            ' Trends\nJan 2001 to Dec 2018')
    else:
        ax.set_title('NSIDC Sea Ice Concentration Trends\nJan 2001 to Dec 2018')
    if(save==True):
        outname = 'modis_trend_gridded_200101_201812'+modis_dict['file_adder']+file_adder+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
    
def plot_grid_time_series(modis_dict,lat_ind,lon_ind,thielsen=False):
    inseason = modis_dict['season_adder'].strip()
    temp_data = modis_dict['grid_modis_conc'][:,lat_ind,lon_ind]
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
    plt.title('Gridded Ice Data: '+str(modis_dict['grid_lat'][lat_ind])+'x'+\
                str(modis_dict['grid_lon'][lon_ind])+'\n'+inseason.title()+\
                ' season of each year)')
    plt.plot(temp_data,label='observations') 
    plt.plot(regress_y,'--',label='trend')
    plt.ylabel('Ice Concentration [%]')
    plt.legend()
    outname = 'nsidc_grid_time_series_'+inseason+'_'+str(int(modis_dict['grid_lat'][lat_ind]))+'x'+\
                str(int(modis_dict['grid_lon'][lon_ind]))+'.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)   
    plt.show()

# tind indicates the reference time to check concentrations and thicknesses before
# plotting time series 
# Syntax to run:
# >>> import sys
# >>> sys.path.append('/home/bsorenson/Research/Ice_analysis')
# >>> from CryoSat2Lib import *
# >>> from IceLib import *
# >>> ice_data = read_ice('all')
# >>> modis_data = read_modis('all','201011','201911')
# >>> albedo_effect_test(modis_data,ice_data,11,1.8)
# This makes the figure "modissat2_melt_time_20120301.png"
def albedo_effect_test(modis_data,ice_data,tind,start_thick,save=False):
#def albedo_effect_test(modis_data,tind,start_thick):

    # Find the time index in the NSIDC ice structure that matches
    # with the time index in the CryoSat2 structure
    ice_ind = 0
    for xi in range(len(ice_data['titles'])):
        if(ice_data['titles'][xi][7:13] == modis_data['dates'][tind]):
            ice_ind = xi

    low_thick   = start_thick - 0.01
    high_thick = start_thick + 0.01

    # Identify regions with the thickness desired by 'start_thick'
    test_con = ice_data['data'][ice_ind,:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
                    (modis_data['ice_thick'][tind,:,:] < high_thick)]
    test_lats = ice_data['lat'][:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
                    (modis_data['ice_thick'][tind,:,:] < high_thick)]
    ##test_con = modis_data['ice_con'][tind,:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
    ##                (modis_data['ice_thick'][tind,:,:] < high_thick)]
    ice_cons = np.zeros((12,test_con.size))

    colors = plt.cm.plasma(np.linspace(0,1,test_con.size))

    # Extract time series for the year (assuming October is the reference
    # month, the next 7 indices
    # -------------------------------------------------------------------
    count = 0
    for ti in range(ice_ind,ice_ind+12):
    #for ti in range(tind,tind+12):
        ice_cons[count,:] = ice_data['data'][ti,:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
                    (modis_data['ice_thick'][tind,:,:] < high_thick)] 
        #ice_cons[count,:] = modis_data['ice_con'][ti,:,:][(modis_data['ice_thick'][tind,:,:] >= low_thick) & \
        #            (modis_data['ice_thick'][tind,:,:] < high_thick)] 
        count+=1

    # Go through and calculate the number of months to melt 50% of the
    # beginning ice concentration. For each grid point, plot starting
    # ice concentration on the x axis and the number of months to melt
    # on the y axis
    # ----------------------------------------------------------------

    # Set up colorbar to color the lines by latitude
    # ----------------------------------------------
    norm = mpl.colors.Normalize(vmin = test_lats.min(), vmax = test_lats.max())
    cmap = mpl.cm.ScalarMappable(norm = norm, cmap = mpl.cm.plasma)
    cmap.set_array([])
    c = np.arange(int(test_lats.min()),int(test_lats.max()))

    # Set up a base datetime object for the x axis
    base_dtm = datetime.strptime(modis_data['dates'][tind],'%Y%m')
    dates = [base_dtm + timedelta(days = int(xval * 31)) for xval in np.arange(12)]

    # Plot the total time series data. Also, calculate
    # the difference in ice concentration over the next four months
    # -------------------------------------------------------------
    diff = np.zeros(ice_cons.shape[1])
    melt_months = np.zeros(ice_cons.shape[1])
    melt_limit = 0  # when zero, means complete melting

    fig, axs = plt.subplots(2,dpi=100)
    fig.set_size_inches(9,5.5)
    for xi in range(ice_cons.shape[1]):
        # Plot the time series for the current grid box
        axs[0].plot(dates,ice_cons[:,xi],c = cmap.to_rgba(test_lats[xi]))

        # Calculate the difference
        diff[xi] = ice_cons[0,xi] - ice_cons[4,xi] 

        # Determine how long it takes 
        try:
            melt_months[xi] = np.argwhere(ice_cons[:,xi] <= melt_limit)[0][0]
        except:
            melt_months[xi] = -9 

    axs[0].xaxis.set_major_formatter(DateFormatter('%b %Y'))    
    
    #fig.colorbar(cmap, ticks = c)
    axs[0].set_ylabel('Ice concentration [%]')
    axs[0].set_title(datetime.strftime(base_dtm,'%B %Y') + '\n' + str(start_thick) + ' m Thickness ')

    # Plot the differences
    # ---------------------
    #fig2,ax2 = plt.subplots(dpi=100)
    p_scat = axs[1].scatter(ice_cons[0,:][melt_months > -9],melt_months[melt_months > -9],cmap = mpl.cm.plasma,c = cmap.to_rgba(test_lats[melt_months > -9]))
    axs[1].set_ylabel('Months needed to melt\n below ' + str(melt_limit) + '% ice concentration')
    #axs[1].set_ylim(bottom = 0)
    #p_scat = axs[1].scatter(ice_cons[0,:],diff[:],cmap = mpl.cm.plasma,s = 12,c = cmap.to_rgba(test_lats[:]))
    #axs[1].set_ylabel('Ice concentration decrease \nfrom '+datetime.strftime(base_dtm,'%B') + \
    #    ' to ' + datetime.strftime(base_dtm + timedelta(days = int(4 * 31)),'%B') + ' [%]')
    cbar_ax = fig.add_axes([0.92,0.11,0.02,0.77])
    fig.colorbar(cmap,cax = cbar_ax,ticks = c[::2],label = 'Latitude [$^{o}$]')
    #fig2.colorbar(cmap, ticks = c)
    axs[1].set_xlabel(datetime.strftime(base_dtm,'%B %Y') + ' Concentration')

    if(save == True):
        outname = 'modissat2_melt_time_' + datetime.strftime(base_dtm,'%Y%m%d') + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()
    return


    # Plot 
    #mask_data = np.ma.masked_where(data < -999., data)

    if(variable[:4] == 'grid'):
        lat_vals = modis_data['grid_lat']
        lon_vals = modis_data['grid_lon']
    else:
        lat_vals = modis_data['lat']
        lon_vals = modis_data['lon']

    plt.close()
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180,180,60,90])
    mesh = ax.pcolormesh(lon_vals,lat_vals,mask_data,\
            transform = datacrs, cmap = colormap,vmin=min_dict[variable],\
            vmax=max_dict[variable])
    #CS = ax.contour(longitude,latitude,smooth_thick,[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],transform = datacrs)
    
    # Adjust and make it look good
    ax.add_feature(cfeature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
    cbar = plt.colorbar(mesh,ticks = tick_dict[variable],orientation='horizontal',pad=0,aspect=50,label=variable)
    cbar.ax.set_xticklabels(tick_label_dict[variable])
    ax.set_xlim(-4170748.535086173,4167222.438879491)
    ax.set_ylim(-2913488.8763307533,2943353.899053069)
    ax.set_title(modis_data['titles'][tind])
    plt.show()

