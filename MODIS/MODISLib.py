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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
from scipy.stats import pearsonr,spearmanr
import subprocess
from scipy import stats
from pyhdf import SD

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
    "2021-08-05": {
        '2125': {
            'Lat': [39.5, 42.5],
            'Lon': [-121.5, -119.5]
        }
    },
    "2021-08-06": {
        '2025': {
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0]
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
 
        print('yay')    

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

def plot_MODIS_channel(filename,channel,zoom=True):

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
    ax = plt.axes(projection = ccrs.LambertConformal())

    mesh = ax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
        MODIS_data['data'],cmap = MODIS_data['colors'], shading='auto', \
        transform = ccrs.PlateCarree()) 

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

def compare_MODIS_3panel(filename,channel1,channel2,channel3,zoom=True,save=False):

    if(channel1== 'red'):
        channel1= 1
    elif(channel1== 'green'):
        channel1= 4
    elif(channel1== 'blue'):
        channel1= 3

    # Step 1: Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    MODIS_data1 = read_MODIS_channel(filename, channel1, zoom = False)
    MODIS_data2 = read_MODIS_channel(filename, channel2, zoom = False)
    MODIS_data3 = read_MODIS_channel(filename, channel3, zoom = False)

    print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",np.min(MODIS_data1['data']))
    print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",np.min(MODIS_data2['data']))
    print("Data3 max = ",np.max(MODIS_data3['data']), "  Data3 min = ",np.min(MODIS_data3['data']))

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
    fig = plt.figure(figsize=(14.5,5))
    ax0 = fig.add_subplot(1,3,1,projection = mapcrs)
    ax1 = fig.add_subplot(1,3,2,projection = mapcrs)
    ax2 = fig.add_subplot(1,3,3,projection = mapcrs)

    # Step 3: Plot the MODIS channel data in the first 2 panels
    # ---------------------------------------------------------
    mesh0 = ax0.pcolormesh(MODIS_data1['lon'],MODIS_data1['lat'],\
        MODIS_data1['data'],cmap = MODIS_data1['colors'], shading='auto', \
        vmin = np.nanmin(cpy_1), vmax = np.nanmax(cpy_1), transform = datacrs) 

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
        vmin = np.nanmin(cpy_2), vmax = np.nanmax(cpy_2), transform = datacrs) 

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

    # Plot channel 3
    mesh2 = ax2.pcolormesh(MODIS_data3['lon'],MODIS_data3['lat'],\
        MODIS_data3['data'],cmap = MODIS_data3['colors'], shading='auto', \
        vmin = np.nanmin(cpy_3), vmax = np.nanmax(cpy_3), transform = datacrs) 

    cbar2 = plt.colorbar(mesh2,ax=ax2,orientation='vertical',\
        pad=0.03,label=MODIS_data3['variable'])
    #cbar1 = plt.colorbar(mesh1,cax = ax1,orientation='horizontal',pad=1,\
    #    aspect=51,shrink = 1.851)
    #cbar1.ax.tick_params(labelsize=14)
    #cbar1.set_label(MODIS_data2['variable'],fontsize=16,weight='bold')
    
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

    cross_date = MODIS_data1['cross_date']
    file_time  = MODIS_data1['file_time']
    if(save):
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel1) + '_vs_ch' + \
            str(channel2) + '_ch' + str(channel3) + '_' + pdate + '_3panel.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
        plt.show()
    

def compare_MODIS_channels(filename,channel1,channel2,zoom=True,save=False):

    if(channel1 == 'red'):
        channel1 = 1
    elif(channel1 == 'green'):
        channel1 = 4
    elif(channel1 == 'blue'):
        channel1 = 3

    # Step 1: Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    MODIS_data1 = read_MODIS_channel(filename, channel1, zoom = zoom)
    MODIS_data2 = read_MODIS_channel(filename, channel2, zoom = zoom)

    print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",np.min(MODIS_data1['data']))
    print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",np.min(MODIS_data2['data']))

    print(MODIS_data1['data'].shape, MODIS_data2['data'].shape)

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

    
    # Step 4: Plot scatter MODIS channel comparison in third panel
    # ------------------------------------------------------------
    #if(MODIS_data1['variable'] == 'Blackbody Temperature [K]'):
    #    max_ch = 350.
    #if(MODIS_data2['variable'] == 'Blackbody Temperature [K]'):
    max_ch = 350.

    tmp_data1 = np.copy(MODIS_data1['data'])
    tmp_data2 = np.copy(MODIS_data2['data'])


    tmp_data1 = np.ma.masked_where( (MODIS_data1['data'] > max_ch) | \
        (MODIS_data2['data'] > max_ch), tmp_data1)
    tmp_data2 = np.ma.masked_where( (MODIS_data1['data'] > max_ch) | \
        (MODIS_data2['data'] > max_ch), tmp_data2)


    plot_data1 = tmp_data1.compressed()
    plot_data2 = tmp_data2.compressed()

    rval_p = pearsonr(plot_data1,plot_data2)[0]
    #rval_s = spearmanr(tmp_data1,tmp_data2)[0]
    print("Pearson:  ",rval_p)
    #print("Spearman: ",rval_s)

    xy = np.vstack([plot_data1,plot_data2])
    z = stats.gaussian_kde(xy)(xy)

    ax2.scatter(plot_data1,plot_data2,c=z,s=6)
    ax2.set_xlabel('Ch. ' + str(MODIS_data1['channel']) +' ' + MODIS_data1['variable'])
    ax2.set_ylabel('Ch. ' + str(MODIS_data2['channel']) +' ' + MODIS_data2['variable'])
    ax2.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    if(save):
        cross_date = MODIS_data1['cross_date']
        file_time  = MODIS_data1['file_time']
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel1) + '_ch' + str(channel2) + '_' + pdate + '_3panel.png'
        plt.savefig(outname,dpi=300)
        print('Saved image',outname)
    else:
        plt.show()

def compare_MODIS_3scatter(filename,channel0,channel1,channel2,channel3,save=False):

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

    plot_data0 = tmp_data0.compressed()
    plot_data1 = tmp_data1.compressed()
    plot_data2 = tmp_data2.compressed()
    plot_data3 = tmp_data3.compressed()

    cross_date = MODIS_data0['cross_date']
    file_time  = MODIS_data0['file_time']

    plt.close('all')
    fig = plt.figure(figsize=(5,13))
    ax0 = fig.add_subplot(3,1,1)
    ax1 = fig.add_subplot(3,1,2)
    ax2 = fig.add_subplot(3,1,3)

    xy = np.vstack([plot_data0,plot_data1])
    z = stats.gaussian_kde(xy)(xy)
    ax0.scatter(plot_data0,plot_data1,c=z,s=6)
    ax0.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
        str(np.average(channel_dict[str(channel0)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data0['variable'])
    ax0.set_ylabel('Ch. ' + str(MODIS_data1['channel']) +' [' + \
        str(np.average(channel_dict[str(channel1)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data1['variable'])
    ax0.set_title('Aqua MODIS ' + cross_date + ' ' + file_time)
    #ax0.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    xy = np.vstack([plot_data0,plot_data2])
    z = stats.gaussian_kde(xy)(xy)
    ax1.scatter(plot_data0,plot_data2,c=z,s=6)
    ax1.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
        str(np.average(channel_dict[str(channel0)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data0['variable'])
    ax1.set_ylabel('Ch. ' + str(MODIS_data2['channel']) +' [' + \
        str(np.average(channel_dict[str(channel2)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data2['variable'])
    #ax0.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    xy = np.vstack([plot_data0,plot_data3])
    z = stats.gaussian_kde(xy)(xy)
    ax2.scatter(plot_data0,plot_data3,c=z,s=6)
    ax2.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
        str(np.average(channel_dict[str(channel0)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data0['variable'])
    ax2.set_ylabel('Ch. ' + str(MODIS_data3['channel']) +' [' + \
        str(np.average(channel_dict[str(channel3)]['Bandwidth'])) \
        + ' μm] ' + MODIS_data3['variable'])
    #ax0.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    if(save):
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel0) + '_vs_ch' + \
            str(channel1) + '_ch' + str(channel2) + '_ch' + \
            str(channel3) + '_' + pdate + '_3scatter.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
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

