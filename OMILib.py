#!/usr/bin/env python
"""
  NAME:

  PURPOSE:
    Plot trends in OMI-measured aerosol index across the northern hemisphere.
    The average AI values are read from files produced by 
    Average_AerosolIndexCalculator.pro (included in ai_trend_codes_YYMMDD.tar.gz)

  PYTHON VERSION:
    2.6.6

  MODULES:
    - Custom module AILib
    - Matplotlib
    - mpl_toolkits.basemap
    - datetime
    - sys
    - numpy
    
"""

import numpy as np
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.colors as color
from matplotlib.colors import rgb2hex,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
from scipy import stats
from scipy.interpolate import interp2d
from netCDF4 import Dataset
import gzip
import h5py
import subprocess
from scipy.stats import pearsonr,spearmanr

# Class MidpointNormalize is used to center the colorbar on zero
class MidpointNormalize(Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

# Function onclick performs an operation at the location of a click on the
# map. 
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    #print ix, iy  
    mapLon, mapLat = m(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = [OMI_data[dictkey][date]['avg'] for date in                           \
        sorted(OMI_data[dictkey].keys())]
    dates = sorted(OMI_data[dictkey].keys())
    x_vals = np.arange(0,len(OMI_data[dictkey].keys()))

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals,avgs)
    print(dictkey, slope*len(x_vals))
    #print slope/len(OMI_data[dictkey].keys())
    regress_y = x_vals*slope+intercept

    #The slope
    S=0
    sm=0
    nx = len(avgs)
    num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
    Sn=np.zeros(num_d)
    for si in range(0,nx-1):
        for sj in range(si+1,nx):
            # Find the slope between the two points
            Sn[sm] = (avgs[si]-avgs[sj])/(si-sj) 
            sm=sm+1
        # Endfor
    # Endfor
    Snsorted=sorted(Sn)
    sm=int(num_d/2.)
    if(2*sm    == num_d):
        tsslope=0.5*(Snsorted[sm]+Snsorted[sm+1])
    if(2*sm+1 == num_d): 
        tsslope=Snsorted[sm+1]
    regress_ts = x_vals*tsslope+intercept

    # Convert the dates to datetime objects
    ddates = [datetime.strptime(date,'%Y%m') for date in dates] 

    label_dates = [ddate.strftime('%b %Y') for ddate in ddates]

    tickNumber = 12

    fig1 = plt.figure()
    plt.plot(avgs)
    plt.plot(x_vals,regress_y,linestyle='--',color='red',label='Control')
    plt.plot(x_vals,regress_ts,linestyle='--',color='blue',label='Thiel-Sen')
    plt.title(dictkey)
    plt.xticks(np.arange(0,len(avgs))[::-int(len(avgs)/tickNumber)],label_dates[::\
        -int(len(avgs)/tickNumber)],rotation=45)
    plt.ylabel('Average Aerosol Index')
    plt.legend()
    plt.show()

def onclick_climo(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = m(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([OMI_data[dictkey][date] for date in sorted(OMI_data[dictkey].keys())])
    ai_avg = np.average(avgs)
    print(dictkey, ai_avg)

def readOMI(inputfile,start_date,end_date,key=None):
    global OMI_data
    OMI_data = {}

    if(key is not None):
        OMI_data[key]={}

    if(inputfile.strip().split('/')[-1].split('.')[-1]=='gz'):
        f = gzip.open(inputfile,'rb')
    else:
        f = open(inputfile,'r')
    #with open(inputfile,'r') as f:
    # Skip the header line
    for line in f:
        templine = line.strip().split()
        if(len(templine)>1):
            if((int(templine[0])>=start_date) & (int(templine[0])<=end_date)):
                if(key is not None):
                    if(templine[1]==key):
                        OMI_data[key][templine[0].decode('utf-8')] = {}
                        OMI_data[key][templine[0].decode('utf-8')]['avg']=float(templine[2])
                        OMI_data[key][templine[0].decode('utf-8')]['#_obs']=int(templine[3])
                else:
                    # If the current lat/lon pair are not found in the dictionary's
                    # keys, then make a new subdictionary for it.
                    if(templine[1].decode('utf-8') not in OMI_data.keys()):
                        OMI_data[templine[1].decode('utf-8')] = {}
                    # If the current lat/lon pair are already in the dictionary's
                    # keys, then add the new data to the subdictionary
                    OMI_data[templine[1].decode('utf-8')][templine[0].decode('utf-8')]={}
                    OMI_data[templine[1].decode('utf-8')][templine[0].decode('utf-8')]['avg']=float(templine[2])
                    OMI_data[templine[1].decode('utf-8')][templine[0].decode('utf-8')]['#_obs']=int(templine[3])
    f.close()    

    return OMI_data

def readOMI_single_swath(plot_time,row_max):
    n_p = 1440
    nl = 720
    lonmin = -180
    lonmax = 180
    latmin = 60.
    latmax = 90.
    # Set up values for gridding the AI data
    lat_gridder = latmin * 4.
    
    lat_ranges = np.arange(latmin,90.1,0.25)
    lon_ranges = np.arange(-180,180.1,0.25)

    g_NRAD_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_NRAD_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_NRAD_500 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_REFL_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_REFL_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_SALB_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    g_UVAI_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_NRAD_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_NRAD_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_NRAD_500 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_REFL_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_REFL_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_SALB_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_UVAI_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

    algae = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

    year = plot_time[:4]
    date = plot_time[4:8]
    if(len(plot_time)==13):
        time = plot_time[9:]
    elif(len(plot_time)==12):
        time = plot_time[8:]
    else:
        time = ''
    base_path = '/home/bsorenson/data/OMI/H5_files/'
    total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+year+'m'+date+'t'+time+'*.he5',\
              shell=True).decode('utf-8').strip().split('\n')

    for fileI in range(len(total_list)):
        # read in data directly from HDF5 files
        data = h5py.File(total_list[fileI],'r')
        print(total_list[fileI])
    
        LAT     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
        LON     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
        NRAD    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/NormRadiance']
        REFL    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
        SALB    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
        UVAI    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
        XTRACK  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
    
        #albedo = ALBEDO[:,:,0]   
        #reflectance = REFLECTANCE[:,:,0]   
        counter = 0
        #AI = AI[:,:,0]   
        # Loop over the values and rows 
        #for i in range(0,int(CBA2)):
        #for i in range(albedo.shape[0]):
        for i in range(NRAD.shape[0]):
            for j in range(0,row_max):
                #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
                if(NRAD[i,j,0]>-2e5):
                #if(plotAI[i,j]>-20):
                    # Only plot if XTrack flag is met
                    if((XTRACK[i,j] == 0) | (XTRACK[i,j] == 4)):
                        # Print values to text file
                        if(LAT[i,j] > latmin):
                            counter+=1
                        #if((plotXTrack[i,j] == 0) | ((plotXTrack[i,j] & 4 == 4))):
                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                            index2 = int(np.floor(LON[i,j]*4 + 720.))
                            #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
                            #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
                        
                            if(index1 < 0): index1 = 0
                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                            if(index2 < 0): index2 = 0                                                                                            
                            if(index2 > 1439): index2 = 1439
                   
                            #if(diff<0.2): 
                            #    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
                            g_NRAD_354[index2, index1] = (g_NRAD_354[index2,index1]*count_NRAD_354[index2,index1] + NRAD[i,j,0])/(count_NRAD_354[index2,index1]+1)
                            g_NRAD_388[index2, index1] = (g_NRAD_388[index2,index1]*count_NRAD_388[index2,index1] + NRAD[i,j,1])/(count_NRAD_388[index2,index1]+1)
                            g_NRAD_500[index2, index1] = (g_NRAD_500[index2,index1]*count_NRAD_500[index2,index1] + NRAD[i,j,2])/(count_NRAD_500[index2,index1]+1)
                            g_REFL_354[index2, index1] = (g_REFL_354[index2,index1]*count_REFL_354[index2,index1] + REFL[i,j,0])/(count_REFL_354[index2,index1]+1)
                            g_REFL_388[index2, index1] = (g_REFL_388[index2,index1]*count_REFL_388[index2,index1] + REFL[i,j,1])/(count_REFL_388[index2,index1]+1)
                            g_SALB_354[index2, index1] = (g_SALB_354[index2,index1]*count_SALB_354[index2,index1] + SALB[i,j,0])/(count_SALB_354[index2,index1]+1)
                            g_UVAI_354[index2, index1] = (g_UVAI_354[index2,index1]*count_UVAI_354[index2,index1] + UVAI[i,j])/(count_UVAI_354[index2,index1]+1)
                            count_NRAD_354[index2,index1] = count_NRAD_354[index2,index1] + 1
                            count_NRAD_388[index2,index1] = count_NRAD_388[index2,index1] + 1
                            count_NRAD_500[index2,index1] = count_NRAD_500[index2,index1] + 1
                            count_REFL_354[index2,index1] = count_REFL_354[index2,index1] + 1
                            count_REFL_388[index2,index1] = count_REFL_388[index2,index1] + 1
                            count_SALB_354[index2,index1] = count_SALB_354[index2,index1] + 1
                            count_UVAI_354[index2,index1] = count_UVAI_354[index2,index1] + 1

   
    # Apply algae screening to 500 nm normalized radiance data 
    print("Applying algae screening to 500 nm normalized radiance data")
    mask_rad500 = np.ma.masked_where(((count_NRAD_500 == 0)), g_NRAD_500)
    mask_rad500 = np.ma.masked_where(((g_SALB_354 > 0.09)), mask_rad500)
    mask_rad500 = np.ma.masked_where(((g_REFL_354 > 0.18)), mask_rad500)
    mask_rad500 = np.ma.masked_where(((g_REFL_354 < 0.09)), mask_rad500)
    mask_NRAD500 = np.ma.masked_where(((g_UVAI_354 < 0.6)), mask_rad500)

    # Apply algae screening to UVAI data 
    print("Applying algae screening to UVAI data")
    mask_uvai = np.ma.masked_where(((count_UVAI_354 == 0)), g_UVAI_354)
    mask_uvai = np.ma.masked_where(((g_SALB_354 > 0.09)), mask_uvai)
    mask_uvai = np.ma.masked_where(((g_REFL_354 > 0.18)), mask_uvai)
    mask_uvai = np.ma.masked_where(((g_REFL_354 < 0.09)), mask_uvai)
    mask_UVAI354 = np.ma.masked_where(((g_UVAI_354 < 0.6)), mask_uvai)

    OMI_single_dict = {}
    OMI_single_dict['lat'] = lat_ranges
    OMI_single_dict['lon'] = lon_ranges
    OMI_single_dict['AI'] = g_UVAI_354
    OMI_single_dict['AI_count'] = count_UVAI_354
    OMI_single_dict['NRAD500'] = g_NRAD_500
    OMI_single_dict['NRAD500_count'] = count_NRAD_500
    OMI_single_dict['AI_algae'] = mask_UVAI354
    OMI_single_dict['NRAD500_algae'] = mask_NRAD500

    return OMI_single_dict

def writeOMI_toNCDF(OMI_data,minlat=30):
    lat_ranges = np.arange(minlat,90.5,1.0)
    lon_ranges = np.arange(-180,180,1.0)
    
    nc = Dataset('/home/bsorenson/Research/OMI/omi_ai_V003_2005_2020.nc','w',format='NETCDF4')
  
    # Dimensions = lat, lon, time
    testdict = OMI_data['50x5']
    testkeys = list(testdict.keys())
    num_lat = len(lat_ranges)
    num_lon = len(lon_ranges)
    num_time = len(testdict.keys())
    times = np.arange(num_time)
    
    n_time = nc.createDimension('nmth',num_time)
    n_lat  = nc.createDimension('dlat',num_lat)
    n_lon  = nc.createDimension('dlon',num_lon)

    MONTH=nc.createVariable('MONTH','i2',('nmth'))
    MONTH.description='Months since January 2005'
    LAT=nc.createVariable('Latitude','i2',('dlat','dlon'))
    LAT.description='Latitude'
    LAT.units='Degrees'
    LON=nc.createVariable('Longitude','i2',('dlat','dlon'))
    LON.description='Longitude'
    LON.units='Degrees'
    AI = nc.createVariable('AI','f4',('nmth','dlat','dlon'))
    AI.description='Monthly Averaged Aerosol Index'
    OB_COUNT=nc.createVariable('OB_COUNT','i2',('nmth','dlat','dlon'))
    OB_COUNT.description='# of OMI AI measurements used in each monthly average'

    # Fill in dimension variables
    for i in range(num_time):
        MONTH[i] = times[i]
    for i in range(num_lat):
        for j in range(num_lon):
            LAT[i,j]=lat_ranges[i]
            LON[i,j]=lon_ranges[j]

    # Fill in actual variables
    for i in range(num_lat):
        print(lat_ranges[i])
        for j in range(num_lon):
            dictkey = (str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j])))
            if(dictkey not in OMI_data):
                # Insert missing values for AI and count
                AI[:,i,j] = [-999.9 for m in range(num_time)]
                OB_COUNT[:,i,j] = [-99 for m in range(num_time)]
            else:
                for m in range(num_time):
                    timekey = testkeys[m]
                    if(timekey not in OMI_data[dictkey]):
                        AI[m,i,j] = -999.9
                        OB_COUNT[m,i,j] = -99
                    else:
                        AI[m,i,j] = OMI_data[dictkey][timekey]['avg']
                        OB_COUNT[m,i,j] = OMI_data[dictkey][timekey]['#_obs']
    nc.close()
   
def readOMI_NCDF(infile='/home/bsorenson/Research/OMI/omi_ai_V003_2005_2020.nc',minlat=50):
    # Read in data to netCDF object
    in_data = Dataset(infile,'r')

    # Set up dictionary to hold data
    OMI_data = {}
    
    OMI_data['AI'] = in_data['AI'][:,:,:]
    OMI_data['LAT'] = in_data['Latitude'][:,:]
    OMI_data['LON'] = in_data['Longitude'][:,:]
    OMI_data['MONTH'] = in_data['MONTH'][:]

    # Set up date strings in the file
    start_date = datetime(year = 2005, month = 1, day = 1)
    OMI_data['DATES'] = \
        [(start_date + relativedelta(months=mi)).strftime('%Y%m') for mi in \
        OMI_data['MONTH']]

    # to add months to datetime object, do
    ###from dateutil.relativedelta import relativedelta
    ###datetime.datetime(year=2004,month=10,day=1) + relativedelta(months=1)
       
    return OMI_data

# This function assumes the data is being read from the netCDF file
# NOTE: Assume user is using new OMI climo file which starts in January
def calcOMI_MonthClimo(OMI_data):

    # Set up arrays to hold monthly climatologies
    month_climo = np.zeros((12,OMI_data['AI'].shape[1],OMI_data['AI'].shape[2]))

    # Mask the monthly averages
    local_data = np.copy(OMI_data['AI'][:,:,:])
    local_mask = np.ma.masked_where(local_data == -999.9, local_data)
 
    # Calculate monthly climatologies
    for m_i in range(12):
        month_climo[m_i,:,:] = np.nanmean(local_mask[m_i::12,:,:],axis=0)
        print("Month: ",OMI_data['DATES'][m_i][4:]) 
    #month_climo[0,:,:]  = np.nanmean(local_mask[0::12,:,:],axis=0)  # January 
    #month_climo[1,:,:]  = np.nanmean(local_mask[1::12,:,:],axis=0)  # February
    #month_climo[2,:,:]  = np.nanmean(local_mask[2::12,:,:],axis=0)  # March 
    #month_climo[3,:,:]  = np.nanmean(local_mask[3::12,:,:],axis=0)  # April 
    #month_climo[4,:,:]  = np.nanmean(local_mask[4::12,:,:],axis=0)  # May 
    #month_climo[5,:,:]  = np.nanmean(local_mask[5::12,:,:],axis=0)  # June 
    #month_climo[6,:,:]  = np.nanmean(local_mask[6::12,:,:],axis=0)  # July 
    #month_climo[7,:,:]  = np.nanmean(local_mask[7::12,:,:],axis=0)  # August 
    #month_climo[8,:,:]  = np.nanmean(local_mask[8::12,:,:],axis=0)  # September 
    #month_climo[9,:,:]  = np.nanmean(local_mask[9::12,:,:],axis=0)  # October 
    #month_climo[10,:,:] = np.nanmean(local_mask[10::12,:,:],axis=0) # November
    #month_climo[11,:,:] = np.nanmean(local_mask[11::12,:,:],axis=0) # December

    # Insert data into dictionary
    OMI_data['MONTH_CLIMO'] = month_climo

    return OMI_data
 

def plotOMI_MK(OMI_data,start_date,end_date,save=False,file_type='XR123',season='',minlat=30.):
    if(file_type=='NXAR'):
        title_flabel = 'No XTrack and All Rows'
        outname_flabel = 'noX_allRows'
    elif(file_type=='XAR'):
        title_flabel = 'XTrack and All Rows'
        outname_flabel = 'xtrack_allRows'
    elif(file_type=='XR123'):
        title_flabel = 'XTrack and Rows 1-23'
        outname_flabel = 'xtrack_rows1to23'
    lat_ranges = np.arange(minlat,90,1.0)
    #lat_ranges = np.arange(lowest_lat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False
    if(season=='spring'):
        spring = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    OMI_data[lkey].pop(tkey)

    # Find the lowest lat in the file
    #lowest_lat = float(sorted(OMI_data.keys())[0].split('x')[0])
     
    # Set up the polar stereographic map
    fig1 = plt.figure()
    ax = plt.subplot(111)
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
   
    # Define the colors for each range of P-Value
    # Ranges: >=0.4,0.4-0.3,0.3-0.2,0.2-0.15,0.15-0.1,0.1-0.05,<0.05
    # Seven colors
    # - White,yellow,orange-yellow,orange,red-orange,red,dark red
    color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
    pRED   = np.array([255.,255.,255.,255.,255.,150.,150.])
    pGREEN = np.array([255.,255.,150., 75.,  0.,  0.,  0.])
    pBLUE  = np.array([255.,  0.,  0.,  0.,  0.,  0.,150.])
    pRED   = pRED/256.
    pGREEN = pGREEN/256.
    pBLUE  = pBLUE/256.

    noRED = 0.
    noGREEN = 0.
    noBLUE = 0.
    
    # - White,blue-green1,orange-yellow,orange,red-orange,red,dark red
    nRED   = np.array([255.,100.,  0.,  0.,  0.,  0.,  0.])
    nGREEN = np.array([255.,255.,255.,255.,200.,100.,  0.])
    nBLUE  = np.array([255.,  0.,  0.,255.,255.,255.,150.])
    nRED   = nRED/256.
    nGREEN = nGREEN/256.
    nBLUE  = nBLUE/256.
 
    ##cmap = plt.get_cmap('bwr')
    ### Center the colorbar on zero
    ##norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    ###norm = Normalize(vmin=vmin,vmax=vmax)
    ##mapper = ScalarMappable(norm=norm,cmap=cmap)

    # Set up initial values for the analysis regions
    canada_avg = 0.
    canada_counts = 0
    siberia_avg = 0.
    siberia_counts = 0 
    eus_avg = 0.
    eus_counts = 0
    wus_avg = 0.
    wus_counts = 0
    ea_avg = 0.
    ea_counts = 0
    wa_avg = 0.
    wa_counts = 0
    europe_avg = 0.
    europe_counts = 0
    
    # Loop over all the keys and print(the regression slopes 
    # Grab the averages for the key
    max_pval = -10.
    min_pval = 10.
    #for i in range(15,20):
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = (str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j])))
            #print(dictkey)
            if(dictkey=='48x-97'):
                #print("sorted(OMI_data[dictkey].keys())=",sorted(OMI_data[dictkey].keys()))
                min_date = sorted(OMI_data[dictkey].keys())[0]
                max_date = sorted(OMI_data[dictkey].keys())[-1]
    
            # If no data are present for the curent lat/lon box, fill it with
            # black
            if(dictkey not in OMI_data.keys()):
                colorRED = 0.
                colorGREEN = 0.
                colorBLUE = 0.
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(OMI_data[dictkey].keys())==0):
                print(dictkey,pval,color_index_value,'NO DATA')
                #print('Here 2')
                colorRED = 0.
                colorGREEN = 0.
                colorBLUE = 0.
            else:
                avgs = [OMI_data[dictkey][date]['avg'] for date in \
                    sorted(OMI_data[dictkey].keys())]
                if(len(avgs)<2):
                    #print('Here 3')
                    colorRED = 0.
                    colorGREEN = 0.
                    colorBLUE = 0.
                else:
                    # Check the current max and min
                    x_vals = np.arange(0,len(OMI_data[dictkey].keys()))
                    avgs = np.ma.masked_array([OMI_data[dictkey][date]['avg'] for date in sorted(OMI_data[dictkey].keys())])
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(OMI_data[dictkey].keys()))
                    dates = temp_dates
                    #dates = temp_dates[np.where(avgs.mask==False)[0]]
                    x_vals = np.arange(0,len(dates))
    
                    nx=len(avgs)
                    num_d=int(nx*(nx-1)/2)  # The number of elements in d
                    Sn=np.zeros(num_d)
                    S=0
                    for ti in range(nx-1):
                        for tj in range(ti+1,nx):
                            S+=np.sign(avgs[tj]-avgs[ti])
                        # Endfor
                    # Endfor

                    # Find the unique values in the data
                    uniq = np.unique(avgs)
                    g = len(uniq)
                    if(nx==g):
                        # No ties
                        Vs = (nx*(nx-1.)*(2.*nx-5.))/18.
                    else:
                        tp = np.zeros(uniq.shape)
                        for si in range(g):
                            tp[si] = sum(avgs==uniq[si])
                        Vs = (nx*(nx-1.)*(2.*nx-5.)-np.sum(tp*(tp-1.)*(2.*tp-5.)))/18.
                    if (S > 0.): 
                        z=(S-1.)/np.sqrt(Vs)
                    elif (S < 0.): 
                        z=(S+1.)/np.sqrt(Vs)
                    else: 
                        z=0.
                    # Calculate the p value of the trend
                    pval=2*(1.-stats.norm.cdf(abs(z)))  # (two-side)
                    alpha=0.05
                    ##h = abs(z) > stats.norm.ppf(1.-alpha/2.)
                    ### Determine the trend type 
                    ##if(z>0) and h:
                    ##    trend='increasing'
                    ##elif(z<0) and h:
                    ##    trend='decreasing'
                    ##else:
                    ##    trend='no trend'

                    if(pval<min_pval):
                        min_pval=pval
                    elif(pval>max_pval):
                        max_pval=pval

                    #color = mapper.to_rgba(slope)
    #color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
                    color_index_value = 0
                    color_index = np.where(pval<color_range)[0]
                    if(len(color_index)>0):
                        color_index_value = color_index[-1]+1
                    if(S>0):
                        colorRED = pRED[color_index_value]
                        colorGREEN = pGREEN[color_index_value]
                        colorBLUE = pBLUE[color_index_value]
                    else:
                        colorRED = nRED[color_index_value]
                        colorGREEN = nGREEN[color_index_value]
                        colorBLUE = nBLUE[color_index_value]
    
                    if(np.isnan(pval) == False):
                        # Add value to analysis regions (if possible)
                        if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                            canada_avg+=pval
                            canada_counts+=1
                        elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                            siberia_avg+=pval
                            siberia_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                            eus_avg+=pval
                            eus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                            wus_avg+=pval
                            wus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                            ea_avg+=pval
                            ea_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                            wa_avg+=pval
                            wa_counts+=1
                        elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                            europe_avg+=pval
                            europe_counts+=1
    
                        print(dictkey,pval,color_index_value)
            # Find the y coordinates of the current LatxLon box
            y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
            # Find the x coordinates of the current LatxLon box
            x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
            # Convert x and y into map coordinates
            mx, my = m(x,y)
            mxy = zip(mx,my)
            pair1 = (mx[0],my[0])
            pair2 = (mx[1],my[1])
            pair3 = (mx[2],my[2])
            pair4 = (mx[3],my[3])
            #mxy = zip(mx,my)
            # Plot the box on the map using color
            colors = [colorRED,colorGREEN,colorBLUE]
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors)
            #poly = Polygon(mxy,facecolor=color2,edgecolor=color2)
            plt.gca().add_patch(poly)
    
    y = [0.,0.,0.,0.]
    # Find the x coordinates of the current LatxLon box
    x = [0.,0.,0.,0.]
    color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
    # Convert x and y into map coordinates
    mx, my = m(x,y)
    pair1 = (mx[0],my[0])
    pair2 = (mx[1],my[1])
    pair3 = (mx[2],my[2])
    pair4 = (mx[3],my[3])
    # Positive
    # Add legend for 0.3-0.2
    colors = [pRED[1],pGREEN[1],pBLUE[1]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.4-0.3')
    plt.gca().add_patch(poly)
    # Add legend for 0.3-0.2
    colors = [nRED[1],nGREEN[1],nBLUE[1]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.4-0.3')
    plt.gca().add_patch(poly)
    # Add legend for 0.2-0.15
    colors = [pRED[2],pGREEN[2],pBLUE[2]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.3-0.2')
    plt.gca().add_patch(poly)
    # Add legend for 0.15-0.1
    colors = [nRED[2],nGREEN[2],nBLUE[2]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.3-0.2')
    plt.gca().add_patch(poly)
    # Add legend for 0.1-0.05
    colors = [pRED[3],pGREEN[3],pBLUE[3]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.2-0.15')
    plt.gca().add_patch(poly)
    # Add legend for <0.05
    colors = [nRED[3],nGREEN[3],nBLUE[3]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.2-0.15')
    plt.gca().add_patch(poly)
    colors = [pRED[4],pGREEN[4],pBLUE[4]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.15-0.1')
    plt.gca().add_patch(poly)
    # Negative
    # Add legend for 0.2-0.15
    colors = [nRED[4],nGREEN[4],nBLUE[4]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.15-0.1')
    plt.gca().add_patch(poly)
    # Add legend for 0.15-0.1
    colors = [pRED[5],pGREEN[5],pBLUE[5]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ 0.1-0.05')
    plt.gca().add_patch(poly)
    # Add legend for 0.1-0.05
    colors = [nRED[5],nGREEN[5],nBLUE[5]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- 0.1-0.05')
    plt.gca().add_patch(poly)
    # Add legend for <0.05
    colors = [pRED[6],pGREEN[6],pBLUE[6]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='+ <0.05')
    plt.gca().add_patch(poly)
    colors = [nRED[6],nGREEN[6],nBLUE[6]]
    poly = Polygon([pair1,pair2,pair3,pair4],facecolor=colors,edgecolor=colors,label='- <0.05')
    plt.gca().add_patch(poly)

    box = ax.get_position()
    ax.set_position([box.x0,box.y0+box.height*0.1,box.width,box.height*0.9])
    ax.legend(loc='lower center',prop={'size':5},bbox_to_anchor=(0.48,-0.1),ncol=6)

    minmax_diff = (max_pval-min_pval) 
    minmax_range = np.arange(min_pval,max_pval,minmax_diff/6.0)
    print("min pval = ",min_pval)
    print("max pval = ",max_pval)
    print("Range = ",minmax_range)

    if(canada_counts>0):
        canada_avg = canada_avg/canada_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
    if(siberia_counts>0):
        siberia_avg = siberia_avg/siberia_counts
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
    if(eus_counts>0):
        eus_avg = eus_avg/eus_counts
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
    if(wus_counts>0):
        wus_avg = wus_avg/wus_counts
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
    if(ea_counts>0):
        ea_avg = ea_avg/ea_counts
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
    if(wa_counts>0):
        wa_avg = wa_avg/wa_counts
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
    if(europe_counts>0):
        europe_avg = europe_avg/europe_counts
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

    #start_date = min_date.decode("utf-8")
    #end_date = max_date.decode("utf-8")
    title_string = 'OMI Average Ultraviolet Aerosol Index Trend P-Values\n'+\
                   start_date+' to '+end_date
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    elif(summer is True):
        title_string = title_string+'\nJune, July, August'
    elif(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    elif(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    # Set up the colorbar
    #cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    #cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    #color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
    
    plt.xticks(rotation=45,fontsize=6)
    #plt.legend(prop={'size':10})
    fig1.canvas.mpl_connect('button_press_event',onclick)
    if(save is True):
        if(spring is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_spring_newData.png'
        elif(summer is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_summer_newData.png'
        elif(autumn is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_autumn_newData.png'
        elif(winter is True):
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_winter_newData.png'
        else:
            filename = 'omi_ai_pvalues_'+start_date+end_date+'_'+str(int(lowest_lat))+'to90_whole_year_newData.png'
        
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

def plotOMI(OMI_data,start_date,end_date,save=False,trend_type='standard',file_type='XR123',season='',minlat=30.):
    if(file_type=='NXAR'):
        title_flabel = 'No XTrack and All Rows'
        outname_flabel = 'noX_allRows'
    elif(file_type=='XAR'):
        title_flabel = 'XTrack and All Rows'
        outname_flabel = 'xtrack_allRows'
    elif(file_type=='XR123'):
        title_flabel = 'XTrack and Rows 1-23'
        outname_flabel = 'xtrack_rows1to23'
    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='_thielSen'
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False

    # If only summer months are being analyzed, remove all data except 
    # in summer
    if(season=='spring'):
        spring = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    OMI_data[lkey].pop(tkey)


    # Find the lowest lat in the file
    #lowest_lat = 50.

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # Set up the polar stereographic map
    fig1 = plt.figure()
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,         \
                resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    

    cmap = plt.get_cmap('bwr')
    #bcmap = plt.cm.set_cmap(cmap)
    if(summer is True):
        v_min = -0.350 # Summer values
        mid_val = 0
        v_max = 0.900
    else:
        v_min = -0.350  # Whole-year values
        mid_val = 0
        v_max = 0.900

    if(minlat>30):
        v_max = 0.6
        mid_val = 0
        v_min = -0.6

    # Center the colorbar on zero
    norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    #norm = Normalize(vmin=vmin,vmax=vmax)
    mapper = ScalarMappable(norm=norm,cmap=cmap)

    # Set up initial values for the analysis regions
    canada_avg = 0.
    canada_counts = 0
    siberia_avg = 0.
    siberia_counts = 0 
    eus_avg = 0.
    eus_counts = 0
    wus_avg = 0.
    wus_counts = 0
    ea_avg = 0.
    ea_counts = 0
    wa_avg = 0.
    wa_counts = 0
    europe_avg = 0.
    europe_counts = 0


    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    max_slope = -10.
    min_slope = 10.
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = (str(int(lat_ranges[i]))+'x'+                            \
                       str(int(lon_ranges[j])))

            # If no data are present for the curent lat/lon box, fill it with
            # black
            keylist = [ky for ky in OMI_data.keys()]
            if(dictkey not in OMI_data.keys()):
                color=(0,0,0,0)
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(OMI_data[dictkey].keys())==0):
                color=(0,0,0,0)
            else:
                avgs = [OMI_data[dictkey][date]['avg'] for date in \
                    sorted(OMI_data[dictkey].keys())]
                # Check the current max and min
                x_vals = np.arange(0,len(OMI_data[dictkey].keys()))
                temp_dates = np.array(sorted(OMI_data[dictkey].keys()))
                # Find the slope of the line of best fit for the time series of
                # average data
                if(trend_type=='standard'): 
                    slope, intercept, r_value, p_value, std_err = \
                        stats.linregress(x_vals,avgs)
                    slope *= len(x_vals)
                else:
                    #The slope
                    S=0
                    sm=0
                    nx = len(avgs)
                    num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
                    Sn=np.zeros(num_d)
                    for si in range(0,nx-1):
                        for sj in range(si+1,nx):
                            # Find the slope between the two points
                            Sn[sm] = (avgs[si]-avgs[sj])/(si-sj) 
                            sm=sm+1
                        # Endfor
                    # Endfor
                    Snsorted=sorted(Sn)
                    sm=int(num_d/2.)
                    print(dictkey,len(Snsorted))
                    if(len(Snsorted)==1):
                        color=(0,0,0,0) 
                    else:
                        if(2*sm    == num_d):
                            slope=0.5*(Snsorted[sm]+Snsorted[sm+1])
                        if((2*sm)+1 == num_d): 
                            slope=Snsorted[sm+1]
                        slope = slope*len(avgs)

                if(slope>max_slope):
                    max_slope=slope
                elif(slope<min_slope):
                    min_slope=slope
        
                color = mapper.to_rgba(slope)

                # Add value to analysis regions (if possible)
                if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                    canada_avg+=slope
                    canada_counts+=1
                elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                    siberia_avg+=slope
                    siberia_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                    eus_avg+=slope
                    eus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                    wus_avg+=slope
                    wus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                    ea_avg+=slope
                    ea_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                    wa_avg+=slope
                    wa_counts+=1
                elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                    europe_avg+=slope
                    europe_counts+=1
                


            # Find the y coordinates of the current LatxLon box
            y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
            # Find the x coordinates of the current LatxLon box
            x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
            # Convert x and y into map coordinates
            mx, my = m(x,y)
            mxy = zip(mx,my)
            pair1 = (mx[0],my[0])
            pair2 = (mx[1],my[1])
            pair3 = (mx[2],my[2])
            pair4 = (mx[3],my[3])
            color2 = rgb2hex(color)
            #color2 = rgb2hex(color[:2]+color[3:])
            # Plot the box on the map using color
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,         \
                edgecolor=color2)
            plt.gca().add_patch(poly)
    
    #start_date = str(temp_dates[0].decode("utf-8"))
    #end_date = str(temp_dates[-1].decode("utf-8"))
    minmax_diff = (max_slope-min_slope) 
    minmax_range = np.arange(min_slope,max_slope,minmax_diff/8.0)
    print("min slope = ",min_slope)
    print("max slope = ",max_slope)
    print("Range = ",minmax_range)

    if(canada_counts>0):
        canada_avg = canada_avg/canada_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
    if(siberia_counts>0):
        siberia_avg = siberia_avg/siberia_counts
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
    if(eus_counts>0):
        eus_avg = eus_avg/eus_counts
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
    if(wus_counts>0):
        wus_avg = wus_avg/wus_counts
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
    if(ea_counts>0):
        ea_avg = ea_avg/ea_counts
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
    if(wa_counts>0):
        wa_avg = wa_avg/wa_counts
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
    if(europe_counts>0):
        europe_avg = europe_avg/europe_counts
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

    title_string = 'Change in OMI Ultraviolet Average Aerosol Index\n'+        \
        title_flabel+'\n'+str(start_date)+' to '+str(end_date)
        ##start_date+' to '+end_date+'\nXTrack and all rows'
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    if(summer is True):
        title_string = title_string+'\nJune, July, August'
    if(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    if(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    # Set up the colorbar
    cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    plt.xticks(rotation=45,fontsize=6)
    fig1.canvas.mpl_connect('button_press_event',onclick)
    if(save is True):
        if(spring is True):
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_spring_'+outname_flabel+'_newData.png'
        elif(summer is True):     
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_summer_'+outname_flabel+'_newData.png'
        elif(autumn is True):     
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_autumn_'+outname_flabel+'_newData.png'
        elif(winter is True):    
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_winter_'+outname_flabel+'_newData.png'
        else:                   
            filename = 'omi_ai'+trend_label+'_trends_'+str(start_date)+str(end_date)+'_'+str(int(minlat))+'to90_whole_year_'+outname_flabel+'_newData.png'
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

def plotOMI_Climo(OMI_data,start_date,end_date,save=False,trend_type='standard',file_type='XR123',season='',minlat=30.):
    if(file_type=='NXAR'):
        title_flabel = 'No XTrack and All Rows'
        outname_flabel = 'noX_allRows'
    elif(file_type=='XAR'):
        title_flabel = 'XTrack and All Rows'
        outname_flabel = 'xtrack_allRows'
    elif(file_type=='XR123'):
        title_flabel = 'XTrack and Rows 1-23'
        outname_flabel = 'xtrack_rows1to23'
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False

    # If only summer months are being analyzed, remove all data except 
    # in summer
    if(season=='spring'):
        spring = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    OMI_data[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(OMI_data.keys()):
            for tkey in sorted(OMI_data[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    OMI_data[lkey].pop(tkey)


    # Find the lowest lat in the file

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # Set up the polar stereographic map
    fig1 = plt.figure()
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,         \
                resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    

    cmap = plt.get_cmap('jet')
    #bcmap = plt.cm.set_cmap(cmap)
    if(summer is True):
        v_min = -0.350 # Summer values
        mid_val = 0
        v_max = 0.900
    else:
        v_min = -0.350  # Whole-year values
        mid_val = 0
        v_max = 0.900
    v_min = -1.000  # Whole-year values
    mid_val = 0
    v_max = 1.500
    if(minlat>30.):
        v_min = 0.0
        mid_val = 0
        v_max = 1.00

    # Center the colorbar on zero
    #norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    norm = Normalize(vmin=v_min,vmax=v_max)
    mapper = ScalarMappable(norm=norm,cmap=cmap)

    # Set up initial values for the analysis regions
    canada_avg = 0.
    canada_counts = 0
    siberia_avg = 0.
    siberia_counts = 0 
    eus_avg = 0.
    eus_counts = 0
    wus_avg = 0.
    wus_counts = 0
    ea_avg = 0.
    ea_counts = 0
    wa_avg = 0.
    wa_counts = 0
    europe_avg = 0.
    europe_counts = 0

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    max_avg_uvai = -10.
    min_avg_uvai = 10.
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = (str(int(lat_ranges[i]))+'x'+                            \
                       str(int(lon_ranges[j])))

            # If no data are present for the curent lat/lon box, fill it with
            # black
            keylist = [ky for ky in OMI_data.keys()]
            if(dictkey not in OMI_data.keys()):
                color=(0,0,0,0)
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(OMI_data[dictkey].keys())==0):
                color=(0,0,0,0)
            else:
                avgs = np.array([OMI_data[dictkey][date]['avg'] for date in \
                    sorted(OMI_data[dictkey].keys())])
                counts = np.array([OMI_data[dictkey][date]['#_obs'] for date in \
                    sorted(OMI_data[dictkey].keys())])
                avg_uvai = sum(avgs*counts)/sum(counts)
                #avg_uvai = np.average(avgs)
                temp_dates = np.array(sorted(OMI_data[dictkey].keys()))
                # Check the current max and min
                #x_vals = np.arange(0,len(OMI_data[dictkey].keys()))
                # Find the slope of the line of best fit for the time series of
                # average data
                #slope, intercept, r_value, p_value, std_err = \
                #    stats.linregress(x_vals,avgs)
                #slope *= len(x_vals)

                if(avg_uvai>max_avg_uvai):
                    max_avg_uvai=avg_uvai
                elif(avg_uvai<min_avg_uvai):
                    min_avg_uvai=avg_uvai
        
                color = mapper.to_rgba(avg_uvai)

                # Add value to analysis regions (if possible)
                if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                    canada_avg+=avg_uvai
                    canada_counts+=1
                elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                    siberia_avg+=avg_uvai
                    siberia_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                    eus_avg+=avg_uvai
                    eus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                    wus_avg+=avg_uvai
                    wus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                    ea_avg+=avg_uvai
                    ea_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                    wa_avg+=avg_uvai
                    wa_counts+=1
                elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                    europe_avg+=avg_uvai
                    europe_counts+=1

            # Find the y coordinates of the current LatxLon box
            y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
            # Find the x coordinates of the current LatxLon box
            x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
            # Convert x and y into map coordinates
            mx, my = m(x,y)
            mxy = zip(mx,my)
            pair1 = (mx[0],my[0])
            pair2 = (mx[1],my[1])
            pair3 = (mx[2],my[2])
            pair4 = (mx[3],my[3])
            color2 = rgb2hex(color)
            #color2 = rgb2hex(color[:2]+color[3:])
            # Plot the box on the map using color
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,         \
                edgecolor=color2)
            plt.gca().add_patch(poly)
    
    #start_date = str(temp_dates[0].decode("utf-8"))
    #end_date = str(temp_dates[-1].decode("utf-8"))
    minmax_diff = (max_avg_uvai-min_avg_uvai) 
    minmax_range = np.arange(min_avg_uvai,max_avg_uvai,minmax_diff/8.0)
    print("min avg_uvai = ",min_avg_uvai)
    print("max avg_uvai = ",max_avg_uvai)
    print("Range = ",minmax_range)

    if(canada_counts>0):
        canada_avg = canada_avg/canada_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
    if(siberia_counts>0):
        siberia_avg = siberia_avg/siberia_counts
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
    if(eus_counts>0):
        eus_avg = eus_avg/eus_counts
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
    if(wus_counts>0):
        wus_avg = wus_avg/wus_counts
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
    if(ea_counts>0):
        ea_avg = ea_avg/ea_counts
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
    if(wa_counts>0):
        wa_avg = wa_avg/wa_counts
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
    if(europe_counts>0):
        europe_avg = europe_avg/europe_counts
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

    start_date = str(start_date)
    end_date = str(end_date)
    title_string = 'OMI Ultraviolet Average Aerosol Index Climatology\n'+      \
        title_flabel+'\n'+start_date+' to '+end_date
        ##start_date+' to '+end_date+'\nXTrack and All Rows'
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    if(summer is True):
        title_string = title_string+'\nJune, July, August'
    if(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    if(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    # Set up the colorbar
    cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    plt.xticks(rotation=45,fontsize=6)
    plt.xlabel('Ultraviolet Aerosol Index',fontsize=6)
    fig1.canvas.mpl_connect('button_press_event',onclick_climo)
    if(save is True):
        if(spring is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_spring_'+outname_flabel+'_newData.png'
        elif(summer is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_summer_'+outname_flabel+'_newData.png'
        elif(autumn is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_autumn_'+outname_flabel+'_newData.png'
        elif(winter is True):
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_winter_'+outname_flabel+'_newData.png'
        else:
            filename = 'omi_ai_climo_'+start_date+end_date+'_'+str(int(minlat))+'to90_whole_year_'+outname_flabel+'_newData.png'
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

def plotOMI_NCDF_Climo(OMI_data,start_idx=0,end_idx=169,season = '',minlat=60,\
                       save=False):

    season_dict = {
        'spring': '\nMAM',\
        'summer': '\nJJA',\
        'autumn': '\nSON',\
        'winter': '\nDJF',\
        '': ''
    }
    
    # Make copy of OMI_data array
    local_data = np.copy(OMI_data['AI'])

    start_date = datetime(year=2004,month=10,day=1)

    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False

    month_objects = []
    keepers = np.arange(1,13)
    if(season=='spring'):
        spring = True 
        keepers = [3,4,5]
    elif(season=='summer'):
        summer = True 
        keepers = [6,7,8]
    elif(season=='autumn'):
        autumn = True 
        keepers = [9,10,11]
    elif(season=='winter'):
        winter = True 
        keepers = [12,1,2]
   
    for m_idx in OMI_data['MONTH']:
        new_date = start_date + relativedelta(months=m_idx)
        if(new_date.month not in keepers):  
            local_data[m_idx,:,:] = -999.9
        else:
            month_objects.append(new_date)

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Mask any missing values
    mask_AI = np.ma.masked_where(local_data == -999.9, local_data)

    # Calculate climatology between desired indices
    OMI_climo = np.nanmean(mask_AI[start_idx:end_idx,:,:],axis=0)

    # Make figure title
    first_date = month_objects[0].strftime("%Y%m")
    last_date = month_objects[-1].strftime("%Y%m")
    print(month_objects[0].strftime("%Y%m"),month_objects[-1].strftime("%Y%m"))
    title = 'OMI AI Climatology\n'+first_date + ' - ' + last_date + season_dict[season]

    # Make figure
    plt.close()
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(OMI_data['LON'], OMI_data['LAT'],OMI_climo,transform = datacrs,cmap = colormap,\
            vmin = -1.0, vmax = 1.5)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='Aerosol Index')
    ax.set_title(title)

    if(save == True):
        season_adder = ''
        if(len(season.strip()) != 0):
            season_adder = '_' + season.strip()
        outname = 'omi_ai_climo_' + first_date + '_' + last_date + season_adder + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# The same as the above function, but plots all 4 seasons in a single 4-panel
# plot.
def plotOMI_NCDF_Climo_FourPanel(OMI_data,start_idx=0,end_idx=169,minlat=60,\
                       save=False):

    season_dict = {
        'spring': '\nMAM',\
        'summer': '\nJJA',\
        'autumn': '\nSON',\
        'winter': '\nDJF',\
        '': ''
    }
    
    # Make copy of OMI_data array
    local_spring = np.copy(OMI_data['AI'])
    local_summer = np.copy(OMI_data['AI'])
    local_autumn = np.copy(OMI_data['AI'])
    local_winter = np.copy(OMI_data['AI'])

    start_date = datetime(year=2004,month=10,day=1)

    month_objects = []
    spring_keepers = [3,4,5]
    summer_keepers = [6,7,8]
    autumn_keepers = [9,10,11]
    winter_keepers = [12,1,2]
   
    for m_idx in OMI_data['MONTH']:
        new_date = start_date + relativedelta(months=m_idx)
        if(new_date.month in spring_keepers):  
            local_spring[m_idx,:,:] = -999.9
        if(new_date.month in spring_keepers):  
            local_spring[m_idx,:,:] = -999.9
        if(new_date.month in spring_keepers):  
            local_spring[m_idx,:,:] = -999.9
        if(new_date.month in spring_keepers):  
            local_spring[m_idx,:,:] = -999.9

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Mask any missing values
    mask_AI = np.ma.masked_where(local_data == -999.9, local_data)

    # Calculate climatology between desired indices
    OMI_climo = np.nanmean(mask_AI[start_idx:end_idx,:,:],axis=0)

    # Make figure title
    first_date = month_objects[0].strftime("%Y%m")
    last_date = month_objects[-1].strftime("%Y%m")
    print(month_objects[0].strftime("%Y%m"),month_objects[-1].strftime("%Y%m"))
    title = 'OMI AI Climatology\n'+first_date + ' - ' + last_date + season_dict[season]

    # Make figure
    plt.close()
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(OMI_data['LON'], OMI_data['LAT'],OMI_climo,transform = datacrs,cmap = colormap,\
            vmin = -1.0, vmax = 1.5)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='Aerosol Index')
    ax.set_title(title)

    if(save == True):
        season_adder = ''
        if(len(season.strip()) != 0):
            season_adder = '_' + season.strip()
        outname = 'omi_ai_climo_' + first_date + '_' + last_date + season_adder + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Plot a monthly climatology 
def plotOMI_MonthClimo(OMI_data,month_idx,minlat = 60,save=False):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(OMI_data['DATES'][0],'%Y%m')
    end_date   = datetime.strptime(OMI_data['DATES'][-1],'%Y%m')

    # Make figure title
    date_month = datetime(year = 1,month = month_idx+1, day = 1).strftime('%B')
    title = 'OMI AI ' + date_month + ' Climatology\n'+\
        start_date.strftime('%b. %Y') + ' - ' + end_date.strftime('%b. %Y')

    # Make figure
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(OMI_data['LON'], OMI_data['LAT'],\
            OMI_data['MONTH_CLIMO'][month_idx,:,:],transform = datacrs,\
            cmap = colormap, vmin = -1.0, vmax = 1.5)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5), \
        orientation='horizontal',pad=0,aspect=50,shrink = 0.905,\
        label='Aerosol Index')
    ax.set_title(title)

    if(save == True):
        out_name = 'omi_ai_month_climo_' + date_month + '.png'
        plt.savefig(out_name,dpi=300)
        print("Saved image",out_name)
    else:
        plt.show()

# Plot a single swath of OMI data with total climatology subtracted
# mask_weakAI: removes AI values below 0.8 when plotting
def single_swath_anomaly_climo(OMI_data,swath_date,month_climo = True,\
        minlat = 60,row_max = 60,mask_weakAI = False, save=False): 
    # - - - - - - - - - - - - - - - -
    # Read in the single swath data
    # - - - - - - - - - - - - - - - -
    latmin = 60 
    
    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Set up values for gridding the AI data
    lat_gridder = latmin * 4.
    
    lat_ranges = np.arange(latmin,90.1,0.25)
    lon_ranges = np.arange(-180,180.1,0.25)
    base_path = '/home/bsorenson/data/OMI/H5_files/'
    year = swath_date[:4]
    date = swath_date[4:8]
    if(len(swath_date)==13):
        time = swath_date[9:]
    elif(len(swath_date)==12):
        time = swath_date[8:]
    else:
        time = ''
    total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+year+'m'+date+'t'+time+'*.he5',\
              shell=True).decode('utf-8').strip().split('\n')

    UVAI = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

    for fileI in range(len(total_list)):
        # read in data directly from HDF5 files
        print(total_list[fileI])
        data = h5py.File(total_list[fileI],'r')
        #ALBEDO = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
        #REFLECTANCE = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
        #CLD = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/CloudFraction']
        AI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
        #PIXEL= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/PixelQualityFlags']
        #MSMNT= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/MeasurementQualityFlags']
        #VZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/ViewingZenithAngle']
        #SZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle']
        #RAZ  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle']
        LAT = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
        LON = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
        #GRND = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags']
    
        #albedo = ALBEDO[:,:,0]   
        #reflectance = REFLECTANCE[:,:,0]   
        counter = 0
        #AI = AI[:,:,0]   
        # Loop over the values and rows 
        #for i in range(0,int(CBA2)):
        #for i in range(albedo.shape[0]):
        for i in range(AI.shape[0]):
            for j in range(0,row_max):
                #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
                if(AI[i,j]>-20):
                #if(plotAI[i,j]>-20):
                    # Only plot if XTrack flag is met
                    if((XTRACK[i,j] == 0) | ((XTRACK[i,j] & 4 == 4))):
                        # Print values to text file
                        if(LAT[i,j] > minlat):
                            counter+=1
    #                        fout.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.6f} {9:.6f} {10:.6f} {11:.6f}\n".format(\
    #                            LAT[i,j],LON[i,j],AI[i,j],0.5,SZA[i,j],VZA[i,j],RAZ[i,j], \
    #                            ALBEDO[i,j,0],ALBEDO[i,j,1],REFLECTANCE[i,j,0],\
    #                            REFLECTANCE[i,j,1],CLD[i,j]))
    
    
                        #if((plotXTrack[i,j] == 0) | ((plotXTrack[i,j] & 4 == 4))):
                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                            index2 = int(np.floor(LON[i,j]*4 + 720.))
                            #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
                            #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
                            
                            if(index1 < 0): index1 = 0
                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                            if(index2 < 0): index2 = 0                                                                                            
                            if(index2 > 1439): index2 = 1439
                       
                            #diff = reflectance[i,j] - albedo[i,j]
                            #if(diff<min_diff):
                            #    min_diff = diff
                            #if(diff>max_diff):
                            #    max_diff = diff
                       
                            #if(diff<0.2): 
                            #    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
                            UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
                            #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                                #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                            count[index2, index1] = count[index2,index1] + 1
    
    # Calculate the row-average AI for the secondary plot
    mask_avgs = np.nanmean(np.ma.masked_where(AI[:,:] < -20, AI[:,:]),axis=0)


    if(month_climo == True):
        # Determine which month to use based on the swath date
        monther = int(date[:2])-1
        local_climo = OMI_data['MONTH_CLIMO'][monther,:,:].T
        title_adder = 'Monthly Anomaly '
        file_adder = 'mnth_anom_'
    else:
        local_data = np.copy(OMI_data['AI'][:,:,:])
        local_mask = np.ma.masked_where(local_data == -999.9, local_data)
        local_climo = np.nanmean(local_mask, axis=0).T
        title_adder = 'Anomaly '
        file_adder = 'anom_'

    # Interpolate 2D data
    # Find where the high-res data fits in to the climo data
    begin_x = np.where(OMI_data['LAT'][:,0] == lat_ranges[0])[0][0]
    high_lats = np.arange(OMI_data['LAT'][begin_x,0], OMI_data['LAT'][-1,0]+0.25, 0.25)
    high_lons = np.arange(OMI_data['LON'][0,0], OMI_data['LON'][0,-1] + 0.25, 0.25)

    high_climo = np.zeros((high_lons.shape[0], high_lats.shape[0]))

    for yi in range(len(OMI_data['LON'][0,:])):
        for xj in range(len(OMI_data['LAT'][begin_x:,0])):
            high_climo[yi*4:yi*4+4, xj*4:xj*4+4] = local_climo[yi, xj+begin_x]

    # Find out where the high res climo data fits in to the
    # single swath data. The single swath data are assumed to be 
    # of a larger size than the climo data
    end_xj = np.where(lat_ranges == high_lats[-1])[0][0]
    end_yi = np.where(lon_ranges == high_lons[-1])[0][0]

    # Calculate anomalies
    UVAI_anomaly = UVAI[0:end_yi+1, 0:end_xj+1] - high_climo  

    # Mask grid boxes where the ob counts are zero
    plot_lat, plot_lon = np.meshgrid(high_lats,high_lons)
    mask_UVAI_anom = np.ma.masked_where(count[0:end_yi+1, 0:end_xj+1] == 0, UVAI_anomaly)
  
    ## Mask any data below AI values of 0.8
    #mask_UVAI_anom = np.ma.masked_where(mask_UVAI_anom < 0.8, mask_UVAI_anom)
 
    plt.close()
    fig1 = plt.figure(figsize=(8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
 
    plt.title('OMI AI '+title_adder + swath_date)
    #plt.title('OMI Reflectivity - Surface Albedo '+swath_date)
    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI_anom,transform = datacrs,cmap = colormap,\
            vmin = -2.0, vmax = 3.0)
    ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    #ax.set_xlim(-4170748.535086173,4167222.438879491)
    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='Aerosol Index Anomaly')
    ##cax = fig.add_axes([0.16,0.075,0.7,0.025])
    ##cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    ##cb.ax.set_xlabel('Aerosol Index')
    #cb.ax.set_xlabel('Reflectivity - Surface Albedo')
    #out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+str(row_max)+'.png'       
    if(save == True):
        out_name = 'omi_single_pass_ai_'+file_adder+swath_date+'_rows_0to'+str(row_max)+'.png'       
        plt.savefig(out_name)
        print('Saved image '+out_name)
    else: 
        plt.show()
    #axs[1].plot(mask_avgs)
    #axs[1].set_xlabel('Sensor Row')
    #axs[1].set_ylabel('Row Average Aerosol Index')
    

# Plot a single swath of OMI data with single-time swath climatology subtracted
# single_swath = the swath that will be corrected (format = YYYYMMDDHHMM)
# climo_date   = the date on the swath climatology directory (format = MMDDHHMM)
def single_swath_anomaly_time(single_swath,climo_date,minlat = 60,row_max = 60): 

    # - - - - - - - - - - - - - - - -
    # Read in the single swath data
    # - - - - - - - - - - - - - - - -
    latmin = 60 
    
    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Set up values for gridding the AI data
    lat_gridder = latmin * 4.
    
    lat_ranges = np.arange(latmin,90.1,0.25)
    lon_ranges = np.arange(-180,180.1,0.25)
    climo_base_path = '/home/bsorenson/data/OMI/swath_anomaly_files/'
    single_base_path = '/home/bsorenson/data/OMI/H5_files/'

    # Grab the path date/time for finding the swath climatology files
    year = single_swath[0:4]
    date = single_swath[4:8]
    if(len(single_swath)==13):
        time = single_swath[9:]
    elif(len(single_swath)==12):
        time = single_swath[8:]
    else:
        time = ''

    total_climo_list = subprocess.check_output('ls '+climo_base_path+climo_date+'/*.he5',\
              shell=True).decode('utf-8').strip().split('\n')

    total_single_list = subprocess.check_output('ls '+single_base_path+'OMI-Aura_L2-OMAERUV_'+\
              year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')

    UVAI_climo = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_climo = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    UVAI_single = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count_single = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

    print("Reading path climatology data")
    for fileI in range(len(total_climo_list)):
        # read in data directly from HDF5 files
        print(total_climo_list[fileI])
        data = h5py.File(total_climo_list[fileI],'r')
        #ALBEDO = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
        #REFLECTANCE = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
        #CLD = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/CloudFraction']
        AI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
        #PIXEL= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/PixelQualityFlags']
        #MSMNT= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/MeasurementQualityFlags']
        #VZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/ViewingZenithAngle']
        #SZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle']
        #RAZ  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle']
        LAT = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
        LON = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
        #GRND = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags']
    
        #albedo = ALBEDO[:,:,0]   
        #reflectance = REFLECTANCE[:,:,0]   
        counter = 0
        #AI = AI[:,:,0]   
        # Loop over the values and rows 
        #for i in range(0,int(CBA2)):
        #for i in range(albedo.shape[0]):
        for i in range(AI.shape[0]):
            for j in range(0,row_max):
                #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
                if(AI[i,j]>-20):
                #if(plotAI[i,j]>-20):
                    # Only plot if XTrack flag is met
                    if((XTRACK[i,j] == 0) | ((XTRACK[i,j] & 4 == 4))):
                        # Print values to text file
                        if(LAT[i,j] > minlat):
                            counter+=1
    #                        fout.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.6f} {9:.6f} {10:.6f} {11:.6f}\n".format(\
    #                            LAT[i,j],LON[i,j],AI[i,j],0.5,SZA[i,j],VZA[i,j],RAZ[i,j], \
    #                            ALBEDO[i,j,0],ALBEDO[i,j,1],REFLECTANCE[i,j,0],\
    #                            REFLECTANCE[i,j,1],CLD[i,j]))
    
    
                        #if((plotXTrack[i,j] == 0) | ((plotXTrack[i,j] & 4 == 4))):
                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                            index2 = int(np.floor(LON[i,j]*4 + 720.))
                            #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
                            #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
                            
                            if(index1 < 0): index1 = 0
                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                            if(index2 < 0): index2 = 0                                                                                            
                            if(index2 > 1439): index2 = 1439
                       
                            #diff = reflectance[i,j] - albedo[i,j]
                            #if(diff<min_diff):
                            #    min_diff = diff
                            #if(diff>max_diff):
                            #    max_diff = diff
                       
                            #if(diff<0.2): 
                            #    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
                            UVAI_climo[index2, index1] = (UVAI_climo[index2,index1]*count_climo[index2,index1] + \
                                AI[i,j])/(count_climo[index2,index1]+1)
                            #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                                #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                            count_climo[index2, index1] = count_climo[index2,index1] + 1
        data.close() 

    # Grab the single-pass data
    print("Reading single-swath data")
    for fileI in range(len(total_single_list)):
        # read in data directly from HDF5 files
        print(total_single_list[fileI])
        data2 = h5py.File(total_single_list[fileI],'r')
        AI  = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
        LAT = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
        LON = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
        XTRACK = data2['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
    
        counter = 0
        # Loop over the values and rows 
        for i in range(AI.shape[0]):
            for j in range(0,row_max):
                if(AI[i,j]>-20):
                #if(plotAI[i,j]>-20):
                    # Only plot if XTrack flag is met
                    if((XTRACK[i,j] == 0) | ((XTRACK[i,j] & 4 == 4))):
                        # Print values to text file
                        if(LAT[i,j] > minlat):
                            counter+=1
    
                            index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                            index2 = int(np.floor(LON[i,j]*4 + 720.))
                            
                            if(index1 < 0): index1 = 0
                            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                            if(index2 < 0): index2 = 0                                                                                            
                            if(index2 > 1439): index2 = 1439
                       
                            UVAI_single[index2, index1] = (UVAI_single[index2,index1]*count_single[index2,index1] + \
                                AI[i,j])/(count_single[index2,index1]+1)
                            count_single[index2, index1] = count_single[index2,index1] + 1
        data2.close() 
    ##!## Calculate the row-average AI for the secondary plot
    ##!#mask_avgs = np.nanmean(np.ma.masked_where(AI[:,:] < -20, AI[:,:]),axis=0)

    # Calculate anomalies
    UVAI_anomaly = UVAI_single - UVAI_climo  

    plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)
    mask_UVAI_climo = np.ma.masked_where(count_climo[:,:] == 0, UVAI_climo)
    mask_UVAI_anom = np.ma.masked_where(count_single[:,:] == 0, UVAI_anomaly)

    # Plot the climatology
    plt.close()
    fig1 = plt.figure(figsize=(8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
 
    plt.title('OMI Aerosol Index Climatology '+single_swath)
    #plt.title('OMI Reflectivity - Surface Albedo '+swath_date)
    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI_climo,transform = datacrs,cmap = colormap,\
            vmin = -2.0, vmax = 3.0)
    ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    #ax.set_xlim(-4170748.535086173,4167222.438879491)
    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='Aerosol Index')
    ##cax = fig.add_axes([0.16,0.075,0.7,0.025])
    ##cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    ##cb.ax.set_xlabel('Aerosol Index')
    #cb.ax.set_xlabel('Reflectivity - Surface Albedo')
    #out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+str(row_max)+'.png'       
    out_name = 'omi_single_pass_ai_single_climo_'+single_swath+'_rows_0to'+str(row_max)+'.png'       
    #out_name = 'omi_single_pass_refl_albedo_diff_'+swath_date+'_rows_0to'+str(row_max)+'.png'       
    plt.savefig(out_name)
    print('Saved image '+out_name)
  
    # Plot the anomalies
    plt.close()
    fig1 = plt.figure(figsize=(8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
 
    plt.title('OMI Aerosol Index Anomaly '+single_swath)
    #plt.title('OMI Reflectivity - Surface Albedo '+swath_date)
    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI_anom,transform = datacrs,cmap = colormap,\
            vmin = -1.0, vmax = 1.5)
    ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    #ax.set_xlim(-4170748.535086173,4167222.438879491)
    #ax.set_ylim(-2913488.8763307533,2943353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='Aerosol Index Anomaly')
    ##cax = fig.add_axes([0.16,0.075,0.7,0.025])
    ##cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    ##cb.ax.set_xlabel('Aerosol Index')
    #cb.ax.set_xlabel('Reflectivity - Surface Albedo')
    #out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+str(row_max)+'.png'       
    out_name = 'omi_single_pass_ai_single_anom_'+single_swath+'_rows_0to'+str(row_max)+'.png'       
    #out_name = 'omi_single_pass_refl_albedo_diff_'+swath_date+'_rows_0to'+str(row_max)+'.png'       
    plt.savefig(out_name)
    print('Saved image '+out_name)
    
    #axs[1].plot(mask_avgs)
    #axs[1].set_xlabel('Sensor Row')
    #axs[1].set_ylabel('Row Average Aerosol Index')
    
    plt.show()

def plot_OMI_v_MODIS(MODIS_data,OMI_single_swath,save=False):

    # Pull out data from dictionaries
    local_modis = np.copy(MODIS_data['data'])[0,:,:].T
    local_omi   = np.ma.MaskedArray.copy(OMI_single_swath['AI_algae'])

    # shape = (1441, 121)
    plot_lat, plot_lon = np.meshgrid(MODIS_data['lat'],MODIS_data['lon'])

    # Mask any data outside the Finland region

    ##!#local_modis[(local_modis == -9) |\
    ##!#            (np.argwhere(np.isnan(local_modis))) |\
    ##!#            ((plot_lon < 17.) | (plot_lon > 54)) |\
    ##!#            ((plot_lat < 68.) | (plot_lat > 76)) |\
    ##!#            (local_modis > 1.2) |\
    ##!#            (local_omi == np.nan)] = np.nan
    ##!#local_omi[(local_modis == -9) |\
    ##!#            (np.argwhere(np.isnan(local_modis))) |\
    ##!#            ((plot_lon < 17.) | (plot_lon > 54)) |\
    ##!#            ((plot_lat < 68.) | (plot_lat > 76)) |\
    ##!#            (local_modis > 1.2) |\
    ##!#            (local_omi == np.nan)] = np.nan

    ##!#local_modis = np.ma.masked_invalid(local_modis)
    ##!#local_omi   = np.ma.masked_invalid(local_omi)

    local_modis = np.ma.masked_where(local_modis == -9, local_modis)
    local_modis = np.ma.masked_invalid(local_modis)
    local_modis = np.ma.masked_where(((plot_lon < 17.) | (plot_lon > 54)), local_modis)
    local_modis = np.ma.masked_where(((plot_lat < 68.) | (plot_lat > 76)), local_modis)
    local_modis = np.ma.masked_where((local_modis > 1.2), local_modis)
    local_omi   = np.ma.masked_where(((plot_lon < 17.) | (plot_lon > 54)), local_omi)
    local_omi   = np.ma.masked_where(((plot_lat < 68.) | (plot_lat > 76)), local_omi)
    #local_modis = np.ma.masked_where(local_omi == np.nan,local_modis)  

    local_modis = local_modis[~local_omi.mask]
    local_omi   = local_omi[~local_omi.mask]

    local_omi   = local_omi[~local_modis.mask]
    local_modis = local_modis[~local_modis.mask]


    print(len(local_modis),len(local_omi))

    splitter = MODIS_data['titles'][0].split('/')[-1]
    plot_date = datetime.strptime(splitter[1:5],'%Y') + \
        relativedelta(days = int(splitter[5:8])-1)
    ##!#datacrs = ccrs.PlateCarree() 
    ##!#colormap = plt.cm.jet
    ##!#ax = plt.axes(projection =ccrs.NorthPolarStereo(central_longitude = 35.) )
    ##!#ax.set_extent([17,54,68,76],datacrs)
    ##!#ax.gridlines()
    ##!#ax.coastlines(resolution='50m')
    ##!#mesh = ax.pcolormesh(plot_lon,plot_lat,local_omi,\
    ##!##mesh = ax.pcolormesh(plot_lon,plot_lat,local_modis,\
    ##!#        transform = datacrs, \
    ##!##        transform = datacrs, norm = cm.LogNorm(vmin = 0.01, vmax = 20.),\
    ##!#        cmap = colormap)
    ##!#ax.set_title('OMI AI\n'+plot_date.strftime('%Y%m%d'))
    ##!##ax.set_title('MODIS Chlorophyll-α\n'+plot_date.strftime('%Y%m%d'))
    ##!#cbar = plt.colorbar(mesh,ticks = [-2.0,-1.0,0.0,1.0,2.0,3.0,],orientation='horizontal',pad=0,\
    ##!##cbar = plt.colorbar(mesh,ticks = [0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0],orientation='horizontal',pad=0,\
    ##!#    aspect=50,shrink = 0.905, label='Chlorophyll Concentration, OCI Algorithm (mg m$^{-3}$)')
    ##!#cbar.ax.set_xticklabels(['-2.0','-1.0','0.0','1.0','2.0','3.0'])
    ##!##cbar.ax.set_xticklabels(['0.01','0.02','0.05','0.1','0.2','0.5','1','2','5','10','20'])
    ##!#plt.show()
    final_modis = local_modis.flatten()
    final_omi   = local_omi.flatten()
    plt.close()
    fig1 = plt.figure()
    plt.scatter(local_modis,local_omi,color='black')
    #plt.scatter(final_modis,final_omi)
    #plt.scatter(MODIS_data['chlor_a'][0,:,:].T,OMI_single_swath['AI_algae'])

    plt.plot(np.unique(final_modis), np.poly1d(np.polyfit(final_modis,final_omi,1))(np.unique(final_modis)),\
        color='red')

    print("Pearson:  ",pearsonr(final_modis,final_omi))
    print("Spearman: ",spearmanr(final_modis,final_omi))

    plt.title(plot_date.strftime('%Y%m%d'))
    plt.xlabel('MODIS ' + MODIS_data['ptitle'])
    plt.ylabel('OMI AI')

    if(save == True):
        outname = 'omi_algae_v_modis_'+MODIS_data['grabber'] + \
            '_' + plot_date.strftime('%Y%m%d') + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()
