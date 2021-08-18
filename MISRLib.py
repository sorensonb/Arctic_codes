#!/home/blake.sorenson/anaconda3/bin/python
"""
  NAME:
    MISRLib.py   

  PURPOSE:

  PYTHON VERSION:
    2.6.6

  MODULES:
    - Matplotlib
    - mpl_toolkits.basemap
    - datetime
    - sys
    - numpy
    - netCDF4
    
"""

# Define the python version
python_version=3

import numpy as np
import sys
import gzip
import importlib
from datetime import datetime
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
from matplotlib.patches import Polygon
import matplotlib.colors as color
from matplotlib.colors import rgb2hex,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
from scipy import stats
from netCDF4 import Dataset
from glob import glob
# The commands module was discontinued for Python 3, so if the user
# is using python 2, import commands instead
if(sys.version_info[0] < 3):
    import commands
    python_version=2
else:
    import subprocess

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

    mapLon, mapLat = m(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([MISR_dict[dictkey][date]['avg'] for date in sorted(MISR_dict[dictkey].keys())])
    #avgs = avgs[np.where(avgs.mask==False)[0]]
    temp_dates = np.array(sorted(MISR_dict[dictkey].keys()))
    dates = temp_dates
    #dates = temp_dates[np.where(avgs.mask==False)[0]]
    x_vals = np.arange(0,len(dates))

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals,avgs)
    print(dictkey, slope*len(x_vals))
    #print(slope/len(test_dict[dictkey].keys())
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
    #print(len(Snsorted),Snsorted)
    if(2*sm    == num_d):
        tsslope=0.5*(Snsorted[sm]+Snsorted[sm+1])
    if(2*sm+1 == num_d): 
        tsslope=Snsorted[sm+1]
    regress_ts = x_vals*tsslope+intercept

    # Convert the dates to datetime objects
    ddates = [datetime.strptime(date,'%Y%m') for date in dates] 
    #ddates = [datetime.strptime(date.decode("utf-8"),'%Y%m') for date in dates] 

    label_dates = [ddate.strftime('%b %Y') for ddate in ddates]

    tickNumber = 12

    fig1 = plt.figure()
    plt.plot(avgs)
    plt.plot(x_vals,regress_y,linestyle='--',color='red',label='Control')
    plt.plot(x_vals,regress_ts,linestyle='--',color='blue',label='Thiel-Sen')
    plt.title(dictkey)
    plt.xticks(np.arange(0,len(avgs))[::-int(len(avgs)/tickNumber)],label_dates[::-int(len(avgs)/tickNumber)],rotation=45)
    plt.ylabel('Aerosol Optical Thickness')
    plt.legend()
    plt.show()

def onclick_climo(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = m(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([MISR_dict[dictkey][date] for date in sorted(MISR_dict[dictkey].keys())])
    aot_avg = np.average(avgs)
    print(dictkey, aot_avg)

def onclick_pval(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = m(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([MISR_dict[dictkey][date] for date in sorted(MISR_dict[dictkey].keys())])
    
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
    print(dictkey,pval)

def readMISR(inputfile,start_date,end_date):
    global MISR_dict 
    MISR_dict = {}
    lowest_lat = 30
    
    lat_ranges = np.arange(lowest_lat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)
    
    if(inputfile.strip().split('/')[-1].split('.')[-1]=='gz'):
        fin = gzip.open(inputfile,'rb')
    else:
        fin = open(inputfile,'r')
    for line in fin:
        templine = line.strip().split()
        if(len(templine)>1):
            # Create the key
            key = str(int(float(templine[1])))+'x'+str(int(float(templine[2])))

            #if((int(templine[0])>=start_date) & (int(templine[0])<=end_date)):
            #if(key is not None):
            #    if(templine[1]==key):
            #        #MISR_dict[key][templine[0]] = {}
            #        MISR_dict[key][templine[0]]=float(templine[2])
            #else:
            # If the current lat/lon pair are not found in the dictionary's
            # keys, then make a new subdictionary for it.
            if(key not in MISR_dict.keys()):
                MISR_dict[key] = {}
            # If the current lat/lon pair are already in the dictionary's
            # keys, then add the new data to the subdictionary
            #MISR_dict[templine[1]][templine[0]]={}
            if((int(templine[0]) >= start_date) & (int(templine[0]) <= end_date)):
                MISR_dict[key][templine[0]]=float(templine[3])
    
    # end for loop (file)
    fin.close()
    return MISR_dict

def readMISR_New(start_date,end_date,minlat=45,filetype='clr100'):
#def readMISR_New(start_date,end_date,dtype='aod_no_flags',minlat=45):
    global MISR_dict 
    MISR_dict = {}
    lowest_lat = minlat 
    
    lat_ranges = np.arange(lowest_lat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)
  
    #if(dtype=='aod_no_flags'): 
    #    misr_files = subprocess.check_output('ls /home/blake.sorenson/MISR/*newC*txt.gz',shell=True).decode('utf-8').strip().split('\n')
    #    #misr_files = subprocess.check_output('ls /home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/30to90/aod_no_flags/*txt.gz',shell=True).decode('utf-8').strip().split('\n')
    #    MISR_dict['data_type'] = 'aod_no_flags'
    #else:
    #    misr_files = subprocess.check_output('ls /home/blake.sorenson/MISR/*newC*txt.gz',shell=True).decode('utf-8').strip().split('\n')
    #    #misr_files = subprocess.check_output('ls /home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/30to90/raw_aod_no_flags/*txt.gz',shell=True).decode('utf-8').strip().split('\n')
    misr_files = subprocess.check_output('ls /home/blake.sorenson/MISR/clr100_V4/*'+\
        filetype+'*txt.gz',shell=True).decode('utf-8').strip().split('\n')

    
    for file in misr_files:
        print(file)
        #if(dtype=='aod_no_flags'):
        #    fdate = file.split('/')[-1].split('_')[-2]
        #else:
        #    fdate = file.split('/')[-1].split('_')[-1][:6]
        fdate = file.split('/')[-1].split('_')[-2]
        #if(file.split('/')[-1].split('.')[-1]=='gz'):
        #    fin = gzip.open(inputfile,'rb')
        #else:
        #    fin = open(inputfile,'r')
        if((int(fdate)>=int(start_date)) & (int(fdate)<=int(end_date))):
            fin = gzip.open(file,'rb')
            for line in fin:
                templine = line.strip().split()
                #if(len(templine)>1):
                # Create the key
                if(float(templine[2])!=-999.):
                    key = str(int(templine[0]))+'x'+str(int(templine[1]))

                    #if((int(templine[0])>=start_date) & (int(templine[0])<=end_date)):
                    #if(key is not None):
                    #    if(templine[1]==key):
                    #        #MISR_dict[key][templine[0]] = {}
                    #        MISR_dict[key][templine[0]]=float(templine[2])
                    #else:
                    # If the current lat/lon pair are not found in the dictionary's
                    # keys, then make a new subdictionary for it.
                    if(key not in MISR_dict.keys()):
                        MISR_dict[key] = {}
                    # If the current lat/lon pair are already in the dictionary's
                    # keys, then add the new data to the subdictionary
                    #MISR_dict[templine[1]][templine[0]]={}
                    #if((int(templine[0]) >= start_date) & (int(templine[0]) <= end_date)):
                    MISR_dict[key][fdate]={}
                    MISR_dict[key][fdate]['avg']=float(templine[2])
                    MISR_dict[key][fdate]['#_obs']=int(templine[3])
            
            # end for loop (file)
            fin.close()
    return MISR_dict

# Data period is of format YYYYMMDDHH
def readMISR_hrly(data_dt,minlat=60.0,season='all'):
    MISR_data_hrly = {}
    
    lat_ranges = np.arange(minlat,90.0,0.25)
    lon_ranges = np.arange(-180.0,180.0,0.25)

    # Grab all available MISR names
    base_path = '/home/bsorenson/data/MISR/'
    total_list = sorted(glob(base_path+'MISR*.nc'))

    # Convert the desired dt to a datetime object to use for finding the file
    # -----------------------------------------------------------------------
    day = False
    if(len(data_dt) == 10):
        str_fmt = "%Y%m%d%H"
        hour_adder = 1
        hour_subtracter = 1
    elif(len(data_dt) == 8):
        print("Daily average")
        day = True
        str_fmt = "%Y%m%d"
        hour_subtracter = 0
        hour_adder = 24
        day_adder = 'daily_'
    else:
        print("INVALID PLOT TIME")
        return    
    
    dt_data_begin = datetime.strptime(data_dt,str_fmt) 
    #    relativedelta(hours = hour_subtracter)
    dt_data_end   = datetime.strptime(data_dt,str_fmt) + \
         relativedelta(hours = hour_adder)

    print(dt_data_begin,dt_data_end)

    # Loop over the list of current MISR data files and find the one that  
    # corresponds to the desired datetime
    # --------------------------------------------------------------------
    good_list = []
    for tfile in total_list:
        data = Dataset(tfile,'r')
        time_str = data['4.4_KM_PRODUCTS/Time'].units
        base_date = datetime.strptime(time_str.split('.')[0],"seconds since %Y-%m-%dT%H:%M:%S")
        begin_date = base_date + relativedelta(seconds = data['4.4_KM_PRODUCTS/Time'][0])
        end_date = base_date + relativedelta(seconds = data['4.4_KM_PRODUCTS/Time'][-1])
        data.close()

        #if(((file_begin >= dt_data_begin) & (file_end <= dt_data_end)) | \
        #((file_begin >= dt_data_begin) & (file_end > dt_data_end)) | \
        #((file_end < dt_data_end) & (file_begin < dt_data_begin))):

        #print(begin_date,end_date)
        if(((end_date >= dt_data_begin) & (begin_date <= dt_data_end)) | \
            ((day) & (begin_date.day == dt_data_begin.day) & \
            (begin_date.month == dt_data_begin.month) & \
            (begin_date.day == dt_data_begin.day))):
        ##!#if(((dt_data_begin >= begin_date) & (dt_data_begin <= end_date)) | \
        ##!#   ((dt_data_end >= begin_date) & (dt_data_end < end_date))):
          # (dt_data_begin == end_date) | \
          # ((day == True) & (dt_data_begin == begin_date)) & \
          #  ((dt_data_end - relativedelta(hours=1)) == end_date)):
            print("Found matching file",tfile)
            good_list.append(tfile)
   
    if(len(good_list) == 0):
        print("ERROR: No valid files found")
        return MISR_data_hrly
 
    # Set up values for gridding the AI data
    lat_gridder = minlat * 4.

    aod_grid = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
 
    for fileI in range(len(good_list)):
        # read in data directly from HDF5 files
        print(good_list[fileI])
        data = Dataset(good_list[fileI],'r')
        lat   = data['4.4_KM_PRODUCTS/Latitude'][:,:].flatten()
        lon   = data['4.4_KM_PRODUCTS/Longitude'][:,:].flatten()
        aod   = data['4.4_KM_PRODUCTS/Aerosol_Optical_Depth'][:,:].flatten()
  
        # Loop over the values and rows 
        for ii in range(aod.shape[0]):
            #for jj in range(aod.shape[1]): 
            #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
            #if((local_time >= dt_data_begin) & (local_time < dt_data_end)):
            if((aod[ii] < 5000) and (aod[ii] > 0) and \
                    (lat[ii] >= minlat)):
                index1 = int(np.floor(lat[ii]*4 - lat_gridder))
                index2 = int(np.floor(lon[ii]*4 + 720.))
                
                if(index1 < 0): index1 = 0
                if(index1 > 719): index1 = 719
                if(index2 < 0): index2 = 0                                                                                            
                if(index2 > 1439): index2 = 1439
            
                try: 
                    #swf_grid[index2, index1] = (swf_grid[index2,index1]*count[index2,index1] + flux[i])/(count[index2,index1]+1)
                    aod_grid[index2, index1] = aod_grid[index2,index1] + aod[ii]
                except IndexError:
                    print(lat[ii,jj],lon[ii,jj],index1,index2)
                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                count[index2, index1] = count[index2,index1] + 1
        data.close()
  
    # Finally, average the values
    # --------------------------- 
    aod_grid[count > 0] = aod_grid[count > 0] / count[count > 0]
     
    MISR_data_hrly['param'] = 'aod'
    MISR_data_hrly['aod']   = aod_grid
    MISR_data_hrly['counts'] = count
    MISR_data_hrly['date'] = data_dt
    MISR_data_hrly['lat'] = lat_ranges
    MISR_data_hrly['lon'] = lon_ranges
     
    #start_date = datetime.strptime(str(start_date),'%Y%m')
    #end_date = datetime.strptime(str(end_date),'%Y%m') +timedelta(days=31)
    #
    #tempdate = data.variables['time'].units
    #orig_date = datetime.strptime(tempdate.split()[2]+' '+tempdate.split()[3],'%Y-%m-%d %H:%M:%S')

    return MISR_data_hrly

def readMISR_NCDF(infile='/home/bsorenson/Research/MISR/misr_aod_clr100_2000_2019.nc',\
                calc_month = True,minlat=65):
    # Read in data to netCDF object
    in_data = Dataset(infile,'r')

    # Set up dictionary to hold data
    MISR_data = {}

    # Set up date strings in the file
    MISR_data['MONTH'] = in_data['MONTH'][:]
    start_date = datetime(year = 2000, month = 3, day = 1)
    MISR_data['DATES'] = \
        [(start_date + relativedelta(months=mi)).strftime('%Y%m') for mi in \
        MISR_data['MONTH']]
    
    MISR_data['AOD'] = in_data['AOD'][:,:,:]
    MISR_data['OB_COUNT'] = in_data['OB_COUNT'][:,:,:]
    MISR_data['LAT'] = in_data['Latitude'][:,:]
    MISR_data['LON'] = in_data['Longitude'][:,:]
    MISR_data['VERSION'] = version = infile.split('/')[-1].split('_')[2]

    if(calc_month == True):
        MISR_data = calcMISR_MonthClimo(MISR_data)

    # to add months to datetime object, do
    ###from dateutil.relativedelta import relativedelta
    ###datetime.datetime(year=2004,month=10,day=1) + relativedelta(months=1)
       
    return MISR_data

def writeMISR_toNCDF(MISR_dict,filename,minlat=-90):
    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)
    
    nc = Dataset(filename,'w',format='NETCDF4')
 
    # Dimensions = lat, lon, time
    testdict = MISR_dict['50x5']
    testkeys = list(testdict.keys())
    num_lat = len(lat_ranges)
    num_lon = len(lon_ranges)
    #num_time = len(testdict.keys())
    #times = np.arange(num_time)

    base_date = datetime(year=2000,month=3,day=1)
    times = np.array([(datetime.strptime(tmpx,'%Y%m').year - base_date.year) \
        * 12 + datetime.strptime(tmpx,'%Y%m').month - base_date.month \
        for tmpx in testkeys])
    num_time = len(times) 
    
    n_time = nc.createDimension('nmth',num_time)
    n_lat  = nc.createDimension('dlat',num_lat)
    n_lon  = nc.createDimension('dlon',num_lon)

    MONTH=nc.createVariable('MONTH','i2',('nmth'))
    MONTH.description='Months since March 2000'
    LAT=nc.createVariable('Latitude','i2',('dlat','dlon'))
    LAT.description='Latitude'
    LAT.units='Degrees'
    LON=nc.createVariable('Longitude','i2',('dlat','dlon'))
    LON.description='Longitude'
    LON.units='Degrees'
    AOD = nc.createVariable('AOD','f4',('nmth','dlat','dlon'))
    AOD.description='Monthly Averaged Aerosol Optical Depth'
    OB_COUNT=nc.createVariable('OB_COUNT','i2',('nmth','dlat','dlon'))
    OB_COUNT.description='# of MISR AOD measurements used in each monthly average'

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
            if(dictkey not in MISR_dict):
                # Insert missing values for AOD and count
                AOD[:,i,j] = [-999.9 for m in range(num_time)]
                OB_COUNT[:,i,j] = [-99 for m in range(num_time)]
            else:
                for m in range(num_time):
                    timekey = testkeys[m]
                    if(timekey not in MISR_dict[dictkey]):
                        AOD[m,i,j] = -999.9
                        OB_COUNT[m,i,j] = -99
                    else:
                        AOD[m,i,j] = MISR_dict[dictkey][timekey]['avg']
                        OB_COUNT[m,i,j] = MISR_dict[dictkey][timekey]['#_obs']
    nc.close()
    print("Saved file",filename)

# This function assumes the data is being read from the netCDF file
# NOTE: Assume user is using new MISR climo file which starts in January
def calcMISR_MonthClimo(MISR_data):

    # Set up arrays to hold monthly climatologies
    month_climo = np.zeros((12,MISR_data['AOD'].shape[1],MISR_data['AOD'].shape[2]))

    # Mask the monthly averages
    local_data = np.copy(MISR_data['AOD'][:,:,:])
    local_mask = np.ma.masked_where(local_data == -999.9, local_data)
 
    # Calculate monthly climatologies
    for m_i in range(12):
        month_climo[m_i,:,:] = np.nanmean(local_mask[m_i::12,:,:],axis=0)
        print("Month: ",MISR_data['DATES'][m_i][4:]) 
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
    MISR_data['MONTH_CLIMO'] = month_climo

    return MISR_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#  Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def plotMISR(MISR_dict,save=False,season='',trend_type='standard',start_date='200003',end_date='201812',minlat=30.):
    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='thielSen'
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
        for lkey in sorted(MISR_dict.keys()):
            try:
                for tkey in sorted(MISR_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        MISR_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(MISR_dict.keys()):
            try:
                for tkey in sorted(MISR_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        MISR_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(MISR_dict.keys()):
            try:
                for tkey in sorted(MISR_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        MISR_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(MISR_dict.keys()):
            try:
                for tkey in sorted(MISR_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        MISR_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    # Find the lowest lat in the file
     
    # Set up the polar stereographic map
    fig1 = plt.figure()
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    if(minlat>10):
        m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
    else:
        m = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
    
    cmap = plt.get_cmap('bwr')
    #bcmap = plt.cm.set_cmap(cmap)
    if(season==''):
        v_min = -0.150 # Summer values
        mid_val = 0
        v_max = 0.150
    else:
        v_min = -0.400  # Whole-year values
        mid_val = 0
        v_max = 0.40
        #v_min = -0.150  # Whole-year total extinction values 
        #mid_val = 0
        #v_max = 0.15
    v_min=-0.2
    v_max = 0.2
    
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
    
    # Loop over all the keys and print(the regression slopes 
    # Grab the averages for the key
    max_slope = -10.
    min_slope = 10.
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
            #if(dictkey=='48x-97'):
            #    #print("sorted(MISR_dict[dictkey].keys())=",sorted(MISR_dict[dictkey].keys()))
            #    min_date = sorted(MISR_dict[dictkey].keys())[0].decode("utf-8")
            #    max_date = sorted(MISR_dict[dictkey].keys())[-1].decode("utf-8")
    
            # If no data are present for the curent lat/lon box, fill it with
            # black
            if(dictkey not in MISR_dict.keys()):
                color=(0,0,0,0)
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(MISR_dict[dictkey].keys())==0):
                color=(0,0,0,0)
            else:
                avgs = [MISR_dict[dictkey][date]['avg'] for date in \
                    sorted(MISR_dict[dictkey].keys())]
                if(len(avgs)<2):
                    color=(0,0,0,0)
                else:
                    # Check the current max and min
                    x_vals = np.arange(0,len(MISR_dict[dictkey].keys()))
                    # DEAL WITH MISSING VALUES
                    # Find the slope of the line of best fit for the time series of
                    # average data
                    avgs = np.ma.masked_array([MISR_dict[dictkey][date]['avg'] for date in sorted(MISR_dict[dictkey].keys())])
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(MISR_dict[dictkey].keys()))
                    dates = temp_dates
                    #dates = temp_dates[np.where(avgs.mask==False)[0]]
                    x_vals = np.arange(0,len(dates))
    
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
            #mxy = zip(mx,my)
            color2 = rgb2hex(color)
            #color2 = rgb2hex(color[:2]+color[3:])
            # Plot the box on the map using color
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,edgecolor=color2)
            #poly = Polygon(mxy,facecolor=color2,edgecolor=color2)
            plt.gca().add_patch(poly)
    
    minmax_diff = (max_slope-min_slope) 
    minmax_range = np.arange(min_slope,max_slope,minmax_diff/8.0)
    print("min slope = ",min_slope)
    print("max slope = ",max_slope)
    print("Range = ",minmax_range )

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

    #start_date = min_date
    #end_date = max_date
    title_string = 'Change in MISR Average Aerosol Optical'+      \
        ' Thickness\n'+start_date+' to '+end_date
    name_adder = '_aod_noflags'
    if(MISR_dict['data_type']!='aod_no_flags'):
        title_string = 'Change in MISR Average Aerosol Optical'+      \
            ' Thickness\n'+start_date+' to '+end_date+'\nRaw AOD No Flags'
        name_adder = ''
    #title_string = 'Change in MISR-2 Average Aerosol Optical Depth\n'+        \
    #    start_date+' to '+end_date
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
    cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    plt.xticks(rotation=45,fontsize=6)
    fig1.canvas.mpl_connect('button_press_event',onclick)
    if(save is True):
        if(spring is True):
            filename = 'misr_aot'+name_adder+'_'+trend_label+'_trends_'+str(int(minlat))+'to90_spring.png'
        elif(summer is True):
            filename = 'misr_aot'+name_adder+'_'+trend_label+'_trends_'+str(int(minlat))+'to90_summer.png'
        elif(autumn is True):
            filename = 'misr_aot'+name_adder+'_'+trend_label+'_trends_'+str(int(minlat))+'to90_autumn.png'
        elif(winter is True):
            filename = 'misr_aot'+name_adder+'_'+trend_label+'_trends_'+str(int(minlat))+'to90_winter.png'
        else:
            filename = 'misr_aot'+name_adder+'_'+trend_label+'_trends_'+str(int(minlat))+'to90_whole_year.png'
        
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

def plotMISR_MK(MISR_dict,save=False,season='',start_date='20003',end_date='201812'):
    lat_ranges = np.arange(30,90,1.0)
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
        for lkey in sorted(MISR_dict.keys()):
            for tkey in sorted(MISR_dict[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    MISR_dict[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(MISR_dict.keys()):
            if(lkey=='data_type'):
                continue
            else:
                for tkey in sorted(MISR_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        MISR_dict[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(MISR_dict.keys()):
            for tkey in sorted(MISR_dict[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    MISR_dict[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(MISR_dict.keys()):
            for tkey in sorted(MISR_dict[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    MISR_dict[lkey].pop(tkey)

    # Find the lowest lat in the file
    lowest_lat = float(sorted(MISR_dict.keys())[0].split('x')[0])
    lowest_lat=30
     
    # Set up the polar stereographic map
    fig1 = plt.figure()
    ax = plt.subplot(111)
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    if(minlat>10):
        m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
    else:
        m = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
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
            dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
            #print(dictkey)
            #if(dictkey=='48x-97'):
            #    #print("sorted(MISR_dict[dictkey].keys())=",sorted(MISR_dict[dictkey].keys()))
            #    min_date = sorted(MISR_dict[dictkey].keys())[0].decode("utf-8")
            #    max_date = sorted(MISR_dict[dictkey].keys())[-1].decode("utf-8")
    
            # If no data are present for the curent lat/lon box, fill it with
            # black
            if(dictkey not in MISR_dict.keys()):
                #print('Here 1')
                colorRED = 0.
                colorGREEN = 0.
                colorBLUE = 0.
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(MISR_dict[dictkey].keys())==0):
                #print('Here 2')
                colorRED = 0.
                colorGREEN = 0.
                colorBLUE = 0.
            else:
                avgs = [MISR_dict[dictkey][date]['avg'] for date in \
                    sorted(MISR_dict[dictkey].keys())]
                if(len(avgs)<2):
                    #print('Here 3')
                    colorRED = 0.
                    colorGREEN = 0.
                    colorBLUE = 0.
                else:
                    # Check the current max and min
                    x_vals = np.arange(0,len(MISR_dict[dictkey].keys()))
                    avgs = np.ma.masked_array([MISR_dict[dictkey][date]['avg'] for date in sorted(MISR_dict[dictkey].keys())])
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(MISR_dict[dictkey].keys()))
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
                    ##alpha=0.05
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
    ##color_range = [0.4,0.3,0.2,0.15,0.1,0.05]
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

    #start_date = min_date
    #end_date = max_date
    title_string = 'MISR Average Aerosol Optical'+      \
        ' Thickness Trend P-Values\n'+start_date+' to '+end_date
    name_adder = '_aod_noflags'
    if(MISR_dict['data_type']!='aod_no_flags'):
        title_string = 'MISR Average Aerosol Optical'+      \
            ' Thickness Trend P-Values\n'+start_date+' to '+end_date+'\nRaw AOD No Flags'
        name_adder = ''
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    elif(summer is True):
        title_string = title_string+'\nJune, July, August'
    elif(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    elif(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)
    fig1.canvas.mpl_connect('button_press_event',onclick)

    if(save is True):
        if(spring is True):
            filename = 'misr_aot'+name_adder+'_pvalues_'+str(int(lowest_lat))+'to90_spring.png'
        elif(summer is True):    
            filename = 'misr_aot'+name_adder+'_pvalues_'+str(int(lowest_lat))+'to90_summer.png'
        elif(autumn is True):   
            filename = 'misr_aot'+name_adder+'_pvalues_'+str(int(lowest_lat))+'to90_autumn.png'
        elif(winter is True):  
            filename = 'misr_aot'+name_adder+'_pvalues_'+str(int(lowest_lat))+'to90_winter.png'
        else:                 
            filename = 'misr_aot'+name_adder+'_pvalues_'+str(int(lowest_lat))+'to90_whole_year.png'
        
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

def plotMISR_Climo(MISR_dict,save=False,trend_type='standard',season='',start_date='200003',end_date='201812',minlat=30.):
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
        for lkey in sorted(MISR_dict.keys()):
            try:
                for tkey in sorted(MISR_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        MISR_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(MISR_dict.keys()):
            try:
                for tkey in sorted(MISR_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        MISR_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(MISR_dict.keys()):
            try:
                for tkey in sorted(MISR_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        MISR_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(MISR_dict.keys()):
            try:
                for tkey in sorted(MISR_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        MISR_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    # Find the lowest lat in the file
     
    # Set up the polar stereographic map
    fig1 = plt.figure()
    #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
    global m
    if(minlat>10):
        m = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
    else:
        m = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
    
    cmap = plt.get_cmap('jet')
    #cmap = plt.get_cmap('nipy_spectral')
    #bcmap = plt.cm.set_cmap(cmap)
    if(season==''):
        v_min = -0.150 # Summer values
        mid_val = 0
        v_max = 0.150
    else:
        v_min = -0.400  # Whole-year values
        mid_val = 0
        v_max = 0.40
        #v_min = -0.150  # Whole-year total extinction values 
        #mid_val = 0
        #v_max = 0.15
    v_min = 0.000  # Whole-year values
    mid_val = 0
    v_max = 0.50
    if(minlat>30.):
        v_min = 0.0
        mid_val = 0
        v_max = 0.3
    
    # Center the colorbar on zero
    #norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
    norm = Normalize(vmin=v_min,vmax=v_max)
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
    
    # Loop over all the keys and print(the regression slopes 
    # Grab the averages for the key
    max_avg = -10.
    min_avg = 10.
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
            #if(dictkey=='48x-97'):
            #    #print("sorted(MISR_dict[dictkey].keys())=",sorted(MISR_dict[dictkey].keys()))
            #    min_date = sorted(MISR_dict[dictkey].keys())[0].decode("utf-8")
            #    max_date = sorted(MISR_dict[dictkey].keys())[-1].decode("utf-8")
    
            # If no data are present for the curent lat/lon box, fill it with
            # black
            if(dictkey not in MISR_dict.keys()):
                color=(0,0,0,0)
            # If, for some reason, the key is made for the current lat/lon box
            # but the dictionary is empty. fill with black.
            elif(len(MISR_dict[dictkey].keys())==0):
                color=(0,0,0,0)
            else:
                avgs = [MISR_dict[dictkey][date]['avg'] for date in \
                    sorted(MISR_dict[dictkey].keys())]
                # Check the current max and min
                #x_vals = np.arange(0,len(MISR_dict[dictkey].keys()))
                # DEAL WITH MISSING VALUES
                # Find the slope of the line of best fit for the time series of
                # average data
                avgs = np.ma.masked_array([MISR_dict[dictkey][date]['avg'] for date in sorted(MISR_dict[dictkey].keys())])
                counts = np.ma.masked_array([MISR_dict[dictkey][date]['#_obs'] for date in sorted(MISR_dict[dictkey].keys())])
                #avg_aot = np.average(avgs)
                avg_aot = sum(avgs*counts)/sum(counts)
                #avgs = avgs[np.where(avgs.mask==False)[0]]
                temp_dates = np.array(sorted(MISR_dict[dictkey].keys()))
                #dates = temp_dates
                #dates = temp_dates[np.where(avgs.mask==False)[0]]
                #x_vals = np.arange(0,len(dates))
    
                #slope, intercept, r_value, p_value, std_err = \
                #    stats.linregress(x_vals,avgs)
                #slope *= len(x_vals)
    
                if(avg_aot>max_avg):
                    max_avg=avg_aot
                elif(avg_aot<min_avg):
                    min_avg=avg_aot
                if(avg_aot>v_max):
                    avg_aot=v_max
        
                color = mapper.to_rgba(avg_aot)

                # Add value to analysis regions (if possible)
                if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                    canada_avg+=avg_aot
                    canada_counts+=1
                elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                    siberia_avg+=avg_aot
                    siberia_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                    eus_avg+=avg_aot
                    eus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                    wus_avg+=avg_aot
                    wus_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                    ea_avg+=avg_aot
                    ea_counts+=1
                elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                    wa_avg+=avg_aot
                    wa_counts+=1
                elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                    europe_avg+=avg_aot
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
            #mxy = zip(mx,my)
            color2 = rgb2hex(color)
            #color2 = rgb2hex(color[:2]+color[3:])
            # Plot the box on the map using color
            poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,edgecolor=color2)
            #poly = Polygon(mxy,facecolor=color2,edgecolor=color2)
            plt.gca().add_patch(poly)
    
    minmax_diff = (max_avg-min_avg) 
    minmax_range = np.arange(min_avg,max_avg,minmax_diff/8.0)
    print("min avg = ",min_avg)
    print("max avg = ",max_avg)
    print("Range = ",minmax_range )

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

    #start_date = min_date
    #end_date = max_date
    title_string = 'MISR Average Aerosol Optical'+      \
        ' Thickness Climatology\n'+start_date+' to '+end_date
    #title_string = 'Change in MISR-2 Average Aerosol Optical Depth\n'+        \
    #    start_date+' to '+end_date
    name_adder = '_aod_noflags'
    if(MISR_dict['data_type']!='aod_no_flags'):
        title_string = 'MISR Average Aerosol Optical'+      \
            ' Thickness Climatology\n'+start_date+' to '+end_date+'\nRaw AOD No Flags'
        name_adder = ''
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
    cax = fig1.add_axes([0.27,0.1,0.5,0.05])
    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
    #cb.ax.set_xlabel('Change in average aerosol index')
    plt.xticks(rotation=45,fontsize=6)
    plt.xlabel('Aerosol Optical Thickness at 558 nm',fontsize=6)
    fig1.canvas.mpl_connect('button_press_event',onclick_climo)
    if(save is True):
        if(spring is True):
            filename = 'misr_aot'+name_adder+'_climo_'+str(int(minlat))+'to90_spring.png'
        elif(summer is True):
            filename = 'misr_aot'+name_adder+'_climo_'+str(int(minlat))+'to90_summer.png'
        elif(autumn is True):
            filename = 'misr_aot'+name_adder+'_climo_'+str(int(minlat))+'to90_autumn.png'
        elif(winter is True):
            filename = 'misr_aot'+name_adder+'_climo_'+str(int(minlat))+'to90_winter.png'
        else:
            filename = 'misr_aot'+name_adder+'_climo_'+str(int(minlat))+'to90_whole_year.png'
        
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()

# Designed to work with the netCDF data
def plotMISR_MonthTrend(MISR_data,month_idx=None,save=False,\
        trend_type='standard',season='',minlat=30.,return_trend=False):
    version = MISR_data['VERSION']
    label_dict = {
        'clr100': '0% Cloud Fraction'
    }
    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='_thielSen'

    if(month_idx == None):
        month_adder = ''
        month_idx = 0
        index_jumper = 1
        do_month = False
        v_max = 0.5
        v_min = -0.5
    else:
        month_adder = '_month'
        if(version == 'V003'):
            index_jumper = 12
        else:
            index_jumper = 6 
        do_month = True
        v_max = 0.7
        v_min = -0.7

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.bwr
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(MISR_data['DATES'][month_idx::index_jumper][0],'%Y%m')
    end_date   = datetime.strptime(MISR_data['DATES'][month_idx::index_jumper][-1],'%Y%m')

    # Make copy of MISR_data array
    print(MISR_data['DATES'][month_idx::index_jumper])
    local_data   = np.copy(MISR_data['AI'][month_idx::index_jumper,:,:])
    local_counts = np.copy(MISR_data['OB_COUNT'][month_idx::index_jumper,:,:])
    local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where(local_mask == -999.9, local_mask)
    ai_trends = np.zeros(local_data.shape[1:])
    print(local_data.shape)

    plt.close('all')
    fig0 = plt.figure()
    # 79 x 152: 14, 332
    latx = 14
    lonx = 332
    ## 75 x -150: 10, 30
    #latx = 10
    #lonx = 30
    plt.plot(MISR_data['DATES'][month_idx::index_jumper],local_mask[:,latx,lonx])
    plt.title(str(int(lat_ranges[latx])) + 'x'+str(int(lon_ranges[lonx]))+\
                '\n'+start_date.strftime('%b') + ' ' + MISR_data['VERSION']) 
    if(save == True):
        month_adder = ''
        if(do_month == True):
            month_adder = '_' + start_date.strftime('%b') 
        outname = 'omi_ai_time_series' + month_adder + '_' + \
            str(int(lat_ranges[latx])) + 'x'+str(int(lon_ranges[lonx])) + \
            '_' + version + '.png'
        plt.savefig(outname)
        print("Saved image",outname)
    #return local_mask[:,latx,lonx]
    

    # Make figure title
    #date_month = datetime(year = 1,month = month_idx+1, day = 1).strftime('%B')
    month_string = ''
    if(do_month == True):
        month_string = start_date.strftime('%B') + ' '

    title = 'MISR AI ' + month_string + 'Trends ('+version+')\n'+\
        start_date.strftime('%b. %Y') + ' - ' + end_date.strftime('%b. %Y')

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    for i in range(0,len(lat_ranges)-1):
        for j in range(0,len(lon_ranges)-1):
            # Check the current max and min
            x_vals = np.arange(0,len(local_mask[:,i,j]))
            # Find the slope of the line of best fit for the time series of
            # average data
            if(trend_type=='standard'): 
                slope, intercept, r_value, p_value, std_err = \
                    stats.linregress(x_vals,local_mask[:,i,j])
                ai_trends[i,j] = slope * len(x_vals)
            else:
                #The slope
                S=0
                sm=0
                nx = len(local_mask[:,i,j])
                num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
                Sn=np.zeros(num_d)
                for si in range(0,nx-1):
                    for sj in range(si+1,nx):
                        # Find the slope between the two points
                        Sn[sm] = (local_mask[si,i,j]-local_mask[sj,i,j])/\
                                 (si-sj) 
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
                    ai_trends[i,j] = slope*len(avgs)

    print(np.min(ai_trends),np.max(ai_trends))

    # Make figure
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(MISR_data['LON'], MISR_data['LAT'],\
            ai_trends,transform = datacrs,\
            cmap = colormap,vmin=v_min,vmax=v_max)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.2),\
        orientation='horizontal',pad=0,aspect=50,shrink = 0.845)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('UV Aerosol Index Trend',fontsize=16,weight='bold')
    ax.set_title(title)

    if(save == True):
        month_adder = ''
        if(do_month == True):
            month_adder = '_' + start_date.strftime('%b') 
        out_name = 'omi_ai_trend' + month_adder + '_' + \
            start_date.strftime('%Y%m') + '_' + end_date.strftime('%Y%m') + \
            '_' + version + '.png'
        plt.savefig(out_name,dpi=300)
        print("Saved image",out_name)
    else:
        plt.show()

    if(return_trend == True):
        return ai_trends

# Plots a single month of MISR climatology data (assumed to be from the 
# netCDF file).
def plotMISR_NCDF_SingleMonth(MISR_data,time_idx,clr=100,minlat=60,save=False):

    if(clr is not None):
        cld = 100 - clr

        comparer = '<'
        if(cld == 0):
            comparer = '='
        titler = "(cloud frac. " + comparer + ' '+str(cld) + '%)' 
        filer = '_clr'+str(clr)
    else:
        clr = ''
        titler = "(No extra cleaning)" 
        filer = ''


    # Make copy of MISR_data array
    local_data = np.copy(MISR_data['AOD'][time_idx,:,:])

    #start_date = datetime(year=2004,month=10,day=1)
    #start_date = datetime.strptime(MISR_data['DATES'][0],"%Y%m")

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Mask any missing values
    mask_AOD = np.ma.masked_where(local_data == -999.9, local_data)
    
    # Mask any data below the threshold latitude
    mask_AOD = np.ma.masked_where(MISR_data['LAT'] < minlat,mask_AOD)

    # Make figure title
    first_date = MISR_data['DATES'][time_idx]
    title = 'MISR AOD ' + titler + '\n'+first_date

    # Make figure
    plt.close()
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(MISR_data['LON'], MISR_data['LAT'],mask_AOD,transform = datacrs,cmap = colormap,\
            vmin = 0.0, vmax = 0.5)
    ax.set_extent([-180,180,minlat,90],datacrs)
    #ax.set_xlim(-3430748.535086173,3430748.438879491)
    #ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(0.0,4.1,0.1),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.845)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('Aerosol Optical Depth',fontsize=16,weight='bold')
    ax.set_title(title)

    if(save == True):
        outname = 'misr_aod_single_month_' + first_date + filer + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plotMISR_NCDF_Climo(MISR_data,start_idx=0,end_idx=190,season = '',minlat=60,\
                       save=False):

    season_dict = {
        'spring': '\nMAM',\
        'summer': '\nJJA',\
        'autumn': '\nSON',\
        'winter': '\nDJF',\
        '': ''
    }
    
    # Make copy of MISR_data array
    local_data = np.copy(MISR_data['AOD'])

    #start_date = datetime(year=2004,month=10,day=1)
    start_date = datetime.strptime(MISR_data['DATES'][0],"%Y%m")

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
   
    for m_idx in MISR_data['MONTH']:
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
    mask_AOD = np.ma.masked_where(local_data == -999.9, local_data)

    # Calculate climatology between desired indices
    MISR_climo = np.nanmean(mask_AOD[start_idx:end_idx,:,:],axis=0)

    # Make figure title
    first_date = month_objects[0].strftime("%Y%m")
    last_date = month_objects[-1].strftime("%Y%m")
    print(month_objects[0].strftime("%Y%m"),month_objects[-1].strftime("%Y%m"))
    title = 'MISR AOD Climatology\n'+first_date + ' - ' + last_date + season_dict[season]

    # Make figure
    plt.close()
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(MISR_data['LON'], MISR_data['LAT'],MISR_climo,transform = datacrs,cmap = colormap,\
            vmin = 0.0, vmax = 0.5)
    ax.set_extent([-180,180,minlat,90],datacrs)
    #ax.set_xlim(-3430748.535086173,3430748.438879491)
    #ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(0.0,1.1,0.1),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label='Aerosol Optical Depth')
    ax.set_title(title)

    if(save == True):
        season_adder = ''
        if(len(season.strip()) != 0):
            season_adder = '_' + season.strip()
        outname = 'misr_aod_climo_' + first_date + '_' + last_date + season_adder + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Plot a monthly climatology 
def plotMISR_MonthClimo(MISR_data,month_idx,minlat = 60,save=False):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(MISR_data['DATES'][month_idx],'%Y%m')
    end_date   = datetime.strptime(MISR_data['DATES'][-1],'%Y%m')

    # Make figure title
    #date_month = datetime(year = 1,month = month_idx+1, day = 1).strftime('%B')
    #date_month = datetime.strftime(start_date,"%Y%m").strptime("%B")
    title = 'MISR AOD ' + start_date.strftime('%B')+ ' Climatology\n'+\
        start_date.strftime('%b. %Y') + ' - ' + end_date.strftime('%b. %Y')

    # Make figure
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(MISR_data['LON'], MISR_data['LAT'],\
            MISR_data['MONTH_CLIMO'][month_idx,:,:],transform = datacrs,\
            cmap = colormap, vmin = 0.0, vmax = 0.5)
    ax.set_extent([-180,180,minlat,90],datacrs)
    #ax.set_xlim(-3430748.535086173,3430748.438879491)
    #ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5), \
        orientation='horizontal',pad=0,aspect=50,shrink = 0.905,\
        label='AOD')
    ax.set_title(title)

    if(save == True):
        out_name = 'misr_aod_month_climo_' + date_month + '.png'
        plt.savefig(out_name,dpi=300)
        print("Saved image",out_name)
    else:
        plt.show()

def plotMISR_hrly(MISR_data_hrly,minlat=60,save=False):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    local_data = np.copy(MISR_data_hrly['aod'])

    plot_lat, plot_lon = np.meshgrid(MISR_data_hrly['lat'],MISR_data_hrly['lon'])
    mask_flux = np.ma.masked_where(MISR_data_hrly['counts'] == 0, local_data)
   
    plt.close('all') 
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution = '50m')
    plt.title('MISR ' + MISR_data_hrly['param'] + ' ' + MISR_data_hrly['date'])
    #plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
    mesh = ax.pcolormesh(plot_lon, plot_lat,mask_flux,transform = datacrs,cmap = colormap,\
            vmin = 0., \
            vmax = 0.5)
    ax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
            #vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])
    cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.845,label=MISR_data_hrly['param'])
    
    if(save == True):
        out_name = 'misr_single_pass_aod_' + MISR_data_hrly['date'] + '.png'
        plt.savefig(out_name,dpi=300)
        print('Saved image '+out_name)
    else:
        plt.show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# The following two function are used to generate box plots that show the
# percent of data north of each latitude band that is removed by screening
# out all MISR scenes with clear pixel fractions below 100%. 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Calculate the percentage of data removed by the cloud screening parameter.
# Calculates the difference 
def count_diffs(MISR_data1,MISR_data2,start_idx = 202,minlat=65):
    # Data assumed to be the same shape
    local_data1 = np.copy(MISR_data1['OB_COUNT'])
    local_data2 = np.copy(MISR_data2['OB_COUNT'])

    mask_data1 = np.ma.masked_where(((local_data1 == -99) | \
        (MISR_data1['LAT'] < minlat)),local_data1)
    mask_data2 = np.ma.masked_where(((local_data2 == -99) | \
        (MISR_data2['LAT'] < minlat)),local_data2)

    ratios = np.zeros(len(range(start_idx,len(MISR_data1['DATES']))))

    for ii in range(start_idx,len(MISR_data1['DATES'])):
        count1 = np.nansum(mask_data1[ii,:,:])
        count2 = np.nansum(mask_data2[ii,:,:])  
        ratios[ii] = 1.0 - count2/count1
    #    print(MISR_data1['DATES'][ii],count2,count1,np.round((count2/count1),3))

    return ratios

# Generate box plots for the percentage of data removed by the
# cloud screening parameter for data north of each latitude band
# between 45 and 85 degrees north. 
# To make the figure, read in the misr_aod_control_200003_201912.nc and
# misr_aod_clr100_200003_201912.nc files into MISR_control and MISR_clr100,
# respectively, using readMISR_NCDF, and run the function using
# >>> count_diffs_box_plot(MISR_control,MISR_clr100)
def count_diffs_box_plot(MISR_data1,MISR_data2,save=False):
    lats = np.arange(45,85,5)
    total_ratios = np.zeros((len(lats),len(MISR_data1['DATES'])))

    plt.close('all')
    fig1,ax = plt.subplots()
    for ii in range(len(lats)):
        total_ratios[ii,:] = count_diffs(MISR_data1,\
            MISR_data2,start_idx=0,minlat=lats[ii])

    ax.violinplot([total_ratios[ii,:][~np.isnan(total_ratios[ii,:])]\
        for ii in range(len(lats))])
    labels = [str(int(lat))+'$^{o}$ N' for lat in np.arange(40,85,5)]
    ax.set_xticklabels(labels)
    yvals = ax.get_yticks()
    ax.set_yticklabels(['{:,.0%}'.format(x) for x in yvals])
    ax.set_title("Percent of MISR AOD data north of each latitude\nremoved"+\
        " by inclusion of only 100% clear scenes")
    ax.grid(axis='y')
    if(save == True):
        outname = "misr_aod_pcnt_removal_ctrl_vs_clr100.png"
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()
