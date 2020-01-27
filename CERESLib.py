#!/usr/bin/env python
"""
  NAME:
    CERESLib.py   

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
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Polygon
import matplotlib.colors as color
from matplotlib.colors import rgb2hex,Normalize
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
from scipy import stats
from netCDF4 import Dataset
from os import system
import glob
# The commands module was discontinued for Python 3, so if the user
# is using python 2, import commands instead

#global sm
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
def onclick_swf(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = sm(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([CERES_dict[dictkey][date]['swf'] for date in sorted(CERES_dict[dictkey].keys())])
    #avgs = avgs[np.where(avgs.mask==False)[0]]
    temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
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
    plt.ylabel('TOA Short Wave Flux')
    plt.legend()
    plt.show()

def onclick_lwf(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = lm(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([CERES_dict[dictkey][date]['lwf'] for date in sorted(CERES_dict[dictkey].keys())])
    #avgs = avgs[np.where(avgs.mask==False)[0]]
    temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
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

    ## Test np.polyfit
    #testz = np.polyfit(x_vals,avgs,1)
    #testp = np.poly1d(testz)



    # Convert the dates to datetime objects
    ddates = [datetime.strptime(date,'%Y%m') for date in dates] 
    #ddates = [datetime.strptime(date.decode("utf-8"),'%Y%m') for date in dates] 

    label_dates = [ddate.strftime('%b %Y') for ddate in ddates]

    tickNumber = 12

    fig1 = plt.figure()
    plt.plot(avgs)
    plt.plot(x_vals,regress_y,linestyle='--',color='red',label='Control')
    plt.plot(x_vals,regress_ts,linestyle='--',color='blue',label='Thiel-Sen')
    #plt.plot(x_vals,testp(x_vals),':',color='green',label='np.polyfit')
    plt.title(dictkey)
    plt.xticks(np.arange(0,len(avgs))[::-int(len(avgs)/tickNumber)],label_dates[::-int(len(avgs)/tickNumber)],rotation=45)
    plt.ylabel('TOA Long Wave Flux')
    plt.legend()
    plt.show()

def onclick_climo_swf(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = sm(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([CERES_dict[dictkey][date]['swf'] for date in sorted(CERES_dict[dictkey].keys())])
    aot_avg = np.average(avgs)
    print(dictkey, aot_avg)

def onclick_climo_lwf(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = lm(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([CERES_dict[dictkey][date]['lwf'] for date in sorted(CERES_dict[dictkey].keys())])
    aot_avg = np.average(avgs)
    print(dictkey, aot_avg)

def onclick_pval_swf(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = sm(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([CERES_dict[dictkey][date] for date in sorted(CERES_dict[dictkey].keys())])
    
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

def onclick_pval_lwf(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    mapLon, mapLat = lm(ix,iy,inverse=True)
    mapLat = int(np.floor(mapLat))
    mapLon = int(np.floor(mapLon))

    dictkey = str(mapLat)+'x'+str(mapLon)
    avgs = np.ma.masked_array([CERES_dict[dictkey][date]['lwf'] for date in sorted(CERES_dict[dictkey].keys())])
    
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

def readCERES_NewGridRaw(start_date,end_date,minlat=30.,dtype='all'):
    global CERES_dict 
    CERES_dict = {}
    
    if(dtype=='clear'):
        longkey = 'toa_lw_clr_mon'
        shortkey= 'toa_sw_clr_mon'
    else:
        longkey = 'toa_lw_all_mon'
        shortkey= 'toa_sw_all_mon'
    CERES_dict['dtype']='gridded'
    CERES_dict['type_flag']=''

    lat_ranges = np.arange(-90,90,1.0)
    lon_ranges = np.arange(0,360,1.0)
    time_ranges = np.arange(0,226)
 
    data = Dataset('/home/bsorenson/HighLatitudeStudy/CERES/CERES_EBAF-TOA_Ed4.0_Subset_200003-201812.nc','r')

    # Use minlat to find the lat indices
    lat_start = np.where(data.variables['lat'][:]-0.5==minlat)
    lat_indices = np.arange(lat_start[0][0],len(lat_ranges),1)
    lon_indices = np.arange(0,len(lon_ranges))
     
    start_date = datetime.strptime(str(start_date),'%Y%m')
    end_date = datetime.strptime(str(end_date),'%Y%m') +timedelta(days=31)
    
    tempdate = data.variables['time'].units
    orig_date = datetime.strptime(tempdate.split()[2]+' '+tempdate.split()[3],'%Y-%m-%d %H:%M:%S')

    for ii in time_ranges:
        curr_date = orig_date+timedelta(days=int(data.variables['time'][ii]))
        if((curr_date>=start_date) & (curr_date<=end_date)):
            #extract data
            timekey = curr_date.strftime('%Y%m')
            print(timekey)
            for latx in lat_indices:
                print(lat_ranges[latx])
                for ij in lon_indices:
                    if(lon_ranges[ij]>179):
                        key = str(int(lat_ranges[latx]))+'x'+str(int(lon_ranges[ij])-360)
                    else:
                        key = str(int(lat_ranges[latx]))+'x'+str(int(lon_ranges[ij]))
                    if(key not in CERES_dict.keys()):
                        CERES_dict[key] = {}
                    CERES_dict[key][timekey] = {}
                    CERES_dict[key][timekey]['swf']=float(data.variables[shortkey][ii,latx,ij])
                    CERES_dict[key][timekey]['swf_cc']=1
                    CERES_dict[key][timekey]['lwf']=float(data.variables[longkey][ii,latx,ij])
                    CERES_dict[key][timekey]['lwf_cc']=1
                    
    # end for loop (file)
    data.close()
    return CERES_dict

def writeCERES_NewGrid(CERES_dict):
    filestart='/home/bsorenson/HighLatitudeStudy/CERES/average_files/Both/ceres_both_gridded_flux_'
    outkeys = CERES_dict['0x0'].keys()
    lat_ranges = np.arange(-90,90,1)
    lon_ranges = np.arange(-180,180,1)
    for okey in outkeys:
        fname = filestart+okey+'.txt'
        outfile = open(fname,'w')
        for latx in lat_ranges:
            for lony in lon_ranges:
                dictkey = str(int(latx))+'x'+str(int(lony))
                checky = lony+180
                if(checky>179):
                    truekey = str(int(latx))+'x'+str(int(checky)-360)
                else:
                    truekey = str(int(latx))+'x'+str(int(checky))
                outfile.write(truekey+' '+str(np.round(CERES_dict[dictkey][okey]['lwf'],3))+' '+\
                              str(np.round(CERES_dict[dictkey][okey]['swf'],3))+'\n')
        outfile.close()
        system('gzip '+fname)
        print('gzip '+fname)

    
def readCERES_EBAF_Raw(start_date,end_date,minlat=30.,sfc_toa='sfc',dtype='total'):
    global CERES_dict 
    CERES_dict = {}
    
    if(dtype=='clear'):
        up_key= sfc_toa+'_sw_up_clr_c_mon'
        down_key= sfc_toa+'_sw_down_clr_c_mon'
        CERES_dict['type_flag']='_clr_clr'
    else:
        up_key= sfc_toa+'_sw_up_clr_t_mon'
        down_key= sfc_toa+'_sw_down_clr_t_mon'
        CERES_dict['type_flag']='_tot_clr'
    CERES_dict['dtype']='EBAF_gridded'

    lat_ranges = np.arange(-90,90,1.0)
    lon_ranges = np.arange(0,360,1.0)
    time_ranges = np.arange(0,226)
 
    data_path = '/home/shared/CERES/EBAF/'
    # Grab all the data files here
    file_list = glob.glob(data_path+'*.nc')
    #data = Dataset('/home/bsorenson/HighLatitudeStudy/CERES/CERES_EBAF-TOA_Ed4.0_Subset_200003-201812.nc','r')

    # Use minlat to find the lat indices
    lon_indices = np.arange(0,len(lon_ranges))
     
    for infile in file_list:
        # Grab the file month from the file name
        file_date = infile[48:54]
        if((int(file_date)>=start_date) & (int(file_date)<=end_date)):
            #extract data
            print(file_date)
            data = Dataset(infile)
            lat_start = np.where(data.variables['lat'][:]-0.5==minlat)
            lat_indices = np.arange(lat_start[0][0],len(lat_ranges),1)
            for latx in lat_indices:
                print(lat_ranges[latx])
                for ij in lon_indices:
                    if(lon_ranges[ij]>179):
                        key = str(int(lat_ranges[latx]))+'x'+str(int(lon_ranges[ij])-360)
                    else:
                        key = str(int(lat_ranges[latx]))+'x'+str(int(lon_ranges[ij]))
                    if(key not in CERES_dict.keys()):
                        CERES_dict[key] = {}
                    CERES_dict[key][file_date] = {}
                    CERES_dict[key][file_date]['clr_alb']=float(data.variables[sfc_toa+'_sw_up_clr_c_mon'][0,latx,ij]/  \
                                                                data.variables[sfc_toa+'_sw_down_clr_c_mon'][0,latx,ij])
                    CERES_dict[key][file_date]['clr_alb_cc']=1
                    CERES_dict[key][file_date]['tot_alb']=float(data.variables[sfc_toa+'_sw_up_clr_t_mon'][0,latx,ij]/  \
                                                                data.variables[sfc_toa+'_sw_down_clr_t_mon'][0,latx,ij])
                    CERES_dict[key][file_date]['tot_alb_cc']=1
                    
    # end for loop (file)
            data.close()
    return CERES_dict

def readCERES(start_date,end_date,lowest_lat,dtype='terra',threshold=False,lat_fix=True,albedo=False,clr=False):
    global CERES_dict 
    CERES_dict = {}
    
    lat_ranges = np.arange(lowest_lat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)
 
    if(dtype=='terra'): 
        CERES_dict['dtype']='terra'
        if(threshold==True):
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/TERRA/gridded/N90to90/threshold_40/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']='_threshold40'
        elif(lat_fix==True):
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/TERRA/gridded/N90to90/lat_fix/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']='_latfix'
        elif(albedo==True):
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/TERRA/temp_files/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']='_albedo'
        elif(clr==True):
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/TERRA/temp_clr_files/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']='_clr'
        else:
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/TERRA/gridded/N90to90/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']=''
    elif(dtype=='aqua'):
        CERES_dict['dtype']='aqua'
        if(threshold==True):
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/AQUA/gridded/N90to90/threshold_40/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']='_threshold40'
        elif(lat_fix==True):
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/AQUA/gridded/N90to90/lat_fix/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']='_latfix'
        elif(albedo==True):
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/AQUA/temp_files/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']='_albedo'
        elif(clr==True):
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/AQUA/temp_clr_files/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']='_clr'
        else:
            ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/AQUA/gridded/N90to90/ceres_avg_gridded_flux*txt.gz')
            CERES_dict['type_flag']=''
    elif(dtype=='both'):
        CERES_dict['dtype']='both'
        CERES_dict['type_flag']=''
        ceres_files = glob.glob('/home/bsorenson/HighLatitudeStudy/CERES/average_files/Both/ceres_both_gridded_flux*txt.gz')

    
    for file in ceres_files:
        fdate = file.split('/')[-1].split('_')[-1][:6]
        print(file,fdate)
        #if(file.split('/')[-1].split('.')[-1]=='gz'):
        #    fin = gzip.open(inputfile,'rb')
        #else:
        #    fin = open(inputfile,'r')
        if((int(fdate)>=int(start_date)) & (int(fdate)<=int(end_date))):
            fin = gzip.open(file,'rb')
            for line in fin:
                templine = line.strip().split()
                if(float(templine[2])!=-999.):
                    # If reading terra-only data, use this method
                    if((dtype=='terra') | (dtype=='aqua')):
                        if(int(templine[1])>179):
                            key = str(int(templine[0]))+'x'+str(int(templine[1])-360)
                        else:
                            key = str(int(templine[0]))+'x'+str(int(templine[1]))
                    else:
                        key = templine[0].decode('utf-8')

                    #if((int(templine[0])>=start_date) & (int(templine[0])<=end_date)):
                    #if(key is not None):
                    #    if(templine[1]==key):
                    #        #CERES_dict[key][templine[0]] = {}
                    #        CERES_dict[key][templine[0]]=float(templine[2])
                    #else:
                    # If the current lat/lon pair are not found in the dictionary's
                    # keys, then make a new subdictionary for it.
                    if(key not in CERES_dict.keys()):
                        CERES_dict[key] = {}
                    # If the current lat/lon pair are already in the dictionary's
                    # keys, then add the new data to the subdictionary
                    #CERES_dict[templine[1]][templine[0]]={}
                    #if((int(templine[0]) >= start_date) & (int(templine[0]) <= end_date)):
                    CERES_dict[key][fdate]= {}
                    if(dtype=='both'):
                        CERES_dict[key][fdate]['lwf'] = float(templine[1])
                        CERES_dict[key][fdate]['lwf_cc'] = 1 
                        CERES_dict[key][fdate]['swf'] = float(templine[2])
                        CERES_dict[key][fdate]['swf_cc'] = 1 
                    else:
                        CERES_dict[key][fdate]['lwf'] = float(templine[2])
                        CERES_dict[key][fdate]['lwf_cc'] = int(templine[3])
                        CERES_dict[key][fdate]['swf'] = float(templine[4])
                        CERES_dict[key][fdate]['swf_cc'] = int(templine[5])
                    if((albedo==True) | (clr==True)):
                        CERES_dict[key][fdate]['alb'] = float(templine[6])
                        CERES_dict[key][fdate]['alb_cc'] = int(templine[7])
            
            # end for loop (file)
            fin.close()
    return CERES_dict

def readCERES_Lat(start_date,end_date,minlat):
    global CERES_dict 
    CERES_dict = {}
    
    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)
 
    if(minlat==30): 
        ceres_files = subprocess.check_output('ls /home/bsorenson/HighLatitudeStudy/CERES/average_files/TERRA/30to90/ceres_avg_flux_*txt',shell=True).decode('utf-8').strip().split('\n')
    elif(minlat==-20):
        ceres_files = subprocess.check_output('ls /home/bsorenson/HighLatitudeStudy/CERES/average_files/TERRA/N20to90/ceres_avg_flux_*txt',shell=True).decode('utf-8').strip().split('\n')

    
    for file in ceres_files:
        print(file)
        fdate = file.split('/')[-1].split('_')[-2]
        fdate = file.split('/')[-1].split('_')[-1][:6]
        #if(file.split('/')[-1].split('.')[-1]=='gz'):
        #    fin = gzip.open(inputfile,'rb')
        #else:
        #    fin = open(inputfile,'r')
        if((int(fdate)>=int(start_date)) & (int(fdate)<=int(end_date))):
            fin = open(file,'r')
            for line in fin:
                templine = line.strip().split()
                #if(len(templine)>1):
                # Create the key
                if(float(templine[2])!=-999.):
                    key = str(int(templine[0]))

                    #if((int(templine[0])>=start_date) & (int(templine[0])<=end_date)):
                    #if(key is not None):
                    #    if(templine[1]==key):
                    #        #CERES_dict[key][templine[0]] = {}
                    #        CERES_dict[key][templine[0]]=float(templine[2])
                    #else:
                    # If the current lat/lon pair are not found in the dictionary's
                    # keys, then make a new subdictionary for it.
                    if(key not in CERES_dict.keys()):
                        CERES_dict[key] = {}
                    # If the current lat/lon pair are already in the dictionary's
                    # keys, then add the new data to the subdictionary
                    #CERES_dict[templine[1]][templine[0]]={}
                    #if((int(templine[0]) >= start_date) & (int(templine[0]) <= end_date)):
                    CERES_dict[key][fdate] = {}
                    CERES_dict[key][fdate]['lwf']=float(templine[1])
                    CERES_dict[key][fdate]['swf']=float(templine[3])
            
            # end for loop (file)
            fin.close()
    return CERES_dict


def plotCERES(CERES_dict,save=False,ptype='swf',season='',trend_type='standard',start_date='200101',end_date='201512',minlat=30.,nozero=False, \
              relative_trend=False):
    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='thielSen'
    if(relative_trend is True):
        trend_label='relative'
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
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    # Find the lowest lat in the file
    #lowest_lat = float(sorted(CERES_dict.keys())[0].split('x')[0])
    #lowest_lat=30
    print("Got here, ptype=",ptype) 
    if((ptype=='swf') | (ptype=='both')):
        # Set up the polar stereographic map
        fig1 = plt.figure()
        global sm
        #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
        if(minlat>10):
            sm = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
        else:
            sm = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
        sm.drawcoastlines()
        sm.drawparallels(np.arange(-80.,81.,20.))
        sm.drawmeridians(np.arange(-180.,181.,20.))    
        
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
        v_min=-30.
        v_max = 30.
        if(relative_trend is True):
            v_min = -20.
            v_max = 20.
        
        # Center the colorbar on zero
        norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
        #norm = Normalize(vmin=vmin,vmax=vmax)
        mapper = ScalarMappable(norm=norm,cmap=cmap)

        # Set up initial values for the analysis regions
        swf_canada_avg = 0.
        swf_canada_counts = 0
        swf_siberia_avg = 0.
        swf_siberia_counts = 0 
        swf_eus_avg = 0.
        swf_eus_counts = 0
        swf_wus_avg = 0.
        swf_wus_counts = 0
        swf_ea_avg = 0.
        swf_ea_counts = 0
        swf_wa_avg = 0.
        swf_wa_counts = 0
        swf_europe_avg = 0.
        swf_europe_counts = 0

        # Loop over all the keys and print(the regression slopes 
        # Grab the averages for the key
        max_slope = -10.
        min_slope = 10.
        for i in range(0,len(lat_ranges)-1):
            for j in range(0,len(lon_ranges)-1):
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                #if(dictkey=='48x-97'):
                #    #print("sorted(CERES_dict[dictkey].keys())=",sorted(CERES_dict[dictkey].keys()))
                #    min_date = sorted(CERES_dict[dictkey].keys())[0].decode("utf-8")
                #    max_date = sorted(CERES_dict[dictkey].keys())[-1].decode("utf-8")
        
                # If no data are present for the curent lat/lon box, fill it with
                # black
                if(dictkey not in CERES_dict.keys()):
                    color=(0,0,0,0)
                # If, for some reason, the key is made for the current lat/lon box
                # but the dictionary is empty. fill with black.
                elif(len(CERES_dict[dictkey].keys())==0):
                    color=(0,0,0,0)
                else:
                    # Deal with short wave fluxes first
                    swf_avgs = [CERES_dict[dictkey][date]['swf'] for date in \
                        sorted(CERES_dict[dictkey].keys())]
                    if(len(swf_avgs)<2):
                        color=(0,0,0,0)
                    else:
                        # Check the current max and min
                        x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                        # DEAL WITH MISSING VALUES
                        # Find the slope of the line of best fit for the time series of
                        # average data
                        swf_avgs = np.ma.masked_array([CERES_dict[dictkey][date]['swf'] for date in sorted(CERES_dict[dictkey].keys())])
                        #avgs = avgs[np.where(avgs.mask==False)[0]]
                        temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
                        dates = temp_dates

                        if(relative_trend is True): 
                            flux_avg = np.average(swf_avgs)

                        if(nozero == True):
                            # Remove zeroes
                            zero_locs = np.where(swf_avgs!=0)
                            swf_avgs = swf_avgs[zero_locs]
                            dates = dates[zero_locs]

                        #dates = temp_dates[np.where(avgs.mask==False)[0]]
                        x_vals = np.arange(0,len(dates))
        
                        if(trend_type=='standard'): 
                            slope, intercept, r_value, p_value, std_err = \
                                stats.linregress(x_vals,swf_avgs)
                            slope *= len(x_vals)
                        else:
                            #The slope
                            S=0
                            sm=0
                            nx = len(swf_avgs)
                            num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
                            Sn=np.zeros(num_d)
                            for si in range(0,nx-1):
                                for sj in range(si+1,nx):
                                    # Find the slope between the two points
                                    Sn[sm] = (swf_avgs[si]-swf_avgs[sj])/(si-sj) 
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
                                slope = slope*len(swf_avgs)
        
                        if(slope>max_slope):
                            max_slope=slope
                        elif(slope<min_slope):
                            min_slope=slope
           
                        if(relative_trend is True):
                            if(lat_ranges[i]==0.):
                                print(str(lon_ranges[j])+" Slope=",str(slope)," Avg=",str(flux_avg)," new=",str((slope/flux_avg)*100.))
                            slope = (slope/flux_avg)*100.
                        color = mapper.to_rgba(slope)

                        # Add value to analysis regions (if possible)
                        if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                            swf_canada_avg+=slope
                            swf_canada_counts+=1
                        elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                            swf_siberia_avg+=slope
                            swf_siberia_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                            swf_eus_avg+=slope
                            swf_eus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                            swf_wus_avg+=slope
                            swf_wus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                            swf_ea_avg+=slope
                            swf_ea_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                            swf_wa_avg+=slope
                            swf_wa_counts+=1
                        elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                            swf_europe_avg+=slope
                            swf_europe_counts+=1
        
                # Find the y coordinates of the current LatxLon box
                y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
                # Find the x coordinates of the current LatxLon box
                x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
                # Convert x and y into map coordinates
                mx, my = sm(x,y)
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

        swf_canada_avg = swf_canada_avg/swf_canada_counts
        swf_siberia_avg = swf_siberia_avg/swf_canada_counts
        swf_eus_avg = swf_eus_avg/swf_eus_counts
        swf_wus_avg = swf_wus_avg/swf_wus_counts
        swf_ea_avg = swf_ea_avg/swf_ea_counts
        swf_wa_avg = swf_wa_avg/swf_wa_counts
        swf_europe_avg = swf_europe_avg/swf_europe_counts
        print("Canada       avg = ",swf_canada_avg,"  counts = ",swf_canada_counts)
        print("Siberia      avg = ",swf_siberia_avg,"  counts = ",swf_siberia_counts)
        print("Eastern US   avg = ",swf_eus_avg,"  counts = ",swf_eus_counts)
        print("Western US   avg = ",swf_wus_avg,"  counts = ",swf_wus_counts)
        print("Eastern Asia avg = ",swf_ea_avg,"  counts = ",swf_ea_counts)
        print("Western Asia avg = ",swf_wa_avg,"  counts = ",swf_wa_counts)
        print("Europe       avg = ",swf_europe_avg,"  counts = ",swf_europe_counts)

        #start_date = min_date
        #end_date = max_date
        title_string = 'Change in CERES Average TOA Short Wave Flux\n'+       \
                       start_date+' to '+end_date
        if(relative_trend is True):
            title_string = 'Relative change in CERES Average TOA Short Wave Flux\n'+       \
                           start_date+' to '+end_date
        #title_string = 'Change in CERES-2 Average Aerosol Optical Depth\n'+        \
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
        fig1.canvas.mpl_connect('button_press_event',onclick_swf)
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))
        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_spring.png'
            elif(summer is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_summer.png'
            elif(autumn is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_autumn.png'
            elif(winter is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_winter.png'
            else:
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            if(ptype!='both'):
                plt.show()
    if((ptype=='lwf') | (ptype=='both')):
        # Set up the polar stereographic map
        fig2 = plt.figure()
        #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
        global lm
        if(minlat>10):
            lm = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
        else:
            lm = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
        lm.drawcoastlines()
        lm.drawparallels(np.arange(-80.,81.,20.))
        lm.drawmeridians(np.arange(-180.,181.,20.))    
        
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
        v_min=-15.
        v_max = 15.0

        if(relative_trend is True):
            v_min = -(1.0)*100
            v_max = (1.0)*100
        
        # Center the colorbar on zero
        norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
        #norm = Normalize(vmin=vmin,vmax=vmax)
        mapper = ScalarMappable(norm=norm,cmap=cmap)

        # Set up initial values for the analysis regions
        lwf_canada_avg = 0.
        lwf_canada_counts = 0
        lwf_siberia_avg = 0.
        lwf_siberia_counts = 0 
        lwf_eus_avg = 0.
        lwf_eus_counts = 0
        lwf_wus_avg = 0.
        lwf_wus_counts = 0
        lwf_ea_avg = 0.
        lwf_ea_counts = 0
        lwf_wa_avg = 0.
        lwf_wa_counts = 0
        lwf_europe_avg = 0.
        lwf_europe_counts = 0
        
        # Loop over all the keys and print(the regression slopes 
        # Grab the averages for the key
        max_slope = -10.
        min_slope = 10.
        for i in range(0,len(lat_ranges)-1):
            for j in range(0,len(lon_ranges)-1):
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                #if(dictkey=='48x-97'):
                #    #print("sorted(CERES_dict[dictkey].keys())=",sorted(CERES_dict[dictkey].keys()))
                #    min_date = sorted(CERES_dict[dictkey].keys())[0].decode("utf-8")
                #    max_date = sorted(CERES_dict[dictkey].keys())[-1].decode("utf-8")
        
                # If no data are present for the curent lat/lon box, fill it with
                # black
                if(dictkey not in CERES_dict.keys()):
                    color=(0,0,0,0)
                # If, for some reason, the key is made for the current lat/lon box
                # but the dictionary is empty. fill with black.
                elif(len(CERES_dict[dictkey].keys())==0):
                    color=(0,0,0,0)
                else:
                    # Deal with short wave fluxes first
                    lwf_avgs = [CERES_dict[dictkey][date]['lwf'] for date in \
                        sorted(CERES_dict[dictkey].keys())]
                    if(len(lwf_avgs)<2):
                        color=(0,0,0,0)
                    else:
                        # Check the current max and min
                        x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                        # DEAL WITH MISSING VALUES
                        # Find the slope of the line of best fit for the time series of
                        # average data
                        lwf_avgs = np.ma.masked_array([CERES_dict[dictkey][date]['lwf'] for date in sorted(CERES_dict[dictkey].keys())])
                        #avgs = avgs[np.where(avgs.mask==False)[0]]
                        temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
                        dates = temp_dates


                        if(nozero == True):
                            # Remove zeroes
                            zero_locs = np.where(avgs!=0)
                            lwf_avgs = lwf_avgs[zero_locs]
                            dates = dates[zero_locs]

                        if(relative_trend is True): 
                            flux_avg = np.average(lwf_avgs)
                            slope = (slope/flux_avg)*100.
                        color = mapper.to_rgba(slope)

                        #dates = temp_dates[np.where(avgs.mask==False)[0]]
                        x_vals = np.arange(0,len(dates))
        
                        if(trend_type=='standard'): 
                            slope, intercept, r_value, p_value, std_err = \
                                stats.linregress(x_vals,lwf_avgs)
                            slope *= len(x_vals)
                        else:
                            #The slope
                            S=0
                            sm=0
                            nx = len(lwf_avgs)
                            num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
                            Sn=np.zeros(num_d)
                            for si in range(0,nx-1):
                                for sj in range(si+1,nx):
                                    # Find the slope between the two points
                                    Sn[sm] = (lwf_avgs[si]-lwf_avgs[sj])/(si-sj) 
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
                                slope = slope*len(lwf_avgs)
        
                        if(slope>max_slope):
                            max_slope=slope
                        elif(slope<min_slope):
                            min_slope=slope
            
                        color = mapper.to_rgba(slope)

                        # Add value to analysis regions (if possible)
                        if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                            lwf_canada_avg+=slope
                            lwf_canada_counts+=1
                        elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                            lwf_siberia_avg+=slope
                            lwf_siberia_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                            lwf_eus_avg+=slope
                            lwf_eus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                            lwf_wus_avg+=slope
                            lwf_wus_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                            lwf_ea_avg+=slope
                            lwf_ea_counts+=1
                        elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                            lwf_wa_avg+=slope
                            lwf_wa_counts+=1
                        elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                            lwf_europe_avg+=slope
                            lwf_europe_counts+=1
        
                # Find the y coordinates of the current LatxLon box
                y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
                # Find the x coordinates of the current LatxLon box
                x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
                # Convert x and y into map coordinates
                mx, my = lm(x,y)
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

        lwf_canada_avg = lwf_canada_avg/lwf_canada_counts
        lwf_siberia_avg = lwf_siberia_avg/lwf_canada_counts
        lwf_eus_avg = lwf_eus_avg/lwf_eus_counts
        lwf_wus_avg = lwf_wus_avg/lwf_wus_counts
        lwf_ea_avg = lwf_ea_avg/lwf_ea_counts
        lwf_wa_avg = lwf_wa_avg/lwf_wa_counts
        lwf_europe_avg = lwf_europe_avg/lwf_europe_counts
        print("Canada       lwf_avg = ",lwf_canada_avg,"  counts = ",lwf_canada_counts)
        print("Siberia      lwf_avg = ",lwf_siberia_avg,"  counts = ",lwf_siberia_counts)
        print("Eastern US   lwf_avg = ",lwf_eus_avg,"  counts = ",lwf_eus_counts)
        print("Western US   lwf_avg = ",lwf_wus_avg,"  counts = ",lwf_wus_counts)
        print("Eastern Asia lwf_avg = ",lwf_ea_avg,"  counts = ",lwf_ea_counts)
        print("Western Asia lwf_avg = ",lwf_wa_avg,"  counts = ",lwf_wa_counts)
        print("Europe       lwf_avg = ",lwf_europe_avg,"  counts = ",lwf_europe_counts)

        #start_date = min_date
        #end_date = max_date
        title_string = 'Change in CERES Average TOA Long Wave Flux\n'+       \
                       start_date+' to '+end_date
        if(relative_trend is True):
            title_string = 'Relative change in CERES Average TOA Short Wave Flux\n'+       \
                           start_date+' to '+end_date
        #title_string = 'Change in CERES-2 Average Aerosol Optical Depth\n'+        \
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
        cax = fig2.add_axes([0.27,0.1,0.5,0.05])
        cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
        #cb.ax.set_xlabel('Change in average aerosol index')
        plt.xticks(rotation=45,fontsize=6)
        fig2.canvas.mpl_connect('button_press_event',onclick_lwf)
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))
        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_spring.png'
            elif(summer is True):                                      
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_summer.png'
            elif(autumn is True):                                      
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_autumn.png'
            elif(winter is True):                                      
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_winter.png'
            else:                                                      
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            plt.show()

def plotCERES_MK(CERES_dict,ptype='swf',save=False,season='',start_date='200003',end_date='201812',minlat=30.):
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
        for lkey in sorted(CERES_dict.keys()):
            for tkey in sorted(CERES_dict[lkey].keys()):
                if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                    CERES_dict[lkey].pop(tkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            if(lkey=='data_type'):
                continue
            else:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            for tkey in sorted(CERES_dict[lkey].keys()):
                if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                    CERES_dict[lkey].pop(tkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            for tkey in sorted(CERES_dict[lkey].keys()):
                if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                    CERES_dict[lkey].pop(tkey)

    # Find the lowest lat in the file
    #lowest_lat = float(sorted(CERES_dict.keys())[0].split('x')[0])
    #lowest_lat=30
     
   
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

    if((ptype=='swf') | (ptype=='both')):
        # Set up the polar stereographic map
        fig1 = plt.figure()
        ax = plt.subplot(111)
        #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
        global sm
        if(minlat>10):
            sm = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
        else:
            sm = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
        sm.drawcoastlines()
        sm.drawparallels(np.arange(-80.,81.,20.))
        sm.drawmeridians(np.arange(-180.,181.,20.))    

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
                #    #print("sorted(CERES_dict[dictkey].keys())=",sorted(CERES_dict[dictkey].keys()))
                #    min_date = sorted(CERES_dict[dictkey].keys())[0].decode("utf-8")
                #    max_date = sorted(CERES_dict[dictkey].keys())[-1].decode("utf-8")
        
                # If no data are present for the curent lat/lon box, fill it with
                # black
                if(dictkey not in CERES_dict.keys()):
                    #print('Here 1')
                    colorRED = 0.
                    colorGREEN = 0.
                    colorBLUE = 0.
                # If, for some reason, the key is made for the current lat/lon box
                # but the dictionary is empty. fill with black.
                elif(len(CERES_dict[dictkey].keys())==0):
                    #print('Here 2')
                    colorRED = 0.
                    colorGREEN = 0.
                    colorBLUE = 0.
                else:
                    avgs = [CERES_dict[dictkey][date]['swf'] for date in \
                        sorted(CERES_dict[dictkey].keys())]
                    if(len(avgs)<2):
                        #print('Here 3')
                        colorRED = 0.
                        colorGREEN = 0.
                        colorBLUE = 0.
                    else:
                        # Check the current max and min
                        x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                        avgs = np.ma.masked_array([CERES_dict[dictkey][date]['swf'] for date in sorted(CERES_dict[dictkey].keys())])
                        #avgs = avgs[np.where(avgs.mask==False)[0]]
                        temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
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
                mx, my = sm(x,y)
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
        mx, my = sm(x,y)
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

        print(canada_avg,canada_counts)
        canada_avg = canada_avg/canada_counts
        siberia_avg = siberia_avg/canada_counts
        eus_avg = eus_avg/eus_counts
        wus_avg = wus_avg/wus_counts
        ea_avg = ea_avg/ea_counts
        wa_avg = wa_avg/wa_counts
        europe_avg = europe_avg/europe_counts
        print("Canada       avg = ",canada_avg,"  counts = ",canada_counts)
        print("Siberia      avg = ",siberia_avg,"  counts = ",siberia_counts)
        print("Eastern US   avg = ",eus_avg,"  counts = ",eus_counts)
        print("Western US   avg = ",wus_avg,"  counts = ",wus_counts)
        print("Eastern Asia avg = ",ea_avg,"  counts = ",ea_counts)
        print("Western Asia avg = ",wa_avg,"  counts = ",wa_counts)
        print("Europe       avg = ",europe_avg,"  counts = ",europe_counts)

        #start_date = min_date
        #end_date = max_date
        title_string = 'CERES Average TOA Short Wave Flux Trend P-Values\n'+start_date+' to '+end_date
        name_adder = '_aod_noflags'
        if(spring is True):
            title_string = title_string+'\nMarch, April, May'
        elif(summer is True):
            title_string = title_string+'\nJune, July, August'
        elif(autumn is True):
            title_string = title_string+'\nSeptember, October, November'
        elif(winter is True):
            title_string = title_string+'\nDecember, January, February'
        plt.title(title_string,fontsize=8)
        fig1.canvas.mpl_connect('button_press_event',onclick_swf)
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))

        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_pvalues_'+minstring+'to90_spring.png'
            elif(summer is True):                                     
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_pvalues_'+minstring+'to90_summer.png'
            elif(autumn is True):                                      
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_pvalues_'+minstring+'to90_autumn.png'
            elif(winter is True):                                     
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_pvalues_'+minstring+'to90_winter.png'
            else:                                                     
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_pvalues_'+minstring+'to90_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            plt.show()
    if((ptype=='lwf') | (ptype=='both')):
        # Set up the polar stereographic map
        fig2 = plt.figure()
        ax = plt.subplot(111)
        #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
        global lm
        if(minlat>10):
            lm = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
        else:
            lm = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
        lm.drawcoastlines()
        lm.drawparallels(np.arange(-80.,81.,20.))
        lm.drawmeridians(np.arange(-180.,181.,20.))    

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
                #    #print("sorted(CERES_dict[dictkey].keys())=",sorted(CERES_dict[dictkey].keys()))
                #    min_date = sorted(CERES_dict[dictkey].keys())[0].decode("utf-8")
                #    max_date = sorted(CERES_dict[dictkey].keys())[-1].decode("utf-8")
        
                # If no data are present for the curent lat/lon box, fill it with
                # black
                if(dictkey not in CERES_dict.keys()):
                    #print('Here 1')
                    colorRED = 0.
                    colorGREEN = 0.
                    colorBLUE = 0.
                # If, for some reason, the key is made for the current lat/lon box
                # but the dictionary is empty. fill with black.
                elif(len(CERES_dict[dictkey].keys())==0):
                    #print('Here 2')
                    colorRED = 0.
                    colorGREEN = 0.
                    colorBLUE = 0.
                else:
                    avgs = [CERES_dict[dictkey][date]['lwf'] for date in \
                        sorted(CERES_dict[dictkey].keys())]
                    if(len(avgs)<2):
                        #print('Here 3')
                        colorRED = 0.
                        colorGREEN = 0.
                        colorBLUE = 0.
                    else:
                        # Check the current max and min
                        x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                        avgs = np.ma.masked_array([CERES_dict[dictkey][date]['lwf'] for date in sorted(CERES_dict[dictkey].keys())])
                        #avgs = avgs[np.where(avgs.mask==False)[0]]
                        temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
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
                mx, my = lm(x,y)
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
        mx, my = lm(x,y)
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

        print(canada_avg,canada_counts)
        canada_avg = canada_avg/canada_counts
        siberia_avg = siberia_avg/canada_counts
        eus_avg = eus_avg/eus_counts
        wus_avg = wus_avg/wus_counts
        ea_avg = ea_avg/ea_counts
        wa_avg = wa_avg/wa_counts
        europe_avg = europe_avg/europe_counts
        print("Canada       lwf avg = ",canada_avg,"  counts = ",canada_counts)
        print("Siberia      lwf avg = ",siberia_avg,"  counts = ",siberia_counts)
        print("Eastern US   lwf avg = ",eus_avg,"  counts = ",eus_counts)
        print("Western US   lwf avg = ",wus_avg,"  counts = ",wus_counts)
        print("Eastern Asia lwf avg = ",ea_avg,"  counts = ",ea_counts)
        print("Western Asia lwf avg = ",wa_avg,"  counts = ",wa_counts)
        print("Europe       lwf avg = ",europe_avg,"  counts = ",europe_counts)

        #start_date = min_date
        #end_date = max_date
        title_string = 'CERES Average TOA Long Wave Flux Trend P-Values\n'+start_date+' to '+end_date
        if(spring is True):
            title_string = title_string+'\nMarch, April, May'
        elif(summer is True):
            title_string = title_string+'\nJune, July, August'
        elif(autumn is True):
            title_string = title_string+'\nSeptember, October, November'
        elif(winter is True):
            title_string = title_string+'\nDecember, January, February'
        plt.title(title_string,fontsize=8)
        fig1.canvas.mpl_connect('button_press_event',onclick_lwf)
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))

        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_pvalues_'+minstring+'to90_spring.png'
            elif(summer is True):                                      
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_pvalues_'+minstring+'to90_summer.png'
            elif(autumn is True):                                      
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_pvalues_'+minstring+'to90_autumn.png'
            elif(winter is True):                                      
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_pvalues_'+minstring+'to90_winter.png'
            else:                                                      
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_pvalues_'+minstring+'to90_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            plt.show()

def plotCERES_Climo(CERES_dict,ptype='swf',save=False,trend_type='standard',season='',\
                    start_date='200101',end_date='201512',minlat=30.,anomaly=False,total_avg_dict=None):
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
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    # Find the lowest lat in the file
    #lowest_lat = float(sorted(CERES_dict.keys())[0].split('x')[0])
    #lowest_lat=30
    
    if((ptype=='swf') | (ptype=='both')):
     
        # Set up the polar stereographic map
        fig1 = plt.figure()
        #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
        global sm
        if(minlat>10):
            sm = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
        else:
            sm = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')

        sm.drawcoastlines()
        sm.drawparallels(np.arange(-80.,81.,20.))
        sm.drawmeridians(np.arange(-180.,181.,20.))    
        
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
        v_min = 50.00  # Whole-year values
        mid_val = 0
        v_max = 200.
        
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
        max_avg = -1000.
        min_avg = 1000.
        for i in range(0,len(lat_ranges)-1):
            for j in range(0,len(lon_ranges)-1):
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                #if(dictkey=='48x-97'):
                #    #print("sorted(CERES_dict[dictkey].keys())=",sorted(CERES_dict[dictkey].keys()))
                #    min_date = sorted(CERES_dict[dictkey].keys())[0].decode("utf-8")
                #    max_date = sorted(CERES_dict[dictkey].keys())[-1].decode("utf-8")
        
                # If no data are present for the curent lat/lon box, fill it with
                # black
                if(dictkey not in CERES_dict.keys()):
                    color=(0,0,0,0)
                # If, for some reason, the key is made for the current lat/lon box
                # but the dictionary is empty. fill with black.
                elif(len(CERES_dict[dictkey].keys())==0):
                    color=(0,0,0,0)
                else:
                    avgs = [CERES_dict[dictkey][date]['swf'] for date in \
                        sorted(CERES_dict[dictkey].keys())]
                    # Check the current max and min
                    #x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                    # DEAL WITH MISSING VALUES
                    # Find the slope of the line of best fit for the time series of
                    # average data
                    avgs = np.ma.masked_array([CERES_dict[dictkey][date]['swf'] for date in sorted(CERES_dict[dictkey].keys())])
                    counts = np.ma.masked_array([CERES_dict[dictkey][date]['swf_cc'] for date in sorted(CERES_dict[dictkey].keys())])
                    avg_aot = sum(avgs*counts)/sum(counts)
                    #avg_aot = np.average(avgs)
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
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
                mx, my = sm(x,y)
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

        canada_avg = canada_avg/canada_counts
        siberia_avg = siberia_avg/canada_counts
        eus_avg = eus_avg/eus_counts
        wus_avg = wus_avg/wus_counts
        ea_avg = ea_avg/ea_counts
        wa_avg = wa_avg/wa_counts
        europe_avg = europe_avg/europe_counts
        print("Canada       swf avg = ",canada_avg,"  counts = ",canada_counts)
        print("Siberia      swf avg = ",siberia_avg,"  counts = ",siberia_counts)
        print("Eastern US   swf avg = ",eus_avg,"  counts = ",eus_counts)
        print("Western US   swf avg = ",wus_avg,"  counts = ",wus_counts)
        print("Eastern Asia swf avg = ",ea_avg,"  counts = ",ea_counts)
        print("Western Asia swf avg = ",wa_avg,"  counts = ",wa_counts)
        print("Europe       swf avg = ",europe_avg,"  counts = ",europe_counts)

        #start_date = min_date
        #end_date = max_date
        title_string = 'CERES Average TOA Short Wave Flux Climatology\n'+start_date+' to '+end_date
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
        plt.xlabel('TOA Short Wave Flux',fontsize=6)
        fig1.canvas.mpl_connect('button_press_event',onclick_climo_swf)
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))
        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_spring.png'
            elif(summer is True):                                                                                
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_summer.png'
            elif(autumn is True):                                                                                
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_autumn.png'
            elif(winter is True):                                                                                
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_winter.png'
            else:                                                                                               
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            if(ptype!='both'):
                plt.show()
    if((ptype=='lwf') | (ptype=='both')):
        # Set up the polar stereographic map
        fig2 = plt.figure()
        #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
        global lm
        if(minlat>10):
            lm = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
        else:
            lm = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
        lm.drawcoastlines()
        lm.drawparallels(np.arange(-80.,81.,20.))
        lm.drawmeridians(np.arange(-180.,181.,20.))    
        
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
        v_min = 150.0  # Whole-year values
        mid_val = 0
        v_max = 300.
        
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
        max_avg = -1000.
        min_avg = 1000.
        for i in range(0,len(lat_ranges)-1):
            for j in range(0,len(lon_ranges)-1):
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                #if(dictkey=='48x-97'):
                #    #print("sorted(CERES_dict[dictkey].keys())=",sorted(CERES_dict[dictkey].keys()))
                #    min_date = sorted(CERES_dict[dictkey].keys())[0].decode("utf-8")
                #    max_date = sorted(CERES_dict[dictkey].keys())[-1].decode("utf-8")
        
                # If no data are present for the curent lat/lon box, fill it with
                # black
                if(dictkey not in CERES_dict.keys()):
                    color=(0,0,0,0)
                # If, for some reason, the key is made for the current lat/lon box
                # but the dictionary is empty. fill with black.
                elif(len(CERES_dict[dictkey].keys())==0):
                    color=(0,0,0,0)
                else:
                    avgs = [CERES_dict[dictkey][date]['lwf'] for date in \
                        sorted(CERES_dict[dictkey].keys())]
                    # Check the current max and min
                    #x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                    # DEAL WITH MISSING VALUES
                    # Find the slope of the line of best fit for the time series of
                    # average data
                    avgs = np.ma.masked_array([CERES_dict[dictkey][date]['lwf'] for date in sorted(CERES_dict[dictkey].keys())])
                    counts = np.ma.masked_array([CERES_dict[dictkey][date]['lwf_cc'] for date in sorted(CERES_dict[dictkey].keys())])
                    avg_aot = sum(avgs*counts)/sum(counts) 
                    #avg_aot = np.average(avgs)
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
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
                mx, my = lm(x,y)
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

        canada_avg = canada_avg/canada_counts
        siberia_avg = siberia_avg/canada_counts
        eus_avg = eus_avg/eus_counts
        wus_avg = wus_avg/wus_counts
        ea_avg = ea_avg/ea_counts
        wa_avg = wa_avg/wa_counts
        europe_avg = europe_avg/europe_counts
        print("Canada       lwf avg = ",canada_avg,"  counts = ",canada_counts)
        print("Siberia      lwf avg = ",siberia_avg,"  counts = ",siberia_counts)
        print("Eastern US   lwf avg = ",eus_avg,"  counts = ",eus_counts)
        print("Western US   lwf avg = ",wus_avg,"  counts = ",wus_counts)
        print("Eastern Asia lwf avg = ",ea_avg,"  counts = ",ea_counts)
        print("Western Asia lwf avg = ",wa_avg,"  counts = ",wa_counts)
        print("Europe       lwf avg = ",europe_avg,"  counts = ",europe_counts)

        #start_date = min_date
        #end_date = max_date
        title_string = 'CERES Average TOA Long Wave Flux Climatology\n'+start_date+' to '+end_date
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
        cax = fig2.add_axes([0.27,0.1,0.5,0.05])
        cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
        #cb.ax.set_xlabel('Change in average aerosol index')
        plt.xticks(rotation=45,fontsize=6)
        plt.xlabel('TOA Long Wave Flux',fontsize=6)
        fig2.canvas.mpl_connect('button_press_event',onclick_climo_lwf)
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))
        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_spring.png'
            elif(summer is True):                                    
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_summer.png'
            elif(autumn is True):                                   
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_autumn.png'
            elif(winter is True):                                  
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_winter.png'
            else:                                                  
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            plt.show()
    if(ptype=='alb'):

        # Terra 2003 - 2015 April average
        average_val = 0.578

        avg_dict = {}

        # Set up the polar stereographic map
        if(anomaly==True):
            fig2 = plt.figure(figsize=(9,9))
        else:
            fig2 = plt.figure(figsize=(8,8))
        fig = plt.gcf()
        #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
        global am
        if(minlat>10):
            am = Basemap(projection='npstere',boundinglat=minlat,lon_0=0,resolution='l')
        else:
            am = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
        am.drawcoastlines()
        am.drawparallels(np.arange(-80.,81.,20.))
        am.drawmeridians(np.arange(-180.,181.,20.))    
        
        cmap = plt.get_cmap('jet')
        if(anomaly==True):
            cmap = plt.get_cmap('bwr')
            
        #cmap = plt.get_cmap('nipy_spectral')
        #bcmap = plt.cm.set_cmap(cmap)
        v_min = 0.0  # Whole-year values
        mid_val = 0
        v_max = 1.
        if(anomaly==True):
            v_min = -0.2
            mid_val = 0
            v_max = 0.2
        
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
        bering_avg = 0.
        bering_counts = 0
       
        total_average = 0
        total_counts  = 0
        
        # Loop over all the keys and print(the regression slopes 
        # Grab the averages for the key
        max_avg = -1000.
        min_avg = 1000.
        for i in range(0,len(lat_ranges)-1):
            for j in range(0,len(lon_ranges)-1):
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                #if(dictkey=='48x-97'):
                #    #print("sorted(CERES_dict[dictkey].keys())=",sorted(CERES_dict[dictkey].keys()))
                #    min_date = sorted(CERES_dict[dictkey].keys())[0].decode("utf-8")
                #    max_date = sorted(CERES_dict[dictkey].keys())[-1].decode("utf-8")
        
                # If no data are present for the curent lat/lon box, fill it with
                # black
                if(dictkey not in CERES_dict.keys()):
                    color=(0,0,0,0)
                    avg_dict[dictkey]=-999.
                # If, for some reason, the key is made for the current lat/lon box
                # but the dictionary is empty. fill with black.
                elif(len(CERES_dict[dictkey].keys())==0):
                    color=(0,0,0,0)
                    avg_dict[dictkey]=-999.
                else:
                    # Check the current max and min
                    #x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                    # DEAL WITH MISSING VALUES
                    # Find the slope of the line of best fit for the time series of
                    # average data
                    val_key = 'alb'
                    if(CERES_dict['dtype']=='EBAF_gridded'):
                        val_key='tot_alb'
                        #val_key='clr_alb'
                    avgs = np.ma.masked_array([CERES_dict[dictkey][date][val_key] for date in sorted(CERES_dict[dictkey].keys())])
                    counts = np.ma.masked_array([CERES_dict[dictkey][date][val_key+'_cc'] for date in sorted(CERES_dict[dictkey].keys())])
                    avg_aot = sum(avgs*counts)/sum(counts) 
                    #avg_aot = np.average(avgs)
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
                    #dates = temp_dates
                    #dates = temp_dates[np.where(avgs.mask==False)[0]]
                    #x_vals = np.arange(0,len(dates))
        
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,avgs)
                    #slope *= len(x_vals)
       
                    avg_dict[dictkey]=avg_aot
                    # Account for anomaly, if needed
                    if(anomaly==True):
                        avg_aot = avg_aot-total_avg_dict[dictkey]
        
                    if(avg_aot>max_avg):
                        max_avg=avg_aot
                    elif(avg_aot<min_avg):
                        min_avg=avg_aot
                    if(avg_aot>v_max):
                        avg_aot=v_max
           
                    #                    ceres_lwf[lat_index][lon_index] = (float)((ceres_lwf[lat_index][lon_index]*
                    #                        (float)ceres_lwf_cc[lat_index][lon_index])+lwf[j])/(float)(ceres_lwf_cc[lat_index][lon_index]+1);
                    #                    //avg_lwf = ((avg_lwf*lwf_count)+lwf[j])/(lwf_count+1);
                    #                    ceres_lwf_cc[lat_index][lon_index]++;

                    if(total_counts==0):
                        total_average = avg_aot
                    else:
                        total_average = (float)(((total_average*total_counts)+avg_aot)/(total_counts+1))
                        total_counts+=1.
            
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
                    elif(((lat_ranges[i] >= 68) & (lat_ranges[i]<73)) & ((lon_ranges[j]>=-179) & (lon_ranges[j]<-165))):
                        bering_avg+=avg_aot
                        bering_counts+=1
        
                # Find the y coordinates of the current LatxLon box
                y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
                # Find the x coordinates of the current LatxLon box
                x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
                # Convert x and y into map coordinates
                mx, my = am(x,y)
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

        print("")
        print("Total average in the range: ",total_average)

        #canada_avg = canada_avg/canada_counts
        #siberia_avg = siberia_avg/canada_counts
        #eus_avg = eus_avg/eus_counts
        #wus_avg = wus_avg/wus_counts
        #ea_avg = ea_avg/ea_counts
        #wa_avg = wa_avg/wa_counts
        #europe_avg = europe_avg/europe_counts
        bering_avg = bering_avg/bering_counts
        #print("Canada       alb avg = ",canada_avg,"  counts = ",canada_counts)
        #print("Siberia      alb avg = ",siberia_avg,"  counts = ",siberia_counts)
        #print("Eastern US   alb avg = ",eus_avg,"  counts = ",eus_counts)
        #print("Western US   alb avg = ",wus_avg,"  counts = ",wus_counts)
        #print("Eastern Asia alb avg = ",ea_avg,"  counts = ",ea_counts)
        #print("Western Asia alb avg = ",wa_avg,"  counts = ",wa_counts)
        #print("Europe       alb avg = ",europe_avg,"  counts = ",europe_counts)
        print("Bering       alb avg = ",bering_avg,"  counts = ",bering_counts)

        #start_date = min_date
        #end_date = max_date
        data_type = 'Broadband Surface Albedo'
        if(CERES_dict['dtype']=='EBAF_gridded'):
            data_type = 'Calculated Albedo'
        if(CERES_dict['type_flag'][:4]=='_clr'):
            data_type = 'Clear-sky Cloud-free '+data_type
        elif(CERES_dict['type_flag'][:4]=='_tot'):
            data_type = 'Clear-sky Total '+data_type
        title_string = 'CERES '+CERES_dict['dtype']+' Average April '+data_type+' Climatology\n'+start_date
        if(anomaly==True):
            title_string = 'CERES '+CERES_dict['dtype']+' Average April '+data_type+' Anomaly\n'+start_date
        if(start_date!=end_date):
            title_string = title_string+' - '+end_date

        if(spring is True):
            title_string = title_string+'\nMarch, April, May'
        elif(summer is True):
            title_string = title_string+'\nJune, July, August'
        elif(autumn is True):
            title_string = title_string+'\nSeptember, October, November'
        elif(winter is True):
            title_string = title_string+'\nDecember, January, February'
        plt.title(title_string)
        # Set up the colorbar
        cax = fig.add_axes([0.16,0.075,0.7,0.025])
        cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
        #cb.ax.set_xlabel('Change in average aerosol index')
        plt.xticks(rotation=45)
        plt.xlabel('Broadband Surface Albedo')
        if(anomaly==True):
            plt.xlabel('Broadband Surface Albedo Anomaly')
        fig2.canvas.mpl_connect('button_press_event',onclick_climo_lwf)
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))
        image_type = 'climo'
        if(anomaly==True):
            image_type = 'anomaly'
        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_alb_'+image_type+'_'+start_date+'_'+end_date+'_'+minstring+'to90_spring.png'
            elif(summer is True):                                    
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_alb_'+image_type+'_'+start_date+'_'+end_date+'_'+minstring+'to90_summer.png'
            elif(autumn is True):                                   
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_alb_'+image_type+'_'+start_date+'_'+end_date+'_'+minstring+'to90_autumn.png'
            elif(winter is True):                                  
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_alb_'+image_type+'_'+start_date+'_'+end_date+'_'+minstring+'to90_winter.png'
            else:                                                  
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_alb_'+image_type+'_'+start_date+'_'+end_date+'_'+minstring+'to90_whole_year.png'
            
            plt.savefig(filename)
            print("Saved image:",filename)
        else:
            plt.show()
        return avg_dict

def plotCERES_GlobalTrend(CERES_dict,save=False,ptype='swf',season='',trend_type='standard',start_date='200101',end_date='201512',minlat=-89.):
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
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    loop_keys = list(CERES_dict['31x0'].keys())
    
    if((ptype=='swf') | (ptype=='both')):
        
        mnth_avgs = np.zeros(len(loop_keys))
        mnth_cnts = np.zeros(len(loop_keys))
        for lki in range(0,len(loop_keys)):
            # Calculate the global average for each month
            for i in range(0,len(lat_ranges)-1):
                for j in range(0,len(lon_ranges)-1):
                    dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                    if(dictkey in CERES_dict):
                        if(CERES_dict[dictkey][loop_keys[lki]]['swf']!=-999.0):
                            mnth_avgs[lki] = ((mnth_avgs[lki]*mnth_cnts[lki])+CERES_dict[dictkey][loop_keys[lki]]['swf'])/(mnth_cnts[lki]+1)
                            mnth_cnts[lki]+=1
           
        x_vals = np.arange(0,len(loop_keys)) 
        z = np.polyfit(x_vals,mnth_avgs,1)
        p = np.poly1d(z)
        #if(trend_type=='standard'): 
        #    slope, intercept, r_value, p_value, std_err = \
        #        stats.linregress(x_vals,mnth_avgs)
        #    slope *= len(x_vals)
        #else:
        #    #The slope
        #    S=0
        #    sm=0
        #    nx = len(mnth_avgs)
        #    num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
        #    Sn=np.zeros(num_d)
        #    for si in range(0,nx-1):
        #        for sj in range(si+1,nx):
        #            # Find the slope between the two points
        #            Sn[sm] = (mnth_avgs[si]-mnth_avgs[sj])/(si-sj) 
        #            sm=sm+1
        #        # Endfor
        #    # Endfor
        #    Snsorted=sorted(Sn)
        #    sm=int(num_d/2.)
        #    print(dictkey,len(Snsorted))
        #    if(len(Snsorted)==1):
        #        color=(0,0,0,0) 
        #    else:
        #        if(2*sm    == num_d):
        #            slope=0.5*(Snsorted[sm]+Snsorted[sm+1])
        #        if((2*sm)+1 == num_d): 
        #            slope=Snsorted[sm+1]
        #        slope = slope*len(swf_avgs)
            
                
        #regress_y = x_vals*slope+intercept

        global_trend = p(x_vals[-1])-p(x_vals[0])
        print("Global SWF trend: ",global_trend)

        fig1 = plt.figure()
        plt.plot(x_vals,mnth_avgs)
        plt.plot(x_vals,p(x_vals),'--')

        title_string = 'Global CERES '+CERES_dict['dtype']+' SWF Averages and Trend\n'+       \
                       start_date+' to '+end_date
        #title_string = 'Change in CERES-2 Average Aerosol Optical Depth\n'+        \
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
        # Convert the dates to datetime objects
        ddates = [datetime.strptime(date,'%Y%m') for date in loop_keys] 
        label_dates = [ddate.strftime('%b %Y') for ddate in ddates]
        tickNumber=12
        plt.xticks(np.arange(0,len(x_vals))[::-int(len(x_vals)/tickNumber)],label_dates[::-int(len(x_vals)/tickNumber)],rotation=45,fontsize=6)
        plt.ylabel('SWF Monthly Average (W/m2)')
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))
        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_spring.png'
            elif(summer is True):                                   
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_summer.png'
            elif(autumn is True):                                  
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_autumn.png'
            elif(winter is True):                                 
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_winter.png'
            else:                                                
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_swf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            if(ptype!='both'):
                plt.show()
    if((ptype=='lwf') | (ptype=='both')):
        mnth_avgs = np.zeros(len(loop_keys))
        mnth_cnts = np.zeros(len(loop_keys))
        for lki in range(0,len(loop_keys)):
            # Calculate the global average for each month
            for i in range(0,len(lat_ranges)-1):
                for j in range(0,len(lon_ranges)-1):
                    dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                    if(dictkey in CERES_dict):
                        if(CERES_dict[dictkey][loop_keys[lki]]['lwf']!=-999.0):
                            mnth_avgs[lki] = ((mnth_avgs[lki]*mnth_cnts[lki])+CERES_dict[dictkey][loop_keys[lki]]['lwf'])/(mnth_cnts[lki]+1)
                            mnth_cnts[lki]+=1
           
        x_vals = np.arange(0,len(loop_keys)) 
        #if(trend_type=='standard'): 
        #    slope, intercept, r_value, p_value, std_err = \
        #        stats.linregress(x_vals,mnth_avgs)
        #    slope *= len(x_vals)
        #else:
        #    #The slope
        #    S=0
        #    sm=0
        #    nx = len(mnth_avgs)
        #    num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
        #    Sn=np.zeros(num_d)
        #    for si in range(0,nx-1):
        #        for sj in range(si+1,nx):
        #            # Find the slope between the two points
        #            Sn[sm] = (mnth_avgs[si]-mnth_avgs[sj])/(si-sj) 
        #            sm=sm+1
        #        # Endfor
        #    # Endfor
        #    Snsorted=sorted(Sn)
        #    sm=int(num_d/2.)
        #    print(dictkey,len(Snsorted))
        #    if(len(Snsorted)==1):
        #        color=(0,0,0,0) 
        #    else:
        #        if(2*sm    == num_d):
        #            slope=0.5*(Snsorted[sm]+Snsorted[sm+1])
        #        if((2*sm)+1 == num_d): 
        #            slope=Snsorted[sm+1]
        #        slope = slope*len(swf_avgs)
            
                
        #regress_y = x_vals*slope+intercept
        z = np.polyfit(x_vals,mnth_avgs,1)
        p = np.poly1d(z)
        global_trend = p(x_vals[-1])-p(x_vals[0])
        print("Global LWF trend: ",global_trend)

        fig1 = plt.figure()
        plt.plot(x_vals,mnth_avgs)
        plt.plot(x_vals,p(x_vals),'--')

        #start_date = min_date
        #end_date = max_date
        title_string = 'Global CERES '+CERES_dict['dtype']+' LWF Averages and Trend\n'+       \
                       start_date+' to '+end_date
        #title_string = 'Change in CERES-2 Average Aerosol Optical Depth\n'+        \
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
       # plt.xticks(rotation=45,fontsize=6)

        # Convert the dates to datetime objects
        ddates = [datetime.strptime(date,'%Y%m') for date in loop_keys] 
        label_dates = [ddate.strftime('%b %Y') for ddate in ddates]
        tickNumber=12
        plt.xticks(np.arange(0,len(x_vals))[::-int(len(x_vals)/tickNumber)],label_dates[::-int(len(x_vals)/tickNumber)],rotation=45,fontsize=6)
        plt.ylabel('LWF Monthly Average (W/m2)')

        tickNumber = 12
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))
        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_spring.png'
            elif(summer is True):                                   
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_summer.png'
            elif(autumn is True):                                  
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_autumn.png'
            elif(winter is True):                                 
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_winter.png'
            else:                                                
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_global_lwf_'+trend_label+'_trends_'+start_date+'_'+end_date+'_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            plt.show()

def plotCERES_LatTrend(CERES_dict,save=False,ptype='swf',season='',trend_type='standard',start_date='200101',end_date='201512',minlat=30.,check_lat=85.,nozero=False):
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
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    loop_keys = list(CERES_dict['31x0'].keys())
   
    check_lats = np.array([])
   
    lat_swf_stdev = np.zeros(len(lat_ranges))
    lat_lwf_stdev = np.zeros(len(lat_ranges))
    
    lat_swf_avgs = np.zeros(len(lat_ranges))
    lat_lwf_avgs = np.zeros(len(lat_ranges))
    lat_swf_cnts = np.zeros(len(lat_ranges))
    lat_lwf_cnts = np.zeros(len(lat_ranges))
    final_xvals = np.arange(0,len(lat_ranges))
    for i in range(0,len(lat_ranges)):
        tmp_lwf_stdev = []
        tmp_swf_stdev = []

        lat_lwf_trends = 0
        lat_lwf_counts = 0
        lat_swf_trends = 0
        lat_swf_counts = 0
        print(lat_ranges[i])
        for j in range(0,len(lon_ranges)):
            dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
            if(dictkey in CERES_dict):
                # Deal with short wave fluxes first
                swf_avgs = [CERES_dict[dictkey][date]['swf'] for date in \
                    sorted(CERES_dict[dictkey].keys())]
                lwf_avgs = [CERES_dict[dictkey][date]['lwf'] for date in \
                    sorted(CERES_dict[dictkey].keys())]
                # Check the current max and min
                x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                # DEAL WITH MISSING VALUES
                # Find the slope of the line of best fit for the time series of
                # average data
                swf_avgs = np.ma.masked_array([CERES_dict[dictkey][date]['swf'] for date in sorted(CERES_dict[dictkey].keys())])
                lwf_avgs = np.ma.masked_array([CERES_dict[dictkey][date]['lwf'] for date in sorted(CERES_dict[dictkey].keys())])


                #avgs = avgs[np.where(avgs.mask==False)[0]]
                temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
                swf_dates = temp_dates
                lwf_dates = temp_dates

                if(nozero == True):
                    # Remove zeroes
                    swf_zero_locs = np.where(swf_avgs!=0)
                    swf_avgs = swf_avgs[swf_zero_locs]
                    swf_dates = swf_dates[swf_zero_locs]
                    lwf_zero_locs = np.where(lwf_avgs!=0)
                    lwf_avgs = lwf_avgs[lwf_zero_locs]
                    lwf_dates = lwf_dates[lwf_zero_locs]

                #dates = temp_dates[np.where(avgs.mask==False)[0]]
                x_vals = np.arange(0,len(swf_dates))
                sz = np.polyfit(x_vals,swf_avgs,1)
                sp = np.poly1d(sz)
                swf_trend = sp(x_vals[-1])-sp(x_vals[0])
                lat_swf_trends+=swf_trend
                tmp_swf_stdev.append(swf_trend)
                lat_swf_counts+=1
                x_vals = np.arange(0,len(lwf_dates))
                lz = np.polyfit(x_vals,lwf_avgs,1)
                lp = np.poly1d(lz)
                lwf_trend = lp(x_vals[-1])-lp(x_vals[0])
                lat_lwf_trends+=lwf_trend
                tmp_lwf_stdev.append(lwf_trend)
                lat_lwf_counts+=1

                # To test what's going on with the trends, look at a single
                # latitude band
                # NOTE: The check_lats array is used as input for 
                #       plotCERES_AllTrendatLat
                if(lat_ranges[i]==check_lat):
                    check_lats = np.append(check_lats[:],swf_trend)
        # end lon loop        

        
        tmp_swf_stdev = np.array(tmp_swf_stdev)
        tmp_lwf_stdev = np.array(tmp_lwf_stdev)
        lat_swf_stdev[i] = np.std(tmp_swf_stdev)
        lat_lwf_stdev[i] = np.std(tmp_lwf_stdev)

        lat_swf_trend=lat_swf_trends/lat_swf_counts
        lat_lwf_trend=lat_lwf_trends/lat_lwf_counts
        lat_swf_avgs[i]=lat_swf_trend
        lat_lwf_avgs[i]=lat_lwf_trend
       
    # Don't plot zeros
    lwf_zero_locs = np.where(lat_lwf_avgs!=0)
    swf_zero_locs = np.where(lat_swf_avgs!=0)

    cmap = plt.get_cmap("tab10")

    fig1 = plt.figure()
    #plt.errorbar(lat_ranges[lwf_zero_locs],lat_lwf_avgs[lwf_zero_locs],yerr=lat_lwf_stdev[lwf_zero_locs],errorevery=3,label='Long-wave Flux')
    #plt.errorbar(lat_ranges[swf_zero_locs],lat_swf_avgs[swf_zero_locs],yerr=lat_swf_stdev[swf_zero_locs],errorevery=4,label='Short-wave Flux')
    test_swf_plus = lat_swf_avgs+lat_swf_stdev
    test_swf_minus = lat_swf_avgs-lat_swf_stdev
    test_lwf_plus = lat_lwf_avgs+lat_lwf_stdev
    test_lwf_minus = lat_lwf_avgs-lat_lwf_stdev
    plt.plot(lat_ranges[lwf_zero_locs],lat_lwf_avgs[lwf_zero_locs],label='Long-wave Flux')
    plt.plot(lat_ranges[swf_zero_locs],lat_swf_avgs[swf_zero_locs],label='Short-wave Flux')
    plt.plot(lat_ranges[lwf_zero_locs],test_lwf_plus[lwf_zero_locs],'--',color=cmap(0))
    plt.plot(lat_ranges[lwf_zero_locs],test_lwf_minus[lwf_zero_locs],'--',color=cmap(0))
    plt.plot(lat_ranges[swf_zero_locs],test_swf_plus[swf_zero_locs],'--',color=cmap(1))
    plt.plot(lat_ranges[swf_zero_locs],test_swf_minus[swf_zero_locs],'--',color=cmap(1))

    title_string = 'CERES '+CERES_dict['dtype']+' Zonal Average TOA Flux Trends\n'+start_date+' to '+end_date

    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    elif(summer is True):
        title_string = title_string+'\nJune, July, August'
    elif(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    elif(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.xlabel('Latitude')
    plt.ylabel('CERES TOA Flux/Study Period')
    #plt.ylim(-15,15)
    plt.title(title_string)
    plt.legend()
    #plt.show()
    # Convert the dates to datetime objects
    #ddates = [datetime.strptime(date,'%Y%m') for date in loop_keys] 
    #label_dates = [ddate.strftime('%b %Y') for ddate in ddates]
    #tickNumber=12
    #plt.xticks(np.arange(0,len(x_vals))[::-int(len(x_vals)/tickNumber)],label_dates[::-int(len(x_vals)/tickNumber)],rotation=45,fontsize=6)
    #plt.ylabel('SWF Monthly Average (W/m2)')
    minstring = str(int(minlat))
    if(minlat<0):
        minstring = 'N'+str(int(abs(minlat)))
    if(save is True):
        if(spring is True):
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_spring.png'
        elif(summer is True):                                   
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_summer.png'
        elif(autumn is True):                                  
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_autumn.png'
        elif(winter is True):                                 
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_winter.png'
        else:                                               
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_whole_year.png'
        
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()
    return lat_swf_avgs,lat_lwf_avgs
    #return lat_swf_avgs,check_lats,check_lat

def plotCERES_LatClimo(CERES_dict,minlat=31.,ptype='swf',save=False,season='',trend_type='standard',sfc_type='all',start_date='200101',end_date='201512'):
    trend_label=''
    lat_ranges = np.arange(minlat,90,1)
    lon_ranges = np.arange(-179,180,1)
    dates = [2005,2007,2008,2009,2010,2011,2012,2013,2014,2015]
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False
    if(season=='spring'):
        spring = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    # Set up the array to hold the latitudinal averages
    lat_total_avgs = np.zeros(len(lat_ranges))
    lat_land_avgs = np.zeros(len(lat_ranges))
    lat_ocean_avgs = np.zeros(len(lat_ranges))

    # Find the lowest lat in the file
    #lowest_lat = 31.


    # Load in the Surface Type data
    bfile = open('/data/landmask/global_18.map','rb')
    img = bytearray(bfile.read())
    img = np.reshape(img,(1080,2160))
    scl = 0.166667
    st_lat = 89.9
    st_lon=-179.9

    
    # Set up the polar stereographic map
    fig1 = plt.figure()
    
    # Loop over all the keys and print(the regression slopes 
    # Grab the averages for the key
    for i in range(0,len(lat_ranges)):
        lat_total_avg = 0. 
        lat_total_count=0.
        lat_land_avg = 0. 
        lat_land_count=0.
        lat_ocean_avg = 0. 
        lat_ocean_count=0.
        for j in range(0,len(lon_ranges)):
    
            #print(type_num)
            # Do surface type checks
            #if(sfc_type=='land'):
            #    si = int(np.round((st_lat-lat_ranges[i])/scl))
            #    sj = int(np.round((lon_ranges[j]-st_lon)/scl))
            #    type_num = img[si,sj]
            #elif(sfc_type=='ocean'):
            #    si = int(np.round((st_lat-lat_ranges[i])/scl))
            #    sj = int(np.round((lon_ranges[j]-st_lon)/scl))
            #    type_num = img[si,sj]

            #if(sfc_flag == True):
            dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
    
            avgs = np.array([CERES_dict[dictkey][date][ptype] for date in \
                sorted(CERES_dict[dictkey].keys())])
            ##avgs = [CERES_dict[dictkey][date]['cloud'+cloud_flag] for date in \
            ##    sorted(CERES_dict[dictkey].keys())]
            # Check the current max and min
            x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
            # DEAL WITH MISSING VALUES
            temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))

            avgs = np.array([CERES_dict[dictkey][date][ptype] for date \
                   in sorted(CERES_dict[dictkey].keys())])
            cnts = np.array([CERES_dict[dictkey][date][ptype+'_cc'] for date \
                   in sorted(CERES_dict[dictkey].keys())])
            temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
            ##avgs = np.ma.masked_array([CERES_dict[dictkey][date]['cloud'+ \
            ##   cloud_flag] for date in sorted(CERES_dict[dictkey].keys())])
            good_avgs = np.where(np.isnan(avgs) == False)
            avgs = avgs[good_avgs]
            cnts = cnts[good_avgs]
            temp_dates = temp_dates[good_avgs] 
            dates = temp_dates
            if(len(avgs)>0):
                avg_temp = np.average(avgs,weights=cnts)

                si = int(np.round((st_lat-lat_ranges[i])/scl))
                sj = int(np.round((lon_ranges[j]-st_lon)/scl))
                type_num = img[si,sj]
                lat_total_avg+=avg_temp
                lat_total_count+=1 
                if(type_num!=17):
                    lat_land_avg+=avg_temp
                    lat_land_count+=1 
                if(type_num==17):
                    lat_ocean_avg+=avg_temp
                    lat_ocean_count+=1 

        if(lat_total_count==0):
            print(lat_ranges[i],'No counts')
        else:
            lat_total_avg=lat_total_avg/lat_total_count
            lat_total_avgs[i] = lat_total_avg
            print(lat_ranges[i],lat_total_avgs[i])
        if(lat_land_count!=0):
            lat_land_avg=lat_land_avg/lat_land_count
            lat_land_avgs[i] = lat_land_avg
        if(lat_ocean_count!=0):
            lat_ocean_avg=lat_ocean_avg/lat_ocean_count
            lat_ocean_avgs[i] = lat_ocean_avg
  
    # Don't plot zeros
    ocean_zero_locs = np.where(lat_ocean_avgs!=0)
    land_zero_locs = np.where(lat_land_avgs!=0)
    total_zero_locs = np.where(lat_total_avgs!=0)
   
    plt.plot(lat_ranges[ocean_zero_locs],lat_ocean_avgs[ocean_zero_locs],label='Ocean')
    plt.plot(lat_ranges[land_zero_locs], lat_land_avgs[land_zero_locs],label='Land')
    plt.plot(lat_ranges[total_zero_locs],lat_total_avgs[total_zero_locs],label='Total')
    plt.legend()
    plt.title('CERES '+CERES_dict['dtype']+' '+ptype+' Climo: '+\
                sfc_type+'\n'+season+ ' '+start_date+ ' to '+end_date)
    plt.xlabel('Latitude')
    plt.ylabel('CERES '+ptype+' [W/m2]')
    minstring = str(int(minlat))
    if(minlat<0):
        minstring = 'N'+str(int(abs(minlat)))
    if(save == True):
        name_end = '.png'
        if(spring==True):
            name_end = '_spring.png'
        elif(summer==True):
            name_end = '_summer.png'
        elif(autumn==True):
            name_end = '_autumn.png'
        elif(winter==True):
            name_end = '_winter.png'
        outname = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_'+ptype+'_lat_climo_'+start_date+'_'+end_date+'_'+minstring+'to90_'+name_end

        plt.savefig(outname,dpi=300)
        print('Saved image '+outname)
    else:
        plt.show()


def plotCERES_LatTrendTrend(CERES_dict,minlat=31.,ptype='swf',save=False,season='',trend_type='standard',sfc_type='all',start_date='200101',end_date='201512',nozero=False):
    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='thielSen_'
    global ext
    ext=False
    lat_ranges = np.arange(minlat,90,1)
    lon_ranges = np.arange(-179,180,1)

    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False
    if(season=='spring'):
        spring = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    # Set up the array to hold the latitudinal averages
    lat_total_avgs = np.zeros(len(lat_ranges))
    lat_ocean_avgs = np.zeros(len(lat_ranges))
    lat_land_avgs = np.zeros(len(lat_ranges))

    # Find the lowest lat in the file
    #lowest_lat = 31.


    # Load in the Surface Type data
    bfile = open('/data/landmask/global_18.map','rb')
    img = bytearray(bfile.read())
    img = np.reshape(img,(1080,2160))
    scl = 0.166667
    st_lat = 89.9
    st_lon=-179.9

    # NOTE: Assume that all GISTEMP data are being analyzed (2000 to 2018)
    years = np.arange(int(start_date[:4]),2016,1)
    months= ['01','02','03','04','05','06','07','08','09','10','11','12']

    test_2d_trends = np.zeros([len(years),len(lat_ranges)])
    yearly_total_trends = []
    yearly_ocean_trends = []
    yearly_land_trends = []
   
    # For testing purposes, only look at 75 N for now
    #lat_ranges = np.array([75])
    
    all_values = []
   
    #Loop over the times first 
    for iy in range(0,len(years)):
        lat_total_avg = 0. 
        lat_total_count=0.
        lat_ocean_avg = 0. 
        lat_ocean_count=0.
        lat_land_avg = 0. 
        lat_land_count=0.

        # Loop over all the keys and print(the regression slopes 
        # Grab the averages for the key
        for i in range(0,len(lat_ranges)):
            lat_total_avg = 0. 
            lat_total_count=0.
            lat_ocean_avg = 0. 
            lat_ocean_count=0.
            lat_land_avg = 0. 
            lat_land_count=0.

            # Loop over all longitudes and find monthly trends at each 1x1 
            # degree box. After finding the trend for the current month, add
            # the trend to a running total of trends. The running total is
            # divided by the number of valid trends at the end of the loop
            # to find the total zonal trend for this year's averages
            for j in range(0,len(lon_ranges)):
                tavgs = []
        
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
       
                for mnth in months:
                    datekey = str(int(years[iy]))+mnth
                    if(datekey in CERES_dict[dictkey].keys()):
                        if(np.isnan(CERES_dict[dictkey][datekey][ptype])==False):
                            tavgs.append(CERES_dict[dictkey][datekey][ptype])
                            if(j==10):
                                all_values.append(CERES_dict[dictkey][datekey][ptype])

                # Convert the basic list of the monthly zonal averages at the 
                # current latitude into a numpy array.
                avgs = np.array(tavgs)

                #if(j==10):
                #    all_values.append(np.average(avgs))
                #good_avgs = np.where(np.isnan(avgs) == False)
                #avgs = avgs[good_avgs]
                #temp_dates = temp_dates[good_avgs] 

                if(nozero == True):
                    # Remove zeroes
                    zero_locs = np.where(avgs!=0)
                    avgs = avgs[zero_locs]

                x_vals = np.arange(0,len(avgs))

                # Calculate the trend of the monthy zonal averages over the current
                # year.
                if(len(avgs)>0):
                    #dates = temp_dates[np.where(avgs.mask==False)[0]]
                    if(trend_type=='standard'): 
                        #slope, intercept, r_value, p_value, std_err = \
                        #    stats.linregress(x_vals,avgs)
                        #slope *= len(x_vals)
                        #regress_y = x_vals*slope+intercept
                        try:
                            z = np.polyfit(x_vals,avgs,1)
                            p = np.poly1d(z)
                            trend = p(x_vals[-1])-p(x_vals[0])
                   #         if(j==10):
                   #             fig1 = plt.figure()
                   #             plt.plot(x_vals,avgs)
                   #             plt.plot(x_vals,p(x_vals))
                   #             #plt.plot(x_vals,regress_y)
                   #             plt.title(str(int(years[iy]))+mnth)
                   #             plt.show()
                        except:
                            print('Alg Error: ',avgs)
                            trend=0
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

                    # Add the 
                    if(np.isnan(trend) == False):
                        #print(lat_ranges[i],lon_ranges[j],slope)
                        si = int(np.round((st_lat-lat_ranges[i])/scl))
                        sj = int(np.round((lon_ranges[j]-st_lon)/scl))
                        type_num = img[si,sj]
                        # Do surface type checks
                        #if(sfc_type=='land'):
                        #    si = int(np.round((st_lat-lat_ranges[i])/scl))
                        #    sj = int(np.round((lon_ranges[j]-st_lon)/scl))
                        #    type_num = img[si,sj]
                        if(type_num!=17):
                            sfc_flag = True
                            lat_land_avg+=trend
                            lat_land_count+=1.
                        #elif(sfc_type=='ocean'):
                        #    si = int(np.round((st_lat-lat_ranges[i])/scl))
                        #    sj = int(np.round((lon_ranges[j]-st_lon)/scl))
                        #    type_num = img[si,sj]
                        if(type_num==17):
                            sfc_flag = True
                            lat_ocean_avg+=trend
                            lat_ocean_count+=1.
                        ###else:
                        ###    sfc_flag = True           

                        lat_total_avg+=trend
                        lat_total_count+=1.
            if(lat_total_count==0):
                print(lat_ranges[i],'No counts')
            else:
                lat_total_avg=lat_total_avg/lat_total_count
                lat_total_avgs[i] = lat_total_avg
                test_2d_trends[iy,i] = lat_total_avg
                print(years[iy],lat_ranges[i],lat_total_avgs[i])
            if(lat_land_count!=0):
                lat_land_avg=lat_land_avg/lat_land_count
                lat_land_avgs[i] = lat_land_avg
            if(lat_ocean_count!=0):
                lat_ocean_avg=lat_ocean_avg/lat_ocean_count
                lat_ocean_avgs[i] = lat_ocean_avg
        yearly_total_trends.append(lat_total_avg)
        yearly_ocean_trends.append(lat_ocean_avg)
        yearly_land_trends.append(lat_land_avg)
 
    print("Yearly total trends",yearly_total_trends)
    print("Yearly ocean trends",yearly_ocean_trends)
    print("Yearly land trends",yearly_land_trends)
  
    # Find trend of trends
    #slope, intercept, r_value, p_value, std_err = \
    #    stats.linregress(x_vals,avgs)
    #slope *= len(x_vals)
    x_vals = np.arange(len(yearly_total_trends))
    #z = np.polyfit(x_vals,yearly_total_trends,1)
    z = np.polyfit(years,yearly_total_trends,1)
    p = np.poly1d(z)
    trend = p(years[-1])-p(years[0])
    #trend = p(x_vals[-1])-p(x_vals[0])

    plotX,plotY = np.meshgrid(years,lat_ranges) 
    fig2 = plt.figure()
    ax = Axes3D(fig2)
    surf = ax.plot_surface(plotX,plotY,test_2d_trends.transpose(),cmap=cm.jet)
    #ax.set_zlim(-6,6)
    fig2.colorbar(surf,shrink=0.5,aspect=5)
    ax.set_xlim(2000,2015)
    ax.set_xlabel('Year')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('CERES '+ptype+ 'trend [W/m2]')
    title_string = 'CERES Yearly Zonal Average TOA '+ptype+' Trends\n'+start_date+' to '+end_date
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    elif(summer is True):
        title_string = title_string+'\nJune, July, August'
    elif(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    elif(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)

    minstring = str(int(minlat))
    if(minlat<0):
        minstring = 'N'+str(int(abs(minlat)))
    name_end = '.png'
    if(spring==True):
        name_end = '_spring.png'
    elif(summer==True):
        name_end = '_summer.png'
    elif(autumn==True):
        name_end = '_autumn.png'
    elif(winter==True):
        name_end = '_winter.png'
    outname = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_3d_'+ptype+'_trend_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_'+name_end
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image ",outname)
    else:
        plt.show()
    return years,lat_ranges,test_2d_trends
  
    ###fig0 = plt.figure()
    ###newvals = np.arange(0,len(all_values))
    ###plt.plot(newvals,all_values)
    ###plt.title('75N')
   
    ####start_date = str(temp_dates[0])
    ####end_date = str(temp_dates[-1])
    ###fig1 = plt.figure()
    ###plt.plot(years,yearly_ocean_trends,label='Ocean')
    ###plt.plot(years,yearly_land_trends,label='Land')
    ###plt.plot(years,yearly_total_trends,label='All')
    ###plt.plot(years,p(years),'--',label='All Trend')
    ###plt.legend()
    ###plt.title('GISS Surface CERESerature Anomaly \n1200 km Smoothing Trend: '+\
    ###            sfc_type+'\n'+season+ ' '+start_date+ ' to '+end_date)
    ###plt.xlabel('Latitude')
    ###plt.ylabel('CERESerature Anomaly Trend Trend [K/year]')
    ###name_end = '.png'
    ###if(save == True):
    ###    if(spring==True):
    ###        name_end = '_spring.png'
    ###    elif(summer==True):
    ###        name_end = '_summer.png'
    ###    elif(autumn==True):
    ###        name_end = '_autumn.png'
    ###    elif(winter==True):
    ###        name_end = '_winter.png'
    ###    if(minlat>=0):
    ###        outname = 'gisstemp_lat_trend_'+start_date+"_"+end_date+"_combined_"+str(int(minlat))+'to90'+name_end
    ###    else:
    ###        outname = 'gisstemp_lat_trend_'+start_date+"_"+end_date+"_combined_N"+str(abs(int(minlat)))+'to90'+name_end
    ###    plt.savefig(outname,dpi=300)
    ###    print('Saved image '+outname)
    ###else:
    ###    plt.show()

    ###return yearly_total_trends,yearly_ocean_trends,yearly_land_trends

def plotCERES_LatClimoTrend(CERES_dict,minlat=31.,ptype='swf',save=False,season='',trend_type='standard',sfc_type='all',start_date='200101',end_date='201512'):
    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='thielSen_'
    global ext
    ext=False
    lat_ranges = np.arange(minlat,90,1)
    lon_ranges = np.arange(-179,180,1)

    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False
    if(season=='spring'):
        spring = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    # Set up the array to hold the latitudinal averages
    lat_total_avgs = np.zeros(len(lat_ranges))
    lat_ocean_avgs = np.zeros(len(lat_ranges))
    lat_land_avgs = np.zeros(len(lat_ranges))

    # Find the lowest lat in the file
    #lowest_lat = 31.


    # Load in the Surface Type data
    bfile = open('/data/landmask/global_18.map','rb')
    img = bytearray(bfile.read())
    img = np.reshape(img,(1080,2160))
    scl = 0.166667
    st_lat = 89.9
    st_lon=-179.9

    # NOTE: Assume that all GISTEMP data are being analyzed (2000 to 2018)
    years = np.arange(int(start_date[:4]),2016,1)
    if(CERES_dict['dtype']=='aqua'):
        years = np.arange(2002,2016,1)
    months= ['01','02','03','04','05','06','07','08','09','10','11','12']

    test_2d_avgs = np.zeros([len(years),len(lat_ranges)])
    yearly_total_avgs = []
    yearly_ocean_avgs = []
    yearly_land_avgs = []
   
    # For testing purposes, only look at 75 N for now
    #lat_ranges = np.array([75])
    
    all_values = []
   
    #Loop over the times first 
    for iy in range(0,len(years)):
        lat_total_avg = 0. 
        lat_total_count=0.
        lat_ocean_avg = 0. 
        lat_ocean_count=0.
        lat_land_avg = 0. 
        lat_land_count=0.

        # Loop over all the keys and print(the regression slopes 
        # Grab the averages for the key
        for i in range(0,len(lat_ranges)):
            lat_total_avg = 0. 
            lat_total_count=0.
            lat_ocean_avg = 0. 
            lat_ocean_count=0.
            lat_land_avg = 0. 
            lat_land_count=0.

            # Loop over all longitudes and find monthly trends at each 1x1 
            # degree box. After finding the trend for the current month, add
            # the trend to a running total of trends. The running total is
            # divided by the number of valid trends at the end of the loop
            # to find the total zonal trend for this year's averages
            for j in range(0,len(lon_ranges)):
                tavgs = []
                tcnts = []
        
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
       
                for mnth in months:
                    datekey = str(int(years[iy]))+mnth
                    if(datekey in CERES_dict[dictkey].keys()):
                        if(np.isnan(CERES_dict[dictkey][datekey][ptype])==False):
                            tavgs.append(CERES_dict[dictkey][datekey][ptype])
                            tcnts.append(CERES_dict[dictkey][datekey][ptype+'_cc'])
                            if(j==10):
                                all_values.append(CERES_dict[dictkey][str(int(years[iy]))+mnth][ptype])
                avgs = np.array(tavgs)
                cnts = np.array(tcnts)
                #if(j==10):
                #    all_values.append(np.average(avgs))
                #good_avgs = np.where(np.isnan(avgs) == False)
                #avgs = avgs[good_avgs]
                #temp_dates = temp_dates[good_avgs] 
                x_vals = np.arange(0,len(avgs))

                if(len(avgs)>0):

                    # Calculate the weighted average of the monthly averages
                    # for the current year at the current lat and lon
                    totavg = np.average(avgs,weights=cnts)
                    #dates = temp_dates[np.where(avgs.mask==False)[0]]
                    si = int(np.round((st_lat-lat_ranges[i])/scl))
                    sj = int(np.round((lon_ranges[j]-st_lon)/scl))
                    type_num = img[si,sj]
                    # Do surface type checks
                    #if(sfc_type=='land'):
                    #    si = int(np.round((st_lat-lat_ranges[i])/scl))
                    #    sj = int(np.round((lon_ranges[j]-st_lon)/scl))
                    #    type_num = img[si,sj]
                    if(type_num!=17):
                        sfc_flag = True
                        lat_land_avg+=totavg
                        lat_land_count+=1.
                    #elif(sfc_type=='ocean'):
                    #    si = int(np.round((st_lat-lat_ranges[i])/scl))
                    #    sj = int(np.round((lon_ranges[j]-st_lon)/scl))
                    #    type_num = img[si,sj]
                    if(type_num==17):
                        sfc_flag = True
                        lat_ocean_avg+=totavg
                        lat_ocean_count+=1.
                    ###else:
                    ###    sfc_flag = True           

                    lat_total_avg+=totavg
                    lat_total_count+=1.
            if(lat_total_count==0):
                print(lat_ranges[i],'No counts')
            else:
                lat_total_avg=lat_total_avg/lat_total_count
                lat_total_avgs[i] = lat_total_avg
                test_2d_avgs[iy,i] = lat_total_avg
                print(years[iy],lat_ranges[i],lat_total_avgs[i])
            if(lat_land_count!=0):
                lat_land_avg=lat_land_avg/lat_land_count
                lat_land_avgs[i] = lat_land_avg
            if(lat_ocean_count!=0):
                lat_ocean_avg=lat_ocean_avg/lat_ocean_count
                lat_ocean_avgs[i] = lat_ocean_avg
        yearly_total_avgs.append(lat_total_avg)
        yearly_ocean_avgs.append(lat_ocean_avg)
        yearly_land_avgs.append(lat_land_avg)
 
  #  print("Yearly total avgs",yearly_total_avgs)
  #  print("Yearly ocean avgs",yearly_ocean_avgs)
  #  print("Yearly land avgs",yearly_land_avgs)
  #
    # Find trend of avgs
    #slope, intercept, r_value, p_value, std_err = \
    #    stats.linregress(x_vals,avgs)
    #slope *= len(x_vals)
    x_vals = np.arange(len(yearly_total_avgs))
    #z = np.polyfit(x_vals,yearly_total_avgs,1)
    z = np.polyfit(years,yearly_total_avgs,1)
    p = np.poly1d(z)
    trend = p(years[-1])-p(years[0])
    #trend = p(x_vals[-1])-p(x_vals[0])

    plotX,plotY = np.meshgrid(years,lat_ranges) 
    fig2 = plt.figure()
    ax = Axes3D(fig2)
    surf = ax.plot_surface(plotX,plotY,test_2d_avgs.transpose(),cmap=cm.jet)
    #ax.set_zlim(-6,6)
    fig2.colorbar(surf,shrink=0.5,aspect=5)
    ax.set_xlim(2000,2015)
    ax.set_xlabel('Year')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('CERES '+ptype+' [W/m2]')
    title_string = 'CERES Yearly Zonal Average TOA '+ptype+' Climatologies\n'+start_date+' to '+end_date
    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    elif(summer is True):
        title_string = title_string+'\nJune, July, August'
    elif(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    elif(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.title(title_string,fontsize=8)

    minstring = str(int(minlat))
    if(minlat<0):
        minstring = 'N'+str(int(abs(minlat)))
    name_end = '.png'
    if(spring==True):
        name_end = '_spring.png'
    elif(summer==True):
        name_end = '_summer.png'
    elif(autumn==True):
        name_end = '_autumn.png'
    elif(winter==True):
        name_end = '_winter.png'
    outname = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_3d_'+ptype+'_climo_trends_'+start_date+'_'+end_date+'_'+minstring+'to90_'+name_end
    if(save == True):
        plt.savefig(outname,dpi=300)
    else:
        plt.show()
    #return years,lat_ranges,test_2d_avgs

def plotCERES_LatTrend_MonthBased(CERES_dict,save=False,ptype='swf',season='',trend_type='standard',start_date='200101',end_date='201512',\
                                  minlat=30.,check_lat=85.):
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
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    loop_keys = list(CERES_dict['31x0'].keys())
   
    check_lats = np.array([])
   
    lat_swf_stdev = np.zeros(len(lat_ranges))
    lat_lwf_stdev = np.zeros(len(lat_ranges))
    
    lat_swf_avgs = np.zeros(len(lat_ranges))
    lat_lwf_avgs = np.zeros(len(lat_ranges))
    lat_swf_cnts = np.zeros(len(lat_ranges))
    lat_lwf_cnts = np.zeros(len(lat_ranges))
    final_xvals = np.arange(0,len(lat_ranges))
    for i in range(0,len(lat_ranges)):
        tmp_lwf_mnth_lst = []
        tmp_swf_mnth_lst = []

        lat_lwf_trends = 0
        lat_lwf_counts = 0
        lat_swf_trends = 0
        lat_swf_counts = 0
        print(lat_ranges[i])
        for mnth in range(0,len(loop_keys)):
            # Make a variable to hold the average of all longitude values 
            # for the current month. Loop over longitudes and grab 
            tmp_lwf_holder = []
            tmp_swf_holder = []


            for j in range(0,len(lon_ranges)):
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                if(dictkey in CERES_dict):
                    tmp_swf_holder.append(CERES_dict[dictkey][loop_keys[mnth]]['swf'])
                    tmp_lwf_holder.append(CERES_dict[dictkey][loop_keys[mnth]]['lwf'])
            # end lon loop        
            
            # Convert holder lists to arrays for easy averaging
            tmp_swf_holder = np.array(tmp_swf_holder)
            tmp_lwf_holder = np.array(tmp_lwf_holder)

            tmp_swf_mnth_lst.append(np.average(tmp_swf_holder))
            tmp_lwf_mnth_lst.append(np.average(tmp_lwf_holder))
        # end mnth loop        

        # Convert list of zonal averages for each month into a np array
        swf_avgs = np.array(tmp_swf_mnth_lst)
        lwf_avgs = np.array(tmp_lwf_mnth_lst)

        # Check the current max and min
        x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
        # DEAL WITH MISSING VALUES
        # Find the slope of the line of best fit for the time series of
        # average data
        #swf_avgs = np.ma.masked_array([CERES_dict[dictkey][date]['swf'] for date in sorted(CERES_dict[dictkey].keys())])
        #lwf_avgs = np.ma.masked_array([CERES_dict[dictkey][date]['lwf'] for date in sorted(CERES_dict[dictkey].keys())])

        #avgs = avgs[np.where(avgs.mask==False)[0]]
        temp_dates = loop_keys
        dates = temp_dates
        #dates = temp_dates[np.where(avgs.mask==False)[0]]
        x_vals = np.arange(0,len(dates))
        sz = np.polyfit(x_vals,swf_avgs,1)
        sp = np.poly1d(sz)
        swf_trend = sp(x_vals[-1])-sp(x_vals[0])
        lat_swf_trends+=swf_trend
        #tmp_swf_stdev.append(swf_trend)
        lat_swf_counts+=1
        lz = np.polyfit(x_vals,lwf_avgs,1)
        lp = np.poly1d(lz)
        lwf_trend = lp(x_vals[-1])-lp(x_vals[0])
        lat_lwf_trends+=lwf_trend
        #tmp_lwf_stdev.append(lwf_trend)
        lat_lwf_counts+=1

        # To test what's going on with the trends, look at a single
        # latitude band
        # NOTE: The check_lats array is used as input for 
        #       plotCERES_AllTrendatLat
        if(lat_ranges[i]==check_lat):
            check_lats = np.append(check_lats[:],swf_trend)
        
        #tmp_swf_stdev = np.array(tmp_swf_stdev)
        #tmp_lwf_stdev = np.array(tmp_lwf_stdev)
        #lat_swf_stdev[i] = np.std(tmp_swf_stdev)
        #lat_lwf_stdev[i] = np.std(tmp_lwf_stdev)

        #lat_swf_trend=lat_swf_trends/lat_swf_counts
        #lat_lwf_trend=lat_lwf_trends/lat_lwf_counts
        lat_swf_avgs[i]=swf_trend
        lat_lwf_avgs[i]=lwf_trend
        
    # end lat loop
       
    # Don't plot zeros
    lwf_zero_locs = np.where(lat_lwf_avgs!=0)
    swf_zero_locs = np.where(lat_swf_avgs!=0)

    cmap = plt.get_cmap("tab10")

    fig1 = plt.figure()
    plt.plot(lat_ranges[lwf_zero_locs],lat_lwf_avgs[lwf_zero_locs],label='Long-wave Flux')
    plt.plot(lat_ranges[swf_zero_locs],lat_swf_avgs[swf_zero_locs],label='Short-wave Flux')

    title_string = 'CERES '+CERES_dict['dtype']+' Zonal Average TOA Flux Trends (Month-first)\n'+start_date+' to '+end_date

    if(spring is True):
        title_string = title_string+'\nMarch, April, May'
    elif(summer is True):
        title_string = title_string+'\nJune, July, August'
    elif(autumn is True):
        title_string = title_string+'\nSeptember, October, November'
    elif(winter is True):
        title_string = title_string+'\nDecember, January, February'
    plt.xlabel('Latitude')
    plt.ylabel('CERES TOA Flux/Study Period')
    plt.title(title_string)
    plt.legend()
    #plt.show()
    # Convert the dates to datetime objects
    #ddates = [datetime.strptime(date,'%Y%m') for date in loop_keys] 
    #label_dates = [ddate.strftime('%b %Y') for ddate in ddates]
    #tickNumber=12
    #plt.xticks(np.arange(0,len(x_vals))[::-int(len(x_vals)/tickNumber)],label_dates[::-int(len(x_vals)/tickNumber)],rotation=45,fontsize=6)
    #plt.ylabel('SWF Monthly Average (W/m2)')
    minstring = str(int(minlat))
    if(minlat<0):
        minstring = 'N'+str(int(abs(minlat)))
    if(save is True):
        if(spring is True):
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_month_first_'+start_date+'_'+end_date+'_'+minstring+'to90_spring.png'
        elif(summer is True):                                   
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_month_first_'+start_date+'_'+end_date+'_'+minstring+'to90_summer.png'
        elif(autumn is True):                                  
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_month_first_'+start_date+'_'+end_date+'_'+minstring+'to90_autumn.png'
        elif(winter is True):                                 
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_month_first_'+start_date+'_'+end_date+'_'+minstring+'to90_winter.png'
        else:                                                
            filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lat_trends_month_first_'+start_date+'_'+end_date+'_'+minstring+'to90_whole_year.png'
        
        plt.savefig(filename,dpi=300)
        print("Saved image:",filename)
    else:
        plt.show()
    return lat_swf_avgs,check_lats,check_lat

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#  Miscellaneous Plot Codes
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# This function plots all of the trends at a given latitude used to find
# the average trend at that latitude. This code was written to plot each
# trend used in validating the summer lat trends.
# 
# INPUT:
#   check_trends: This is an array of all 360 trends calculated at each
#                 longitude at 'check_lat'. The average of this array is 
#                 the average trend seen in the lat_trend plots.
#
#   check_lat: The latitude that is being investigated.
def plotCERES_AllTrendatLat(check_trends,check_lat,season=''):
    avg_trend = np.average(check_trends)
    sat = 'Terra'
    sat_name = 'terra'

    dtype = 'SWF'
    dtype_name = 'swf'

    fig1 = plt.figure()
    plt.plot(np.arange(len(check_trends)),check_trends)
    plt.plot(np.arange(len(check_trends)),np.full(len(check_trends),avg_trend))
    plt.title(sat+' CERES SWF '+season+' Trends at '+str(check_lat)+' N\nEach Total Trend at Each Longitude')
    plt.ylabel(dtype+' Trend [W/m2/Study Period]')
    plt.xlabel('Longitude [Degrees East]')
    minstring = str(int(check_lat))
    if(int(check_lat)<0):
        minstring = 'N'+str(abs(int(check_lat)))
    if(season!=''):
        minstring = minstring+'_'+season
    outname = 'ceres_'+sat_name+'_'+dtype_name+'_all_trends_'+minstring+'.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)

# This plots the ob counts on the full grid.
# Based off of plotCERES_Climo, but with ob counts instead of flux measurements
def plotCERES_ObCount(CERES_dict,ptype='swf',save=False,trend_type='standard',season='',start_date='200101',end_date='201512',minlat=30.):
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
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    # Find the lowest lat in the file
    #lowest_lat = float(sorted(CERES_dict.keys())[0].split('x')[0])
    #lowest_lat=30
    
    if((ptype=='swf') | (ptype=='both')):
     
        # Set up the polar stereographic map
        fig1 = plt.figure()
        #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
        global sm
        if(minlat>10):
            sm = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
        else:
            sm = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')

        sm.drawcoastlines()
        sm.drawparallels(np.arange(-80.,81.,20.))
        sm.drawmeridians(np.arange(-180.,181.,20.))    
        
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
        v_min = 0.00  # Whole-year values
        mid_val = 0
        v_max = 2000.
        
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
        max_avg = -2000.
        min_avg = 4000.
        for i in range(0,len(lat_ranges)-1):
            for j in range(0,len(lon_ranges)-1):
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                #if(dictkey=='48x-97'):
                #    #print("sorted(CERES_dict[dictkey].keys())=",sorted(CERES_dict[dictkey].keys()))
                #    min_date = sorted(CERES_dict[dictkey].keys())[0].decode("utf-8")
                #    max_date = sorted(CERES_dict[dictkey].keys())[-1].decode("utf-8")
        
                # If no data are present for the curent lat/lon box, fill it with
                # black
                if(dictkey not in CERES_dict.keys()):
                    color=(0,0,0,0)
                # If, for some reason, the key is made for the current lat/lon box
                # but the dictionary is empty. fill with black.
                elif(len(CERES_dict[dictkey].keys())==0):
                    color=(0,0,0,0)
                else:
                    avgs = [CERES_dict[dictkey][date]['swf'] for date in \
                        sorted(CERES_dict[dictkey].keys())]
                    # Check the current max and min
                    #x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                    # DEAL WITH MISSING VALUES
                    # Find the slope of the line of best fit for the time series of
                    # average data
                    #avgs = np.ma.masked_array([CERES_dict[dictkey][date]['swf'] for date in sorted(CERES_dict[dictkey].keys())])
                    counts = np.ma.masked_array([CERES_dict[dictkey][date]['swf_cc'] for date in sorted(CERES_dict[dictkey].keys())])
                    avg_cnt = np.average(counts)
                    #avg_aot = sum(avgs*counts)/sum(counts)
                    #avg_aot = np.average(avgs)
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
                    #dates = temp_dates
                    #dates = temp_dates[np.where(avgs.mask==False)[0]]
                    #x_vals = np.arange(0,len(dates))
        
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,avgs)
                    #slope *= len(x_vals)
        
                    if(avg_cnt>max_avg):
                        max_avg=avg_cnt
                    elif(avg_cnt<min_avg):
                        min_avg=avg_cnt
                    if(avg_cnt>v_max):
                        avg_cnt=v_max
            
                    color = mapper.to_rgba(avg_cnt)

                    # Add value to analysis regions (if possible)
                    if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                        canada_avg+=avg_cnt
                        canada_counts+=1
                    elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                        siberia_avg+=avg_cnt
                        siberia_counts+=1
                    elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                        eus_avg+=avg_cnt
                        eus_counts+=1
                    elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                        wus_avg+=avg_cnt
                        wus_counts+=1
                    elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                        ea_avg+=avg_cnt
                        ea_counts+=1
                    elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                        wa_avg+=avg_cnt
                        wa_counts+=1
                    elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                        europe_avg+=avg_cnt
                        europe_counts+=1
        
                # Find the y coordinates of the current LatxLon box
                y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
                # Find the x coordinates of the current LatxLon box
                x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
                # Convert x and y into map coordinates
                mx, my = sm(x,y)
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

        try:
            canada_avg = canada_avg/canada_counts
            print("Canada       swf avg = ",canada_avg,"  counts = ",canada_counts)
        except (ZeroDivisionError):
            print("No Canada Ob Counts") 
        try:
            siberia_avg = siberia_avg/canada_counts
            print("Siberia      swf avg = ",siberia_avg,"  counts = ",siberia_counts)
        except (ZeroDivisionError):
            print("No Siberia Ob Counts") 
        try:
            eus_avg = eus_avg/eus_counts
            print("Eastern US   swf avg = ",eus_avg,"  counts = ",eus_counts)
        except (ZeroDivisionError):
            print("No Eastern US Ob Counts") 
        try:
            wus_avg = wus_avg/wus_counts
            print("Western US   swf avg = ",wus_avg,"  counts = ",wus_counts)
        except (ZeroDivisionError):
            print("No Western US Ob Counts") 
        try:
            ea_avg = ea_avg/ea_counts
            print("Eastern Asia swf avg = ",ea_avg,"  counts = ",ea_counts)
        except (ZeroDivisionError):
            print("No Eastern Asia Ob Counts") 
        try:
            wa_avg = wa_avg/wa_counts
            print("Western Asia swf avg = ",wa_avg,"  counts = ",wa_counts)
        except (ZeroDivisionError):
            print("No Western Asia Ob Counts") 
        try:
            europe_avg = europe_avg/europe_counts
            print("Europe       swf avg = ",europe_avg,"  counts = ",europe_counts)
        except (ZeroDivisionError):
            print("No Europe Ob Counts") 


        #start_date = min_date
        #end_date = max_date
        title_string = 'CERES Average TOA Short Wave Flux Counts\n'+start_date+' to '+end_date
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
        plt.xlabel('TOA Short Wave Flux Ob Counts',fontsize=6)
        fig1.canvas.mpl_connect('button_press_event',onclick_climo_swf)
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))
        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_spring.png'
            elif(summer is True):                                   
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_summer.png'
            elif(autumn is True):                                  
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_autumn.png'
            elif(winter is True):                                 
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_winter.png'
            else:                                                
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_swf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            if(ptype!='both'):
                plt.show()
    if((ptype=='lwf') | (ptype=='both')):
        # Set up the polar stereographic map
        fig2 = plt.figure()
        #m = Basemap(projection='ortho',lon_0=0,lat_0=40,resolution='l')
        global lm
        if(minlat>10):
            lm = Basemap(projection='npstere',boundinglat=minlat-5,lon_0=0,resolution='l')
        else:
            lm = Basemap(projection='mill',boundinglat=0,lon_0=0,resolution='l')
        lm.drawcoastlines()
        lm.drawparallels(np.arange(-80.,81.,20.))
        lm.drawmeridians(np.arange(-180.,181.,20.))    
        
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
        v_min = 0.0  # Whole-year values
        mid_val = 0
        v_max = 2000.
        
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
        max_avg = -2000.
        min_avg = 4000.
        for i in range(0,len(lat_ranges)-1):
            for j in range(0,len(lon_ranges)-1):
                dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
                #if(dictkey=='48x-97'):
                #    #print("sorted(CERES_dict[dictkey].keys())=",sorted(CERES_dict[dictkey].keys()))
                #    min_date = sorted(CERES_dict[dictkey].keys())[0].decode("utf-8")
                #    max_date = sorted(CERES_dict[dictkey].keys())[-1].decode("utf-8")
        
                # If no data are present for the curent lat/lon box, fill it with
                # black
                if(dictkey not in CERES_dict.keys()):
                    color=(0,0,0,0)
                # If, for some reason, the key is made for the current lat/lon box
                # but the dictionary is empty. fill with black.
                elif(len(CERES_dict[dictkey].keys())==0):
                    color=(0,0,0,0)
                else:
                    avgs = [CERES_dict[dictkey][date]['lwf'] for date in \
                        sorted(CERES_dict[dictkey].keys())]
                    # Check the current max and min
                    #x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
                    # DEAL WITH MISSING VALUES
                    # Find the slope of the line of best fit for the time series of
                    # average data
                    #avgs = np.ma.masked_array([CERES_dict[dictkey][date]['lwf'] for date in sorted(CERES_dict[dictkey].keys())])
                    counts = np.ma.masked_array([CERES_dict[dictkey][date]['lwf_cc'] for date in sorted(CERES_dict[dictkey].keys())])
                    avg_cnt = np.average(counts) 
                    #avg_aot = sum(avgs*counts)/sum(counts) 
                    #avg_aot = np.average(avgs)
                    #avgs = avgs[np.where(avgs.mask==False)[0]]
                    temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
                    #dates = temp_dates
                    #dates = temp_dates[np.where(avgs.mask==False)[0]]
                    #x_vals = np.arange(0,len(dates))
        
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,avgs)
                    #slope *= len(x_vals)
        
                    if(avg_cnt>max_avg):
                        max_avg=avg_cnt
                    elif(avg_cnt<min_avg):
                        min_avg=avg_cnt
                    if(avg_cnt>v_max):
                        avg_cnt=v_max
            
                    color = mapper.to_rgba(avg_cnt)

                    # Add value to analysis regions (if possible)
                    if(((lat_ranges[i] >= 50) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=-130) & (lon_ranges[j]<-70))):
                        canada_avg+=avg_cnt
                        canada_counts+=1
                    elif(((lat_ranges[i] >= 51) & (lat_ranges[i]<75)) & ((lon_ranges[j]>=60) & (lon_ranges[j]<150))):
                        siberia_avg+=avg_cnt
                        siberia_counts+=1
                    elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-100) & (lon_ranges[j]<-70))):
                        eus_avg+=avg_cnt
                        eus_counts+=1
                    elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<49)) & ((lon_ranges[j]>=-125) & (lon_ranges[j]<-101))):
                        wus_avg+=avg_cnt
                        wus_counts+=1
                    elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=100) & (lon_ranges[j]<140))):
                        ea_avg+=avg_cnt
                        ea_counts+=1
                    elif(((lat_ranges[i] >= 30) & (lat_ranges[i]<50)) & ((lon_ranges[j]>=50) & (lon_ranges[j]<99))):
                        wa_avg+=avg_cnt
                        wa_counts+=1
                    elif(((lat_ranges[i] >= 40) & (lat_ranges[i]<60)) & ((lon_ranges[j]>=-10) & (lon_ranges[j]<40))):
                        europe_avg+=avg_cnt
                        europe_counts+=1
        
                # Find the y coordinates of the current LatxLon box
                y = [lat_ranges[i],lat_ranges[i+1],lat_ranges[i+1],lat_ranges[i]]
                # Find the x coordinates of the current LatxLon box
                x = [lon_ranges[j],lon_ranges[j],lon_ranges[j+1],lon_ranges[j+1]]
                # Convert x and y into map coordinates
                mx, my = lm(x,y)
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

        try:
            canada_avg = canada_avg/canada_counts
            print("Canada       swf avg = ",canada_avg,"  counts = ",canada_counts)
        except (ZeroDivisionError):
            print("No Canada Ob Counts") 
        try:
            siberia_avg = siberia_avg/canada_counts
            print("Siberia      swf avg = ",siberia_avg,"  counts = ",siberia_counts)
        except (ZeroDivisionError):
            print("No Siberia Ob Counts") 
        try:
            eus_avg = eus_avg/eus_counts
            print("Eastern US   swf avg = ",eus_avg,"  counts = ",eus_counts)
        except (ZeroDivisionError):
            print("No Eastern US Ob Counts") 
        try:
            wus_avg = wus_avg/wus_counts
            print("Western US   swf avg = ",wus_avg,"  counts = ",wus_counts)
        except (ZeroDivisionError):
            print("No Western US Ob Counts") 
        try:
            ea_avg = ea_avg/ea_counts
            print("Eastern Asia swf avg = ",ea_avg,"  counts = ",ea_counts)
        except (ZeroDivisionError):
            print("No Eastern Asia Ob Counts") 
        try:
            wa_avg = wa_avg/wa_counts
            print("Western Asia swf avg = ",wa_avg,"  counts = ",wa_counts)
        except (ZeroDivisionError):
            print("No Western Asia Ob Counts") 
        try:
            europe_avg = europe_avg/europe_counts
            print("Europe       swf avg = ",europe_avg,"  counts = ",europe_counts)
        except (ZeroDivisionError):
            print("No Europe Ob Counts") 

        #start_date = min_date
        #end_date = max_date
        title_string = 'CERES Average TOA Long Wave Flux Ob Counts\n'+start_date+' to '+end_date
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
        cax = fig2.add_axes([0.27,0.1,0.5,0.05])
        cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
        #cb.ax.set_xlabel('Change in average aerosol index')
        plt.xticks(rotation=45,fontsize=6)
        plt.xlabel('TOA Long Wave Flux Ob Counts',fontsize=6)
        fig2.canvas.mpl_connect('button_press_event',onclick_climo_lwf)
        minstring = str(int(minlat))
        if(minlat<0):
            minstring = 'N'+str(int(abs(minlat)))
        if(save is True):
            if(spring is True):
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_spring.png'
            elif(summer is True):                                   
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_summer.png'
            elif(autumn is True):                                  
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_autumn.png'
            elif(winter is True):                                 
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_winter.png'
            else:                                                
                filename = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_lwf_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_whole_year.png'
            
            plt.savefig(filename,dpi=300)
            print("Saved image:",filename)
        else:
            plt.show()

def plotCERES_LatObCount(CERES_dict,minlat=31.,ptype='swf',save=False,season='',trend_type='standard',sfc_type='all',start_date='200101',end_date='201512'):
    trend_label=''
    lat_ranges = np.arange(minlat,90,1)
    lon_ranges = np.arange(-179,180,1)
    dates = [2005,2007,2008,2009,2010,2011,2012,2013,2014,2015]
    # If only summer months are being analyzed, remove all data except 
    # in summer
    spring = False
    summer = False
    autumn = False
    winter = False
    if(season=='spring'):
        spring = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<3) or (int(tkey[4:])>5)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='summer'):
        summer = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<6) or (int(tkey[4:])>8)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='autumn'):
        autumn = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])<9) or (int(tkey[4:])>11)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)
    elif(season=='winter'):
        winter = True 
        for lkey in sorted(CERES_dict.keys()):
            try:
                for tkey in sorted(CERES_dict[lkey].keys()):
                    if((int(tkey[4:])>2) and (int(tkey[4:])!=12)):
                        CERES_dict[lkey].pop(tkey)
            except AttributeError:
                print('Bad key',lkey)

    # Set up the array to hold the latitudinal averages
    lat_total_avgs = np.zeros(len(lat_ranges))
    lat_land_avgs = np.zeros(len(lat_ranges))
    lat_ocean_avgs = np.zeros(len(lat_ranges))

    # Find the lowest lat in the file
    #lowest_lat = 31.


    # Load in the Surface Type data
    bfile = open('/data/landmask/global_18.map','rb')
    img = bytearray(bfile.read())
    img = np.reshape(img,(1080,2160))
    scl = 0.166667
    st_lat = 89.9
    st_lon=-179.9

    
    # Set up the polar stereographic map
    fig1 = plt.figure()
    
    # Loop over all the keys and print(the regression slopes 
    # Grab the averages for the key
    for i in range(0,len(lat_ranges)):
        lat_total_avg = 0. 
        lat_total_count=0.
        lat_land_avg = 0. 
        lat_land_count=0.
        lat_ocean_avg = 0. 
        lat_ocean_count=0.
        for j in range(0,len(lon_ranges)):
    
            #print(type_num)
            # Do surface type checks
            #if(sfc_type=='land'):
            #    si = int(np.round((st_lat-lat_ranges[i])/scl))
            #    sj = int(np.round((lon_ranges[j]-st_lon)/scl))
            #    type_num = img[si,sj]
            #elif(sfc_type=='ocean'):
            #    si = int(np.round((st_lat-lat_ranges[i])/scl))
            #    sj = int(np.round((lon_ranges[j]-st_lon)/scl))
            #    type_num = img[si,sj]

            #if(sfc_flag == True):
            dictkey = str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j]))
    
            avgs = np.array([CERES_dict[dictkey][date][ptype] for date in \
                sorted(CERES_dict[dictkey].keys())])
            ##avgs = [CERES_dict[dictkey][date]['cloud'+cloud_flag] for date in \
            ##    sorted(CERES_dict[dictkey].keys())]
            # Check the current max and min
            x_vals = np.arange(0,len(CERES_dict[dictkey].keys()))
            # DEAL WITH MISSING VALUES
            temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))

            #avgs = np.array([CERES_dict[dictkey][date][ptype] for date \
            #       in sorted(CERES_dict[dictkey].keys())])
            cnts = np.array([CERES_dict[dictkey][date][ptype+'_cc'] for date \
                   in sorted(CERES_dict[dictkey].keys())])
            temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
            ##avgs = np.ma.masked_array([CERES_dict[dictkey][date]['cloud'+ \
            ##   cloud_flag] for date in sorted(CERES_dict[dictkey].keys())])
            #good_avgs = np.where(np.isnan(avgs) == False)
            #avgs = avgs[good_avgs]
            #cnts = cnts[good_avgs]
            #temp_dates = temp_dates[good_avgs] 
            dates = temp_dates
            if(len(cnts)>0):
                avg_cnts = np.average(cnts)

                si = int(np.round((st_lat-lat_ranges[i])/scl))
                sj = int(np.round((lon_ranges[j]-st_lon)/scl))
                type_num = img[si,sj]
                lat_total_avg+=avg_cnts
                lat_total_count+=1 
                if(type_num!=17):
                    lat_land_avg+=avg_cnts
                    lat_land_count+=1 
                if(type_num==17):
                    lat_ocean_avg+=avg_cnts
                    lat_ocean_count+=1 

        if(lat_total_count==0):
            print(lat_ranges[i],'No counts')
        else:
            lat_total_avg=lat_total_avg/lat_total_count
            lat_total_avgs[i] = lat_total_avg
            print(lat_ranges[i],lat_total_avgs[i])
        if(lat_land_count!=0):
            lat_land_avg=lat_land_avg/lat_land_count
            lat_land_avgs[i] = lat_land_avg
        if(lat_ocean_count!=0):
            lat_ocean_avg=lat_ocean_avg/lat_ocean_count
            lat_ocean_avgs[i] = lat_ocean_avg
  
    # Don't plot zeros
    ocean_zero_locs = np.where(lat_ocean_avgs!=0)
    land_zero_locs = np.where(lat_land_avgs!=0)
    total_zero_locs = np.where(lat_total_avgs!=0)
   
    plt.plot(lat_ranges[ocean_zero_locs],lat_ocean_avgs[ocean_zero_locs],label='Ocean')
    plt.plot(lat_ranges[land_zero_locs], lat_land_avgs[land_zero_locs],label='Land')
    plt.plot(lat_ranges[total_zero_locs],lat_total_avgs[total_zero_locs],label='Total')
    plt.legend()
    plt.title('CERES '+CERES_dict['dtype']+' '+ptype+' Ob Counts: '+\
                sfc_type+'\n'+season+ ' '+start_date+ ' to '+end_date)
    plt.xlabel('Latitude')
    plt.ylabel('CERES '+ptype+' Ob Counts')
    minstring = str(int(minlat))
    if(minlat<0):
        minstring = 'N'+str(int(abs(minlat)))
    if(save == True):
        name_end = '.png'
        if(spring==True):
            name_end = '_spring.png'
        elif(summer==True):
            name_end = '_summer.png'
        elif(autumn==True):
            name_end = '_autumn.png'
        elif(winter==True):
            name_end = '_winter.png'
        outname = 'ceres_'+CERES_dict['dtype']+CERES_dict['type_flag']+'_'+ptype+'_lat_counts_'+start_date+'_'+end_date+'_'+minstring+'to90_'+name_end

        plt.savefig(outname,dpi=300)
        print('Saved image '+outname)
    else:
        plt.show()

def plotSWF_TimeSeries(CERES_dict,lat,lon,save=False,season='',start_date='200101',end_date='201512',nozero=False):

    dictkey = str(lat)+'x'+str(lon)
    avgs = np.ma.masked_array([CERES_dict[dictkey][date]['swf'] for date in sorted(CERES_dict[dictkey].keys())])

    temp_dates = np.array(sorted(CERES_dict[dictkey].keys()))
    dates = temp_dates
    if(nozero == True):
        # Remove zeroes
        zero_locs = np.where(avgs!=0)
        avgs = avgs[zero_locs]
        dates = temp_dates[zero_locs]

    #avgs = avgs[np.where(avgs.mask==False)[0]]
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
    plt.ylabel('TOA Short Wave Flux')
    plt.legend()
    plt.show()

# This function automatically regenerates all known figures 
def plotCERES_AllPlots(CERES_dict,CERES_dict_summer,CERES_dict_winter,minlat=31., \
                       season='',trend_type='standard',sfc_type='all', \
                       start_date='200101',end_date='201512'):
    min_lat = minlat
    startdate = start_date
    enddate = end_date

    print("plotCERES whole_year all")
    plotCERES(CERES_dict,save=True,ptype='both',season='',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_Climo whole_year all")
    plotCERES_Climo(CERES_dict,save=True,ptype='both',season='',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_GlobalTrend whole_year all")
    plotCERES_GlobalTrend(CERES_dict,save=True,ptype='both',season='',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_LatTrend whole_year")
    plotCERES_LatTrend(CERES_dict,save=True,ptype='swf',season='',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat,check_lat=85.)
    print("plotCERES_LatClimo whole_year swf")
    plotCERES_LatClimo(CERES_dict,minlat=min_lat,ptype='swf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimo whole_year lwf")
    plotCERES_LatClimo(CERES_dict,minlat=min_lat,ptype='lwf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend whole_year swf")
    plotCERES_LatTrendTrend(CERES_dict,minlat=min_lat,ptype='swf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend whole_year lwf")
    plotCERES_LatTrendTrend(CERES_dict,minlat=min_lat,ptype='lwf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend whole_year swf")
    plotCERES_LatClimoTrend(CERES_dict,minlat=min_lat,ptype='swf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend whole_year lwf")
    plotCERES_LatClimoTrend(CERES_dict,minlat=min_lat,ptype='lwf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    if(CERES_dict['dtype']!='both'):
        print("plotCERES_ObCount whole_year swf")
        plotCERES_ObCount(CERES_dict,save=True,ptype='swf',season='',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_ObCount whole_year lwf")
        plotCERES_ObCount(CERES_dict,save=True,ptype='lwf',season='',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_LatObCount whole_year swf")
        plotCERES_LatObCount(CERES_dict,minlat=min_lat,ptype='swf',save=True,season='',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
        print("plotCERES_LatObCount whole_year lwf")
        plotCERES_LatObCount(CERES_dict,minlat=min_lat,ptype='lwf',save=True,season='',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)

    print("plotCERES summer all")
    plotCERES(CERES_dict_summer,save=True,ptype='both',season='summer',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_Climo summer all")
    plotCERES_Climo(CERES_dict_summer,save=True,ptype='both',season='summer',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_GlobalTrend summer all")
    plotCERES_GlobalTrend(CERES_dict_summer,save=True,ptype='both',season='summer',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_LatTrend summer")
    plotCERES_LatTrend(CERES_dict_summer,save=True,ptype='swf',season='summer',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat,check_lat=85.)
    print("plotCERES_LatClimo summer swf")
    plotCERES_LatClimo(CERES_dict_summer,minlat=min_lat,ptype='swf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimo summer lwf")
    plotCERES_LatClimo(CERES_dict_summer,minlat=min_lat,ptype='lwf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend summer swf")
    plotCERES_LatTrendTrend(CERES_dict_summer,minlat=min_lat,ptype='swf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend summer lwf")
    plotCERES_LatTrendTrend(CERES_dict_summer,minlat=min_lat,ptype='lwf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend summer swf")
    plotCERES_LatClimoTrend(CERES_dict_summer,minlat=min_lat,ptype='swf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend summer lwf")
    plotCERES_LatClimoTrend(CERES_dict_summer,minlat=min_lat,ptype='lwf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    if(CERES_dict_summer['dtype']!='both'):
        print("plotCERES_ObCount summer swf")
        plotCERES_ObCount(CERES_dict_summer,save=True,ptype='swf',season='summer',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_ObCount summer lwf")
        plotCERES_ObCount(CERES_dict_summer,save=True,ptype='lwf',season='summer',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_LatObCount summer swf")
        plotCERES_LatObCount(CERES_dict_summer,minlat=min_lat,ptype='swf',save=True,season='summer',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
        print("plotCERES_LatObCount summer lwf")
        plotCERES_LatObCount(CERES_dict_summer,minlat=min_lat,ptype='lwf',save=True,season='summer',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)


    print("plotCERES winter both")
    plotCERES(CERES_dict_winter,save=True,ptype='both',season='winter',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_Climo winter both")
    plotCERES_Climo(CERES_dict_winter,save=True,ptype='both',season='winter',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_GlobalTrend winter both")
    plotCERES_GlobalTrend(CERES_dict_winter,save=True,ptype='both',season='winter',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_LatTrend winter")
    plotCERES_LatTrend(CERES_dict_winter,save=True,ptype='swf',season='winter',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat,check_lat=85.)
    print("plotCERES_LatClimo winter swf")
    plotCERES_LatClimo(CERES_dict_winter,minlat=min_lat,ptype='swf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimo winter lwf")
    plotCERES_LatClimo(CERES_dict_winter,minlat=min_lat,ptype='lwf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend winter swf")
    plotCERES_LatTrendTrend(CERES_dict_winter,minlat=min_lat,ptype='swf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend winter lwf")
    plotCERES_LatTrendTrend(CERES_dict_winter,minlat=min_lat,ptype='lwf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend winter swf")
    plotCERES_LatClimoTrend(CERES_dict_winter,minlat=min_lat,ptype='swf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend winter lwf")
    plotCERES_LatClimoTrend(CERES_dict_winter,minlat=min_lat,ptype='lwf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    if(CERES_dict_winter['dtype']!='both'):
        print("plotCERES_ObCount winter swf")
        plotCERES_ObCount(CERES_dict_winter,save=True,ptype='swf',season='winter',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_ObCount winter lwf")
        plotCERES_ObCount(CERES_dict_winter,save=True,ptype='lwf',season='winter',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_LatObCount winter swf")
        plotCERES_LatObCount(CERES_dict_winter,minlat=min_lat,ptype='swf',save=True,season='winter',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
        print("plotCERES_LatObCount winter lwf")
        plotCERES_LatObCount(CERES_dict_winter,minlat=min_lat,ptype='lwf',save=True,season='winter',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)

