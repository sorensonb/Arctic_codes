#!/usr/bin/env python
"""
  NAME:
    CERESLib.py   

  PURPOSE:

  PYTHON VERSION:
    2.6.6

  MODULES:
    - Matplotlib
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
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Polygon
#import matplotlib.colors as color
from matplotlib.colors import rgb2hex,Normalize
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs
from os import system
import glob
# The commands module was discontinued for Python 3, so if the user
# is using python 2, import commands instead

def covariance(x,y):
    avg_x = np.average(x)
    avg_y = np.average(y)
    N = len(x)
    if(len(x)!=len(y)):
        print("ERROR: Arrays are not the same size.\nArray x has len=",len(x),\
              "\nArray y has len=",len(y))
    else:
        cov = (np.sum((x-avg_x)*(y-avg_y)))/(N-1)
        return cov

def correlation(x,y):
    avg_x = np.average(x)
    avg_y = np.average(y)
    N = len(x)
    if(len(x)!=len(y)):
        print("ERROR: Arrays are not the same size.\nArray x has len=",len(x),\
              "\nArray y has len=",len(y))
    else:
        cov = (np.sum((x-avg_x)*(y-avg_y)))/(N-1)
        std_x = np.std(x)
        std_y = np.std(y)
        corr = cov/(std_x*std_y)
        return(corr)

def ice_trend_calc(years,months,ice,avg_ice,index,str_month):
    interpx = years[np.where(months==index)]
    interper = np.poly1d(np.polyfit(interpx,ice,1)) 
    # Normalize trend by dividing by number of years
    trend = (interper(interpx[-1])-interper(interpx[0]))
    pcnt_change = (trend/avg_ice)*100.
    print(str_month+" trend: ",np.round(trend,3)," % month mean: ",np.round(pcnt_change,3))
    return trend

def trend_calc(years,months,ice,avg_ice,index,str_month):
    interpx = years[np.where(months==index)]
    # Find the number of decades being analyzed
    num_dec = 216./120.
    if(len(interpx)!=0):
        #interpx = years
        interper = np.poly1d(np.polyfit(interpx,ice,1)) 
        # Normalize trend by dividing by number of years
        trend = (interper(interpx[-1])-interper(interpx[0]))
        pcnt_change = (trend/avg_ice)*100.
        # Find the percent change per decade
        pcnt_chg_dec = pcnt_change/num_dec
        print(str_month+" trend:\t",np.round(trend,3),"\t% month mean: ",np.round(pcnt_change,3),"\t%/decade: ",np.round(pcnt_chg_dec,3))
        return trend
    else:
        return -99.

def readgridCERES(start_date,end_date,param,minlat=70.5,season='all'):
    global CERES_dict 
    CERES_dict = {}
  
    spring=False
    summer=False
    autumn=False
    winter=False
    sunlight=True
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
   
    lat_ranges = np.arange(-89.5,90.5,1.0)
    lon_ranges = np.arange(0.5,360.5,1.0)

    # Grab all the files
    base_path = '/home/bsorenson/data/CERES/SSF_1Deg/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    total_list = sorted(glob.glob(base_path+'*'))

    # Loop over all files and find the ones that match with the desired times
    final_list = []
    for f in total_list:
        fdate = f.split('_')[-1][:6]
        if((int(fdate)>=int(start_date)) & (int(fdate)<=int(end_date))):
            if(spring==True):
                if((int(fdate[-2:])>=3) & (int(fdate[-2:])<6)):
                    final_list.append(f)
            elif(summer==True):
                if((int(fdate[-2:])>=6) & (int(fdate[-2:])<9)):
                    final_list.append(f)
            elif(autumn==True):
                if((int(fdate[-2:])>=9) & (int(fdate[-2:])<12)):
                    final_list.append(f)
            elif(winter==True):
                if((int(fdate[-2:])==12) | (int(fdate[-2:])<3)):
                    final_list.append(f)
            elif(sunlight==True):
                if((int(fdate[-2:])>3) & (int(fdate[-2:])<10)):
                    final_list.append(f)
            else:
                final_list.append(f)
    time_dim = len(final_list)

    lat_indices = np.where(lat_ranges>=minlat)[0]

    CERES_dict['param'] = param
    CERES_dict['data']   = np.zeros((time_dim,len(lat_indices),len(lon_ranges)))
    CERES_dict['trends'] = np.zeros((len(lat_indices),len(lon_ranges)))
    CERES_dict['dates'] = [] 
    data = Dataset(final_list[0],'r')
    CERES_dict['parm_name'] = data.variables[param].standard_name 
    CERES_dict['parm_unit'] = data.variables[param].units 
    CERES_dict['lat'] = lat_ranges[lat_indices]
    CERES_dict['lon'] = lon_ranges
    CERES_dict['month_fix'] = ''
    CERES_dict['season']=season
    data.close()

    # Loop over good files and insert data into dictionary
    i_count = 0
    for ff in final_list:
        print(ff)
        data = Dataset(ff,'r')
        CERES_dict['data'][i_count,:,:] = data.variables[param][0,lat_indices,:]
        CERES_dict['dates'].append(ff.split('_')[-1][:6])
    #    print(data.variables[param][0,lat_indices,:])
        data.close() 
        i_count+=1
     
    #start_date = datetime.strptime(str(start_date),'%Y%m')
    #end_date = datetime.strptime(str(end_date),'%Y%m') +timedelta(days=31)
    #
    #tempdate = data.variables['time'].units
    #orig_date = datetime.strptime(tempdate.split()[2]+' '+tempdate.split()[3],'%Y-%m-%d %H:%M:%S')

    return CERES_dict

def calc_avgCERES(CERES_dict,save=False,start_date='200101',end_date='201812',minlat=70.5,month_fix=False):
    
    print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
    print("\nAlbedo Data\n")

    if(month_fix==True):
        CERES_dict['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = CERES_dict['lat']
    lon_ranges = CERES_dict['lon']
    # Create array to hold monthly averages over region
    initial_avgs = np.full(len(CERES_dict['dates']),-9999.)
    initial_years  = np.zeros(len(initial_avgs))
    initial_months = np.zeros(len(initial_avgs))
    
    # Loop over times and calculate average at each time
    for time in range(len(CERES_dict['dates'])):
        #temp = np.copy(CERES_dict['data'][time,:,:])
        #temp[temp==-999] = np.nan
        initial_years[time] = int(CERES_dict['dates'][time][:4])
        initial_months[time] = int(CERES_dict['dates'][time][4:])

        good_indices = np.where(CERES_dict['data'][time,:,:]!=-999)
        if(len(good_indices[0])!=0):
            temp = np.copy(CERES_dict['data'][time,:,:])
            temp[temp==-999] = np.nan
            initial_avgs[time] = np.nanmean(temp)
        else:
            initial_avgs[time] = -420
        # Remove October and November if desired to match up with Pistore et al 2014
        if((month_fix==True) & ((initial_months[time]==10) | (initial_months[time]==11) | (initial_months[time]==2))):
            initial_avgs[time] = -420

    # Remove missing values
    good_indices = np.where(initial_avgs!=-420)
    years = initial_years[good_indices]
    months = initial_months[good_indices]
    regional_avgs = initial_avgs[good_indices]

    #tttt[tttt==-999]=np.nan
    #np.nanmean(tttt)
    jan_alb = initial_avgs[good_indices][np.where(months==1)]
    feb_alb = initial_avgs[good_indices][np.where(months==2)]
    mar_alb = initial_avgs[good_indices][np.where(months==3)]
    apr_alb = initial_avgs[good_indices][np.where(months==4)]
    may_alb = initial_avgs[good_indices][np.where(months==5)]
    jun_alb = initial_avgs[good_indices][np.where(months==6)]
    jul_alb = initial_avgs[good_indices][np.where(months==7)]
    aug_alb = initial_avgs[good_indices][np.where(months==8)]
    sep_alb = initial_avgs[good_indices][np.where(months==9)]
    oct_alb = initial_avgs[good_indices][np.where(months==10)]
    nov_alb = initial_avgs[good_indices][np.where(months==11)]
    dec_alb = initial_avgs[good_indices][np.where(months==12)]

    # Calculate the mean albedo for each month
    avg_jan = np.average(jan_alb)
    avg_feb = np.average(feb_alb)
    avg_mar = np.average(mar_alb)
    avg_apr = np.average(apr_alb)
    avg_may = np.average(may_alb)
    avg_jun = np.average(jun_alb)
    avg_jul = np.average(jul_alb)
    avg_aug = np.average(aug_alb)
    avg_sep = np.average(sep_alb)
    avg_oct = np.average(oct_alb)
    avg_nov = np.average(nov_alb)
    avg_dec = np.average(dec_alb)
  
    # Calculate the average albedo value over the period
    avg_alb = np.average(initial_avgs[good_indices]) 
    print("Average albedo (N of ",str(int(CERES_dict['lat'][0])),") over period: ",avg_alb) 

    # Deseasonalize the data
    deseasonal_alb = np.copy(initial_avgs)
    #deseasonal_alb = np.copy(regional_avgs)
    
    deseasonal_alb[np.where(initial_months==1)]  = deseasonal_alb[np.where(initial_months==1)]  - avg_jan
    deseasonal_alb[np.where(initial_months==2)]  = deseasonal_alb[np.where(initial_months==2)]  - avg_feb
    deseasonal_alb[np.where(initial_months==3)]  = deseasonal_alb[np.where(initial_months==3)]  - avg_mar
    deseasonal_alb[np.where(initial_months==4)]  = deseasonal_alb[np.where(initial_months==4)]  - avg_apr
    deseasonal_alb[np.where(initial_months==5)]  = deseasonal_alb[np.where(initial_months==5)]  - avg_may
    deseasonal_alb[np.where(initial_months==6)]  = deseasonal_alb[np.where(initial_months==6)]  - avg_jun
    deseasonal_alb[np.where(initial_months==7)]  = deseasonal_alb[np.where(initial_months==7)]  - avg_jul
    deseasonal_alb[np.where(initial_months==8)]  = deseasonal_alb[np.where(initial_months==8)]  - avg_aug
    deseasonal_alb[np.where(initial_months==9)]  = deseasonal_alb[np.where(initial_months==9)]  - avg_sep
    deseasonal_alb[np.where(initial_months==10)] = deseasonal_alb[np.where(initial_months==10)] - avg_oct
    deseasonal_alb[np.where(initial_months==11)] = deseasonal_alb[np.where(initial_months==11)] - avg_nov
    deseasonal_alb[np.where(initial_months==12)] = deseasonal_alb[np.where(initial_months==12)] - avg_dec
    avg_deseasonal = np.average(deseasonal_alb[good_indices]) 
    print("Average deseasonal albedo over period: ",avg_deseasonal) 
    
    # Calculate the trend of the total data
    # ignore missing values
    interpx = np.arange(len(initial_years))
    total_interper = np.poly1d(np.polyfit(interpx[good_indices],initial_avgs[good_indices],1)) 
    # Normalize trend by dividing by number of years
    total_trend = (total_interper(interpx[good_indices][-1])-total_interper(interpx[good_indices][0]))
    print("Total albedo trend (200101 - 201812): ",np.round(total_trend,3))
    pcnt_change = (total_trend/avg_alb)*100.
    print("     % of average: ",pcnt_change)

    ##slope,intercept,r_value,p_value,std_err = scipy.stats.linregress(deseasonal_ext[good_indices],deseasonal_alb[good_indices])
    ##x_range = np.arange(min(deseasonal_ext[good_indices]),max(deseasonal_ext[good_indices]),0.01)
    ##newy = x_range*slope+intercept
    ##trend = newy[-1]-newy[0]

    # Calculate the trend of the deseasonalized data
    de_interpx = np.arange(len(initial_years))
    total_de_interper = np.poly1d(np.polyfit(de_interpx[good_indices],deseasonal_alb[good_indices],1)) 
    # Normalize trend by dividing by number of years
    total_de_trend = (total_de_interper(de_interpx[good_indices][-1])-total_de_interper(de_interpx[good_indices][0]))

    # Calculate the trend again using scipy
    # Also finds r_value and p_value
    slope,intercept,r_value,p_value,std_err = stats.linregress(de_interpx[good_indices],deseasonal_alb[good_indices])
    newy = de_interpx[good_indices]*slope+intercept
    test_total_de_trend = newy[-1]-newy[0] 

    print("Total deseasonal albedo trend (200101 - 201812): ",np.round(total_de_trend,5))
    print("            r_value                            : ",np.round(r_value,5))
    print("            p_value                            : ",np.round(p_value,5))
    pcnt_change = (total_de_trend/avg_alb)*100.
    print("     % of average: ",pcnt_change)
    
    # Calculate trend for each month data
    print("\nMonthly trends")
    trends = []
    trends.append(trend_calc(years,months,jan_alb,avg_jan,1,'January'))
    trends.append(trend_calc(years,months,feb_alb,avg_feb,2,'February'))
    trends.append(trend_calc(years,months,mar_alb,avg_mar,3,'March'))
    trends.append(trend_calc(years,months,apr_alb,avg_apr,4,'April'))
    trends.append(trend_calc(years,months,may_alb,avg_may,5,'May'))
    trends.append(trend_calc(years,months,jun_alb,avg_jun,6,'June'))
    trends.append(trend_calc(years,months,jul_alb,avg_jul,7,'July'))
    trends.append(trend_calc(years,months,aug_alb,avg_aug,8,'August'))
    trends.append(trend_calc(years,months,sep_alb,avg_sep,9,'September'))
    trends.append(trend_calc(years,months,oct_alb,avg_oct,10,'October'))
    trends.append(trend_calc(years,months,nov_alb,avg_nov,11,'November'))
    trends.append(trend_calc(years,months,dec_alb,avg_dec,12,'December'))
    
    
    fig0 = plt.figure(figsize=(9,9))
    ax1 = fig0.add_subplot(211)
    ax1.plot(interpx[good_indices],initial_avgs[good_indices],label='Average Albedo')
    ax1.plot(interpx[good_indices],total_interper(interpx[good_indices]),'--',label='Trend ('+str(np.round(total_trend,3))+'/Study Period)')
    ax1.set_title('Arctic (North of '+str(int(CERES_dict['lat'][0]))+' N) Average Albedo')
    #ax1.set_xlabel('Months After January 2001')
    ax1.set_xticks(np.arange(len(initial_years)+1)[::24])
    ax1.set_xticklabels(np.arange(2001,2020)[::2])
    ax1.set_ylabel('Albedo')
    ax1.legend()
    
    ax2 = fig0.add_subplot(212)
    ax2.plot(de_interpx[good_indices],deseasonal_alb[good_indices],label='Deseasonalized Average Albedo')
    ax2.plot(de_interpx[good_indices],total_de_interper(de_interpx[good_indices]),'--',label='Trend ('+str(np.round(total_de_trend,3))+'/Study Period)')
    ax2.text(de_interpx[good_indices][2],-0.02,"y = "+str(np.round(slope,4))+"x")
    ax2.text(de_interpx[good_indices][2],-0.03,"p = "+str(np.round(p_value,4)))
    ax2.set_title('Arctic (North of '+str(int(CERES_dict['lat'][0]))+' N) Deseasonalized Average Albedo')
    #ax2.set_xlabel('Months After January 2001')
    ax2.set_xticks(np.arange(len(initial_years)+1)[::24])
    ax2.set_xticklabels(np.arange(2001,2020)[::2])
    ax2.set_ylabel('Albedo Anomaly')
    ax2.legend()
   
    # Removing the monthly time series 
    ###ax2 = fig0.add_subplot(212)
    #### Plot the change for each month
    ####fig1 = plt.figure()
    ###ax2.plot(years[np.where(months==1)], jan_alb,label='January')
    ###ax2.plot(years[np.where(months==2)], feb_alb,label='February')
    ###ax2.plot(years[np.where(months==3)], mar_alb,label='March')
    ###ax2.plot(years[np.where(months==4)], apr_alb,label='April')
    ###ax2.plot(years[np.where(months==5)], may_alb,label='May')
    ###ax2.plot(years[np.where(months==6)], jun_alb,label='June')
    ###ax2.plot(years[np.where(months==7)], jul_alb,'--',label='July')
    ###ax2.plot(years[np.where(months==8)], aug_alb,'--',label='August')
    ###ax2.plot(years[np.where(months==9)], sep_alb,'--',label='September')
    ###ax2.plot(years[np.where(months==10)],oct_alb,'--',label='October')
    ###ax2.plot(years[np.where(months==11)],nov_alb,'--',label='November')
    ###ax2.plot(years[np.where(months==12)],dec_alb,'--',label='December')
    ###ax2.set_ylabel('Albedo')
    ###ax2.set_xticks(np.arange(2001,2019)[::2])
    ###ax2.legend(loc='upper center',bbox_to_anchor=(0.5,-0.05), ncol=6)
    if(save==True):
        outname = 'albedo_changes_with_deseason'+CERES_dict['month_fix']+'_'+str(int(CERES_dict['lat'][0]))+'.png'
        fig0.savefig(outname,dpi=300)
        print("Saved image ",outname)
   
    trends = np.array(trends)
    good_trends = trends[trends!=-99.] 
    fig3,ax = plt.subplots()
    labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    ax.plot(np.arange(2,len(good_trends)+2),good_trends)
    ax.set_xticks(np.arange(1,13))
    ax.set_xlim(1,13)
    ax.set_xticklabels(labels)
    ax.set_title('Trend in Monthly Arctic Albedo\n2001 - 2018')
    ax.set_ylabel('Albedo Trend (/Study Period)')
    plt.grid()
    if(save==True):
        outname = 'monthly_albedo_trends'+CERES_dict['month_fix']+'_'+str(int(CERES_dict['lat'][0]))+'.png'
        fig3.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
        plt.show()
    return initial_avgs,trends,deseasonal_alb

def calc_CERES_trend(CERES_dict,save=False,start_date='200012',end_date='201812',minlat=45.5,month_fix=False,\
                     adjusted=False):
    
    print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
    print("\nAlbedo Data\n")

    if(month_fix==True):
        CERES_dict['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = CERES_dict['lat']
    lon_ranges = CERES_dict['lon']
    ## Create array to hold monthly averages over region
    #initial_avgs = np.full(len(CERES_dict['dates']),-9999.)
    #initial_years  = np.zeros(len(initial_avgs))
    #initial_months = np.zeros(len(initial_avgs))
   
    # Loop over the lat and lon dimensions
    for xi in range(len(CERES_dict['lat'])):
        max_trend = -999
        min_trend = 999
        for yj in range(len(CERES_dict['lon'])):
            # Calculate the trend at the current box
            interp_data = CERES_dict['data'][:,xi,yj]
            good_indices = np.where(interp_data!=-9.99e+02)
            if(len(good_indices[0])==0):
                total_trend = -9999.
            else:
                interpx = np.arange(len(interp_data))
                #print(CERES_dict['data'][:,xi,yj][good_indices])
                total_interper = np.poly1d(np.polyfit( \
                    interpx[good_indices],\
                    CERES_dict['data'][:,xi,yj][good_indices],1)) 
                # Normalize trend by dividing by number of years
                total_trend = (total_interper(interpx[-1])-\
                               total_interper(interpx[0]))
                if(total_trend>max_trend):
                    max_trend = total_trend
                if(total_trend<min_trend):
                    min_trend = total_trend
            CERES_dict['trends'][xi,yj] = total_trend
        print(xi)
        #print("max trend = ",max_trend)
        #print("min trend = ",min_trend)

    #plt.pcolormesh(plot_good_data,cmap=plt.cm.bwr,vmin=-50,vmax=50)
    #plt.colorbar(label='Albedo Trend')

    if(CERES_dict['param']=='toa_sw_clr_mon'):
        min_val=-40.
        max_val = 40.
        title_adder = 'TOA Clear-Sky SWF'
    elif(CERES_dict['param']=='toa_sw_all_mon'):
        min_val=-25.
        max_val = 25.
        title_adder = 'TOA All-Sky SWF'
    elif(CERES_dict['param']=='toa_sw_cld_mon'):
        min_val=-25.
        max_val = 25.
        title_adder = 'TOA Cloudy-Sky SWF'
    elif(CERES_dict['param']=='toa_lw_clr_mon'):
        min_val=-15.
        max_val = 15.
        title_adder = 'TOA Clear-Sky LWF'
    elif(CERES_dict['param']=='toa_lw_all_mon'):
        min_val=-10.
        max_val = 10.
        title_adder = 'TOA All-Sky LWF'
    elif(CERES_dict['param']=='toa_lw_cld_mon'):
        min_val=-10.
        max_val = 10.
        title_adder = 'TOA Cloudy-Sky LWF'
    elif(CERES_dict['param']=='toa_net_clr_mon'):
        min_val=-60.
        max_val = 60.
        title_adder = 'TOA Clear-Sky Net Flux'
    elif(CERES_dict['param']=='toa_net_all_mon'):
        min_val=-25.
        max_val = 25.
        title_adder = 'TOA All-Sky Net Flux'
    elif(CERES_dict['param']=='toa_net_cld_mon'):
        min_val=-25.
        max_val = 25.
        title_adder = 'TOA Cloudy-Sky Net Flux'
    elif(CERES_dict['param']=='toa_alb_clr_mon'):
        min_val=-0.15
        max_val = 0.15
        title_adder = 'TOA Clear-Sky Albedo'
    elif(CERES_dict['param']=='toa_alb_all_mon'):
        min_val=-0.06
        max_val = 0.06
        title_adder = 'TOA All-Sky Albedo'
    elif(CERES_dict['param']=='toa_alb_cld_mon'):
        min_val=-0.06
        max_val = 0.06
        title_adder = 'TOA Cloudy-Sky Albedo'
        
    colormap = plt.cm.bwr
    #coolwarm = plt.cm.coolwarm
    #ax = plt.axes(projection=ccrs.NorthPolarStereo())
    if(adjusted==True):
        fig1 = plt.figure(figsize=(8,5))
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=45.))
    else:
        fig1 = plt.figure()
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
    ax.gridlines()
    #ax = plt.axes(projection=ccrs.Miller())
    file_season = '_'+CERES_dict['season']
    title_season = ' '+CERES_dict['season']
    if(CERES_dict['season']=='all'):
        file_season=''
        title_season=''
        
    mesh = plt.pcolormesh(CERES_dict['lon'],CERES_dict['lat'],CERES_dict['trends'],\
            transform=ccrs.PlateCarree(),vmin=min_val,vmax=max_val,cmap=colormap)
    ax.set_title('Terra CERES '+title_adder+title_season.title()+' Trend\n'+\
                 CERES_dict['dates'][0]+' - '+CERES_dict['dates'][-1])
    ender = '.png'
    if(adjusted==True):
        ax.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        axins = inset_axes(ax,width="5%",height="100%",loc='lower left',\
                    bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax.transAxes,\
                    borderpad=0)
        cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_dict['parm_unit'])
        ax.set_xlim(-5170748.535086173,5167222.438879491)
        ax.set_ylim(-3913488.8763307533,3943353.899053069)
        ender = '_adjusted'+ender
    else:
        cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_dict['parm_unit'])
    #cbar.set_label(parm_unit)
    ax.coastlines()
    if(save==True):
        outname = 'ceres_grid_terra_trend_'+CERES_dict['param']+'_'+\
                  CERES_dict['dates'][0]+'_'+CERES_dict['dates'][-1]+\
                  file_season+ender
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show() 
   
    #if(save==True):
    #    outname = 'monthly_albedo_trends'+CERES_dict['month_fix']+'_'+str(int(CERES_dict['lat'][0]))+'.png'
    #    fig3.savefig(outname,dpi=300)
    #    print("Saved image",outname)
    #else: 
    #    plt.show()
    #return initial_avgs,trends,deseasonal_alb

def make_scatter(CERES_dict_alb_clr,CERES_dict_sw_clr,\
                 CERES_dict_lw_clr,CERES_dict_net_clr):
    dot_size=8
    # Generate scatter plots
    title_text = 'December 2000 - December 2018'
    fig1 = plt.figure()
    plt.scatter(CERES_dict_alb_clr['trends'].flatten(),\
        CERES_dict_sw_clr['trends'].flatten(),s=dot_size)
    plt.title(title_text)
    plt.xlabel('Terra CERES Clear-Sky Albedo Trend')
    plt.ylabel('Terra CERES Clear-Sky SWF Trend')
    outname = 'ceres_clralb_clrsw_trend_scatter.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)
    
    fig2 = plt.figure()
    plt.scatter(CERES_dict_alb_clr['trends'].flatten(),\
        CERES_dict_lw_clr['trends'].flatten(),s=dot_size)
    plt.title(title_text)
    plt.xlabel('Terra CERES Clear-Sky Albedo Trend')
    plt.ylabel('Terra CERES Clear-Sky LWF Trend')
    outname = 'ceres_clralb_clrlw_trend_scatter.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)
    
    fig3 = plt.figure()
    plt.scatter(CERES_dict_alb_clr['trends'].flatten(),\
        CERES_dict_net_clr['trends'].flatten(),s=dot_size)
    plt.title(title_text)
    plt.xlabel('Terra CERES Clear-Sky Albedo Trend')
    plt.ylabel('Terra CERES Clear-Sky Net Flux Trend')
    outname = 'ceres_clralb_clrnet_trend_scatter.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)
def icePlots(infile,monthfix=False):

    print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
    print("\nIce Data\n")
    
    # Read data into arrays
    with open(infile,'r') as f:
        flines = f.readlines()
        initial_years  = np.zeros(len(flines)-1)
        initial_months = np.zeros(len(flines)-1)
        initial_extent = np.zeros(len(flines)-1)
        initial_area   = np.zeros(len(flines)-1)
        for i, line in enumerate(flines[1:]):
            templine = line.strip().split(',')
            initial_years[i] = templine[0]
            initial_months[i] = templine[1]
            initial_extent[i] = templine[4]
            initial_area[i] = templine[5]

    ### Remove winter months out of solidarity with albedo
    ##good_indices = np.where((initial_months<=2.)&(initial_months>=10.))
    ##years = initial_years[good_indices]
    ##months = initial_months[good_indices]
    ##extent = initial_extent[good_indices]
    years = initial_years
    months = initial_months
    extent = initial_extent
    
    jan_ext = extent[np.where(months==1)]
    feb_ext = extent[np.where(months==2)]
    mar_ext = extent[np.where(months==3)]
    apr_ext = extent[np.where(months==4)]
    may_ext = extent[np.where(months==5)]
    jun_ext = extent[np.where(months==6)]
    jul_ext = extent[np.where(months==7)]
    aug_ext = extent[np.where(months==8)]
    sep_ext = extent[np.where(months==9)]
    oct_ext = extent[np.where(months==10)]
    nov_ext = extent[np.where(months==11)]
    dec_ext = extent[np.where(months==12)]

    # Calculate the mean extent for each month
    avg_jan = np.average(jan_ext)
    avg_feb = np.average(feb_ext)
    avg_mar = np.average(mar_ext)
    avg_apr = np.average(apr_ext)
    avg_may = np.average(may_ext)
    avg_jun = np.average(jun_ext)
    avg_jul = np.average(jul_ext)
    avg_aug = np.average(aug_ext)
    avg_sep = np.average(sep_ext)
    avg_oct = np.average(oct_ext)
    avg_nov = np.average(nov_ext)
    avg_dec = np.average(dec_ext)
  
###    feby = feb_ext-np.full(len(feb_ext),avg_feb) 
###    sepy = sep_ext-np.full(len(sep_ext),avg_sep) 
###    figt = plt.figure()
###    plt.plot(feby,label='Feb anomaly')
###    plt.plot(sepy,label='Sep anomaly')
######    plt.plot(feb_ext,label='Feb')
######    plt.plot(np.full(len(feb_ext),avg_feb),label='Feb avg')
######    plt.plot(sep_ext,label='Sep')
######    plt.plot(np.full(len(sep_ext),avg_sep),label='Sep avg')
###    plt.legend()
###    plt.grid()
###    plt.show() 
    
    # Calculate the average extent value over the period
    avg_ext = np.average(extent) 
    #avg_ext = np.average(extent) 
    print("Average extent over period: ",avg_ext) 
   
    # Deseasonalize the data
    deseasonal_ext = np.copy(extent)
    
    deseasonal_ext[np.where(months==1)]  = deseasonal_ext[np.where(months==1)]  - avg_jan
    deseasonal_ext[np.where(months==2)]  = deseasonal_ext[np.where(months==2)]  - avg_feb
    deseasonal_ext[np.where(months==3)]  = deseasonal_ext[np.where(months==3)]  - avg_mar
    deseasonal_ext[np.where(months==4)]  = deseasonal_ext[np.where(months==4)]  - avg_apr
    deseasonal_ext[np.where(months==5)]  = deseasonal_ext[np.where(months==5)]  - avg_may
    deseasonal_ext[np.where(months==6)]  = deseasonal_ext[np.where(months==6)]  - avg_jun
    deseasonal_ext[np.where(months==7)]  = deseasonal_ext[np.where(months==7)]  - avg_jul
    deseasonal_ext[np.where(months==8)]  = deseasonal_ext[np.where(months==8)]  - avg_aug
    deseasonal_ext[np.where(months==9)]  = deseasonal_ext[np.where(months==9)]  - avg_sep
    deseasonal_ext[np.where(months==10)] = deseasonal_ext[np.where(months==10)] - avg_oct
    deseasonal_ext[np.where(months==11)] = deseasonal_ext[np.where(months==11)] - avg_nov
    deseasonal_ext[np.where(months==12)] = deseasonal_ext[np.where(months==12)] - avg_dec
    avg_deseasonal = np.average(deseasonal_ext) 
    print("Average deseasonal extent over period: ",avg_deseasonal) 

    # Calculate the trend of the total data
    interpx = np.arange(len(years))
    total_interper = np.poly1d(np.polyfit(interpx,extent,1)) 
    # Normalize trend by dividing by number of years
    total_trend = (total_interper(interpx[-1])-total_interper(interpx[0]))
    print("Total extent trend (200012 - 201812): ",np.round(total_trend,3))
    pcnt_change = (total_trend/avg_ext)*100.
    print("     % of average: ",pcnt_change)

    # Calculate the trend of the deseasonalized data
    interpx = np.arange(len(years))
    total_de_interper = np.poly1d(np.polyfit(interpx,deseasonal_ext,1)) 
    # Normalize trend by dividing by number of years
    total_de_trend = (total_de_interper(interpx[-1])-total_de_interper(interpx[0]))

    # Calculate the trend again using scipy
    # Also finds r_value and p_value
    slope,intercept,r_value,p_value,std_err = stats.linregress(interpx,deseasonal_ext)
    newy = interpx*slope+intercept
    test_total_de_trend = newy[-1]-newy[0] 

    print("Total deseasonal extent trend (200012 - 201812): ",np.round(total_de_trend,5))
    print("            r_value                            : ",np.round(r_value,5))
    print("            p_value                            : ",np.round(p_value,5))
    pcnt_change = (total_de_trend/avg_ext)*100.
    print("     % of average: ",pcnt_change)
    
    # Calculate trend for each month data
    print("\nMonthly trends")
    trends = []
    trends.append(trend_calc(years,months,jan_ext,avg_jan,1,'January'))
    trends.append(trend_calc(years,months,feb_ext,avg_feb,2,'February'))
    trends.append(trend_calc(years,months,mar_ext,avg_mar,3,'March'))
    trends.append(trend_calc(years,months,apr_ext,avg_apr,4,'April'))
    trends.append(trend_calc(years,months,may_ext,avg_may,5,'May'))
    trends.append(trend_calc(years,months,jun_ext,avg_jun,6,'June'))
    trends.append(trend_calc(years,months,jul_ext,avg_jul,7,'July'))
    trends.append(trend_calc(years,months,aug_ext,avg_aug,8,'August'))
    trends.append(trend_calc(years,months,sep_ext,avg_sep,9,'September'))
    trends.append(trend_calc(years,months,oct_ext,avg_oct,10,'October'))
    trends.append(trend_calc(years,months,nov_ext,avg_nov,11,'November'))
    trends.append(trend_calc(years,months,dec_ext,avg_dec,12,'December'))
   
    fig0 = plt.figure(figsize=(9,9))
    ax1 = fig0.add_subplot(211)
    ax1.plot(extent,label='Extent')
    ax1.plot(total_interper(interpx),'--',label='Trend ('+str(np.round(total_trend,3))+'/Study Period)')
    ax1.set_title('Arctic Sea Ice Extent')
    #ax1.set_xlabel('Months After January 2001')
    ax1.set_xticks(np.arange(len(years)+1)[::24])
    ax1.set_xticklabels(np.arange(2001,2020)[::2])
    ax1.set_ylabel('Extent (Mkm$^{2}$)')
    ax1.legend()

    ax2 = fig0.add_subplot(212)
    ax2.plot(deseasonal_ext,label='Deseasonalized Extent')
    ax2.plot(total_de_interper(interpx),'--',label='Trend ('+str(np.round(total_de_trend,3))+'/Study Period)')
    ax2.text(interpx[2],-1.00,"y = "+str(np.round(slope,4))+"x")
    ax2.text(interpx[2],-1.20,"p = "+str(np.round(p_value,4)))
    ax2.set_title('Arctic Deseasonalized Ice Extent')
    #ax2.set_xlabel('Months After January 2001')
    ax2.set_xticks(np.arange(len(years)+1)[::24])
    ax2.set_xticklabels(np.arange(2001,2020)[::2])
    ax2.set_ylabel('Extent Anomaly (Mkm$^{2}$)')
    ax2.legend()

    ###### Removing the monthly time series 
    ###ax2 = fig0.add_subplot(212)
    #### Plot the change for each month
    ####fig1 = plt.figure()
    ###ax2.plot(years[np.where(months==1)], jan_ext,label='January')
    ###ax2.plot(years[np.where(months==2)], feb_ext,label='February')
    ###ax2.plot(years[np.where(months==3)], mar_ext,label='March')
    ###ax2.plot(years[np.where(months==4)], apr_ext,label='April')
    ###ax2.plot(years[np.where(months==5)], may_ext,label='May')
    ###ax2.plot(years[np.where(months==6)], jun_ext,label='June')
    ###ax2.plot(years[np.where(months==7)], jul_ext,'--',label='July')
    ###ax2.plot(years[np.where(months==8)], aug_ext,'--',label='August')
    ###ax2.plot(years[np.where(months==9)], sep_ext,'--',label='September')
    ###ax2.plot(years[np.where(months==10)],oct_ext,'--',label='October')
    ###ax2.plot(years[np.where(months==11)],nov_ext,'--',label='November')
    ###ax2.plot(years[np.where(months==12)],dec_ext,'--',label='December')
    ###ax2.set_ylabel('Extent (Mkm$^{2}$)')
    ###ax2.set_xticks(np.arange(2001,2019)[::2])
    ###ax2.legend(loc='upper center',bbox_to_anchor=(0.5,-0.05), ncol=6)
    fig0.savefig('extent_changes_with_deseason.png',dpi=300)
    print("Saved image extent_changes_with_deseason.png") 
    ### Make a figure of the deseasonalized trends
    ##fig6 = plt.figure()
    ##plt.plot(deseasonal_ext)
    ##plt.ylabel('Extent Anomaly (Mkm$^{2}$)')
    ##plt.title('Deseasonalized Arctic Ice Extent')
    ##plt.show()
    
    
    fig3,ax = plt.subplots()
    labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    ax.plot(np.arange(1,13),trends)
    ax.set_xticks(np.arange(1,13))
    ax.set_ylim(-2.0,-0.5)
    ax.set_xticklabels(labels)
    ax.set_title('Trend in Monthly Sea Ice Extent\n2001 - 2018')
    ax.set_ylabel('Extent Trend (Mkm$^{2}$/Study Period)')
    plt.grid()
    fig3.savefig('monthly_extent_trends.png',dpi=300)
    print("Saved image monthly_extent_trends.png") 
    #plt.show()

    return extent,trends,months,years,deseasonal_ext

def all_years(ice_avgs,alb_avgs,years,monthfix,minlat):
    colors = np.arange(2001,2019) 
 
    cmap = cm.bwr
    #norm = Normalize(vmin=vmin,vmax=vmax)
    norm = matplotlib.colors.Normalize(vmin=2001,vmax=2018) 
    mapper = ScalarMappable(norm=norm,cmap=cmap)
    ccolors = ['red','darkorange','gold',\
               'yellow','greenyellow','green',\
               'springgreen','aquamarine','aqua',\
               'deepskyblue','cornflowerblue','blue',\
               'mediumslateblue','blueviolet','darkviolet',\
               'purple','fuchsia','deeppink']
    fig    = plt.figure(figsize=(8,7))
    ax = fig.add_subplot(111)
    #fig.set_figheight(8)
    #fig.set_figwidth(10)
    ax.set_xlim(3,17)
    ax.set_ylim(0.3,0.65) 
    patches = []
    for xi in range(18):
        pairs = []
        for yi in range(12):
            index = xi*12+yi
            if(alb_avgs[index]!=-420):
                pairs.append((ice_avgs[index],alb_avgs[index]))
        color = mapper.to_rgba(2001+xi)
        color2 = rgb2hex(color)
        poly = Polygon(pairs,fill=False,edgecolor=color2) 
        #poly = Polygon(pairs,fill=False,edgecolor=ccolors[xi]) 
        plt.gca().add_patch(poly)
        patches.append(poly) 
    cax = fig.add_axes([0.91,0.11,0.025,0.77])
    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='vertical')
    #p = PatchCollection(patches,alpha=0.4)
    #p.set_array(np.array(colors))
    #ax.add_collection(p)
    #fig.colorbar(p,ax=ax)
    ax.set_xlabel('Sea Ice Extent (Mkm$^{2}$)')
    ax.set_ylabel('Clear-sky Albedo')
    outname = 'yearly_scatter_ice_vs_alb'+CERES_dict['month_fix']+'_'+str(int(CERES_dict['lat'][0]))+'.png'
    fig.savefig(outname,dpi=300)
    print('Saved image ',outname)
    #plt.show()

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

