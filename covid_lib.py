#!/usr/bin/env python

"""


"""

import numpy as np
from subprocess import check_output
import datetime
from metpy.units import units
from metpy.calc import wind_components
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import imageio

def split_metar_file(metar_file):
    # This code reads a METAR file and calculates
    # statistics for temp and wind data. For temp,
    # the average, min, and max temp for each day
    # is calculated. For wind, the average u and v
    # components of the wind are found.

    # Open metar file
    with open(metar_file,'r') as mfile:

        # Grab location from input file
        site_name = metar_file.split('_')[0]
        # Open output file
        # Outfile will have format:
        # Station,Valid,AvgTemp,MinTemp,MaxTemp,AvgU,AvgV
        outname = site_name+'_area_temp_wind.csv'
        outfile = open(outname,'w')
        outfile.write('Site,Valid,AvgTemp,MinTemp,MaxTemp,AvgU,AvgV\n')

        temps = []
        us    = []
        vs    = []
        old_date = ''
        for i, line in enumerate(mfile):
            if(i>0):
                templine = line.strip().split(',')
                #print(i,templine)
                # Grab variables
                ICAO  = templine[0]     # Site ID
                valid = templine[1]     # Valid ime
                lon   = templine[2]     # Lon
                lat   = templine[3]     # Lat
                tmp   = templine[4]     # Temperature
                wdir  = templine[5]     # Direction
                sknt  = templine[6]     # Speed in knots

                # Pull the date out of the valid time
                date = valid.split()[0]
                # Calculate averages for the day if a new day is found
                if((date!=old_date) & (old_date!='')):
                    np_temps = np.array(temps)
                    np_us = np.array(us)
                    np_vs = np.array(vs)

                    mean_temp = np.average(np_temps)
                    min_temp = np.min(np_temps)
                    max_temp = np.max(np_temps)

                    mean_u = np.average(np_us)
                    mean_v = np.average(np_vs)

                    # Print out the data to the files
                    print("Calculated averages: ",old_date)
                    outfile.write('{:},{:},{:.3},{:.3},{:.3},{:.3},{:.3}\n'.format(ICAO,\
                        old_date,mean_temp,min_temp,max_temp,\
                        mean_u,mean_v))

                    # Reset the temporary arrays and old date
                    old_date = date
                    temps = []
                    us = []
                    vs = []
                # Insert values into arrays if not missing
                if(tmp!='M'):
                    temps.append(float(tmp))
                if((wdir!='M') & (sknt!='M')):
                    # convert sknt to m/s
                    wspd = (float(sknt)*units.knots).to('m/s')
                    # Calculate u and v comps
                    u,v = wind_components(wspd,float(wdir)*units.degrees)
                    us.append(u.magnitude)
                    vs.append(v.magnitude)
                if(old_date==''):
                    old_date = date
    outfile.close()
    print("Close filed:",outname)

def read_meteo_data(meteo_file,meteo_dict=None):
    # Grab the daily averaged temp and wind data outputted from
    # split_metar_file()
    site_name = meteo_file.split('_')[0]

    if(meteo_dict==None):
        meteo_dict = {}
    meteo_dict[site_name] = {}

    # Loop over the meteorology file and insert the data
    with open(meteo_file,'r') as mfile:
        mlines = mfile.readlines()
        meteo_dict[site_name]['avg_temp'] = np.zeros((len(mlines)-1))
        meteo_dict[site_name]['min_temp'] = np.zeros((len(mlines)-1))
        meteo_dict[site_name]['max_temp'] = np.zeros((len(mlines)-1))
        meteo_dict[site_name]['avg_u'] = np.zeros((len(mlines)-1))
        meteo_dict[site_name]['avg_v'] = np.zeros((len(mlines)-1))
        meteo_dict[site_name]['dates'] = []
        # Skip the header line and insert stuff
        print("Reading meteorology file")
        for i, line in enumerate(mlines[1:]):
            templine = line.strip().split(',')
            print(templine[1])
            meteo_dict[site_name]['dates'].append(templine[1])
            meteo_dict[site_name]['avg_temp'][i] = float(templine[2])
            meteo_dict[site_name]['min_temp'][i] = float(templine[3])
            meteo_dict[site_name]['max_temp'][i] = float(templine[4])
            meteo_dict[site_name]['avg_u'][i] = float(templine[5])
            meteo_dict[site_name]['avg_v'][i] = float(templine[6])

    return meteo_dict
# num_mets is the number of days of met observations in the
# meteo_dict. Assumes that the same start date is used for 
# the meteo data and the covid data.
def read_covid_data(num_mets):
    # Grab the confirmed cases, deaths, and population csv files
    # from the directory.
    case_file = 'covid_confirmed_usafacts.csv'
    death_file = 'covid_deaths_usafacts.csv'
    pop_file = 'covid_county_population_usafacts.csv'

    # Create dictionary
    covid_dict = {}
    covid_dict['data'] = {}
    # The dictionary will have keys:
    # covid_dict['dates']
    # covid_dict['temp']
    # covid_dict['u']
    # covid_dict['v']
    # covid_dict['data']
    #   > This will be a second dictionary holding all the
    #       counties.
    #   covid_dict['data']['Queens County NY']
    #       > Below this dictionary will be all the daily 
    #         death counts and confirmed case counts
    #       covid_dict['data']['Queens County NY']['deaths']
    #       covid_dict['data']['Queens County NY']['cases']
    #       covid_dict['data']['Queens County NY']['population']

    # Loop over cases file and insert 
    # Assume that the cases and meteorology data have the same
    # start date.
    start_ind = 4
    with open(case_file,'r') as cfile:
        clines = cfile.readlines()
        # Skip the header line and insert stuff
        print("Reading case file")
        for i, line in enumerate(clines[1:]):
            templine = line.strip().split(',')
            county_key = ' '.join([templine[1],templine[2]])
            covid_dict['data'][county_key] = {}
            covid_dict['data'][county_key]['cases'] = templine[start_ind:start_ind+num_mets]
    # Read the death data
    with open(death_file,'r') as dfile:
        dlines = dfile.readlines()
        # Skip the header line and insert stuff
        print("Reading death file")
        for i, line in enumerate(dlines[1:]):
            templine = line.strip().split(',')
            county_key = ' '.join([templine[1],templine[2]])
            try:
                covid_dict['data'][county_key]['deaths'] = templine[start_ind:start_ind+num_mets]
            except:
                county_key = ' '.join([templine[1],'and','City',templine[2]])
                if(county_key in covid_dict['data'].keys()):
                    covid_dict['data'][county_key]['deaths'] = templine[start_ind:start_ind+num_mets]

    # Read the population data
    with open(pop_file,'r') as pfile:
        plines = pfile.readlines()
        # Skip the header line and insert stuff
        print("Reading population file")
        for i, line in enumerate(plines[1:]):
            templine = line.strip().split(',')
            county_key = ' '.join([templine[1],templine[2]])
            try:
                covid_dict['data'][county_key]['population'] = templine[3]
            except:
                county_key = ' '.join([templine[1],'and','City',templine[2]])
                if(county_key in covid_dict['data'].keys()):
                    covid_dict['data'][county_key]['population'] = templine[3]

    return covid_dict
def get_daily_data(case_str,var='olr'):
    latmax = 180  # 19.5 N
    latmin = 0   # 19.5 S

    if(var=='olr'):
        file_var = 'toa_lw_all_daily'
        label = 'OLR Anomaly [W/m2]'
        anom_max = 60
        anom_min = -60
    elif(var=='press'):
        file_var = 'aux_surf_press_daily'
        label = 'SLP Anomaly [hPa]'
        anom_max = 5 
        anom_min = -5
    else:
        print("Invalid variable option. Either 'olr' or 'press'")
        return

    # Path to data files
    data_path = '/data/CERES/SSF_1Deg/daily/Terra/'

    if(case_str=='case1'):
        # Grab the 2000 - 2001 file
        filename = 'CERES_SSF1deg-Day_Terra-MODIS_Ed4A_Subset_20000601-20010630.nc'
    elif(case_str=='case2'):
        # Grab the 2000 - 2001 file
        filename = 'CERES_SSF1deg-Day_Terra-MODIS_Ed4A_Subset_20170501-20180531.nc'
    else:
        print("Invalid case. Try again")
        return  

    # Grab the daily data
    daily_dict = {}
    data = Dataset(data_path+filename,"r")
    #data = Dataset(data_path+"CERES_SSF1deg-Day_Terra-MODIS_Ed4A_Subset_20171101-20180228.nc","r")

    # Grab the times and create datetime objects
    start_date = datetime.date(2000,3,1)
    dates = [start_date + datetime.timedelta(days= \
            int(data.variables['time'][ti].data))  \
            for ti in range(data.dimensions['time'].size)]

    daily_dict['var'] = var
    daily_dict['time'] = data.variables['time'][:]
    daily_dict['lat']  = data.variables['lat'][latmin:latmax]
    daily_dict['lon']  = data.variables['lon'][:]
    print("Accessing dictionary with file_var = ",file_var)
    daily_dict['daily_vals'] = data.variables[file_var][:,latmin:latmax,:]
    daily_dict['dates'] = dates
    daily_dict['anom_max'] = anom_max
    daily_dict['anom_min'] = anom_min
    daily_dict['label'] = label
    
    return daily_dict

def plot_anomaly(daily_dict,climo_dict,tidx):
    # Calculate anomalies
    anomalies = daily_dict['daily_vals'][tidx,:,:] - climo_dict['climo'][:,:]

    data_proj = ccrs.PlateCarree()
    plot_proj = ccrs.Miller()
    colormap = plt.cm.bwr
   
    fig1 = plt.figure() 
    ax = plt.subplot(projection = plot_proj)
    ax.coastlines()
    ax.set_extent((40,170,-40,40),crs=data_proj)
    ax.set_yticks([-40,-30,-20,-10,0,10,20,30,40],crs=data_proj)
    ax.set_xticks([40,60,80,100,120,140,160],crs=data_proj)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    mesh = plt.pcolormesh(daily_dict['lon'], daily_dict['lat'], \
             anomalies, transform=data_proj,cmap=colormap, \
             vmin=daily_dict['anom_min'],vmax=daily_dict['anom_max'])
    cbar = plt.colorbar(mesh,shrink=0.62,label= daily_dict['label'])
    ax.set_title(datetime.datetime.strftime(daily_dict['dates'][tidx],'%Y-%m-%d'))
    outname = 'mjo_'+daily_dict['var']+'_anom_'+datetime.datetime.strftime(daily_dict['dates'][tidx],'%Y%m%d')+'.png'
    print("Saving image ",outname)
    plt.savefig(outname,dpi=100)    
    plt.close('all')

def make_gif(var):
    # Grab all the image files
    filenames = check_output('ls mjo_'+var+'_anom_*.png',shell=True).decode('utf-8').strip().split('\n')
    images = []
    filedate = filenames[0][13:19]
    filedate_end = filenames[-1][13:19]
    for filename in filenames:
        print(filename)
        images.append(imageio.imread(filename))
    gifname = 'mjo_'+var+'_loop_'+filedate+'_'+filedate_end+'.gif'
    imageio.mimsave(gifname,images)
    print("Saved gif "+gifname)

def plot_lon_time(daily_dict,climo_dict):

    # Make array to hold lat averaged values
    latmin = 80   # 10 S
    latmax = 100  # 10 N
    lonmin = 0 
    lonmax = 360

    longitudes = daily_dict['lon'][lonmin:lonmax]
    date_ranges = np.arange(len(daily_dict['dates']))
    lat_avgs = np.zeros((daily_dict['daily_vals'].shape[0],lonmax-lonmin))

    for ti in range(daily_dict['daily_vals'].shape[0]):
        print(daily_dict['dates'][ti])
        # Calculate anomalies
        anomalies = daily_dict['daily_vals'][ti,:,:] - climo_dict['climo'][:,:]

        # Calculate lat average
        lat_avg = np.ma.average(anomalies[latmin:latmax,:].data,axis=0)

        # Insert the desired lat averages over the desired longitude
        # range into the average array
        lat_avgs[ti,:] = lat_avg[lonmin:lonmax]

    if(daily_dict['daily_vals'].shape[0]>180):
        figsize=(5,11)
        label_int = 32
        format_str = '%b %Y'
    else:
        figsize=(8,7)
        label_int = 14
        format_str = '%Y-%m-%d'
    # Get the dates to put in list
    # Go every 2 weeks
    label_text = [datetime.datetime.strftime(dx,format_str) \
                 for dx in daily_dict['dates'][::label_int]]
    myFmt = mdates.DateFormatter(format_str)
    plevels = np.arange(daily_dict['anom_min'],daily_dict['anom_max'])
    #plevels = np.linspace(-62,60,100)
    lat_avgs[lat_avgs < daily_dict['anom_min']] = daily_dict['anom_min']
    # Make the figure
    fig1 = plt.figure(figsize=figsize)
    plt.contourf(longitudes,date_ranges,lat_avgs,\
                 cmap = plt.cm.bwr,levels=plevels,\
                 extend='min')
    plt.colorbar(label=daily_dict['label'])
    plt.yticks(np.arange(daily_dict['daily_vals'].shape[0])[::label_int],\
               label_text,fontsize=8)
    plt.xlabel('Longitude')
    plt.xticks(np.arange(len(daily_dict['lon']))[::60],\
               daily_dict['lon'][::60],fontsize=8)
#    plt.ylabel('Days after 2017/11/01')
    outname = 'mjo_lon_time_'+daily_dict['var']+'_'+datetime.datetime.strftime(daily_dict['dates'][0],'%Y%m%d')+'.png'
    plt.savefig(outname,dpi=300)
    print("Saved image "+outname)
    plt.show()
    ##data_proj = ccrs.PlateCarree()
    ##plot_proj = ccrs.Miller()
    ##colormap = plt.cm.bwr
   
    ##fig1 = plt.figure() 
    ##ax = plt.subplot(projection = plot_proj)
    ##ax.coastlines()
    ##ax.set_extent((40,170,-40,40),crs=data_proj)
    ##ax.set_yticks([-40,-30,-20,-10,0,10,20,30,40],crs=data_proj)
    ##ax.set_xticks([40,60,80,100,120,140,160],crs=data_proj)
    ##lon_formatter = LongitudeFormatter(zero_direction_label=True)
    ##lat_formatter = LatitudeFormatter()
    ##ax.xaxis.set_major_formatter(lon_formatter)
    ##ax.yaxis.set_major_formatter(lat_formatter)
    ##mesh = plt.pcolormesh(daily_dict['lon'], daily_dict['lat'], \
    ##         anomalies, transform=data_proj,cmap=colormap, \
    ##         vmin=-150,vmax=150)
    ##cbar = plt.colorbar(mesh,shrink=0.62,label= 'OLR Anomaly [W/m2]')
    ##ax.set_title(datetime.datetime.strftime(daily_dict['dates'][tidx],'%Y-%m-%d'))
    ##outname = 'mjo_olr_anom_'+datetime.datetime.strftime(daily_dict['dates'][tidx],'%y%m%d')+'.png'
    ##print("Saving image ",outname)
    ##plt.savefig(outname,dpi=100)    
    ##plt.close('all')

