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

county_dict = {
    'nyc': 'Queens County NY',              # Populated
    'newOrleans': 'Orleans Parish LA',      # Populated 
    'denton': 'Denton County TX' ,          # Populated
    'lac': 'Los Angeles County CA',         # Populated
    'chi': 'Cook County IL',                # Populated
    'gfk': 'Grand Forks County ND',         # Rural
    'eagleCoCo': 'Eagle County CO',         # Rural
    'duluth': 'St. Louis County MN',        # Rural
    'abbevilleCoSC': 'Abbeville County SC', # Rural
    'medford': 'Jackson County OR',         # Rural
}

##!## County areas, in square miles
##!#area_dict = {
##!#    'nyc': 108.1,
##!#    'newOrleans': 349.8,
##!#    'denton': 953.,
##!#    'lac': 4751.,
##!#    'chi': 1635.,
##!#    'gfk': 1440.,
##!#    'eagleCoCo': 1692.,
##!#    'duluth': 6860.,
##!#    'abbevilleCoSC': 511.,
##!#    'medford': 2802.,
##!#}

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

                    if(len(np_temps)==0):
                        mean_temp=-99.
                        min_temp=-99.
                        max_temp=-99.
                    else:
                        mean_temp = np.average(np_temps)
                        min_temp = np.min(np_temps)
                        max_temp = np.max(np_temps)

                    if(len(np_us)==0):
                        mean_u=-99.
                        mean_v=-99.
                    else:
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
        print("Reading meteorology file",site_name)
        for i, line in enumerate(mlines[1:]):
            templine = line.strip().split(',')
            meteo_dict[site_name]['dates'].append(templine[1])
            # Average temperature
            if(templine[2]=='-99.0'):
                meteo_dict[site_name]['avg_temp'][i] = np.nan 
            else:
                meteo_dict[site_name]['avg_temp'][i] = float(templine[2])
            # Min temperature 
            if(templine[3]=='-99.0'):
                meteo_dict[site_name]['min_temp'][i] = np.nan
            else:
                meteo_dict[site_name]['min_temp'][i] = float(templine[3])
            # Max temperature
            if(templine[4]=='-99.0'):
                meteo_dict[site_name]['max_temp'][i] = np.nan
            else:
                meteo_dict[site_name]['max_temp'][i] = float(templine[4])
            # U wind
            if(templine[5]=='-99.0'):
                meteo_dict[site_name]['avg_u'][i] = np.nan
            else:
                meteo_dict[site_name]['avg_u'][i] = float(templine[5])
            # V wind
            if(templine[6]=='-99.0'):
                meteo_dict[site_name]['avg_v'][i] = np.nan
            else:
                meteo_dict[site_name]['avg_v'][i] = float(templine[6])

    return meteo_dict
# num_mets is the number of days of met observations in the
# meteo_dict. Assumes that the same start date is used for 
# the meteo data and the covid data.
def read_covid_data():
    # Grab the confirmed cases, deaths, and population csv files
    # from the directory.
    case_file = 'covid_confirmed_usafacts.csv'
    death_file = 'covid_deaths_usafacts.csv'
    pop_file = 'covid_county_population_usafacts.csv'
    area_file = 'county_areas.csv'

    # Create dictionary
    covid_dict = {}
    covid_dict['data'] = {}
    # The dictionary will have keys:
    # covid_dict['dates']
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
            if(i==0):
                covid_dict['dates'] = templine[4:]
            county_key = ' '.join([templine[1],templine[2]])
            covid_dict['data'][county_key] = {}
            covid_dict['data'][county_key]['cases'] = np.array([int(num) for num in templine[start_ind:]])
            # Calculate cases per 
    # Read the death data
    with open(death_file,'r') as dfile:
        dlines = dfile.readlines()
        # Skip the header line and insert stuff
        print("Reading death file")
        for i, line in enumerate(dlines[1:]):
            templine = line.strip().split(',')
            county_key = ' '.join([templine[1],templine[2]])
            try:
                covid_dict['data'][county_key]['deaths'] = np.array([int(num) for num in templine[start_ind:]])
            except:
                county_key = ' '.join([templine[1],'and','City',templine[2]])
                if(county_key in covid_dict['data'].keys()):
                    covid_dict['data'][county_key]['deaths'] = np.array([int(num) for num in templine[start_ind:]])

    # Read the population data
    with open(pop_file,'r') as pfile:
        plines = pfile.readlines()
        # Skip the header line and insert stuff
        print("Reading population file")
        for i, line in enumerate(plines[1:]):
            templine = line.strip().split(',')
            county_key = ' '.join([templine[1],templine[2]])
            try:
                covid_dict['data'][county_key]['population'] = int(templine[3])
            except:
                county_key = ' '.join([templine[1],'and','City',templine[2]])
                if(county_key in covid_dict['data'].keys()):
                    covid_dict['data'][county_key]['population'] = int(templine[3])

    # Read in Matt's county areas file
    with open(area_file,'r') as afile:
        alines = afile.readlines() 
        # Skip the header line and insert stuff
        print("Reading area file")
        for i, line in enumerate(alines[1:]):
            templine = line.strip().split(',')
            county_key = templine[0]
            if(county_key in covid_dict['data'].keys()):
                covid_dict['data'][county_key]['area'] = float(templine[2])
                # Calculate population density
                covid_dict['data'][county_key]['pop_density'] = \
                    covid_dict['data'][county_key]['population']/\
                    covid_dict['data'][county_key]['area']

    return covid_dict

def plot_single_county(meteo_dict,covid_dict,county,meteo_var,covid_var,save=False):

    meteo_key = county
    # Use the user-defined county variable to grab the county
    # key for the covid dict
    county_key = county_dict[county]

    # Plot stuff
    fig, ax1 = plt.subplots()
    plt.title(county_key)
    color='tab:red'
    ax1.plot(meteo_dict[meteo_key][meteo_var],color=color,label='meteo_var')
    ax1.set_ylabel(meteo_var.title(),color=color)
    
    ax2 = ax1.twinx()
    
    color='tab:blue'
    ax2.plot(covid_dict['data'][county_key][covid_var],color=color,label='Cases')
    ax2.set_ylabel('# '+covid_var.title(),color=color)
    
    image_name = county+'_'+meteo_var+'_'+covid_var+'_single_county.png'
    if(save==True):
        plt.savefig(image_name,dpi=300)
        print("Saved image",image_name)
    else:
        plt.show()
    plt.close()

def plot_counties_covid(covid_dict,covid_var,save=False):
    fig1, ax1 = plt.subplots()
    ax1.set_title(covid_var.title())
    for key in county_dict.keys():
        ax1.plot(covid_dict['data'][county_dict[key]][covid_var],label=county_dict[key])
    ax1.set_ylabel(covid_var.title())
    ax1.set_xlabel('Days after 2020-01-22')
    plt.legend()

    image_name = covid_var+'_all_counties.png'
    if(save==True):
        plt.savefig(image_name,dpi=300)
        print("Saved image",image_name)
    else:
        plt.show()
    plt.close()

def plot_counties_meteo(meteo_dict,meteo_var,save=False):
    fig1, ax1 = plt.subplots()
    ax1.set_title(meteo_var.title())
    for key in meteo_dict.keys():
        #print(key)
        ax1.plot(meteo_dict[key][meteo_var],label=key)
    ax1.set_ylabel(meteo_var.title())
    ax1.set_xlabel('Days after 2020-01-22')
    plt.legend()

    image_name = meteo_var+'_all_counties.png'
    if(save==True):
        plt.savefig(image_name,dpi=300)
        print("Saved image",image_name)
    else:
        plt.show()
    plt.close()


# Plot scatter plots of avg temp vs cases and avg u and v vs cases
def plot_cases_scatter(meteo_dict,covid_dict,county,meteo_var,covid_var,save=False):

    # Assume that the starting date for both the COVID and 
    # meteorology data are the same, so use the smaller of the 
    # two to determine the size

    # Meteorology data are on the x axis
    x1_vals = np.copy(meteo_dict[county][meteo_var])
    #x2_vals = np.copy(meteo_dict[county]['avg_u'])
    #x3_vals = np.copy(meteo_dict[county]['avg_v'])

    # COVID data are on the y axis
    y_vals = np.copy(covid_dict['data'][county_dict[county]][covid_var])

    if(len(x1_vals)>len(y_vals)):
        x1_vals = x1_vals[:len(y_vals)]
        #x2_vals = x2_vals[:len(y_vals)]
        #x3_vals = x3_vals[:len(y_vals)]
    elif(len(x1_vals)<len(y_vals)):
        y_vals = y_vals[:len(x1_vals)]

    # Eliminate data where case or death counts are zero
    x1_vals = x1_vals[y_vals!=0] 
    y_vals = y_vals[y_vals!=0] 

    fig1, axs = plt.subplots()
    #fig1, axs = plt.subplots(3)
    markersize=8
    # Plot avg temp
    axs.scatter(x1_vals,y_vals,s=markersize)
    axs.set_title(county_dict[county]+' '+meteo_var.title())
    axs.set_ylabel(covid_var.title())
    axs.set_xlabel(meteo_var.title())

    ### Plot avg zonal wind
    ##axs[1].scatter(x2_vals,y_vals,s=markersize)
    ##axs[1].set_title('Avg Zonal Wind')
    ##axs[1].set_ylabel(covid_var.title())
    ##axs[1].set_xlabel('Avg Zonal Wind [m/s]')

    ### Plot avg zonal wind
    ##axs[2].scatter(x3_vals,y_vals,s=markersize)
    ##axs[2].set_title('Avg Meridional Wind')
    ##axs[2].set_ylabel(covid_var.title())
    ##axs[2].set_xlabel('Avg Meridional Wind [m/s]')
    if(save==True):
        image_name = county+'_'+meteo_var+'_'+covid_var+'_scatter.png'
        plt.savefig(image_name,dpi=300)
        print("Saved image",image_name)
    else:
        plt.show()
    plt.close()

def plot_pdense_cases_scatter(covid_dict):
    case_key='cases_p100k'
    fig1,ax = plt.subplots()
    for ckey in covid_dict['data'].keys():
        if('pop_density' in covid_dict['data'][ckey].keys()):
            cases = covid_dict['data'][ckey][case_key][-1]
            pop_dense = covid_dict['data'][ckey]['pop_density']
            ax.scatter(pop_dense,cases,color='black')
    ax.set_xlabel('Population density [#/sq mi]') 
    ax.set_ylabel('Cases per 100k as of 2020/04/23') 
    plt.show()
    plt.close() 

def plot_pdense_temp_scatter(meteo_dict,covid_dict):
    # To find county cases for a specific day:
    #[covid_dict['data'][county_dict[ckey]]['cases'][40] for \
    #    ckey in county_dict.keys()]

    date_index = 40

    fig1, ax = plt.subplots()

    # Calculate average temperature after the date index

    # Determine coloring based on cases for each county
    # In future, color by exponential fit parameters
    
    # Plot scatter, with population density on x-axis, avg
    # temperature after date index on y-axis, and colored by
    # cases.
    for ckey in county_dict.keys():
        tavg = np.average(meteo_dict[ckey]['avg_temp'][date_index:])
        pop_dense = covid_dict['data'][county_dict[ckey]]['pop_density']
        ax.scatter(pop_dense,tavg,label=county_dict[ckey])
    
    ax.set_xlabel('Population Density [#/mi2]')
    ax.set_ylabel('Average Temperature During Cases [degC]')
    plt.show()
    plt.close()
