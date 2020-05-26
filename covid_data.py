#!/usr/bin/env python

"""


"""

import numpy as np
import matplotlib.pyplot as plt
from covid_lib import *

all_locs = True
if(all_locs == True):
    for i, key in enumerate(county_dict):
        if(i==0):
            meteo_dict = read_meteo_data(key+'_area_temp_wind.csv')
        else:
            meteo_dict = read_meteo_data(key+'_area_temp_wind.csv',meteo_dict=meteo_dict)
            
## Read NYC and newOrleans meteo data
#meteo_dict = read_meteo_data('nyc_area_temp_wind.csv')
#meteo_dict = read_meteo_data('newOrleans_area_temp_wind.csv',meteo_dict=meteo_dict)
#meteo_dict = read_meteo_data('gfk_area_temp_wind.csv')
#meteo_dict = read_meteo_data('dto_area_temp_wind.csv',meteo_dict=meteo_dict)

# Read the COVID data
covid_dict = read_covid_data()

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

#covid_dict = process_covid_data(covid_dict)

for key in county_dict.keys():
    # Calculate cases / 100000 for the desired counties
    covid_key = county_dict[key]
    covid_dict['data'][covid_key]['cases_p100k'] = 100000*\
        covid_dict['data'][covid_key]['cases']/covid_dict['data'][covid_key]['population']
    covid_dict['data'][covid_key]['deaths_p100k'] = 100000*\
        covid_dict['data'][covid_key]['deaths']/covid_dict['data'][covid_key]['population']

    # Calculate daily cases and deaths
    daily_cases  = np.copy(covid_dict['data'][covid_key]['cases'])
    daily_deaths = np.copy(covid_dict['data'][covid_key]['deaths'])
    daily_cases[0] = 0
    daily_deaths[0] = 0
    for i in range(1,len(daily_cases)):
        daily_cases[i] = covid_dict['data'][covid_key]['cases'][i]-\
                         covid_dict['data'][covid_key]['cases'][i-1]
        daily_deaths[i] = covid_dict['data'][covid_key]['deaths'][i]-\
                         covid_dict['data'][covid_key]['deaths'][i-1]
    covid_dict['data'][covid_key]['daily_cases'] = daily_cases 
    covid_dict['data'][covid_key]['daily_deaths'] = daily_deaths 

    ##!## Calculate population density
    ##!#covid_dict['data'][covid_key]['pop_density'] = \
    ##!#    covid_dict['data'][covid_key]['population']/area_dict[key]

#county='gfk'
#mvar = 'avg_temp'
#cvar = 'cases'
#plot_single_county(meteo_dict,covid_dict,county,mvar,cvar)

##!#covid_variables = ['cases','deaths','cases_p100k','deaths_p100k','daily_cases','daily_deaths']
##!#meteo_variables = ['avg_temp','min_temp','max_temp','avg_u','avg_v']
##!## Make all the images
##!#for lkey in county_dict.keys():
##!#    print("Location = ",lkey)
##!#    covid_key = county_dict[lkey]
##!#    for covid_var in covid_variables:
##!#        print("  COVID variable = ",covid_var)
##!#        plot_counties_covid(covid_dict,covid_var)
##!#        for meteo_var in meteo_variables:
##!#            print("  meteo variable = ",meteo_var)
##!#            plot_counties_meteo(meteo_dict,meteo_var)
##!#            plot_single_county(meteo_dict,covid_dict,lkey,meteo_var,covid_var)
##!#            plot_cases_scatter(meteo_dict,covid_dict,lkey,meteo_var,covid_var)


