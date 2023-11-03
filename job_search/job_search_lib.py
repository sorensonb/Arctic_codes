"""


"""

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from geopy.geocoders import Nominatim
import sys
import os
import json
import shapely
import cartopy.geodesic as cgeodesic

home_dir = os.environ['HOME']

sys.path.append(home_dir)
from python_lib import *

us_state_to_abbrev = {
    "Alabama": "AL",
    "Alaska": "AK",
    "Arizona": "AZ",
    "Arkansas": "AR",
    "California": "CA",
    "Colorado": "CO",
    "Connecticut": "CT",
    "Delaware": "DE",
    "Florida": "FL",
    "Georgia": "GA",
    "Hawaii": "HI",
    "Idaho": "ID",
    "Illinois": "IL",
    "Indiana": "IN",
    "Iowa": "IA",
    "Kansas": "KS",
    "Kentucky": "KY",
    "Louisiana": "LA",
    "Maine": "ME",
    "Maryland": "MD",
    "Massachusetts": "MA",
    "Michigan": "MI",
    "Minnesota": "MN",
    "Mississippi": "MS",
    "Missouri": "MO",
    "Montana": "MT",
    "Nebraska": "NE",
    "Nevada": "NV",
    "New Hampshire": "NH",
    "New Jersey": "NJ",
    "New Mexico": "NM",
    "New York": "NY",
    "North Carolina": "NC",
    "North Dakota": "ND",
    "Ohio": "OH",
    "Oklahoma": "OK",
    "Oregon": "OR",
    "Pennsylvania": "PA",
    "Rhode Island": "RI",
    "South Carolina": "SC",
    "South Dakota": "SD",
    "Tennessee": "TN",
    "Texas": "TX",
    "Utah": "UT",
    "Vermont": "VT",
    "Virginia": "VA",
    "Washington": "WA",
    "West Virginia": "WV",
    "Wisconsin": "WI",
    "Wyoming": "WY",
    "District of Columbia": "DC",
    "American Samoa": "AS",
    "Guam": "GU",
    "Northern Mariana Islands": "MP",
    "Puerto Rico": "PR",
    "United States Minor Outlying Islands": "UM",
    "U.S. Virgin Islands": "VI",
}

mapcrs = ccrs.LambertConformal(central_longitude = -95)
datacrs = ccrs.PlateCarree()
lat_lon_dict_file = 'lat_lon_database.txt'

def read_job_file(infile):
    data = pd.read_csv(infile)

    return data

def calc_distances(data1, data2, data3, data4, threshold = 60, lat_lon_dict = None):

    # Loop over each of the UNIV/WFOs in data1
    if('University' in data1.keys()):
        base_var = 'University'
        base_data = data1
        loop_var = 'WFO'
        loop_data = data2
        
    else:
        base_var = 'WFO'
        base_data = data1
        loop_var = 'University'
        loop_data = data2


    base_list = list(base_data[base_var])
   
    count_dict = {}
 
    # Grab all the lats and lons for the base data
    # --------------------------------------------
    base_lats = np.full(len(base_list), np.nan)
    base_lons = np.full(len(base_list), np.nan)

    for ii in range(len(base_list)):
        count_dict[base_list[ii]] = 0
        xx = base_data['City'][ii] 
        yy = base_data['State'][ii]
        if(len(yy) > 2):
            ll_yy = us_state_to_abbrev[yy]
        else:
            ll_yy = yy
        city_name = xx + ' ' + ll_yy

        if(lat_lon_dict is not None):
            if(city_name in lat_lon_dict.keys()):
                base_lats[ii] = lat_lon_dict[city_name]['Lat']
                base_lons[ii] = lat_lon_dict[city_name]['Lon']
            else:
                base_lats[ii], base_lons[ii] = extract_lat_lon(city_name)
        else:
            base_lats[ii], base_lons[ii] = extract_lat_lon(city_name)
       
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # ATSCI
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

    def check_nws_distance(data, threshold, count_dict, \
            lat_lon_dict = None):

        if('University' in data.keys()):
            loop_var = 'University'
        elif('Lab/Center' in data.keys()):
            loop_var = 'Lab/Center'
        elif('Company' in data.keys()):
            loop_var = 'Company'
            
 
        loop_list = list(data[loop_var])
        # Grab all the lats and lons for the loop data
        # --------------------------------------------
        loop_lats = np.full(len(loop_list), np.nan)
        loop_lons = np.full(len(loop_list), np.nan)

        for ii in range(len(loop_list)):

            xx = loop_data['City'][ii] 
            yy = loop_data['State'][ii]
            if(len(yy) > 2):
                ll_yy = us_state_to_abbrev[yy]
            else:
                ll_yy = yy
            city_name = xx + ' ' + ll_yy

            if(lat_lon_dict is not None):
                if(city_name in lat_lon_dict.keys()):
                    loop_lats[ii] = lat_lon_dict[city_name]['Lat']
                    loop_lons[ii] = lat_lon_dict[city_name]['Lon']
                    if(loop_list[ii] == 'NASA Goddard'):
                        print("GODDARD LAT/LON",loop_lats[ii], loop_lons[ii])
                    if(loop_list[ii] == 'Riverside Technology'):
                        print("RIVERSIDE LAT/LON",loop_lats[ii], loop_lons[ii])
                else:
                    loop_lats[ii], loop_lons[ii] = extract_lat_lon(city_name)
            else:
                loop_lats[ii], loop_lons[ii] = extract_lat_lon(city_name)

        for ii in range(len(base_list)):
            # Allocate an array to hold the distances for each of the other
            # -------------------------------------------------------------
        
            # Check the ATSCI distances
            # -------------------------
            dist_array = np.full(len(loop_list), np.nan)
            dists = np.array([find_distance_between_points(base_lats[ii], base_lons[ii], 
                llat, llon) for llat, llon in zip(loop_lats, loop_lons)])
            num_within = len(np.where(dists < threshold)[0])
            count_dict[base_list[ii]] += num_within
            if(base_data['City'][ii] == 'Indianapolis'):
                if(num_within > 0):
                    print('Indianapolis is close to',np.array(loop_list)[np.where(dists < threshold)[0]])

        return count_dict

    count_dict = check_nws_distance(data2, threshold, count_dict, \
        lat_lon_dict = lat_lon_dict)
    count_dict = check_nws_distance(data3, threshold, count_dict, \
        lat_lon_dict = lat_lon_dict)
    count_dict = check_nws_distance(data4, threshold, count_dict, \
        lat_lon_dict = lat_lon_dict)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # LABS
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 
    ##!#loop_list = list(data3['Lab/Center'])
    ##!## Grab all the lats and lons for the loop data
    ##!## --------------------------------------------
    ##!#loop_lats = np.full(len(loop_list), np.nan)
    ##!#loop_lons = np.full(len(loop_list), np.nan)

    ##!#for ii in range(len(loop_list)):

    ##!#    xx = loop_data['City'][ii] 
    ##!#    yy = loop_data['State'][ii]
    ##!#    if(len(yy) > 2):
    ##!#        ll_yy = us_state_to_abbrev[yy]
    ##!#    else:
    ##!#        ll_yy = yy
    ##!#    city_name = xx + ' ' + ll_yy

    ##!#    if(lat_lon_dict is not None):
    ##!#        if(city_name in lat_lon_dict.keys()):
    ##!#            loop_lats[ii] = lat_lon_dict[city_name]['Lat']
    ##!#            loop_lons[ii] = lat_lon_dict[city_name]['Lon']
    ##!#        else:
    ##!#            loop_lats[ii], loop_lons[ii] = extract_lat_lon(city_name)
    ##!#    else:
    ##!#        loop_lats[ii], loop_lons[ii] = extract_lat_lon(city_name)
 
    ##!#for ii in range(len(base_list)):
    ##!#    # Allocate an array to hold the distances for each of the other
    ##!#    # -------------------------------------------------------------
    ##!#
    ##!#    # Check the ATSCI distances
    ##!#    # -------------------------
    ##!#    dist_array = np.full(len(loop_list), np.nan)
    ##!#    dists = np.array([find_distance_between_points(base_lats[ii], base_lons[ii], 
    ##!#        llat, llon) for llat, llon in zip(loop_lats, loop_lons)])
    ##!#    num_within = len(np.where(dists < threshold)[0])
    ##!#    count_dict[base_list[ii]] += num_within
    ##!#    #if(base_data['City'][ii] == 'Goodland'):
    ##!#    #    if(num_within > 0):
    ##!#    #        print('Goodland is close to',np.array(loop_list)[np.where(dists < threshold)[0]])

    ##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    ##!##
    ##!## PRIVATE SECTOR
    ##!##
    ##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 
    ##!#loop_list = list(data4['Company'])
    ##!## Grab all the lats and lons for the loop data
    ##!## --------------------------------------------
    ##!#loop_lats = np.full(len(loop_list), np.nan)
    ##!#loop_lons = np.full(len(loop_list), np.nan)

    ##!#for ii in range(len(loop_list)):

    ##!#    xx = loop_data['City'][ii] 
    ##!#    yy = loop_data['State'][ii]
    ##!#    if(len(yy) > 2):
    ##!#        ll_yy = us_state_to_abbrev[yy]
    ##!#    else:
    ##!#        ll_yy = yy
    ##!#    city_name = xx + ' ' + ll_yy

    ##!#    if(lat_lon_dict is not None):
    ##!#        if(city_name in lat_lon_dict.keys()):
    ##!#            loop_lats[ii] = lat_lon_dict[city_name]['Lat']
    ##!#            loop_lons[ii] = lat_lon_dict[city_name]['Lon']
    ##!#        else:
    ##!#            loop_lats[ii], loop_lons[ii] = extract_lat_lon(city_name)
    ##!#    else:
    ##!#        loop_lats[ii], loop_lons[ii] = extract_lat_lon(city_name)
 
    ##!#for ii in range(len(base_list)):
    ##!#    # Allocate an array to hold the distances for each of the other
    ##!#    # -------------------------------------------------------------
    ##!#
    ##!#    # Check the ATSCI distances
    ##!#    # -------------------------
    ##!#    dist_array = np.full(len(loop_list), np.nan)
    ##!#    dists = np.array([find_distance_between_points(base_lats[ii], base_lons[ii], 
    ##!#        llat, llon) for llat, llon in zip(loop_lats, loop_lons)])
    ##!#    num_within = len(np.where(dists < threshold)[0])
    ##!#    count_dict[base_list[ii]] += num_within
    ##!#    #if(base_data['City'][ii] == 'Goodland'):
    ##!#    #    if(num_within > 0):
    ##!#    #        print('Goodland is close to',np.array(loop_list)[np.where(dists < threshold)[0]])


    for ii in range(len(base_list)):
        print(count_dict[base_list[ii]], base_list[ii])
 
    #find_distance_between_points(lat1, lon1, lat2, lon2)

def extract_lat_lon(city_name):

    geolocator = Nominatim(user_agent = 'myapplication')
    location = geolocator.geocode(city_name)
    llat = location.latitude
    llon = location.longitude

    return llat, llon

def convert_lat_lon_to_dict(atsci, nws, labs, priv, \
        outname = lat_lon_dict_file):

    print("Converting lats and lons")
   
    def convert_local_data(data, ldict):
        list_cities = list(data['City'])
        list_states = list(data['State'])
        for xx, yy in zip(list_cities, list_states):
            if(len(yy) > 2):
                ll_yy = us_state_to_abbrev[yy]
            else:
                ll_yy = yy
            city_name = xx + ' ' + ll_yy
            llat, llon = extract_lat_lon(city_name)

            if(city_name not in lat_lon_dict.keys()):
                print(city_name,'ADDED')
                lat_lon_dict[city_name] = {}
                lat_lon_dict[city_name]['Lat'] = llat
                lat_lon_dict[city_name]['Lon'] = llon
            else:
                print(city_name,'ALREADY IN DICT')
    
        return ldict
 
    lat_lon_dict = {}
    
    # ATSCI
    # -----
    lat_lon_dict = convert_local_data(atsci, lat_lon_dict)
    lat_lon_dict = convert_local_data(nws, lat_lon_dict)
    lat_lon_dict = convert_local_data(labs, lat_lon_dict)
    lat_lon_dict = convert_local_data(priv, lat_lon_dict)
    
    return lat_lon_dict
    
    
    
def write_lat_lon_dict_to_file(lat_lon_dict, \
        outname = lat_lon_dict_file):

    with open(outname,'w') as fout:
        json.dump(lat_lon_dict, fout, indent = 4, \
            sort_keys = True)
        

def read_lat_lon_from_csv(infile = lat_lon_dict_file):
    print("Reading lats and lons")

    with open(infile,'r') as fin:
        lat_lon_dict = json.load(fin)

    return lat_lon_dict

color_dict = {
    'University': 'tab:blue',
    'WFO': 'tab:red',
    'Lab/Center': 'tab:orange',
    'Company': 'tab:cyan',
}   

size_dict = {
    'University': 12,
    'WFO': 9,
    'Lab/Center': 6,
    'Company': 4,
}

def plot_job_sites_on_map(ax, data, size = None, color = None, \
        lat_lon_dict = None, draw_circles = False, radius = 60):

    list_cities = list(data['City'])
    list_states = list(data['State'])
    dtype =  data.keys()[0]
    print("======= Plotting",dtype,"========")

    if(size is None):
        point_size = size_dict[dtype]
    else:
        point_size = size

    if(color is None):
        color_val  = color_dict[dtype]
    else:
        color_val  = color

    if(draw_circles): radius = radius * 1e3

    for xx, yy in zip(list_cities, list_states):
        if(len(yy) > 2):
            ll_yy = us_state_to_abbrev[yy]
        else:
            ll_yy = yy
        city_name = xx + ' ' + ll_yy
        #print(city_name)

        if(lat_lon_dict is not None):
            if(city_name in lat_lon_dict.keys()):
                llat = lat_lon_dict[city_name]['Lat']
                llon = lat_lon_dict[city_name]['Lon']
            else:
                llat, llon = extract_lat_lon(city_name)
        else:
            llat, llon = extract_lat_lon(city_name)

        plot_point_on_map(ax, llat, llon, markersize = point_size, \
            color = color_val, alpha = 1.0)
        if(draw_circles):
            circle_points = cgeodesic.Geodesic().circle(lon = llon, lat = llat, \
                radius = radius, endpoint = False)
            geom = shapely.geometry.Polygon(circle_points)
            ax.add_geometries((geom,), crs = datacrs, facecolor = 'red', \
                edgecolor = 'red', linewidth = 1, alpha = 0.25)
            
    ax.set_extent([-125, -70, 23, 50], datacrs)
    ax.coastlines()
    ax.add_feature(cfeature.STATES)

def plot_all_sites(atsci, nws, labs, priv, ax = None, lat_lon_dict = None, \
        draw_nws_circles = True, radius = 60):

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    plot_job_sites_on_map(ax, atsci, lat_lon_dict = lat_lon_dict)
    plot_job_sites_on_map(ax, nws, lat_lon_dict = lat_lon_dict, draw_circles = draw_nws_circles, radius = radius)
    plot_job_sites_on_map(ax, labs, lat_lon_dict = lat_lon_dict)
    plot_job_sites_on_map(ax, priv, lat_lon_dict = lat_lon_dict)

    ##!#print("Plotting university sites")
    ##!#list_univs = list(atsci['University'])
    ##!#list_cities = list(atsci['City'])
    ##!#list_states = list(atsci['State'])

    ##!#for xx, yy in zip(list_cities, list_states):
    ##!#    city_name = xx + ' ' + yy
    ##!#    print(city_name)
    ##!#    llat, llon = extract_lat_lon(city_name)

    ##!#    plot_point_on_map(ax, llat, llon, markersize = 7, color = 'tab:blue', alpha = 0.5)

    ##!#print("Plotting NWS sites")
    ##!#list_univs = list(nws['WFO'])
    ##!#list_cities = list(nws['City'])
    ##!#list_states = list(nws['State'])

    ##!#for xx, yy in zip(list_cities, list_states):
    ##!#    city_name = xx + ' ' + yy
    ##!#    print(city_name)
    ##!#    llat, llon = extract_lat_lon(city_name)

    ##!#    plot_point_on_map(ax, llat, llon, markersize = 7, color = 'tab:red', alpha = 0.5)

    ##!#print("Plotting Lab/Center sites")
    ##!#list_labs  = list(labs['Lab/Center'])
    ##!#list_cities = list(labs['City'])
    ##!#list_states = list(labs['State'])

    ##!#for xx, yy in zip(list_cities, list_states):
    ##!#    city_name = xx + ' ' + yy
    ##!#    print(city_name)
    ##!#    llat, llon = extract_lat_lon(city_name)

    ##!#    plot_point_on_map(ax, llat, llon, markersize = 7, color = 'tab:orange', alpha = 0.5)

    ##!#print("Plotting Companies sites")
    ##!#list_labs  = list(priv['Company'])
    ##!#list_cities = list(priv['City'])
    ##!#list_states = list(priv['State'])

    ##!#for xx, yy in zip(list_cities, list_states):
    ##!#    city_name = xx + ' ' + yy
    ##!#    print(city_name)
    ##!#    llat, llon = extract_lat_lon(city_name)

    ##!#    plot_point_on_map(ax, llat, llon, markersize = 4, color = 'tab:green', alpha = 0.5)



    ax.set_extent([-130, -65, 23, 50], datacrs)
    ax.coastlines()
    ax.add_feature(cfeature.STATES)
    if(not in_ax):
        fig.tight_layout()
        plt.show()
