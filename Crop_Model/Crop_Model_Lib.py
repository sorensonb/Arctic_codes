#!/usr/bin/env python

"""


"""

import os
home_dir = os.environ['HOME']
import sys
sys.path.append(home_dir)
import numpy as np
from glob import glob
from netCDF4 import Dataset
from pyproj.transformer import TransformerGroup
import fiona
import geopandas as gpd
import pandas as pd
from python_lib import *


def convert_projections(proj1, proj2, coord1, coord2):
    trans_group = TransformerGroup(proj1, proj2)
    new_coords = trans_group.transformers[0].transform(coord1, coord2)
    return new_coords

def read_layer_names_from_gdb(gdb_file):
    layers = fiona.listlayers(gdb_file)
    return layers

def read_layer_from_gdb(gdb_file, layer_name):
    return gpd.read_file(gdb_file, layer = layer_name)

def find_mukey_from_latlon(inlat, inlon, matchfile = 'test_mukey.csv'):

    # Open the match file
    # -------------------
    df = pd.read_csv(matchfile)

    # Loop over the file, looking for the mukey that most closely matches
    # the lat/lon inputs
    # -------------------------------------------------------------------
    f_lat = 89.999
    f_lon = 179.999
    f_dist = 99999.
    f_mukey = 0

    for ii in range(len(df)):
        slat = df['Lat'][ii]
        slon = df['Lon'][ii]

        dist1 =  ((slat - f_lat) **2. + (slon - f_lon) **2.) ** 0.5
        dist2 = find_distance_between_points(slat, slon, inlat, inlon)

        if(dist2 < f_dist):
            #print("REPLACING VALUE", ii, slat, slon)
            f_lat = slat
            f_lon = slon
            f_mukey = df['Mukey'][ii]
            f_dist = dist2

    print("Lat = ", f_lat, " Lon = ", f_lon, " Mukey = ", f_mukey, " Dist = ", f_dist)
        
