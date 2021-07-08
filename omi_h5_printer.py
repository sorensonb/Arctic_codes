#!/usr/bin/env python

"""
Prints off the lat, lon, and GQPF flag values for an OMI HDF5 file

"""

import sys
import h5py

def get_whole_flags(value):
    return format(value,"016b")

def get_ice_flags(value):
    return int(format(value,"016b")[-15:-8],2)

if(len(sys.argv) != 2):
    print("SYNTAX: omi_h5_printer.py omi_h5_file")
    sys.exit()

data = h5py.File(sys.argv[1])

LAT = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:].flatten()
LON = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:].flatten()
GPQF = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:].flatten()
#GPQF = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags'][:,:].flatten()

for x, y, z in zip(LAT,LON,GPQF):
    print(int(x*1000),'\t',z)
    #print(int(x*1000),'\t',z,'\t',get_whole_flags(z),'\t',get_ice_flags(z),'\t',format(z,"016b")[-15:-8])

data.close()
