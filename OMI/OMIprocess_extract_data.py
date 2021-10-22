"""
 *	Program: OMIprocess-Climo Extract Data
 *
 *	Programmer: Shawn L. Jaker
 *				University of North Dakota
 *				Atmospheric Science
 *
 *	Purpose:
 *
 *			File handling (supporting) module for OMIprocess-Climo
 *			
 *			SEE OMIprocess.py FOR APPLICATION README
 *
 *	Modifications:
 *		SLJ - 2021/06 Add removal of identified bad lines.
 *
"""
# he5 Filename format
## OMI-Aura_L2-OMAERUV_2020m0603t2333-o84502_v003-2020m0605t042557.he5

# csv filename format
## OMI-Aura_L2-OMAERUV_2020m0603t2333-o84502_v003-2020m0605t042557.csv

ROW_ANOMALY_FILE = r"row_anomaly_dates_short.txt"

fn_format = r"OMI-Aura_L2-OMAERUV_????m????t????-o?????_v003-????m????t??????"
re_format = r"OMI-Aura_L2-OMAERUV_([0-9]{4})m([0-9]{2})([0-9]{2})t([0-9]{2})([0-9]{2})-o([0-9]{5})_v003-([0-9]{4})m([0-9]{2})([0-9]{2})t([0-9]{2})([0-9]{2})([0-9]{2})"

import h5py as h5
import pandas as pd
import datetime as dt
import fnmatch as fm
import numpy as np
import os
import re
import subprocess
import random
from pathlib import Path
import csv

from os import listdir
from os.path import isfile, join

# Using a 21600, 43200 grid
#from global_land_mask import globe

def check_data(date_range, location, fn_ext):
	# find all files
	
	dates = dict()
	
	try:
		
		for file in sorted(os.listdir(location)):
			if fm.fnmatch(file, fn_format + fn_ext):
				
				m = re.search(re_format, file + fn_ext)
				
				dt_start = dt.datetime( year   = int(m.group(1)),
										month  = int(m.group(2)),
										day    = int(m.group(3)),
										hour   = int(m.group(4)),
										minute = int(m.group(5)))
				
				orbt = int(m.group(6))
				
				dt_end = dt.datetime(   year   = int(m.group(7)),
										month  = int(m.group(8)),
										day    = int(m.group(9)),
										hour   = int(m.group(10)),
										minute = int(m.group(11)),
										second = int(m.group(12)))
				
				if date_range[0] <= dt_start and dt_start <= date_range[1]:
					dates[orbt] = (dt_start, dt_end)
	except:
		
		print("problem checking data for " + fn_ext + " in " + location)
		pass
			
	return dates

def check_he5_data(date_range):
	# find all he5 files
	
	return check_data(date_range, LOCATION_HE5, ".he5")

def check_csv_data(date_range):
	# find all csv files
	
	return check_data(date_range, LOCATION_CSV, ".csv")

def check_orbits(dates_he5, dates_csv):
	
	for orbit in sorted(dates_csv):
		del dates_he5[orbit]
	
	return dates_he5

def build_filename(location, orbit, date_entry, fn_ext):

	dt_s = date_entry[0]
	dt_e = date_entry[1]

	filename = location \
		+ r"OMI-Aura_L2-OMAERUV_" \
		+ str("%04d" % (dt_s.year)) \
		+ "m" \
		+ str("%02d" % (dt_s.month)) \
		+ str("%02d" % (dt_s.day)) \
		+ "t" \
		+ str("%02d" % (dt_s.hour)) \
		+ str("%02d" % (dt_s.minute)) \
		+ "-o" \
		+ str("%05d" % (orbit)) \
		+ "_v003-" \
		+ str("%04d" % (dt_e.year)) \
		+ "m" \
		+ str("%02d" % (dt_e.month)) \
		+ str("%02d" % (dt_e.day)) \
		+ "t" \
		+ str("%02d" % (dt_e.hour)) \
		+ str("%02d" % (dt_e.minute)) \
		+ str("%02d" % (dt_e.second)) \
		+ fn_ext
			 
	return filename


def extract_variable(path_he5, filename_he5, data_path):
	print("extracting " + path_he5 + filename_he5)
	
	f = h5.File(path_he5 + filename_he5, 'r')
	
	data_field = f[data_path][:]
	data_xfla = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:]
	data_uvai = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:]
	data_lati = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:]
	
	data_field = data_field[(data_xfla == 0) & (data_uvai > -10) & (data_uvai < 10) & (data_lati > 55)]
	
	elems = np.prod(data_field.shape)
	
	data_field = data_field.reshape((elems))
	
	return data_field

def load_bad_tracks(filename = ROW_ANOMALY_FILE):
	bad_tracks = dict()
	
	with open(filename) as f:
		file_lines = f.readlines()
		
		for fline in file_lines:
			file_line = fline.strip()
			m = re.match("([0-9]{8})(.*)", file_line)
			if m and len(m.groups()) > 1 and m[2] != '':
				date = int(m[1])
				tracks = str(m[2]).split()
				tracks = list([int(a) for a in tracks])
				bad_tracks[date] = tracks
				#print(date, bad_tracks[date])
	
	return bad_tracks
	
def mark_bad_tracks(data_trak, data_xfla, bad_tracks, date):
	
	if date in bad_tracks:
		for track in bad_tracks[date]:
			#print(date, track, np.count_nonzero(data_trak == track))
			data_xfla = np.where(data_trak == track - 1, -1, data_xfla)

		print("removed", np.count_nonzero(data_xfla == -1), "from", date)

	return data_xfla

def extract_data(dates):

	bad_tracks = load_bad_tracks()
	
	for orbit in sorted(dates):
	
		date_entry = dates[orbit]
		# build the h5 and csv filenames
		filename_he5 = build_filename(LOCATION_HE5, orbit, date_entry, ".he5")
		filename_csv = build_filename(LOCATION_CSV, orbit, date_entry, ".csv")
		
		str_date = str("%04d%02d%02d" % (date_entry[0].year, date_entry[0].month, date_entry[0].day))
		int_date = int(str_date)
		
		print("extracting " + filename_he5 + " " + str(int_date))
		
		# open the h5 file
		f = h5.File(filename_he5, 'r')

		# get locations
		data_lati = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:]
		data_long = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:]
		
		# get the needed dependent parameter
		data_uvai = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:]
		
		# get angles
		data_szea = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle'][:]
		data_vzea = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/ViewingZenithAngle'][:]
		data_raza = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle'][:]

		# get time (using seconds in day)
		data_time = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SecondsInDay'][:]
		
		# get height (terrain pressure)
		data_pres = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/TerrainPressure'][:]
		
		# get surface albedo, norm radiance, and final AOD in 3 groups
		data_sal0 = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo'][:,:,0]
		data_sal1 = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo'][:,:,1]

		data_nra0 = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/NormRadiance'][:,:,0]
		data_nra1 = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/NormRadiance'][:,:,1]
		
		#data_fao0 = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/FinalAerosolOpticalDepth'][:,:,0]
		#data_fao1 = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/FinalAerosolOpticalDepth'][:,:,1]
		
		#data_sal2 = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo'][:,:,2]
		#data_nra2 = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/NormRadiance'][:,:,2]
		#data_fao2 = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/FinalAerosolOpticalDepth'][:,:,2]
		
		# get cloud fraction		
		data_clou = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/CloudFraction'][:]
		#data_refl = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity'][:]

		# get the needed data quality flags
		data_xfla = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:]
		data_gfla = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags'][:]
		#data_ffla = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/FinalAlgorithmFlags'][:]
		#data_afla = f['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/AlgorithmFlags_AerosolIndex'][:]

		#data_is_land  = np.full((data_lati.shape), False)
		#data_is_ocean = np.full((data_lati.shape), False)
		#data_is_coast = np.full((data_long.shape), False)
		"""
		for i in range(data_lati.shape[0]):
			#print (i)
			for j in range(data_lati.shape[1]):
				#print("\t", j)
				if data_lati[i,j].astype(float) <= 90  and data_lati[i,j].astype(float) >= -90 and \
				   data_long[i,j].astype(float) <= 180 and data_long[i,j].astype(float) >= -180:
					if globe.is_land (data_lati[i,j].astype(float), data_long[i,j].astype(float)):
						data_is_land [i,j] = True
						data_is_ocean[i,j] = False
					else:
						data_is_land [i,j] = False
						data_is_ocean[i,j] = True
				else:
					print("Bad data at ", i, j, data_lati[i,j].astype(float), data_long[i,j].astype(float))
					

		for i in range(data_lati.shape[0]):
			#print (i)
			for j in range(data_lati.shape[1]):
				#print("\t", j)
				
				for k in range (-2, 3):
					for l in range (-1, 2):
						if (k == 0 and l == 0) or (i + k < 0) or (i + k >= data_lati.shape[0]) or (j + l < 0) or (j + l >= data_lati.shape[1]): continue

						if data_is_land[i,j] and data_is_ocean[i+k,j+l]:
							data_is_coast[i,j] = True
							data_is_coast[i+k,j+l] = True
		"""
		# reshape the data to 1-D and columize
		#print(data_uvai.shape)
		elems = np.prod(data_uvai.shape)
		
		#print(elems)
		data_clas = np.full((elems), 1)			# V1

		data_uvai = data_uvai.reshape((elems))	# V2
		
		data_lati = data_lati.reshape((elems))	# V3
		data_long = data_long.reshape((elems))	# V4
	
		data_szea = data_szea.reshape((elems))	# V5
		data_vzea = data_vzea.reshape((elems))	# V6
		data_raza = data_raza.reshape((elems))	# V7
	
		data_sal0 = data_sal0.reshape((elems))	# V8
		data_sal1 = data_sal1.reshape((elems))	# V9

		data_nra0 = data_nra0.reshape((elems))	# V10
		data_nra1 = data_nra1.reshape((elems))	# V11

		#data_fao0 = data_fao0.reshape((elems))	# 
		#data_fao1 = data_fao1.reshape((elems))	# 

		data_clou = data_clou.reshape((elems))	# V12
		
		data_jday = np.full((elems), int(date_entry[0].strftime('%j')))			# V13
		
		data_time = np.repeat(data_time, 60)	# V14
		
		data_orbt = np.full((elems), orbit)		# V15
		
		data_indx = np.arange(0, elems)			# V16
		
		data_date = np.full((elems), int_date)	# V17
		
		data_trak = data_indx % 60				# V18
		
		data_pres = data_pres.reshape((elems))	# V19
		
		data_gfla = data_gfla.reshape((elems))	# V20
		
		data_xfla = data_xfla.reshape((elems))	# V21
		
		#data_is_land = data_is_land.reshape((elems))			# V22
		
		#data_is_ocean = data_is_ocean.reshape((elems))			# V23
		
		#data_is_coast = data_is_coast.reshape((elems))			# V24
		
		#print("data_indx")
		#print(data_indx)
		
		#data_clas = data_lati.astype(int)
		
		# Remove negative clouds before stacking
		data_clou = np.where(data_clou > 0, data_clou, 0)	

		data_xfla = mark_bad_tracks(data_trak, data_xfla, bad_tracks, int_date)
		
		data = np.stack((data_clas, \
						 data_uvai, \
						 data_lati, \
						 data_long, \
						 data_szea, \
						 data_vzea, \
						 data_raza, \
						 data_sal0, \
						 data_sal1, \
						 data_nra0, \
						 data_nra1, \
						 data_clou, \
						 data_jday, \
						 data_time, \
						 data_orbt, \
						 data_indx, \
						 data_date, \
						 data_trak, \
						 #data_is_land, \
						 #data_is_ocean, \
						 #data_is_coast, \
						 data_pres, \
						 data_gfla,
						 data_xfla), axis=1)
						 
		#print (data)
		#print (data.shape)
		
		# Screen data

		#		(data_xfla == 0) & 
		
		data = data[(data_uvai > -10) & (data_uvai < 10) & (data_lati > 55)] #  & (data_pres > 600) & (data_sal0 > 0.07)
		
		# Screen any missing values
		data = np.where(data < -1.267E+030, np.nan, data)
		
		#print(data)
		#print(data.shape)
		
		# save data to the csv file
		if len(data) > 0:
			np.savetxt(filename_csv, data, delimiter=',', fmt=r' %.8f')
		else:
			print("no data for " + filename_csv)
		
		# close all the files
		f.close()
	
def load_csv(orbit, date_entry):
	
	# biuld the csv filename
	filename_csv = build_filename(LOCATION_CSV, orbit, date_entry, ".csv")
	print("loading " + filename_csv);

	result = pd.read_csv(filename_csv, header=None)
	#print(result)
	
	return result

def load_data(date_range, settings):

	global IMAGE_DIR 
	global LOCATION_HE5
	global LOCATION_CSV
	global LOCATION_OUT

	IMAGE_DIR = settings["IMAGE_DIR"]
	LOCATION_HE5 = settings["LOCATION_HE5"]
	LOCATION_CSV = settings["LOCATION_CSV"]
	LOCATION_OUT = settings["LOCATION_OUT"]
	
	# ****
	# Check dates of the he5 files that are in the range
	# ****
	dates_he5 = check_he5_data(date_range)
	#print("dates_he5")
	#print(dates_he5)

	# ****
	# Check dates of the csv files that are in the range
	# ****
	dates_csv = check_csv_data(date_range)
	#print("dates_csv")
	#print(dates_csv)
	
	# ****
	# Compare the lists and remove csv's from he5's list (using orbit #)
	# ****
	dates_missing = check_orbits(dates_he5, dates_csv)
	#prinst("dates_missing")
	#print(dates_missing)

	# ****
	# Extract training data that is not available as csv (yet)
	# ****
	extract_data(dates_missing)
	
	# ****
	# Recheck dates of the csv files that are in the range
	# ****
	dates_csv = check_csv_data(date_range)
	#print("dates_csv")
	#print(dates_csv)
	
	# ****
	# Load and return the csv data
	# ****

	dt_s = date_range[0]
	dt_e = date_range[1]
	
	temp_dir = "/dev/shm/tmp." \
		+ str("%04d" % (dt_s.year)) \
		+ r"m" \
		+ str("%02d" % (dt_s.month)) \
		+ str("%02d" % (dt_s.day)) \
		+ r"_" \
		+ str("%04d" % (dt_e.year)) \
		+ r"m" \
		+ str("%02d" % (dt_e.month)) \
		+ str("%02d" % (dt_e.day)) \
		+ r"." \
		+ str("%06d/" % random.randint(0,100000))
	
	Path(temp_dir).mkdir(parents=True, exist_ok=True)
	
	print("creating all.csv");
	subprocess.call(["rm -f " + temp_dir + "all.csv"], shell=True)
	for orbit in sorted(dates_csv):
		filename_csv = build_filename(LOCATION_CSV, orbit, dates_csv[orbit], ".csv")
		print("concatenating " + filename_csv + " to " + temp_dir + "all.csv")
		subprocess.call(["cat " + filename_csv + " >> " + temp_dir + "all.csv"], shell=True)
		
	print("loading all.csv");
	data = pd.read_csv(temp_dir + "all.csv", header=None)
		
	print("done loading data, removing temp file")

	Path(temp_dir + "all.csv").unlink()

	Path(temp_dir).rmdir()
		
	print("done cleaning up.")
	
	#print(data.shape)
	return data
	
def save_data(orb, data_type, data, data_set, settings):

	global IMAGE_DIR 
	global LOCATION_HE5
	global LOCATION_CSV
	global LOCATION_OUT

	IMAGE_DIR = settings["IMAGE_DIR"]
	LOCATION_HE5 = settings["LOCATION_HE5"]
	LOCATION_CSV = settings["LOCATION_CSV"]
	LOCATION_OUT = settings["LOCATION_OUT"]
	
	data_filter = data_set.V2 - data
	
	data_out = pd.DataFrame({"V3":data_set.V3, \
	                         "V4":data_set.V4, \
	                         "V2":data_set.V2, \
	                         "FIL":data_filter, \
	                         data_type:data, \
	                         "V14":data_set.V14, \
	                         "V5":data_set.V5, \
	                         "V6":data_set.V6, \
	                         "V7":data_set.V7, \
	                         "V8":data_set.V8, \
	                         "V9":data_set.V9, \
	                         "V10":data_set.V10, \
	                         "V11":data_set.V11, \
	                         "V12":data_set.V12, \
	                         "V21":data_set.Gflag})
	                         
	#print("data_out:\n", data_out)
	                         
	filename = str(LOCATION_OUT + str(orb) + r'_' + data_type + r'.txt')
	
	np.savetxt(filename, data_out.values, fmt=r'%16.8f')
	
	return filename
