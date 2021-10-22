"""
 *	Program: OMIprocess-Climo
 *
 *	Programmer: Shawn L. Jaker
 *				University of North Dakota
 *				Atmospheric Science
 *
 *	Purpose:
 *
 *			This program is used to compute and remove the average
 * 			OMI UV aerosol index climatology from northern artic regions
 *			using the viewing zenith angle, solar zenith angle, relative
 *			azimuth angle, and solar albedo.
 *			
 *			Two methods are employed within this code to produce 
 *			comparable results. The method is to bin the
 *			data based on 2.5 degree bins for zenith angles, 1°
 *			degree bins for the azimuth angle, and 0.05 albedo.
 *
 *	Usage command *(SEE NOTE AT END OF README SECTION):
 *
 *			The program uses a single, or multiple training period(s) 
 *			and a reconstruction period, both of which are specified on the
 *			command line, as follows:
 *
 *			Running WITHOUT the dates will use HARD CODED dates of:
 *				date_range_learn = (dt.datetime(2007, 4, 7, 11, 4), dt.datetime(2007, 4, 7, 11, 4))	# 2007/04/07 @ 11:04
 *				date_range_recon = (dt.datetime(2008, 4, 7, 11, 9), dt.datetime(2008, 4, 7, 11, 9))	# 2008/04/07 @ 11:09
 *			
 *			dates in format: YYYYMMDD
 *
 *		Using multiple training periods
 *			$ python3 OMIprocss.py {{{{...} date_start_T, date_end_T} date_start_T, date_end_T}, date_start_R, date_end_R}
 *
 *		Using single training periods
 *			$ python3 OMIprocss.py {date_start_T, date_end_T, date_start_R, date_end_R}
 *
 *		Using hard-coded dates (for ease of debugging)
 *			$ python3 OMIprocss.py
 *
 *	Example command:
 *
 *			$ python3 OMIprocess.py 20050401 20050430 20060401 20060430 20070401 20070430 20080401 20080430 20090401 20090430 20080401 20080430
 *			$ python3 OMIprocess.py 20050401 20050430 20060401 20060430 20070401 20070430 20080401 20080430 20090401 20090430 20080401 20080430
 *			$ python3 OMIprocess.py 20070401 20070430 20080401 20080430
 *			$ python3 OMIprocess.py
 *
 *
 *	Returns:
 *			
 *			Based on configuration, the program will
 *
 *			> create appropriate output directories (if needed)
 *			> check for previsously extracted data
 *			> extract missing training and reconstruction data to CSV and load
 *			> compute climatology from training data and output to climo.sj
 *			> output reconstruction data with the climatology removed
 *			> output reconstruction data using linear regression
 *			> output plots showing the original data, filter values, and reconstructed data
 *
 *	POST PROCESSING FOR IDL FIGURES:
 *
 *		run idl
 *		$ idl
 *
 *		load plotmaps module
 *		> .r plotmap
 *
 *		run script file (to avoid reloading idl for each figure)
 *		> @plotmap.list
 *
 *		run individual file
 *		> plotmap,"/home/user/OMIprocess-Climo/out_files/19837_CLI.txt"
 *
 *		* note that the about command uses the orbit number, use actual file name, of course
 *
 *	__BETTER__ Post processing code (Any one of these will work, depending on desired output)
 *		Each orbit and the mean of all orbits will be output to the 'figs' directory.
 *		
 *		$ python3 plot_single_OMI.py 202004010038 UVAerosolIndex 60
 *		$ python3 plot_single_OMI.py 2020040112 UVAerosolIndex 60
 *		$ python3 plot_single_OMI.py 20200401 UVAerosolIndex 60
 *		$ python3 plot_single_OMI.py 202004 UVAerosolIndex 60
 *
 *	Note:
 *			
 *			If the settings.dct file does not exist, the first time the
 *			program is run, it will be created. The settings file can then be 
 *			changed to suit the user's needs. On the subsequent run, the
 *			needed directories will be created, if needed (except for H5).
 *
 *	Modifications:
 *		SLJ - 2021/06 Remove superfluous code segments and large debug sections.
 *		SLJ - 2021/06 Add removal of identified bad lines.
 *		SLJ - 2021/06 Remove modificaiton above and put into extraction code.
 *		SLJ - 2021/07 Modified to save and load climo to avoid re-binning
 *		SLJ - 2021/07 Modified to reduce run-time
 *
"""


"""
	Configurable Options
"""

USE_EXISTING_CLIMO_FILE = True

# Screen pre-identified rows based on date
SCREEN_XFLAG = True

# Use veiwing and solar geometry to create bins
USE_GEOMETRY = True

# Use ground quality (type) flag to create bins
USE_GROUNDFLAG = True

# Output each orbit when calculations are complete
DO_ORBIT_OUTPUT = True

"""
	End of Configurable Options
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.path as mpath
import numpy.ma as ma

import pickle

import os
import json
from pathlib import Path

import sys
import datetime as dt

import subprocess

import pandas as pd
import cartopy.crs as ccrs
"""
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from scipy import stats
import statsmodels.api as sm
"""
import OMIprocess_extract_data as ed

import regex as re

ROWS = 60 # Actually, these are columns, but the method describes them as 'rows'

szea_bins=np.linspace(0,90,37)		# 37 bin edges >> 36 bins >>  90 / 36 =  2.5 degree bins
vzea_bins=np.linspace(0,90,37)		# 37 bin edges >> 36 bins >>  90 / 36 =  2.5 degree bins
raza_bins=np.linspace(-1,181,92)	# 92 bin edges >> 91 bins >> 182 / 91 =  2.0 degree bins
salb_bins=np.linspace(0.0,1.0,21)	# 21 bin edges >> 20 bins >> 1.0 / 20 = 0.05 albedo bins
flag_bins=np.linspace(0,7,8)		#  8 bin edges >>  7 bins

def extract_flags(byte):
	mask_land_water = 0x000f
	mask_other_flag = 0x00f0
	mask_snow_ice   = 0x7f00
	mask_neighbor_f = 0x8000
	
	# print(byte & mask_land_water, (byte & mask_other_flag) >> 4, (byte & mask_snow_ice) >> 8, (byte & mask_neighbor_f) >> 15)
	
	return byte & mask_land_water, (byte & mask_other_flag) >> 4, (byte & mask_snow_ice) >> 8, (byte & mask_neighbor_f) >> 15
	
def get_class(byte):
	
	flags = extract_flags(byte)
	
	#is_ocean   = 1 if (flags[0] == 0 or flags[0] == 6 or flags[0] == 7) else 0
	#is_land    = 1 if (flags[0] == 1) else 0
	#is_coastal = 1 if (flags[0] == 3) else 0
	
	is_sfland  = 1 if (flags[2] == 0) else 0                         # snow-free land
	is_seaice  = 1 if (flags[2] >= 1 and flags[2] <= 100) else 0     # sea ice concentration
	is_permice = 1 if (flags[2] == 101) else 0                       # permanent ice
	is_drysnow = 1 if (flags[2] == 103) else 0                       # dry snow
	is_ifocean = 1 if (flags[2] == 104) else 0                       # ocean
	is_mixpixl = 1 if (flags[2] == 124) else 0                       # mixed pixel at coast line
	is_other   = 1 if (flags[2] >= 125 and flags[2] <= 127) else 0   # suspect ice, corners, error (other?)
	
	# is_other   = 1 if not (is_ocean or is_land or is_coastal or is_seaice or is_permice or is_drysnow) else 0
	
	#return is_ocean, is_land, is_coastal, is_sfland, is_seaice, is_permice, is_drysnow, is_ifocean, is_mixpixl, is_other
	
	if is_sfland:  return 1
	if is_seaice:  return 2
	if is_permice: return 3
	if is_drysnow: return 4
	if is_ifocean: return 5
	if is_mixpixl: return 6
	if is_other:   return 7

def prepare_data(data):
	print("\tPreparing variables")
	
	data.columns = ["V"+str(i) for i in range(1, len(data.columns)+1)]  # rename column names to be similar to R naming convention
	
	data.V1 = data.V2.astype(int)
	
	elems = np.prod(data.V2.shape)

	# print(data.shape)
	#print(data)
	data_lat = data.V3
	data_lon = data.V4
	data_sza = data.V5.to_numpy(copy=True)
	data_vza = data.V6.to_numpy(copy=True)
	data.V14 = data.V14 / 86400.
	data_idx = data.V16
	data_dat = [int(i) for i in data.V17]
	data_trk = [int(i) for i in data.V18]
	data_prs = data.V19
	data_gfl = [get_class(int(i)) for i in data.V20]
	data_xfl = [int(i) for i in data.V21]
	
	data_cos_sza = data.V5.to_numpy()
	data_cos_vza = data.V6.to_numpy()

	#print (elems)

	data_cos_sza = np.cos(data_sza / 180. * math.pi)
	data_cos_vza = np.cos(data_vza / 180. * math.pi)
		
	data = pd.DataFrame({'V1':data.V1, \
	                     'V2':data.V2, \
	                     'V3':data.V3, \
	                     'V4':data.V4, \
	                     'V5':data_cos_sza, \
	                     'V6':data_cos_vza, \
	                     'V7':data.V7, \
	                     'V8':data.V8, \
	                     'V9':data.V9, \
	                     'V10':data.V10, \
	                     'V11':data.V11, \
	                     'V12':data.V12, \
	                     'V14':data.V14, \
	                     'V15':data.V15, \
	                     'IDX':data_idx, \
	                     #'is_land':data.V17, \
	                     #'is_ocean':data.V18, \
	                     #'is_coast':data.V19, \
	                     'DATE':data_dat, \
	                     'ROW':data_trk, \
	                     'PRS':data_prs, \
	                     'Vszea':data_sza, \
	                     'Vvzea':data_vza, \
	                     'Gflag':data_gfl, \
	                     'Xflag':data_xfl})
	                     
	#data.set_index('V2')
	
	# print(data.shape)
	#print(data)
	
	return data
	
	
def train_(data):
	print("	Training bins")
	
	#print(szea_bins)
	#print(vzea_bins)
	#print(raza_bins)
	#print(salb_bins)
	#print(flag_bins)
	
	len_szea_bins = len(szea_bins)
	len_vzea_bins = len(vzea_bins)
	len_raza_bins = len(raza_bins)
	len_salb_bins = len(salb_bins)
	len_flag_bins = len(flag_bins)
	
	if USE_GEOMETRY and USE_GROUNDFLAG:
		bins_t = data.assign(
			bin_szea = pd.cut(data.Vszea, bins=szea_bins, labels=range(len_szea_bins - 1)),
			bin_vzea = pd.cut(data.Vvzea, bins=vzea_bins, labels=range(len_vzea_bins - 1)),
			bin_raza = pd.cut(data.V7, bins=raza_bins, labels=range(len_raza_bins - 1)),
			bin_salb = pd.cut(data.V9, bins=salb_bins, labels=range(len_salb_bins - 1)),			# Albedo channel 2
			bin_flag = pd.cut(data.Gflag, bins=flag_bins, labels=range(len_flag_bins - 1))
		)
	elif USE_GEOMETRY:
		bins_t = data.assign(
			bin_szea = pd.cut(data.Vszea, bins=szea_bins, labels=range(len_szea_bins - 1)),
			bin_vzea = pd.cut(data.Vvzea, bins=vzea_bins, labels=range(len_vzea_bins - 1)),
			bin_raza = pd.cut(data.V7, bins=raza_bins, labels=range(len_raza_bins - 1)),
			bin_salb = pd.cut(data.V9, bins=salb_bins, labels=range(len_salb_bins - 1))			# Albedo channel 2
		)
	elif USE_GROUNDFLAG:
		bins_t = data.assign(
			bin_flag = pd.cut(data.Gflag, bins=flag_bins, labels=range(len_flag_bins - 1))
		)
	
	"""
	if USE_GEOMETRY:
		print(bins_t.bin_szea)
		print(bins_t.bin_vzea)
		print(bins_t.bin_raza)
		print(bins_t.bin_salb)
	
	if USE_GROUNDFLAG:
		print(bins_t.bin_flag)
	"""
	
	print("	Calculating bin IDs")
	
	if USE_GEOMETRY:
		szea = bins_t.bin_szea.to_numpy()
		vzea = bins_t.bin_vzea.to_numpy()
		raza = bins_t.bin_raza.to_numpy()
		salb = bins_t.bin_salb.to_numpy()

	if USE_GROUNDFLAG:
		flag = bins_t.bin_flag.to_numpy()
	
	if USE_GEOMETRY and USE_GROUNDFLAG:
		dim1 = len_salb_bins
		dim2 = len_salb_bins * len_raza_bins
		dim3 = len_salb_bins * len_raza_bins * len_vzea_bins
		dim4 = len_salb_bins * len_raza_bins * len_vzea_bins * len_szea_bins
		
		bin_id =  \
			   + salb \
			   + (raza * dim1) \
			   + (vzea * dim2) \
			   + (szea * dim3) \
			   + (flag * dim4)

	elif USE_GEOMETRY:
		bin_id =  \
			   + salb \
			   + (raza * len_salb_bins) \
			   + (vzea * len_salb_bins * len_raza_bins) \
			   + (szea * len_salb_bins * len_raza_bins * len_vzea_bins)

	elif USE_GROUNDFLAG:
		bin_id = flag
	
	#print(bin_id)
	
	print("	Assigning bin IDs")
	
	bins = bins_t.assign(
		bin_tuple = pd.Categorical(bins_t.filter(regex='bin_').apply(tuple, 1)),
		bin_ids = bin_id
	)
	
	#print(type(bins))
	#print(bins)
		
	return bins

def train_climate(data, limit_training=False):

	bins = train_(data)
	#print("BINS ARE: ")
	#print(bins.describe())
	#bin_tups = bins.bin_tuple
	
	#print(bins)
	#print(type(bin_tups))
	#print(bin_tups)

	print("\tCalcualting bin means")
	
	if USE_GEOMETRY and USE_GROUNDFLAG:
		bin_mean = np.ndarray((len(szea_bins), len(vzea_bins), len(raza_bins), len(salb_bins), len(flag_bins)), dtype=object)
	
	elif USE_GEOMETRY:
		bin_mean = np.ndarray((len(szea_bins), len(vzea_bins), len(raza_bins), len(salb_bins)), dtype=object)
	
	elif USE_GROUNDFLAG:
		bin_mean = np.ndarray((len(flag_bins)), dtype=object)
	
	
	#print(len(bin_tups))
	#print(len(bins.bin_tuple))
	
	bin_tups = bins.bin_tuple.to_numpy()
	bin_vals = bins.V2.to_numpy()
	
	
	if USE_GEOMETRY:
		bin_szea = bins.Vszea.to_numpy()
		bin_vzea = bins.Vvzea.to_numpy()
		bin_raza = bins.V7.to_numpy()
		bin_salb = bins.V9.to_numpy()
	
	if USE_GROUNDFLAG:
		bin_flag = bins.Gflag.to_numpy()
	"""
	print("	Creating climo bins")
	
	if USE_GEOMETRY and USE_GROUNDFLAG:
		for i in range(len(szea_bins)):
			for j in range(len(vzea_bins)):
				for k in range(len(raza_bins)):
					for l in range(len(salb_bins)):
						for m in range(len(flag_bins)):
							bin_mean[i,j,k,l,m] = list()
	
	elif USE_GEOMETRY:
		for i in range(len(szea_bins)):
			for j in range(len(vzea_bins)):
				for k in range(len(raza_bins)):
					for l in range(len(salb_bins)):
						bin_mean[i,j,k,l] = list()
	
	elif USE_GROUNDFLAG:
		for m in range(len(flag_bins)):
			bin_mean[m] = list()
	"""
	
	print("	Accumulating bin values")
	
	for n in range(len(bins.bin_tuple)):
		if USE_GEOMETRY and USE_GROUNDFLAG:
			(i,j,k,l,m) = bin_tups[n]
			if ((0 <= i < len(szea_bins)) and (0 <= j < len(vzea_bins)) and (0 <= k < len(raza_bins)) and (0 <= l < len(salb_bins)) and (0 <= m < len(flag_bins))):
				if bin_mean[i,j,k,l,m] is None: bin_mean[i,j,k,l,m] = list()
				bin_mean[i,j,k,l,m].append((bin_vals[n], bins.index.values[n]))
			else:
				print(i, j, k, l, m, bin_vals[n], bin_szea[n], bin_vzea[n], bin_raza[n], bin_salb[n], bin_flag[n])
		
		elif USE_GEOMETRY:
			(i,j,k,l) = bin_tups[n]
			if ((0 <= i < len(szea_bins)) and (0 <= j < len(vzea_bins)) and (0 <= k < len(raza_bins)) and (0 <= l < len(salb_bins))):
				bin_mean[i,j,k,l].append((bin_vals[n], bins.index.values[n]))
			else:
				print(i, j, k, l, bin_vals[n], bin_szea[n], bin_vzea[n], bin_raza[n], bin_salb[n])
		
		elif USE_GROUNDFLAG:
			m = bin_tups[n]
			if ((0 <= m < len(flag_bins))):
				bin_mean[m].append((bin_vals[n], bins.index.values[n]))
			else:
				print(m, bin_flag[n])
	
	print("	Calculating")
	
	bin_mean = bin_mean.reshape((len(szea_bins) * len(vzea_bins) * len(raza_bins) * len(salb_bins) * len(flag_bins)))

	all_indxs = np.zeros((len(szea_bins) * len(vzea_bins) * len(raza_bins) * len(salb_bins) * len(flag_bins)))
	
	def f(bin_i):
		if bin_i is None or len(bin_i) <= 2: return np.nan
		le = len(bin_i)
					
		bin_tups_list = sorted(bin_i, key=lambda x:x[0])
			
		if limit_training:
			bin_tups_list = bin_tups_list[int(le/2):]
		
		bin_val_list = [i[0] for i in bin_tups_list]
		bin_idx_list = [i[1] for i in bin_tups_list]
		
		#print(bin_idx_list)
		
		for x in bin_idx_list:
			all_indxs[x] = 1
		
		return np.mean(np.array(bin_val_list))
		#print(i, j, k, l, bin_mean[i,j,k,l,m])
		
	bin_mean = np.array([f(x) for x in bin_mean], dtype=object)
	
	bin_mean = bin_mean.reshape((len(szea_bins), len(vzea_bins), len(raza_bins), len(salb_bins), len(flag_bins)))
	
	count = range(len(szea_bins) * len(vzea_bins) * len(raza_bins) * len(salb_bins) * len(flag_bins))
	
	all_indxs = [i for i in count if all_indxs[i] == 1]

	"""
	elif USE_GEOMETRY and USE_GROUNDFLAG:
		for i in range(len(szea_bins)):
			for j in range(len(vzea_bins)):
				for k in range(len(raza_bins)):
					for l in range(len(salb_bins)):
						for m in range(len(flag_bins)):
							if not bin_mean[i,j,k,l,m] is None:
								le = len(bin_mean[i,j,k,l,m])
								
								if le > 2:
									bin_tups_list = sorted(bin_mean[i,j,k,l,m], key=lambda x:x[0])
										
									if limit_training:
										bin_tups_list = bin_tups_list[int(le/2):]
									
									bin_val_list = [i[0] for i in bin_tups_list]
									bin_idx_list = [i[1] for i in bin_tups_list]
									
									all_indxs.extend(bin_idx_list)
									
									bin_mean[i,j,k,l,m] = np.mean(np.array(bin_val_list))
								else:
									bin_mean[i,j,k,l,m] = np.nan
									#print(i, j, k, l, bin_mean[i,j,k,l,m])
	
	elif USE_GEOMETRY:
		for i in range(len(szea_bins)):
			for j in range(len(vzea_bins)):
				for k in range(len(raza_bins)):
					for l in range(len(salb_bins)):
					
						le = len(bin_mean[i,j,k,l])
						
						if le > 2:
							bin_tups_list = sorted(bin_mean[i,j,k,l], key=lambda x:x[0])
							
							if limit_training:
								bin_tups_list = bin_tups_list[int(le/2):]
							
							bin_val_list = [i[0] for i in bin_tups_list]
							bin_idx_list = [i[1] for i in bin_tups_list]
							
							all_indxs.extend(bin_idx_list)
							
							bin_mean[i,j,k,l] = np.mean(np.array(bin_val_list))
						else:
							bin_mean[i,j,k,l] = np.nan
						#print(i, j, k, l, bin_mean[i,j,k,l])
	
	elif USE_GROUNDFLAG:
		for m in range(len(flag_bins)):
					
			le = len(bin_mean[m])
			
			if le > 2:
				bin_tups_list = sorted(bin_mean[m], key=lambda x:x[0])
				
				if limit_training:
					bin_tups_list = bin_tups_list[int(le/2):]
				
				bin_val_list = [i[0] for i in bin_tups_list]
				bin_idx_list = [i[1] for i in bin_tups_list]
				
				all_indxs.extend(bin_idx_list)
				
				bin_mean[m] = np.mean(np.array(bin_val_list))
			else:
				bin_mean[m] = np.nan
				#print(m, bin_mean[m])
	"""
	trained_data = bins.loc[all_indxs]
	
	print("DONE TRAINING")
	
	return bin_mean, trained_data
	
def remove_climate(data, bin_mean):

	bins = train_(data)
	
	print("\tRemoving climatology")
	#print(len(bins.bin_tuple))
	
	bin_tups = bins.bin_tuple.to_numpy()
	bin_vals = bins.V2.to_numpy()
	
	for n in range(len(bins.bin_tuple)):
		if USE_GEOMETRY and USE_GROUNDFLAG:
			(i,j,k,l,m) = bin_tups[n]
			if not bin_mean[i,j,k,l,m] is None:
				bin_vals[n] = bin_vals[n] - bin_mean[i,j,k,l,m]
			
		elif USE_GEOMETRY:
			(i,j,k,l)   = bin_tups[n]
			bin_vals[n] = bin_vals[n] - bin_mean[i,j,k,l]
		
		elif USE_GROUNDFLAG:
			m = bin_tups[n]
			bin_vals[n] = bin_vals[n] - bin_mean[m]
	
	bins.V2 = bin_vals
	
	return bins

def run_data(data_t, data):
	
	#print("-----\ndata_t")
	#print(data_t.shape)
	#print(data_t)
		
	#print("-----\ndata")
	#print(data.shape)
	#print(data)
	
	#data_l = 
	
	if len(data) < 2 or len(data_t) < 2: return pd.DataFrame(columns=data.columns)
	
	dataY   = data.V2
	dataY_t = data_t.V2
	
	dataX = pd.DataFrame({'V9':data.V9})
	dataX = dataX.loc[:,'V9']
	
	dataX_t = pd.DataFrame({'V9':data_t.V9})
	dataX_t = dataX_t.loc[:,'V9']
	
	if len(dataY) < 2 or len(dataY_t) < 2: return pd.DataFrame(columns=data.columns)
	
	return pd.DataFrame(columns=data.columns)

def load_settings():
	global IMAGE_DIR 
	global LOCATION_HE5
	global LOCATION_CSV
	global LOCATION_OUT
	global LIMIT_TRAINING
	
	if os.path.isfile("settings.dct"):
		with open("settings.dct", 'r') as file:
			settings = json.loads(file.read())
			print ("Loaded settings from file")
			
			IMAGE_DIR = settings["IMAGE_DIR"]
			LOCATION_HE5 = settings["LOCATION_HE5"]
			LOCATION_CSV = settings["LOCATION_CSV"]
			LOCATION_OUT = settings["LOCATION_OUT"]
			LIMIT_TRAINING = settings["LIMIT_TRAINING"]
			
			Path(IMAGE_DIR).mkdir(parents=True, exist_ok=True)
			Path(LOCATION_CSV).mkdir(parents=True, exist_ok=True)
			Path(LOCATION_OUT).mkdir(parents=True, exist_ok=True)
			
	else:
		with open("settings.dct", 'w') as file:
			settings = {"IMAGE_DIR":"img_files/", 
						"LOCATION_HE5":"H5_files/",
						"LOCATION_CSV":"csv_files/",
						"LOCATION_OUT":"out_files/",
						"LIMIT_TRAINING":True}
			file.write(json.dumps(settings))
			print ("Wrote default settings to file, please alter as needed,")
			print ("The directories will be created (except HE5) on the next run.")
			
		sys.exit()
	return settings
	
def load_(date_range, settings):

	data_  = ed.load_data(date_range, settings)
	data_  = prepare_data(data_)
	
	return data_
	
def save_climo(bin_cl, climo_filename):

	print ("Saving climo to", climo_filename)
	
	file = open(climo_filename, "wb")
	
	pickle.dump(bin_cl, file)
	
	file.close()

	return
	
def load_climo(climo_filename):

	print ("Loading climo from", climo_filename)
	
	file = open(climo_filename, "rb")
	
	bin_cl = pickle.load(file)
	
	file.close()
	
	return bin_cl
	
if __name__ == "__main__":
	# Turn interactive mode off so that matplotlib doesn't fail when run in background
	plt.ioff()
	
	# ****
	# Check for settings file and load it.
	# ****

	settings = load_settings()
	
	climo_filename = r"climo.bin"
	
	if USE_EXISTING_CLIMO_FILE and os.path.isfile(climo_filename):
		no_existing_climo_file = False
	else:
		no_existing_climo_file = True
		
	# ****
	# Get date ranges for training from commane line
	# ****

	print("-----\nGet date ranges for training from commane line")	
	
	date_range_learn = list()
	
	# ****
	# Load training data from csv to np.array
	# ****
	
	#date_range_learn.append((dt.datetime(2007, 4, 7, 11, 4), dt.datetime(2007, 4, 7, 11, 4)))
	
	for i in range(1,len(sys.argv)-2,2):
		date_1 = dt.datetime.strptime(sys.argv[i], "%Y%m%d")
		date_2 = dt.datetime.strptime(sys.argv[i + 1], "%Y%m%d") + dt.timedelta(hours = 23, minutes = 59, seconds = 59)
		date_range_learn.append((date_1, date_2))

	print("date_range_learn ", date_range_learn)	

	if no_existing_climo_file:
		print("-----\nLoad training data from csv to np.array")
		
		data_t = pd.DataFrame()
		
		for i in range(0,len(date_range_learn)):
			print("Loading ", date_range_learn[i])
			data_t = pd.concat([data_t, load_(date_range_learn[i], settings)])
			
		#print("data_t:\n", data_t)
			
		#date_range_recon = (dt.datetime(2008, 4, 7, 11, 9), dt.datetime(2008, 4, 7, 11, 9))

		# ****
		# Alter and screen data_t as needed
		# ****
		
		print("-----\nScreening data")	

		data_t = data_t[(data_t.V2.astype(float) >= -10.) & (data_t.V2.astype(float) <=  10.)]
		data_t = data_t[data_t.PRS > 650]
		data_t = data_t[data_t.V3 > 60]
		
		if SCREEN_XFLAG:
			print("Cleaning training xflag")
			data_t = data_t[(data_t.Xflag == 0)]
		
		limit_training=LIMIT_TRAINING
		
		bin_cl, trained_data = train_climate(data_t, limit_training)
		data_t = trained_data
				
		if USE_EXISTING_CLIMO_FILE:
			save_climo(bin_cl, climo_filename)

		del trained_data

	else:
		bin_cl = load_climo(climo_filename)

	
	date_1 = dt.datetime.strptime(sys.argv[len(sys.argv)-2], "%Y%m%d")
	date_2 = dt.datetime.strptime(sys.argv[len(sys.argv)-1], "%Y%m%d") + dt.timedelta(hours = 23, minutes = 59, seconds = 59)
	date_range_recon = (date_1, date_2)

	print("date_range_recon ", date_range_recon)
	
	# ****
	# Load reconstruction data
	# ****	
	
	print("-----\nLoad reconstruction data from csv to np.array")
	
	data = load_(date_range_recon, settings)
	
	#print("data:\n", data)
	
	date_range_l = date_range_learn[0]
	
	"""
	date_str_learn = str("%4dm%02d%02d_%4dm%02d%02d" % (date_range_l[0].year, \
														date_range_l[0].month, \
														date_range_l[0].day, \
														date_range_l[1].year, \
														date_range_l[1].month, \
														date_range_l[1].day))
	
	date_str_recon = str("%4dm%02d%02d_%4dm%02d%02d" % (date_range_recon[0].year, \
														date_range_recon[0].month, \
														date_range_recon[0].day, \
														date_range_recon[1].year, \
														date_range_recon[1].month, \
														date_range_recon[1].day))
	
	"""
	
	data = data[(data.V2.astype(float) >= -10.) & (data.V2.astype(float) <=  10.)]
	data = data[data.PRS > 650]
	data = data[data.V3 > 60]
	
	if SCREEN_XFLAG:
		print("Cleaning data xflag")
		data = data[(data.Xflag == 0)]

	data_c = pd.DataFrame()
		
	#output_climate(bin_cl, date_range_learn)
	data_c = remove_climate(data, bin_cl)
	
	del bin_cl
	
	bin_data, trained_data = train_climate(data)
	data = trained_data

	del trained_data
	
	#sys.exit()
	
	#print (bin_space)
	
	#bin_space = range(len(salb_bins) * len(raza_bins) * len(vzea_bins) * len(szea_bins))
	
	"""
	bin_space = np.unique(np.concatenate([(data_t.bin_ids), (data.bin_ids)]))
	
	print(bin_space)
	#data_l = pd.DataFrame(columns=data.columns)
	
	for bin_i in bin_space:
		
		print(bin_i)
		
		data_t_bin = data_t[data_t.bin_ids == bin_i]
		data_bin = data[data.bin_ids == bin_i]
		
		#print(data1_dgg.shape)
		#print(data2_dgg.shape)
		
		#print(data1_dgg)
		#print(data2_dgg)

		#data_l_dgg = run_data(data_t_bin, data_bin)
		
		#print(data_l1_dgg.shape)
		#print(data_a1_dgg.shape)
		#print(data_r1_dgg.shape)
		
		#print(data_l1_dgg)
		#print(data_a1_dgg)
		#print(data_r1_dgg)
		
		#data_l = data_l.append(pd.concat([data_l_dgg]))

	print(data_l)
	"""
	#if len(data_l) > 0: data_l = data_l.sort_values(by=['IDX'])

	#print(data_l.shape)
	#print(data_a.shape)
	#print(data_r.shape)
	
	data_orbits = pd.DataFrame({'V15':data.V15}).astype(int).to_numpy()
	orbits = np.unique(sorted(data_orbits))
	
	if DO_ORBIT_OUTPUT:
		
		# ****
		# Orbital file output
		# ****	
		
		for orb in orbits:
		
			data_set = data[data.V15.astype(int) == orb]
			data_cli = data_c[data_c.V15.astype(int) == orb]
			#data_lin = data_l[data_l.V15.astype(int) == orb]

			# Construct parallel dataframe with all predictions and (original) independent variables
			data_orb = pd.DataFrame({'V1' :data_set.V1, \
									 'V2' :data_set.V2, \
									 'CLI':data_cli.V2, \
									 #'LIN':data_lin.V2, \
									 'V3' :data_set.V3, \
									 'V4' :data_set.V4, \
									 'V5' :data_set.V5, \
									 'V6' :data_set.V6, \
									 'V7' :data_set.V7, \
									 'V8' :data_set.V8, \
									 'V9' :data_set.V9, \
									 'V10':data_set.V10, \
									 'V11':data_set.V11, \
									 'V12':data_set.V12, \
									 'V14':data_set.V14, \
									 'IDX':data_set.IDX, \
	                        		 'Gflag':data_set.Gflag})

			data_orb = data_orb.sort_values(by=['IDX'])

			# OUTPUT TEXT FILES HERE
			#print("data out for ", str(orb))
			
			filename = ed.save_data(orb, "CLI", data_cli.V2, data_orb, settings)
			print("saved to ", filename)
			
	print ("done")
