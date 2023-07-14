#!/usr/bin/env python

"""
  NAME:
    auto_TROPOMI_download.py

  PURPOSE:
    Downloads all TROPOMI data files for a given month.

  SYNTAX:
    $ python auto_TROPOMI_download.py

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2023/07/12:
      Written

"""

import numpy as np
from bs4 import BeautifulSoup
from datetime import datetime, timedelta
import requests
import subprocess
import os

# listFD grabs all files and directories that end with the provided ext at the
# provided url.
# ----------------------------------------------------------------------------
def listFD(url, ext=''):
    page = requests.get(url).text
    #print page
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') \
        if node.get('href').endswith(ext)]

# Download
def download_TROPOMI_single_day(date_str, dest_dir = './', \
        run_commands = False):

    base_url = 'https://measures.gesdisc.eosdis.nasa.gov/data/AER/TROPOMAER.1'

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = listFD(dt_date_str.strftime(base_url + '/%Y/%j'), ext = '.nc')
    except subprocess.CalledProcessError:
        print("ERROR: No TROPOMI files for the input DTG",date_str)
        return -2

    if(len(files) == 0):
        print("ERROR: No TROPOMI files returned from the request. Exiting")
        return -1
  
    total_files = files[::2]
    files_only = [tfile.strip().split('/')[-1] for tfile in total_files]
 
    # Loop over each filename and download
    # ------------------------------------
    for ii, tfile in enumerate(total_files):
        base_cmnd = "wget --load-cookies ~/.urs_cookies --save-cookies "+\
            "~/.urs_cookies --keep-session-cookies --content-disposition "
        cmnd = base_cmnd + tfile
        print(cmnd)
        if(run_commands): os.system(cmnd)

        # Move the file to the specified destination
        # ------------------------------------------
        if(dest_dir != './'):
            cmnd = 'mv ' + files_only[ii] + ' ' + dest_dir
            print(cmnd)
            if(run_commands): os.system(cmnd)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# BEGIN MAIN
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Change this variable to specify the month for which data is to be downloaded
# ----------------------------------------------------------------------------
date_str = '20220908'
download_TROPOMI_single_day(date_str, \
    run_commands = True)
#begin_dt_date = datetime.strptime(date_str, '%Y%m')
#local_dt_date = begin_dt_date
#
## Begin loop to download all TROPOMI files from all days in the desired
## month
## ---------------------------------------------------------------------
#while(local_dt_date.month == begin_dt_date.month):
#
#    local_str = local_dt_date.strftime('%Y%m%d')
#    
#    # Call the function to download the TROPOMI data for the current
#    # day. Note that changing the 'dest_dir' variable makes the 
#    # function automatically move the downloaded files to that
#    # specified directory. By default, the files are downloaded
#    # to the working directory.
#    # --------------------------------------------------------------
#    download_TROPOMI_single_day(local_str, dest_dir = '/home/bsorenson/data/', \
#        run_commands = True)
#   
#    # Increment the date object by 1 day and continue
#    # ----------------------------------------------- 
#    local_dt_date = local_dt_date + timedelta(days = 1)
