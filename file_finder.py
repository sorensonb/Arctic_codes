#!/usr/bin/env python
"""


"""
import os
import sys
from datetime import datetime
import subprocess

home_dir = os.environ['HOME']

base_dir = home_dir + "/Research/"
dest_dir = home_dir + "/Arctic_codes/"

# First, list all the python files
# ---------------------------------
cmnd = 'find ' + base_dir + ' -type f ' + \
    '-name \'*.py\''
print(cmnd)
#os.system(cmnd)

python_files = subprocess.check_output(cmnd, \
    shell = True).decode('utf-8').strip().split('\n')

# Then, grab any fortran files
# ---------------------------------
cmnd = 'find ' + base_dir + ' -type f ' + \
    '-name \'*.f90\''
print(cmnd)
#os.system(cmnd)

fortran_files = subprocess.check_output(cmnd, \
    shell = True).decode('utf-8').strip().split('\n')

python_files = python_files + fortran_files

#file_info = 

today = datetime.today()

data_dict = {}

for pfile in python_files:
    cmnd = 'ls -l ' + pfile
    info = subprocess.check_output(cmnd, shell = True\
        ).decode('utf-8').strip().split()

    # Make a local datetime object
    # ----------------------------
    if(len(info[7].split(':')) == 2):
        # Modified the current year
        date_str = str(today.year)+info[5]+info[6].zfill(2) + \
            ' ' + info[7]
        dt_date_str = datetime.strptime(date_str, '%Y%b%d %H:%M')
        #print(dt_date_str)
    else:
        # Modified a previous year
        date_str = info[7]+info[5]+info[6].zfill(2)
        dt_date_str = datetime.strptime(date_str, '%Y%b%d')
        #print(dt_date_str)

    data_dict[pfile] = dt_date_str

    #print(pfile, info[5:8])
    #os.system(cmnd)

sorted_dict = sorted(data_dict.items(), key = lambda x:x[1])
converted_dict = dict(sorted_dict)

# Print the sorted dictionary, with the newest files
# at the end of the list
# ---------------------------------------------------
for key in converted_dict.keys():
    print(converted_dict[key], key.strip().split('/')[-1])

