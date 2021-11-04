#!/usr/bin/env python

"""
  NAME:

  PURPOSE:

  SYNTAX:
    ./auto_AIRS_download.py air_file_list.txt

    NOTE: it is assumed that the user has a urs_cookies file set up

  MODIFICATIONS:
    Blake Sorenson

"""
import glob, os, sys

if(len(sys.argv) != 2):
    print("SYNTAX: ./auto_AIRS_donwload.py airs_file_list.txt")
    sys.exit()

# Download the data using wget
cmnd = "wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies -i " + sys.argv[1]
print(cmnd)
os.system(cmnd)


airs_data_path = "/home/bsorenson/data/AIRS/Aqua/"
filenames = glob.glob('*TTP*')
for name in filenames:
    new_name = name.split('&')[2][6:]
    cmnd = "mv '" + name + "' " + airs_data_path + new_name
    print(cmnd)
    os.system(cmnd)
