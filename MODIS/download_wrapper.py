#!/usr/bin/env python

"""


"""

from auto_asos_downloader import *
import os
import sys
import numpy as np
import requests

#files = main()

##!#def get_elevation(lat, lon):
##!#    query = ('https://nationalmap.gov/epqs/pqs.php'f'?x={lon}&y={lat}&units=Meters&output=json')
##!#    r = requests.get(query).json()  # json object, various ways you can extract value
##!#    elevation = r['USGS_Elevation_Point_Query_Service']['Elevation_Query']['Elevation']
##!#    return elevation

args = sys.argv
if(len(sys.argv) < 5):
    print("SYNTAX: ./asos_download_wrapper.py case_date start_date end_date [at least one station]")
    print("     case_date, start_date and end_date are both formatted YYYYMMDD")
    print("     stations must be 3-letter ICAO IDs (GFK, FAR, CKN...)")
    sys.exit()

case_date  = sys.argv[1]
start_date = sys.argv[2]
end_date   = sys.argv[3]
stns = sys.argv[4:]

# Create dictionary to hold the elevations for each station
# ---------------------------------------------------------
##stn_dict = {}

with open("./stations.txt","w") as fout:
    for stn in stns:
        fout.write(stn + "\n")

# Extract year, month, and day information from the start and end dates
# ---------------------------------------------------------------------
s_year  = int(start_date[:4])
s_month = int(start_date[4:6])
s_day   = int(start_date[6:8])

e_year  = int(end_date[:4])
e_month = int(end_date[4:6])
e_day   = int(end_date[6:8])

print(start_date, end_date)

outfile = 'asos_data_' + case_date + '.csv'
if(os.path.exists(outfile)):
    add_num = 2
    outfile = 'asos_data_' + case_date + '_' + str(add_num) + '.csv'
    while(os.path.exists(outfile)):
        print(outfile)
        add_num += 1
        outfile = 'asos_data_' + case_date + '_' + str(add_num) + '.csv'

print(outfile)

# Call the main function from auto_asos_downloader to retrieve the METARS
# for the provided dates and the provided locations. 'files' contains a
# list of all the auto-generated files created by the downloader, so these
# files will be concatenated into one large file.
# ------------------------------------------------------------------------
files = main(s_year,s_month,s_day,e_year,e_month,e_day)
#files = main(year,month,day)

convert_temp = False
with open (outfile,'w') as fout:
    for ii, ffile in enumerate(files):
        print(ffile)
        with open(ffile,'r') as fin:
            flines = fin.readlines()
            if(ii == 0):
                idx = 5
            else:
                idx = 6
            for jj, line in enumerate(flines[idx:]):
                tmpline = line.strip()
                splitline = tmpline.split(',')
                if((idx == 5) & (jj == 0)):
                    ##!#splitline.append('elevation')
                    if('tmpc' not in splitline):
                        convert_temp = True
                        splitline.append('tmpc')
                    tmpline = ','.join(splitline)
                    tmpf_idx = splitline.index('tmpf')
                else:
                    if(convert_temp):
                        if(jj > 0):
                            # See if the elevation has already been retreived
                            # -----------------------------------------------
                            ##!#if(splitline[0] not in stn_dict.keys()):
                            ##!#    # Find elevation for current station
                            ##!#    stn_dict[splitline[0]] = get_elevation(float(splitline[3]), float(splitline[2]))

                            # Insert elevation data to current line
                            # -------------------------------------
                            ##!#splitline.append(str(stn_dict[splitline[0]]))
                     
                            # Insert  tmpc data
                            # -----------------
                            if(splitline[tmpf_idx] == 'M'):
                                tmpc = 'M'
                            else:
                                tmpc = np.round((float(splitline[tmpf_idx]) - 32.) * 5./9., 1)
                            splitline.append(str(tmpc))

                            tmpline = ','.join(splitline)
                    
                fout.write(tmpline+'\n')
        #print("rm "+ffile) 
        os.system("rm " + ffile)
    print("saved file",outfile)
