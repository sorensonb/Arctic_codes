#!/usr/bin/env python

"""


"""

from auto_asos_downloader import *
import os
import sys
import numpy as np

print("yay")

#files = main()

args = sys.argv
if(len(sys.argv) < 3):
    print("SYNTAX: ./asos_download_wrapper.py YYYYMMDD [at least one station]")
    print("     stations must be 3-letter ICAO IDs (GFK, FAR, CKN...)")
    sys.exit()

date = sys.argv[1]
stns = sys.argv[2:]

with open("./stations.txt","w") as fout:
    for stn in sys.argv[2:]:
        fout.write(stn + "\n")

year  = int(date[:4])
month = int(date[4:6])
day   = int(date[6:8])

print(date)
files = main(year,month,day)

outfile = 'asos_data_' + date + '.csv'
if(os.path.exists(outfile)):
    outfile = 'asos_data_' + date + '_2.csv'

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
                    if('tmpc' not in splitline):
                        convert_temp = True
                        splitline.append('tmpc')
                        tmpline = ','.join(splitline)
                else:
                    if(convert_temp):
                        if(jj > 0):
                            if(splitline[4] == 'M'):
                                tmpc = 'M'
                            else:
                                tmpc = np.round((float(splitline[4]) - 32.) * 5./9., 1)

                            splitline.append(str(tmpc))
                            tmpline = ','.join(splitline)
                    
                fout.write(tmpline+'\n')
        #print("rm "+ffile) 
        os.system("rm " + ffile)
    print("saved file",outfile)
