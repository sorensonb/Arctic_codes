#!/usr/bin/env python

"""
  NAME:
    monthly_temperature_analyzer.py
  
  PURPOSE:
    Plot any number of refractive index structure parameter profiles on a single
    graph. 
  
  MODIFICATIONS:
    Blake Sorenson  - 180218???: Written 
  
"""

# Import functions
import numpy
import sys
import os
import commands
import warnings
import matplotlib.pyplot as plt
from time import strptime
from calendar import month_name
from readsounding import *

if(len(sys.argv)<3):
    print "SYNTAX: ./monthly_temperature_analyzer.py time monthly_sounding_archive_file [meter]"
    print "        time - 00 or 12"
    print ""
    print "        OPTIONAL ARGUMENTS:"
    print "        meter - If meter is added as the last argument, the code"
    print "                will calculate average and standard deviations over"
    print "                500 meter windows instead of pressure levels"
    sys.exit()

METER=False
if(len(sys.argv)==4):
    if(sys.argv[3]=='meter'):
        METER=True

work_path='/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/'
# Add check for compare_temperature.py, plot_temperature.py, and readsounding.py
#if (os.path.isfile(work_path+'sounding_splitter.py') is False):
#    print "ERROR: "+work_path+"sounding_splitter.py not found"
#    print "       Please ensure the necessary program is in the correct location"
#    sys.exit(1)

# Extract the individual soundings from the archive file
infile = sys.argv[2]
date = infile.split('/')[-1].split('_')
month = date[0]
month_no = strptime(month[:3],'%b').tm_mon
str_monthno = str(month_no)
if(month_no<10):
    str_monthno = '0'+str_monthno
tmonth = month_name[month_no]
year = date[1]
location = date[2]
file_path = '/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/sounding_archive/'+year+'/'+location+'/soundings/'
archive_path = '/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/sounding_archive/'+year+'/'+location+'/'
# If the files are already extracted and located in 
#   /nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/2017/***/soundings/,
# ignore this step
status,files_there = commands.getstatusoutput('ls '+file_path+year[2:]+str_monthno+'*.txt')
files_there = files_there.split('\n')
if(len(files_there)==1):
    print "Extracting soundings from archive file"
    os.system(work_path+'sounding_splitter.py '+archive_path+infile.strip().split('/')[-1])
else:
    print "Soundings already extracted from archive file. Proceeding"
stime = sys.argv[1]

# Get the file names from the recently split archive file
status,output = commands.getstatusoutput('ls '+file_path+year[2:]+str_monthno+'*_'+stime+'0000_'+location+'_UWYO.txt') 
files = output.split('\n')
list.sort(files)

# Calculate Cn2 for each sounding and calculate the average and standard deviation
# for each sounding
avgs = []
stdevs = []
temperature_sfcto925 = []
temperature_925to850 = []
temperature_850to700 = []
temperature_700to500 = []
temperature_500to400 = []
temperature_400to300 = []
temperature_300to200 = []
temperature_200to150 = []
temperature_150to100 = []
temperature_100to50 = []
temperature_50totop = []
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for f in files:
        data = readsounding(f)
        temperature_sfcto925.append(numpy.average(data['TEMP'][numpy.where(data['PRESS']>925)]))
        temperature_925to850.append(numpy.average(data['TEMP'][numpy.where((data['PRESS']<=925) & (data['PRESS']>850))]))
        temperature_850to700.append(numpy.average(data['TEMP'][numpy.where((data['PRESS']<=850) & (data['PRESS']>700))]))
        temperature_700to500.append(numpy.average(data['TEMP'][numpy.where((data['PRESS']<=700) & (data['PRESS']>500))]))
        temperature_500to400.append(numpy.average(data['TEMP'][numpy.where((data['PRESS']<=500) & (data['PRESS']>400))]))
        temperature_400to300.append(numpy.average(data['TEMP'][numpy.where((data['PRESS']<=400) & (data['PRESS']>300))]))
        temperature_300to200.append(numpy.average(data['TEMP'][numpy.where((data['PRESS']<=300) & (data['PRESS']>200))]))
        temperature_200to150.append(numpy.average(data['TEMP'][numpy.where((data['PRESS']<=200) & (data['PRESS']>150))]))
        temperature_150to100.append(numpy.average(data['TEMP'][numpy.where((data['PRESS']<=150) & (data['PRESS']>100))]))
        temperature_100to50.append(numpy.average(data['TEMP'][numpy.where((data['PRESS']<=100) & (data['PRESS']>50))]))
        temperature_50totop.append(numpy.average(data['TEMP'][numpy.where(data['PRESS']<=50)]))

log_data_list = [numpy.array(temperature_sfcto925),numpy.array(temperature_925to850),\
    numpy.array(temperature_850to700),numpy.array(temperature_700to500),\
    numpy.array(temperature_500to400),numpy.array(temperature_400to300),\
    numpy.array(temperature_300to200),numpy.array(temperature_200to150),\
    numpy.array(temperature_150to100),numpy.array(temperature_100to50),\
    numpy.array(temperature_50totop)]

# Get rid of missing values
for i,s in enumerate(log_data_list):
    log_data_list[i] = log_data_list[i][numpy.logical_not(numpy.isnan(log_data_list[i]))]

plot_title='Daily '+stime+'Z Averaged Temperature Profiles: '+location+' '+tmonth+' '+year

fig1 = plt.figure()
ax = fig1.gca()
bx = ax.boxplot(log_data_list,vert=0)

ax.set_yticklabels(['Sfc-925 mb','925-850 mb','850-700 mb','700-500 mb','500-400 mb',\
    '400-300 mb','300-200 mb','200-150 mb','150-100 mb','100-50 mb','50 mb-top'])
plt.yticks(rotation=45,fontsize=6)
plt.title(plot_title)
plt.xlabel('Temperature [degC]')
plt.xlim(-80,40)

if(METER is True):
    image_name = month+'_'+year+'_'+location+'_'+stime+'Z_temp_meter_boxplots.png'
else:
    image_name = month+'_'+year+'_'+location+'_'+stime+'Z_temp_boxplots.png'
fig1.savefig(image_name,dpi=200)
print "Saved image: "+image_name
