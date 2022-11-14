#!/usr/bin/env python

"""
  NAME:
   monthly_cn2_analyzer.py
  
  PURPOSE:
    Plot any number of refractive index structure parameter profiles on a single
    graph. 
  
  MODIFICATIONS:
    Blake Sorenson  - 170714: Written 
  
"""

# Import functions
import numpy
import sys
import os
import commands
import warnings
from time import strptime
from calendar import month_name
import matplotlib.pyplot as plt
from readsounding import *
from compare_cn2 import *

if(len(sys.argv)<3):
    print "SYNTAX: ./monthly_analyzer.py time monthly_sounding_archive_file [meter]"
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
# Add check for compare_cn2.py, plot_cn2.py, and readsounding.py
#if (os.path.isfile(work_path+'sounding_splitter.py') is False):
#    print "ERROR: "+work_path+"sounding_splitter.py not found"
#    print "       Please ensure the necessary program is in the correct location"
#    sys.exit(1)

# Extract the individual soundings from the archive file
infile = sys.argv[2]
# Grab the month information from the monthly archive file name
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
cn2_sfcto925 = []
cn2_925to850 = []
cn2_850to700 = []
cn2_700to500 = []
cn2_500to400 = []
cn2_400to300 = []
cn2_300to200 = []
cn2_200to150 = []
cn2_150to100 = []
cn2_100to50 = []
cn2_50totop = []
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for f in files:
        data = readsounding(f)
        data = cn2_calc(data) 
    #    log_cn2 = numpy.log10(data['CN2'])
        cn2_sfcto925.append(numpy.average(data['CN2'][numpy.where(data['PRESS']>925)]))
        cn2_925to850.append(numpy.average(data['CN2'][numpy.where((data['PRESS']<=925) & (data['PRESS']>850))]))
        cn2_850to700.append(numpy.average(data['CN2'][numpy.where((data['PRESS']<=850) & (data['PRESS']>700))]))
        cn2_700to500.append(numpy.average(data['CN2'][numpy.where((data['PRESS']<=700) & (data['PRESS']>500))]))
        cn2_500to400.append(numpy.average(data['CN2'][numpy.where((data['PRESS']<=500) & (data['PRESS']>400))]))
        cn2_400to300.append(numpy.average(data['CN2'][numpy.where((data['PRESS']<=400) & (data['PRESS']>300))]))
        cn2_300to200.append(numpy.average(data['CN2'][numpy.where((data['PRESS']<=300) & (data['PRESS']>200))]))
        cn2_200to150.append(numpy.average(data['CN2'][numpy.where((data['PRESS']<=200) & (data['PRESS']>150))]))
        cn2_150to100.append(numpy.average(data['CN2'][numpy.where((data['PRESS']<=150) & (data['PRESS']>100))]))
        cn2_100to50.append(numpy.average(data['CN2'][numpy.where((data['PRESS']<=100) & (data['PRESS']>50))]))
        cn2_50totop.append(numpy.average(data['CN2'][numpy.where(data['PRESS']<=50)]))

log_data_list = [numpy.log10(cn2_sfcto925),numpy.log10(cn2_925to850),\
    numpy.log10(cn2_850to700),numpy.log10(cn2_700to500),\
    numpy.log10(cn2_500to400),numpy.log10(cn2_400to300),\
    numpy.log10(cn2_300to200),numpy.log10(cn2_200to150),\
    numpy.log10(cn2_150to100),numpy.log10(cn2_100to50),\
    numpy.log10(cn2_50totop)]

#avg_data_list = [np.average(x) for x in log_data_list]
#stdev_data_list = [np.std(y) for y in log_data_list]

# Get rid of missing values
for i,s in enumerate(log_data_list):
    log_data_list[i] = log_data_list[i][numpy.logical_not(numpy.isnan(log_data_list[i]))]

# Generate a plot title from the input file name
date = infile.split('/')[-1].split('.')[0].split('_')
month = date[0]
month_no = strptime(month[:3],'%b').tm_mon
tmonth = month_name[month_no]
year = date[1]
location = date[2]
plot_title='Daily '+stime+'Z Averaged $C_n^2$ Profiles: '+location+' '+tmonth+' '+year

fig1 = plt.figure(figsize=(8,6))
ax = fig1.gca()
bx = ax.boxplot(log_data_list, vert=0)

ax.set_yticklabels(['Sfc-925 mb','925-850 mb','850-700 mb','700-500 mb','500-400 mb',\
    '400-300 mb','300-200 mb','200-150 mb','150-100 mb','100-50 mb','50 mb-top'])
plt.yticks(rotation=45,fontsize=9)
plt.title(plot_title)
plt.xlabel('$log_{10}$ $C_n^2$')
plt.ylabel('Pressure levels [mb]',fontsize=10)
plt.xlim(-20,-12)

if(METER is True):
    image_name = month+'_'+year+'_'+location+'_'+stime+'Z_cn2_meter_boxplots.png'
else:
    image_name = month+'_'+year+'_'+location+'_'+stime+'Z_cn2_boxplots.png'
fig1.savefig(image_name,dpi=200)
print "Saved image: "+image_name
