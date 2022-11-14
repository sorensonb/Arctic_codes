#!/usr/bin/env python

"""
  NAME:
  
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
import matplotlib.pyplot as plt
from readsounding import *
from compare_cn2 import *

if(len(sys.argv)==1):
    print "SYNTAX: ./monthly_analyzer.py monthly_sounding_archive_file"
    sys.exit()

work_path='/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/'
# Add check for compare_cn2.py, plot_cn2.py, and readsounding.py
#if (os.path.isfile(work_path+'sounding_splitter.py') is False):
#    print "ERROR: "+work_path+"sounding_splitter.py not found"
#    print "       Please ensure the necessary program is in the correct location"
#    sys.exit(1)

# Extract the individual soundings from the archive file
print "Extracting soundings from archive file"
os.system(work_path+'sounding_splitter.py '+sys.argv[1])

# Get the file names from the recently split archive file
status,output = commands.getstatusoutput('ls 1*_UWYO.txt') 
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
for f in files:
    data = readsounding(f)
    data = cn2_calc(data) 
    log_cn2 = numpy.log10(data['CN2'])
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
    #avgs.append(numpy.average(log_cn2))
    #stdevs.append(numpy.std(log_cn2))

log_data_list = [numpy.log10(cn2_sfcto925),numpy.log10(cn2_925to850),\
    numpy.log10(cn2_850to700),numpy.log10(cn2_700to500),\
    numpy.log10(cn2_500to400),numpy.log10(cn2_400to300),\
    numpy.log10(cn2_300to200),numpy.log10(cn2_200to150),\
    numpy.log10(cn2_150to100),numpy.log10(cn2_100to50),\
    numpy.log10(cn2_50totop)]

# Get rid of missing values
for i,s in enumerate(log_data_list):
    log_data_list[i] = log_data_list[i][numpy.logical_not(numpy.isnan(log_data_list[i]))]

date = sys.argv[1].split('/')[-1].split('.')[0].split('_')
month = date[0]
year = date[1]
location = date[2]
plot_title='Daily Averaged $C_n^2$ Profiles: '+location+' '+month+' '+year

fig1 = plt.figure()
ax = fig1.gca()
bx = ax.boxplot(log_data_list)

ax.set_xticklabels(['Sfc-925','925-850','850-700','700-500','500-400',\
    '400-300','300-200','200-150','150-100','100-50','50-top'])
plt.xticks(rotation=45)
plt.title(plot_title)
plt.ylabel('$log_{10}$ $C_n^2$')
plt.ylim(-20,-13)

image_name = month+'_'+year+'_'+location+'_cn2_boxplots.png'
fig1.savefig(image_name,dpi=200)
print "Saved image: "+image_name
