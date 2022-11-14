#!/usr/bin/env python

"""
NAME:
  sounding_splitter.py

PURPOSE:
  Read UWYO sounding data from an large text file containing many individual 
  UWyo soundings copied from the UWyo online database and creates individual
  files for each sounding.

MODIFICATIONS:
  Blake Sorenson  - 170628: Written 

"""

# Import functions
import numpy
from time import strptime
import sys

if(len(sys.argv)!=2):
    print "SYNTAX: ./data_splitter.py UWYO_archive_file"
    sys.exit()


start_indices = []
end_indices = []

infile = sys.argv[1]
# Get the year from the archive file's name
year = infile.split('/')[-1].split('_')[1]
# Read in all the lines from the archive file with readlines
lines = open(infile,'r').readlines()

total_lines = 0
with open(infile) as g:
    num_lines = sum(1 for line in g)
    total_lines = num_lines

# Find the indices of the beginning line of each sounding in the large file
with open(infile) as f:
    for i, line in enumerate(f):
        if(line[0]=='7'):
            start_indices.append(i)

# Find the indices of the end of each sounding
for i in range(1,len(start_indices)):
    end_indices.append(start_indices[i]-1)

# Add the last index in the file to the end list 
end_indices.append(total_lines-2)

# Convert the lists into numpy arrays
start_indices = numpy.array(start_indices)
end_indices = numpy.array(end_indices)

file_dump = '/nas/home/bsorenson/prg/thermosonde1617/monthly_analysis/sounding_archive/'+year+'/'
# Loop over all the files and make extra files for each sounding
for i in range(0, len(start_indices)):
    # Grab the first line of each sounding that contains date and time info
    first_line = lines[start_indices[i]:end_indices[i]][0].split(' ')
    # Get the sounding 
    location = first_line[1]
    str_year  = first_line[-1].strip()[2:]
    int_month = strptime(first_line[-2],'%b').tm_mon
    if(int_month<10):
        str_month = '0'+str(int_month)
    else:
        str_month = str(int_month)
    str_day   = first_line[-3]
    str_time  = first_line[-4][:2]
    filename = file_dump+location+'/soundings/'+str_year+str_month+str_day+'_'+str_time+'0000_'+location+'_UWYO.txt'
    
    outfile = open(filename,'w')
    outfile.writelines(lines[start_indices[i]:end_indices[i]])
    outfile.close()    
    print "Saved file: "+filename
