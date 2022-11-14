#!/usr/bin/env python

"""
NAME:
  gif_maker.py

PURPOSE:
  Creates a gif out of any number of .png images provided to the program

SYNTAX:
  ./gif_maker.py gif_name any_number_of_png_images

MODIFICATIONS:
  Blake Sorenson  - 180119: Written 

"""

# Import functions
import glob
import os
import sys

if(len(sys.argv) < 3):
    print "SYNTAX: ./gif_maker.py gif_name any_number_of_png_images [clean]"
    print "         clean - removes any .png files from the current directory."
    print "                 Must be at the end of the command"
    sys.exit(1)

# Check if the user forgot to put a gif name
if sys.argv[1].split('.')[-1] == 'png':
    print "ERROR: No gif name added. Exiting"
    sys.exit(1)

gif_name = sys.argv[1]

CLEAN = False
# Check if user wants to remove extra .png files after
for arg in sys.argv:
    if arg == "clean":
        CLEAN = True

if CLEAN is True:
    sys.argv = numpy.delete(sys.argv,-1)

# Make the gif for each of the image files
file_list = sys.argv[2:]
list.sort(file_list)

with open('image_list.txt','w') as file:
    for item in file_list:
        file.write("%s\n" % item)

print "Creating .gif file"
os.system('convert @image_list.txt {}.gif'.format(gif_name))
print "Saved gif: "+gif_name+".gif"
if CLEAN is True:
    print "Removing .png image files"
    os.system('rm *.png')
