#!/usr/bin/env python
"""


"""
import sys

def heightFinder(h):
    if(h==0):
        return 0
    if(h==1):
        return 1
    elif(h==2):
        return 2
    return 1+heightFinder(h-1)+heightFinder(h-2)

if(len(sys.argv)!=2):
    print "SYNTAX: ./recursive_test.py number"
    sys.exit(1)

height = heightFinder(int(sys.argv[1]))
print "Minimum number of nodes in tree of height "+sys.argv[1]+" is "+\
            str(height)
