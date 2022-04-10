#!/usr/bin/env python

"""
NAME:
  constants.py

PURPOSE:
  Parse an aircraft constants XML file and return desired aircraft constants as
  a dictionary.

IMPORTING:
  from constants import constants

USAGE:
  Define desired constants array by
    des_cons = ['these', 'are', 'desired', 'constants']

  Parse the 'constants.xml' file, the default aircraft constants file, by
    const_dict = constants(des_cons)
  
  Parse an aircraft constants file for time-independent constants by
    const_dict = constants(des_cons, fname='constants_file.xml')

  Parse an aircraft constants file for time-dependent constants by
    const_dict = constants(des_cons, fname='constants_file.xml', jul_date=12.34)

EXAMPLES:
  Define desired costants array by
    >>> des_cons = ['Length', 'R_const', 'datasystem_timeoffset']

  Parse the 'constants.xml' file, the default aircraft constants file, by
    >>> constants(des_cons)
    {'Length': 8.1285, 'R_const': 287.0, 'datasystem_timeoffset': 0.0}

  Parse an aircraft constants file for time-independent constants by
    >>> constants(des_cons, fname='N555DS_constants.xml')
    {'Length': 9.3477, 'R_const': None, 'datasystem_timeoffset': None}

  Parse an aircraft constants file for time-dependent constants by
    >>> constants(des_cons, fname='N555DS_constants.xml', jul_date=2455692.0)
    {'Length': 9.3477, 'datasystem_timeoffset': '32400', 'R_const': None}

PARAMETERS:
  desired_constants   - an array of strings that represent the names of the
                        desired constants to be found in the XML file
  fname (optional)    - the XML file to be parsed for the constants
                      - default value is 'constants.xml'
  jul_date (optional) - the Julian Date of the file being processed as a float
                      - default value is 0.0

CALLS:
  Python modules numpy, objectify within lxml, and os

NOTES:
  - All Python subroutines that need physical, calibration, or aircraft-specific
    constants should call this method near the beginning of the subroutine to
    have the constants defined. Any constants needed by the calling routine
    should be placed in the parameter list in the file header.
  - The final time period is used to define constants in the event of multiple
    exact time periods being defined.

MODIFICATIONS:
  Nicholas Gapp - 161123: Written
  Nicholas Gapp - 161128: Added absolute path to locate the constants XML file.
  Nicholas Gapp - 161206: Updated path to constants XML files.


Copyright 2016 David Delene
This program is distributed under the terms of the GNU General Public License

This file is part of Airborne Data Processing and Analysis (ADPAA).

ADPAA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ADPAA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ADPAA.  If not, see <http://www.gnu.org/licenses/>.
"""

from lxml import objectify
import numpy as np
import os

def constants(desired_constants, fname='constants.xml', jul_date=0.0):

    # Get absolute path for the constants file.
    # Constants XML files are currently in /usr/local/ADPAA/share directory.
    adpaa_dir = os.environ['ADPAA_DIR']
    abs_fname = adpaa_dir + '/share/' + fname

    # Try to open the file.
    try:
        fhand = open(abs_fname, 'r')
    except:
        print "ERROR: cannot open '"+abs_fname+"' for reading.  Exiting..."
        return -1

    # Read the contents of the XML file to the 'data' variable.
    data = fhand.read()

    # Create the output dictionary.
    out = {}

    # Create an XML tree from the data within the XML file using the fromstring
    # routine within the objectify module.
    c_obj = objectify.fromstring(data)

    # Find the desired arrays within the XML tree, declare them as empty numpy
    # arrays according to the array declaration, and add them to the 'out'
    # dictionary.
    for ad in c_obj.iter(tag='array_dec'):

        if ad.text.strip() in desired_constants:

            # Get the right data type for the array
            if ad.attrib['data'] == 'struct':
                data_type = 'object'
            else:
                data_type = ad.attrib['data']

            # Set the name of the array as the key to a new item in the
            # out dictionary and set the array as the value of the item
            # in the out dictionary.
            key = ad.text.strip()
            if ad.attrib['rows'] == '0':
                value = np.empty(int(ad.attrib['cols']), dtype=data_type)
            elif ad.attrib['cols'] == '0':
                value = np.empty(int(ad.attrib['rows']), dtype=data_type)
            else:
                value = np.empty([int(ad.attrib['rows']), \
                                  int(ad.attrib['cols'])],\
                                 dtype=data_type)

            # Add new item to out dictionary
            out.update({key:value})

    # Find the desired arrays within the XML tree and fill them inside the
    # 'out' dictionary with the correct values.
    for arr in c_obj.iter(tag='array'):
        if arr.text.strip() in desired_constants:

            # Skip desired array if it has not been declared.
            if arr.text.strip() not in out:
                print "ERROR: no array declaration statement found for '"+\
                       arr.text.strip()+"' and therefore cannot fill the array."
                continue

            # Loop through all elements in the array.
            elements = arr.iter(tag='element')
            for e in elements:

                # Fill the multi-dimensional array.
                if ',' in e.attrib['pos']:

                    # Get position from the element declaration line.
                    pos = str.split(e.attrib['pos'], ',')

                    # If the element is designated for a constant, fill this
                    # element with the constant stored in the text content.
                    if e.attrib['type'] == 'const':
                        print e.text.strip()
                        out[arr.text.strip()][int(pos[0])][int(pos[1])] = \
                                                                  e.text.strip()

                    # If the element is designated for a structure, fill this
                    # element with the structure stored within the element.
                    else:
                        out[arr.text.strip()][int(pos[0])][int(pos[1])] = \
                                                                 e.struct.attrib

                # Fill the single-dimensional array.
                else:

                    # If the element is designated for a constant, fill this
                    # element with the constant stored in the text content.
                    if e.attrib['type'] == 'const':
                        print e.text.strip()
                        out[arr.text.strip()][int(e.attrib['pos'])] = \
                                                                  e.text.strip()

                    # If the element is designated for a structure, fill this
                    # element with the structure stored within the element.
                    else:
                        out[arr.text.strip()][int(e.attrib['pos'])] = \
                                                                 e.struct.attrib

    # Find the desired dictionaries (structures) within the XML tree and add
    # them to the 'out' dictionary.
    for struct in c_obj.iter(tag='struct'):
        if struct.text.strip() in desired_constants:
            out.update({struct.text.strip():struct.attrib})

    # Search through all time periods within the XML tree to find the desired
    # time-dependent constants that are valid for the specified julian date.
    for tp in c_obj.iter(tag='time_period'):

        # Set the 'valid_tp' flag to 1 if the time period is valid for the
        # specified julian date.
        valid_tp = 0
        if (float(tp.attrib['start']) == 0.0) & \
           (jul_date <= float(tp.attrib['end'])):
            # There is no starting time of the time period and the specified
            # julian date is less than or equal to the ending time.
            valid_tp = 1
        elif (jul_date >= float(tp.attrib['start'])) & \
             (float(tp.attrib['end']) == 0.0):
            # There is no ending time of the time period and the specified
            # julian date is greater than or equal to the starting time.
            valid_tp = 1
        elif (jul_date >= float(tp.attrib['start'])) & \
             (jul_date <= float(tp.attrib['end'])):
            # The specified julian date is within the starting and ending times
            # of the time period.
            valid_tp = 1

        if valid_tp == 1:
            
            # Find the desired time-dependent arrays within the XML tree,
            # declare them as empty numpy arrays according to the array
            # declaration, and add them to the 'out' dictionary.
            if tp.attrib['type'] == 'arr_dec':
                if tp.array_dec.text.strip() in desired_constants:

                    # Get the right data type for the array
                    if tp.array_dec.attrib['data'] == 'struct':
                        data_type = 'object'
                    else:
                        data_type = tp.array_dec.attrib['data']

                    # Set the name of the array as the key to a new item in the
                    # out dictionary and set the array as the value of the item
                    # in the out dictionary.
                    key = tp.array_dec.text.strip()
                    if tp.array_dec.attrib['rows'] == '0':
                        value = np.empty(int(tp.array_dec.attrib['cols']), \
                                         dtype=data_type)
                    elif ad.attrib['cols'] == '0':
                        value = np.empty(int(tp.array_dec.attrib['rows']), \
                                         dtype=data_type)
                    else:
                        value = np.empty([int(tp.array_dec.attrib['rows']), \
                                          int(tp.array_dec.attrib['cols'])],\
                                         dtype=data_type)

                    # Add new item to out dictionary
                    out.update({key:value})

            # Find the desired time-dependent arrays within the XML tree and
            # fill them inside the 'out' dictionary with the correct values.
            elif tp.attrib['type'] == 'arr':
                if tp.array.text.strip() in desired_constants:

                    # Skip desired array if it has not been declared.
                    if tp.array.text.strip() not in out:
                        print "ERROR: no array declaration statement found" + \
                              " for '"+arr.text.strip()+",' therefore the" +  \
                              " array cannot be filled."
                        continue

                    # Loop through all elements in the array.
                    elements = tp.array.iter(tag='element')
                    for e in elements:

                        # Fill the multi-dimensional array.
                        if ',' in e.attrib['pos']:

                            # Get position from the element declaration line.
                            pos = str.split(e.attrib['pos'], ',')

                            # If the element is designated for a constant, fill
                            # this element with the constant stored in the text
                            # content.
                            if e.attrib['type'] == 'const':
                                out[tp.array.text.strip()]\
                                   [int(pos[0])][int(pos[1])] = e.text.strip()

                            # If the element is designated for a structure, fill
                            # this element with the structure stored within the
                            # element.
                            else:
                                out[tp.array.text.strip()]\
                                   [int(pos[0])][int(pos[1])] = e.struct.attrib

                        # Fill the single-dimensional array.
                        else:

                            # If the element is designated for a constant, fill
                            # this element with the constant stored in the text
                            # content.
                            if e.attrib['type'] == 'const':
                                out[tp.array.text.strip()]\
                                   [int(e.attrib['pos'])] = e.text.strip()

                            # If the element is designated for a structure, fill
                            # this element with the structure stored within the
                            # element.
                            else:
                                out[tp.array.text.strip()]\
                                   [int(e.attrib['pos'])] = e.struct.attrib

            # Find the desired time-dependent dictionaries (structures) within
            # the XML tree and add them to the 'out' dictionary.
            elif tp.attrib['type'] == 'struct':
                if tp.struct.text.strip() in desired_constants:
                    out.update({tp.struct.text.strip():tp.struct.attrib})

            # Find the desired time-dependent single-value constants within the
            # attributes of time period and add them to the 'out' dictionary.
            elif tp.attrib['type'] == 'const':
                for c in desired_constants:
                    if c in tp.attrib:
                        out.update({c:tp.attrib[c]})

    # Find the desired single-value constants within the XML tree and add them
    # to the 'out' dictionary.
    # Need to check for desired constants not in the 'out' dictionary because
    # the find method will replace constants not found with NoneType which will
    # erase everything else in the 'out' dictionary.
    for c in desired_constants:
        if c not in out:
            out[c] = c_obj.find(c)

    # Close file and return the output dictionary.
    fhand.close()
    return out
