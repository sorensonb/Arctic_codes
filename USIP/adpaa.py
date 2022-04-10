#!/usr/bin/env python
"""
  NAME:
    adpaa.py
 
  PURPOSE:
    Creates an object of the ADPAA class. Built-in methods allow data within
    a NASA-formatted ASCII file to be extracted and stored in an object of the
    ADPAA class and allow information within an object of the ADPAA class to be
    written to a file in proper NASA format.

  CALLS:
    Built-in Python modules numpy/numpy.ma, math, re, and sys
    ADPAA Python modules readfile, writefile, julian_date, and constants

  MODIFICATIONS:
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/09/13: Written
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/09/19:
        ReadFile/WriteFile methods now use name class attribute to read from or
        write to a file, added start_time and end_time class attributes, added
        methods to convert time formats, and added/updated comments.
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/09/23:
        ReadFile/WriteFile methods have been moved to separate files and input
        arguments are now optional for both methods, revised comments. 
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/09/26: Imported re module
    Jared Marquis - 161013:
        Fixed bug in sfm2hms function where values less than 10 had no
        leading zeros.
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/11/23:
        Added constants method to parse constants from aircraft constants
        XML files, added calculate_jul_date method in the ReadFile method to
        automatically calculate the Julain Date based on the 'name' class
        attribute, fixed object constructor to include the 'start_time' and
        'end_time' class attributes, added the __check_versions method to the
        class constructor to check for compatibility, and added/updated/revised
        comments. 
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/11/30: Added comments.
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/12/05: 
        Fixed calculate_jul_date method to parse date and time from 'name'
        class attribute no matter where that information is within the string,
        added functionality for any Julian Date to be calculated with
        calculate_jul_date method, and added comments.
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2018/07/26:
        Updated default class attributes: NNCOML=4, FFI=1001, IVOL=1,
        ONAME="Delene, David", ORG="University of North Dakota", VVOL=1, DX=1.0,
        XNAME="Time [seconds]; UT seconds from midnight on day aircraft flight
        started", RDATE=(current date), and DTYPE="Preliminary Data"; 
        automatically creates a VSCAL array of the proper length based on length
        of VNAME; added help and __check_data methods; and updated, revised, and
        added comments.

  IMPORTING:
    from adpaa import ADPAA

  USAGE:
    Create a default APDAA object by
      A = ADPAA()

    Read data from a NASA-formatted ASCII file by
      Specifying filename:
        A.ReadFile('input_filename')
      Using filename set by class attribute 'name':
        A.ReadFile()

    Write data contained within an ADPAA object to a properly formatted file by
      Specifying filename:
        A.WriteFile('output_filename')
      Using filename set by class attribute 'name':
        A.WriteFile()

    Create a dictionary filled with desired aircraft constants if
      Julian Date can be calculated from the 'name' class attribute:
        A.constants(['string', 'array', 'of', 'constants'])
      Julian Date cannot be calculated from the 'name' class attribute:
        A.constants(['string', 'array', 'of', 'constants'], jd=12345.6)

    Mask the data contents with missing values by
      A.mask_MVC()

    Create a string array of output data formats by
      formats = A.create_format_arr(mvc_arr)

    Convert time from seconds from midnight (sfm) to hour-minute-second
    (HH:MM:SS) format by
      hms = A.sfm2hms(sfm)

    Convert time from hour-minute-second (HH:MM:SS) format to seconds from
    midnight (sfm) format by
      sfm = A.hms2sfm(hh, mm, ss)

    Calculate the Julian Date by
      Using 'name' class attribute:
        jd = A.calculate_jul_date()
      Using array filled with date and time information:
        jd = A.calculate_jul_date(dt=['YY','MM','DD','hh','mm','ss'])

    Get help with a specific ADPAA class attribute by
      A.help('att')

    Get help using a list of all ADPAA class attributes and their values by
      A.help()

    Get help with the methods within the ADPAA class using Python's built-in
    help method by
      help(A)

  EXAMPLES:
    Instantiate an ADPAA object
      >>> air_ascent_1Hz = ADPAA()

    Read in data from the '15_08_08_14_57_48.air.ascent.1Hz' file
      >>> air_ascent_1Hz.ReadFile('15_08_08_14_57_48.air.ascent.1Hz')

      --OR--
      >>> air_ascent_1Hz.name = '15_08_08_14_57_48.air.ascent.1Hz'
      >>> air_ascent_1Hz.ReadFile()

    Write data within air_ascent_1Hz object to new file called
    '15_08_08_14_57_48.air.ascent.new.1Hz'
      >>> air_ascent_1Hz.WriteFile('15_08_08_14_57_48.air.ascent.new.1Hz')

      --OR--
      >>> air_ascent_1Hz.name = '15_08_08_14_57_48.air.ascent.new.1Hz'
      >>> air_ascent_1Hz.WriteFile()

    Create a dictionary filled with desired aircraft constants if
      Julian Date can be calculated from the 'name' class attribute:
        >>> air_ascent_1Hz.constants(['R_const', 'c_slope_temp'])
        {'c_slope_temp': '28.207588', 'R_const': 287.0}
      
      Julian Date cannot be calculated from the 'name' class attribute:
        >>> air_ascent_1Hz.constants(['R_const', 'c_slope_temp'], jd=2457243.12)
        {'c_slope_temp': '28.207588', 'R_const': 287.0}

    Mask the data contents with missing values
      >>> air_ascent_1Hz.mask_MVC()

    Convert time from seconds from midnight (sfm) to hour-minute-second
    (HH:MM:SS) format
      >>> air_ascent_1Hz.sfm2hms(55160.36)
      '15:19:20.3600'

    Convert time from hour-minute-second (HH:MM:SS) format to seconds from
    midnight (sfm) format
      >>> air_ascent_1Hz.hms2sfm(14, 57, 48)
      53868.0
    
    Calculate the Julian Date of the file by
      Using 'name' class attribute:
        >>> air_ascent_1Hz.calculate_jul_date()
        2457243.123472

      Using an array filled with date and time information:
        >>> air_ascent_1Hz.calculate_jul_date(dt=['15','08','08','14','57','48'])
        2457243.123472

        --OR--
        >>> air_ascent_1Hz.calculate_jul_date(dt=[15, 8, 8, 14, 57, 48])
        2457243.123472

    Get help with a specific ADPAA class attribute by
      >>> A.help('VFREQ')
      VFREQ:
        Time frequency of the data.
        Current value: 1 Hz Data

    Get help using a list of all ADPAA class attributes and their values by
      >>> A.help()
      Opens a scrollable list (Linux only) with all ADPAA class attributes and
      their values.  Instructions to navigate list are included at the top
      of the screen when list opens.

    Get help with the methods within the ADPAA class using Python's built-in
    help method by
      >>> help(A)
      Opens a scrollable list with all the methods within the ADPAA class.

  NOTES:
    A detailed description of each class attribute is given below.  Uppercase
    attributes come directly from the file header and lower case attributes are
    derived from the information within the file.

      NLHEAD: number of lines (integer) composing the file header. NLHEAD is the
      first recorded value on the first line of an exchange file.
  
      FFI: ASCII file format number. For the UND Citation aircraft data this
      will always be 1001.
  
      ONAME: a character string specifying the name(s) of the originator(s) of
      the exchange file, last name first. On one line and not exceeding 132
      characters.
  
      ORG: character string specifying the organization or affiliation of the
      originator of the exchange file. Can include address, phone number, email
      address, etc. On one line and not exceeding 132 characters.
  
      SNAME: a character string specifying the source of the measurements or
      model results which compose the primary variables, on one line and not
      exceeding 132 characters. Can include instrument name, measurement
      platform, etc.
      
      MNAME: A character string specifying the name of the field project that
      the data were obtained from. 
  
      IVOL: volume number (integer) of the total number of volumes required to
      store a complete dataset, assuming only one file per volume. To be used in
      conjunction with VVOL to allow data  exchange of large datasets requiring
      more than one volume of the exchange medium (diskette, etc.).
  
      VVOL: total number of volumes (integer) required to store the complete
      dataset, assuming one file per volume. If VVOL>1 then each volume must
      contain a file header with an incremented value for IVOL, and continue the
      data records with monotonic independent variable marks.
  
      DATE: UTC date at which the data within the exchange file begins. For
      aircraft data files DATE is the UTC date of takeoff. DATE is in the form
      YYYY MM DD (year, month, day) with each integer value separated by at least
      one space. For example: 2009 09 21.
  
      RDATE: date of data reduction or revision, in the same form as DATE.
  
      DX(s): interval (real) between values of the s-th independent variable,
      X(i,s), i=1,NX(s); in the same units as specified in XNAME(s). DX(s) is
      zero for a non-uniform interval. DX(s) is non-zero for a constant interval.
      If DX(s) is non-zero then it is required that
      NX(s) = (X(NX(s),s)-X(1,s)) / DX(s) + 1. For some file formats the value of
      DX also depends on the unbounded independent variable and is expressed as
      DX(m,s).
  
      XNAME(s): a character string giving the name and/or description of the
      s-th independent variable, on one line and not exceeding 132 characters.
      Include units of measure and order the independent variable names such
      that, when reading primary variables from the data records, the most
      rapidly varying independent variable is listed first and the most slowly
      varying independent variable is listed last. Currently this is
      Time [Seconds] from midnight on day aircraft flight started for all UND
      exchange files.
  
      NV: number of primary variables in the exchange file (integer). This number
      plus one (for the time value) gives the number of parameters in the data
      file.
  
      VSCAL(n): scale factor by which one multiplies recorded values of the n-th
      primary variable to convert them to the units specified in VNAME(n).
      Currently this is 1 for all UND Citation Aircraft recorded values.
  
      VMISS(n): a quantity indicating missing or erroneous data values for the
      n-th primary variable. VMISS(n) must be larger than any "good" data value,
      of the n-th primary variable, recorded in the file. The value of VMISS(n)
      defined in the file header is the same value that appears in the data
      records for missing/bad values of V(X,n).
  
      VNAME(n): a character string giving the name and/or description of the
      n-th primary variable, on one line and not exceeding 132 characters.
      Include units of measure the data will have after multiplying by the n-th
      scale factor, VSCAL(n). The order in which the primary variable names are
      listed in the file header is the same order in which the primary variables
      are read from the data records, and the same order in which scale factors
      and missing values for the primary variables are read from the file header
      records.
  
      NSCOML: number of special comment lines (integer) within the file header.
      Special comments are reserved to note special problems or circumstances
      concerning the data within a specific exchange file so they may easily be
      found and flagged by those reading the file. If NSCOML=0 then there are
      no special comment lines.
  
      NNCOML: number of normal comment lines (integer) within the file header,
      including blank lines and data column headers, etc. Normal comments are
      those which apply to all of a particular kind of dataset, and can be used
      to more completely describe the contents of the file. If NNCOML=0 then
      there are no normal comment lines.
  
      DTYPE: version description of the data. Typically either "Preliminary Data"
      or "Final Data".
  
      VFREQ: time frequency of the data in the format "## Hz Data," where ## is
      the integer frequency of the data.
  
      VDESC: a character string on a single line containing a short description
      of each variable in the exchange file. No spaces are allowed in each short
      variable description.
  
      VUNITS: a character string on a single line containing the units of each
      variable in the exchange file. No spaces are allowed in each unit's
      description.

      name: the name of the file to read from (for ReadFile method) or write to
      (for WriteFile method).

      data: a dictionary containing the data within the file. The keys of the 
      dictionary are the VDESC values and the values of the dictionary are arrays 
      of floats that correspond to the columns under each VDESC value.

      mvc: a dictionary containing the missing value codes for each parameter
      within the file.  The keys of the dictionary are the VDESC values and the
      values of the dictionary are the corresponding VMISS values.

      format: an string array containing the output data formats depending on
      the format of the corresponding missing value codes.

      start_time: the start time of the data in Universal Coordinated Time (UTC).

      end_time: the end time of the data in Universal Coordinated Time (UTC).

      duration: a string representation of the duration of the file in
      hour-minute-second (HH:MM:SS) format.

      jul_date: the Julian Date based 'name' if it is in the correct format.

    A detailed description of the methods within the class are described below.

      __init__(): constructor for the ADPAA class

      ReadFile(name='input_filename'): calls the external ReadFile method with
      the optional name argument to read a NASA-formatted ASCII file into the
      ADPAA class. 

      WriteFile(name='output_filename'): calls the external WriteFile method
      with the optional name argument write information contained in an ADPAA
      object to a properly-formatted NASA ASCII file.

      constants(constants_arr, jd=0.0): calls the external constants method with
      the optional jd argument to parse aircraft constants XML files and return
      desired aircraft constants in a dictionary.

      mask_MVC(): masks the data contents with missing values.

      create_format_arr(mvc_arr): creates a string array containing the output
      data formats depending on the format of the corresponding missing value
      codes within the "mvc_arr" string array.

      sfm2hms(sfm): converts time from seconds from midnight (sfm) to
      hour-minute-second (HH:MM:SS) form. 'sfm' can be a string, float, or
      integer.

      hms2sfm(hh, mm, ss): converts time from hour-minute-second (HH:MM:SS) form
      to seconds from midnight (sfm). 'hh,' 'mm,' and 'ss' can be strings,
      floats, or integers.

      calculate_jul_date(dt=[]): calculate the Julian Date based on the 'name'
      class attribute or elements in the optional 'dt' array. Elements in the
      'dt' array can be strings, integers, or floats.

      help('att'): gives help on ADPAA class attribute(s) using the optional
      'att' class attribute. Note that help(ADPAA_object) will give help on
      available methods to the currently defined ADPAA object.

    When creating a new data file from scratch, only the SNAME, MNAME, VMISS
    VNAME, VDESC, VUNITS, and data class attributes are necessary to set in
    order to print data to file.  Depending on the data, VFREQ may also have
    to be set manually in order to be correct.

  COPYRIGHT:
    2016, 2017, 2018 David Delene <delene@aero.und.edu>

    This program is distributed under terms of the GNU General Public License
 
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

import math
import numpy
import numpy.ma as ma
import re
import sys
from os import path
from datetime import date

from readfile import ReadFile
from writefile import WriteFile
from julian_date import julian_date
from adpaa_help import adpaa_help

class ADPAA(object):
    def __init__(self):
        """
        Syntax: A=ADPAA()
        Purpose: create default object of the ADPAA class. 
        Input: nothing
        Output: a default ADPAA object.
        Notes: Uppercase class attributes are included in the NASA-formatted
               header and lower case attributes are derived from the information
               within the file.
        """
        # Define class attributes.
        self.NLHEAD = 0
        self.FFI = 1001
        self.ONAME = 'Delene, David'
        self.ORG = 'University of North Dakota'
        self.SNAME = ''
        self.MNAME = ''
        self.IVOL = 1
        self.VVOL = 1
        self.DATE = ''
        self.RDATE = date.today().strftime('%Y %m %d') # today's date as YYYY MM DD
        self.DX = 1.0
        self.XNAME = 'Time [seconds]; UT seconds from midnight on day aircraft flight started'
        self.NV = 0
        self.VSCAL = [1.0]
        self.VMISS = ['']
        self.VNAME = ['']
        self.NSCOML = 0
        self.NNCOML = 4
        self.DTYPE = 'Preliminary Data'
        self.VFREQ = ''
        self.VDESC = ['']
        self.VUNITS = ['']
        self.data = {}
        self.duration = ''
        self.mvc = {}
        self.format = ['']
        self.name = ''
        self.start_time = ''
        self.end_time = ''
        self.jul_date = 0

        # Check versions of Python and its packages for compatibility.
        self.__check_versions()


    def ReadFile(self, name=''):
        """
        Syntax: ReadFile([name= ])
        Purpose: Read a NASA-formatted ASCII file into the ADPAA class.
        Input: 'name' - optional argument to set the filename of the file to
                        input (uses ADPAA class attribute 'name' otherwise).
        Output: nothing
        """
        if (name == '') & (self.name == ''):
            print("WARNING (adpaa.py): class attribute 'name' not set--returning"+\
                  " default ADPAA object.")
            return
        if (name == '') & (self.name != ''):
            name = self.name

        # Call the ReadFile method.
        ReadFile(self, name)

        # Set the Julian Date of the file based on filename.
        self.jul_date = self.calculate_jul_date()


    def WriteFile(self, name=''):
        """
        Syntax: WriteFile([name= ])
        Purpose: Write information contained in an ADPAA object to a properly-
                 formatted NASA ASCII file. 
        Input: 'name' - optional argument to set the filename of the file being
                        written to (uses ADPAA class attribute 'name' otherwise)
        Output: nothing
        """
        if (name == '') & (self.name == ''):
            print("ERROR: Cannot write to file because class attribute 'name'"+\
                  " is not set.")
            return
        if (name == '') & (self.name != ''):
            name = self.name

        # Set number of variables (self.NV) according to the length of the
        # long name array (self.VNAME).
        self.NV = len(self.VNAME)

        # Check that all data has proper data types and is formatted properly.
        self.__check_data(name)

        # Call the WriteFile method.
        WriteFile(self, name)


    def constants(self, constants_arr, jd=0.0):
        """
        Syntax: constants(constants_arr, [jd= ])
        Purpose: Parse aircraft constants XML files and return desired aircraft
                 constants in a dictionary.
        Input: 'constants_arr' - an array of desired constants as strings
               'jd' - optional argument sets Julian Date for time-dependent constants
        Output: a dictionary containing the desired constants that are correct for
                the given Julian Date.
        Notes: If no Julian Date is given, the class attribute 'jul_date' is used.
        """
        # Import the constants method.
        from constants import constants

        # Parse the default constants XML file for the desired constants.
        constants_dict = constants(constants_arr)

        # Get calibration data based on the ADPAA class attribute 'SNAME'.
        if 'Convair 580 (C-FNRC)' in self.SNAME:

            # Create a dictionary 'N_consts' with the desired constants from the
            # N555DS XML constants file using the Julian Date passed into
            # the method.
            if jd > 0.0:
                N_consts = constants(constants_arr, \
                                     fname='C-FNRC_constants.xml', jul_date=jd)

            # Create a dictionary 'N_consts' with the desired constants from the
            # N555DS XML constants file using the Julian Date calculated from
            # the filename.
            elif self.jul_date > 0.0:
                N_consts = constants(constants_arr,                \
                                     fname='C-FNRC_constants.xml', \
                                     jul_date=self.jul_date)

            else:
                print("WARNING (adpaa.py): Julian Date not set--returning only desired"+\
                      " default constants. Use the 'jd' argument to manually"+\
                      " set the Julain Date.")
                return constants_dict

            # Update the 'constant_dict' dictionary with the desired constants
            # only if a value exists for that item in the 'N_const' dictionary
            # as to not override the values already in 'const_dict'.
            for key,val in N_consts.items():
                if val != None:
                    constants_dict.update({key:val})

        if 'UND Citation II' in self.SNAME:

            # Create a dictionary 'N_consts' with the desired constants from the
            # N555DS XML constants file using the Julian Date passed into
            # the method.
            if jd > 0.0:
                N_consts = constants(constants_arr, \
                                     fname='N555DS_constants.xml', jul_date=jd)

            # Create a dictionary 'N_consts' with the desired constants from the
            # N555DS XML constants file using the Julian Date calculated from
            # the filename.
            elif self.jul_date > 0.0:
                N_consts = constants(constants_arr,                \
                                     fname='N555DS_constants.xml', \
                                     jul_date=self.jul_date)

            else:
                print("WARNING (adpaa.py): Julian Date not set--returning only desired"+\
                      " default constants. Use the 'jd' argument to manually"+\
                      " set the Julain Date.")
                return constants_dict

            # Update the 'constant_dict' dictionary with the desired constants
            # only if a value exists for that item in the 'N_const' dictionary
            # as to not override the values already in 'const_dict'.
            for key,val in N_consts.items():
                if val != None:
                    constants_dict.update({key:val})

        else:
            print("WARNING (adpaa.py): invalid ADPAA class attribute 'SNAME'--only "+\
                  "desired default constants will be returned.")

        # Remove constants that have not been found in any constants file.
        for k in constants_dict.keys():
            if constants_dict[k] == None:
                print("WARNING (adpaa.py): constant '"+k+"' does not exist and will not"+\
                      " be returned.")
                del constants_dict[k]

        return constants_dict


    def mask_MVC(self):
        """
        Syntax: mask_MVC()
        Purpose: Mask the data contents with missing values.
        Input: nothing
        Output: nothing
        Notes: Resetting into the initial view needs reinstatiation of the object.

               numpy.ma reference says masked_values() should be used for
               floating-point missing values, however ASCII-files can contain
               integer MVCs (i.e for counts - 9999999). Although masks correctly
               applied, this situation might be further considered.

               With this configuration, integer MVCs are converted to floats.
               So far it's working fine.
        """
        _i = 0
        for keys in self.VDESC[1:]:
            # Mask the missing values.
            self.data[keys] = ma.masked_values(self.data[keys],\
                                               float(self.VMISS[_i]))
            _i += 1


    def create_format_arr(self, mvc_arr):
        """
        Syntax: create_format_arr(mvc_arr)
        Purpose: Creates an array of output data format strings based on the
                 corresponding missing value codes (MVCs).
        Input: 'mvc_arr' - an array of format strings of the MVCs.
        Output: an array of strings that represent the data output format.
        Note: Underscored variables are private to the class.
        """
        # Define output data format for the time data.
        _outfmt = [' {0:11.4f}']

        # Find the output data format for all other data by appending a format
        # string to the _outfmt array depending on the format of the MVC.
        for i in range(len(mvc_arr)):
            if ('E' in mvc_arr[i]) | ('e' in mvc_arr[i]):
                # If the MVC is in scientific notation, keep corresponding data
                # in scientific notation.
                _suffix = len(re.split('\.|E|e',mvc_arr[i])[1])
                _outstr = ' {0:11.'+str(_suffix)+'e}'
                _outfmt.append(_outstr)
            elif '.' in mvc_arr[i]:
                # If the MVC is a float, keep corresponding data in same format
                # (same number of digits behind the decimal).
                _suffix = len(str.split(mvc_arr[i],'.')[1])
                _outstr = ' {0:11.'+str(_suffix)+'f}'
                _outfmt.append(_outstr)
            else:
                # If the MVC is an integer, pad corresponding data with spaces
                # in front.
                _strlen = len(mvc_arr[i])
                _padnum = 11-_strlen
                _padding = ''
                for i in range(_padnum):
                    _padding += ' ' 
                _outstr = ' '+_padding+'{0:'+str(_strlen)+'.0f}'
                _outfmt.append(_outstr)

        return _outfmt


    def sfm2hms(self, sfm):
        """
        Syntax: sfm2hms(sfm)
        Purpose: Converts time from seconds from midnight (sfm) to
                 hour-minute-second (HH:MM:SS) form.
        Input: 'sfm' - a string, float, or integer representing the time in sfm.
        Output: the HH:MM:SS representation of the input time.
        """
        h  = float(sfm) / 3600.0
        hh = math.floor(h)
        m  = (h - hh) * 60.0
        mm = math.floor(m)
        ss = (m - mm) * 60.0
        
        # Add leading zeros to single digit times if needed.
        str_hh='%02.0f'%hh
        str_mm='%02.0f'%mm
        str_ss='%07.4f'%ss

        return '%s:%s:%s'%(str_hh,str_mm,str_ss)


    def hms2sfm(self, hh, mm, ss):
        """
        Syntax: hms2sfm(hh, mm, ss)
        Purpose: Converts time from hour-minute-second (HH:MM:SS) form to
                 seconds from midnight (sfm).
        Input: 'hh' - a string, float, or integer that represents the hour.
               'mm' - a string, float, or integer that represents the minute.
               'ss' - a string, float, or integer that represents the second.
        Output: the sfm representation of the input time.
        """
        return (float(hh)*3600.0) + (float(mm)*60.0) + float(ss)


    def calculate_jul_date(self, dt=[], quiet=0):
        """
        Syntax: calculate_jul_date([dt=[YY,MM,DD,hh,mm,ss],] [quiet=0|1])
        Purpose: Calculates the Julian Date based on the 'name' class attribute
                 or elements in a given array.
        Input: 'dt'    - optional array of strings, integers, or floats representing
                         the date and time in which to calculate the Julian Date.
               'quiet' - optional keyword (0 or 1) to suppress error output.
        Output: the Julian Date or 0.0 if Julian Date cannot be calculated.
        Notes: 'dt' array needs six elements in this order: year (YY), month (MM),
               day (DD), hour (hh), minute (mm), and second (ss). All digits cannot
               be more than than two digits in length.
        """
        # Set default Julian Date.
        jul_date = 0.0

        if len(dt) == 0:
            # Try to calculate the Julian Date from the 'name' class attribute
            # using the julian_date module. If the 'name' class attribute is not
            # in *YY_MM_DD_hh_mm_ss.* format, print a warning message and return
            # a Julian Date of 0.0.
            try:
                # Write *YY_MM_DD_HH_MM_SS.* format using regular expressions so
                # the wanted date/time can be searched for in the 'name'
                # class attribute.
                #            Y    Y       M    M       D    D
                fmt='(.|^)([0-9][0-9][_][0-1][0-9][_][0-3][0-9][_]' + \
                          '[0-5][0-9][_][0-5][0-9][_][0-5][0-9])'
                #            h    h       m    m       s    s

                # Search for and parse the date/time in the 'name' class
                # attribute using the re module.
                srch = re.search(fmt, self.name)
                HMS_str = srch.group(2)

                # Split date/time into parts and calculate the Julian Date using
                # the julian_date module.
                file_date = str.split(HMS_str, '_')
                jul_date = julian_date(int(file_date[0]), int(file_date[1]), \
                                       int(file_date[2]), int(file_date[3]), \
                                       int(file_date[4]), int(file_date[5]))
            except:
                if quiet == 0:
                    print("WARNING (adpaa.py): ADPAA class attribute 'name' not in "+\
                          "*YY_MM_DD_hh_mm_ss.* format so Julian Date cannot be "+\
                          "calculated automatically.")

        elif (len(dt) == 6) & (int(dt[0]) < 100):
            # Calculate the Julian Date from the given 'dt' argument as long as
            # all date and time information is provided and the year is not more
            # than two digits in length.
            jul_date = julian_date(int(dt[0]), int(dt[1]), int(dt[2]), \
                                   int(dt[3]), int(dt[4]), int(dt[5]))
        else:
            if quiet == 0:
                print("ERROR: Cannot calculate Julian Date. The six elements in"+\
                      " 'dt' can be strings, floats, or integers and cannot"+\
                      " exceed two digits in length.")

        return jul_date


    def help(self, att=''):
        """
        Syntax: help([att])
        Purpose: Displays help information about one or all ADPAA class attributes.
        Input: 'att' - optional ADPAA class attribute
        Output: Scrollable help page (Linux only) with information about all ADPAA class
                attributes or a help message with information about the specified
                ADPAA class attribute.
        """
        default_adpaa = ADPAA()
        adpaa_help(self, default_adpaa, att)


    def __check_versions(self):
        """
        Syntax: private to class (can only be called from within the class)
        Purpose: Checks the versions of Python and the necessary packages.
        Input: nothing
        Output: error message(s) and exits program if compatibility is not met.
        """
        # Set error flag.
        error = 0

        # Check if Python is version 2.6 or greater
        if (sys.version_info[0] + sys.version_info[1]) < 8:
            print("ERROR: The ADPAA class needs Python version 2.6 or greater."+\
                  " Current Python version is "+str(sys.version_info[0])+"."+\
                   str(sys.version_info[1])+".")
            error = 1

        # Check if numpy is version 1.4 or greater
        np_ver = numpy.version.version.split('.')
        if (int(np_ver[0]) + int(np_ver[1])) < 5:
            print("ERROR: The ADPAA class needs numpy version 1.4 or greater."+\
                  " Current numpy version is "+numpy.version.version+".")
            error += 1

        # If an error message is found, exit program.
        if error >= 1:
            sys.exit(1)


    def __check_data(self, name):
        """
        Syntax: private to class (can only be called from within the class)
        Purpose: checks data types and formatting prior to writing data to file.
        Input: 'name' - the name of the output file.
        Output: error message(s) and exits program if conditions are not met.
        """
        # -------------------------------------------------------------------- #
        # Check for correct data types.                                        #
        # -------------------------------------------------------------------- #
        dtype_out = ['\nData type errors:']

        if type(self.ONAME) is not str:
            dtype_out.append('  ONAME must be a string.')
        if type(self.ORG) is not str:
            dtype_out.append('  ORG must be a string.')
        if type(self.SNAME) is not str:
            dtype_out.append('  SNAME must be a string.')
        if type(self.MNAME) is not str:
            dtype_out.append('  MNAME must be a string.')
        if type(self.IVOL) is not int:
            dtype_out.append('  IVOL must be an integer.')
        if type(self.VVOL) is not int:
            dtype_out.append('  VVOL must be an integer.')
        if type(self.DATE) is not str:
            dtype_out.append('  DATE must be a string.')
        if type(self.DX) is not float:
            dtype_out.append('  DX must be a float.')
        if type(self.XNAME) is not str:
            dtype_out.append('  XNAME must be a string.')
        if type(self.VSCAL) is not list:
            dtype_out.append('  VSCAL must be an array (list) of strings or floats.')
        if type(self.VNAME) is not list:
            dtype_out.append('  VNAME must be an array (list) of strings.')
        if type(self.NSCOML) is not int:
            dtype_out.append('  NSCOML must be an integer.')
        if type(self.NNCOML) is not int:
            dtype_out.append('  NNCOML must be an integer.')
        if type(self.DTYPE) is not str:
            dtype_out.append('  DTYPE must be a string.')
        if type(self.VFREQ) is not str:
            dtype_out.append('  VFREQ must be a string.')
        if type(self.VDESC) is not list:
            dtype_out.append('  VDESC must be an array (list) of strings.')
        if type(self.VUNITS) is not list:
            dtype_out.append('  VUNITS must be an array (list) of strings.')
        if type(self.data) is not dict:
            dtype_out.append('  data must be a dictionary.')
        else:
            for key in self.data.keys():
                if type(self.data[key]) is not numpy.ndarray:
                    dtype_out.append("  data['"+key+"'] must be a numpy array "+\
                                     "(numpy.array()).")
                else:
                    try:
                        if self.data[key].dtype is not numpy.dtype('float64'):
                            dtype_out.append("  data['"+key+"'] must be a numpy "+\
                                             "array (numpy.array()) of floats "+\
                                             "(numpy.dtype('float64')).")
                    except:
                        dtype_out.append("  data['"+key+"'] must be a numpy "+\
                                         "array (numpy.array()) of floats "+\
                                         "(numpy.dtype('float64')).")
        if type(self.VMISS) is not list:
            dtype_out.append('  VMISS must be an array (list) of strings.')
        else:
            for i in range(len(self.VMISS)):
                if type(self.VMISS[i]) is not str:
                    dtype_out.append('  VMISS['+str(i)+'] must be a string.')

        if len(dtype_out) == 1:
            del dtype_out[0]

        # -------------------------------------------------------------------- #
        # Check formatting.
        # -------------------------------------------------------------------- #
        format_out = ['\nFormatting errors:']

        # Set VSCAL to proper length (length of VNAME) if still at default value.
        if (type(self.VSCAL) is list) & (type(self.VNAME) is list):
            if len(self.VSCAL) == 1:
                self.VSCAL = [1.0] * len(self.VNAME)
            
        # Check for proper lengths of VSCAL, VMISS, VNAME, VDESC, VUNITS, and data.
        if (type(self.VSCAL) is list) & (type(self.VMISS) is list) & \
           (type(self.VNAME) is list) & (type(self.VDESC) is list) & \
           (type(self.VUNITS) is list) & (type(self.data) is dict):
            if len(self.VSCAL) == len(self.VMISS) == len(self.VNAME) == \
               len(self.VDESC)-1 == len(self.VUNITS)-1 == len(self.data)-1:
               foo=1
            else:
                format_out.append('  Array lengths are not equal. Note that '+\
                                  'the length of VNAME may be the issue.')
                if len(self.VSCAL) != self.NV:
                    format_out.append('    len(VSCAL)  = '+str(len(self.VSCAL))+\
                                      ' -- must be equal to len(VNAME) or NV: '+\
                                      str(self.NV))
                else:
                    format_out.append('    len(VSCAL)  = '+str(len(self.VSCAL))+\
                                      ' -- correct (may have been set by default)')
                if len(self.VMISS) != self.NV:
                    format_out.append('    len(VMISS)  = '+str(len(self.VMISS))+\
                                      ' -- must be equal to len(VNAME) or NV: '+\
                                      str(self.NV))
                else:
                    format_out.append('    len(VMISS)  = '+str(len(self.VMISS))+\
                                      ' -- correct')
                if len(self.VDESC) != self.NV+1:
                    format_out.append('    len(VDESC)  = '+str(len(self.VDESC))+\
                                      ' -- must be equal to len(VNAME)+1 or NV+1: '+\
                                      str(self.NV+1))
                else:
                    format_out.append('    len(VDESC)  = '+str(len(self.VDESC))+\
                                      ' -- correct')
                if len(self.VUNITS) != self.NV+1:
                    format_out.append('    len(VUNITS) = '+str(len(self.VUNITS))+\
                                      ' -- must be equal to len(VNAME)+1 or NV+1: '+\
                                      str(self.NV+1))
                else:
                    format_out.append('    len(VUNITS) = '+str(len(self.VUNITS))+\
                                      ' -- correct')
                if len(self.data) != self.NV+1:
                    format_out.append('    len(data)   = '+str(len(self.data))+\
                                      ' -- must be equal to len(VNAME)+1 or NV+1: '+\
                                      str(self.NV+1))
                else:
                    format_out.append('    len(data)   = '+str(len(self.data))+\
                                      ' -- correct')

        # Check that all keys in VDESC are in data and check that all values within
        # data are numpy arrays of floats.
        if (type(self.VDESC) is list) & (type(self.data) is dict):
            for key in self.VDESC:
                if key not in self.data:
                    format_out.append("  data['"+key+"'] does not exist. Check "+\
                                      "that all items in VDESC match the keys "+\
                                      "in data dictionary.")
                elif (type(self.data[key]) is numpy.ndarray) & \
                     (self.VDESC[0] in self.data):
                    if len(self.data[key]) != len(self.data[self.VDESC[0]]):
                        format_out.append("  data['"+key+"'] must be same length "+\
                                          "as first data column: len(data['"+key+\
                                          "'])="+str(len(self.data[key]))+\
                                          ", len(data['"+self.VDESC[0]+"'])="+\
                                          str(len(self.data[self.VDESC[0]]))+".")

        # Check format of DATE if not default value, otherwise attempt to
        # automatically set DATE based on filename of output file.
        if type(self.DATE) is str:
            if self.DATE != '':
                try:
                    #           Y    Y    Y    Y       M    M       D    D
                    fmt   = '([1-2][0-9][0-9][0-9][ ][0-1][0-9][ ][0-3][0-9])'
                    match = re.match(fmt, self.DATE)
                    date  = match.group(0)
                except:
                    format_out.append('  DATE must be formatted as follows: YYYY MM DD')
            else:
                self.name = name
                self.jul_date = self.calculate_jul_date(quiet=1) # check filename format
                if self.jul_date > 0.0:
                    name_splt = str.split(name, '.')
                    date_splt = str.split(name_splt[0], '_')

                    if int(date_splt[0]) < 50:
                        self.DATE = '20'+date_splt[0]+' '+date_splt[1]+' '+date_splt[2]
                    else:
                        self.DATE = '19'+date_splt[0]+' '+date_splt[1]+' '+date_splt[2]
                else:
                    format_out.append('  DATE cannot be set because filename is not '+\
                                      'in YY_MM_DD_hh_mm_ss.* format. Please set DATE '+\
                                      'in the format as follows: YYYY MM DD')

        # Check for proper formats of FFI, SNAME, and MNAME.
        if self.FFI != 1001:
            format_out.append('  FFI must equal the integer 1001.')
        if self.SNAME == '':
            format_out.append('  SNAME cannot be an empty string.')
        if self.MNAME == '':
            format_out.append('  MNAME cannot be an empty string.')

        # Attempt to automatically set VFREQ if still default value and Time is
        # in the data dictionary, otherwise check formatting.
        # Need at least 1 second of data for VFREQ to be set automatically.
        # Asynchronous data may not throw an error if attempting to set VFREQ automatically.
        if (type(self.VFREQ) is str) & (type(self.data) is dict):
            if (self.VFREQ == '') & ('Time' in self.data):
                dur  = self.data['Time'][-1] - self.data['Time'][0]
                freq = len(self.data['Time']) / dur
                rem  = freq % 1
                if (dur > 1.0) & (freq > 1.0) & (rem >= 0.0):
                    self.VFREQ = '%d Hz Data' % (freq - rem)
                else:
                    format_out.append('  VFREQ must not be an empty string and '+\
                                      'cannot be set automatically. Please set '+\
                                      'VFREQ in the format as follows: ## Hz Data '+\
                                      '(where ## is the frequency of the data).')

        # Check that all elements in VDESC do not exceed 11 characters in length.
        if type(self.VDESC) is list:
            for vdesc in self.VDESC:
                if len(vdesc) > 11:
                    format_out.append("  VDESC element '"+vdesc+"' must not exceed "+\
                                      "11 characters in length.")

        # Check that all elements in VUNITS do not exceed 11 characters in length.
        if type(self.VUNITS) is list:
            for vunits in self.VUNITS:
                if len(vunits) > 11:
                    format_out.append("  VUNITS element '"+vunits+"' must not exceed "+\
                                      "11 characters in length.")

        # Check that all VDESC elements are unique.
        if type(self.VDESC) is list:
            if len(self.VDESC) != len(set(self.VDESC)):
                format_out.append('  VDESC must contain unique shortnames. Check that '+\
                                  'no VDESC elements are alike.')

        # Delete formatting error header if no formatting errors exist.
        if len(format_out) == 1:
            del format_out[0]

        # Output errors and exit program if errors exist.
        if (len(dtype_out) > 1) | (len(format_out) > 1):
            print('ERROR (writefile.py/adpaa.py):')
            print('Wrong data type(s) and/or formatting detected when writing'+\
                  ' data to file. Please review list of issues and try again.')
            for do in dtype_out:
                print(do)
            for fo in format_out:
                print(fo)

            sys.exit(1)
