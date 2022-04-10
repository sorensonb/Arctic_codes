"""
  NAME:
    adpaa_help.py
 
  PURPOSE:
    Displays the defintion of a single or all ADPAA class attribute(s) and the
    currently assigned value(s).
    
  CALLS:
    Built-in Python modules numpy, sys, and subprocess

  MODIFICATIONS:
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2018/07/26: Written

  USAGE:
    Create a default ADPAA object by
      A = ADPAA()

    Display the defintion of a single ADPAA class attribute and its currently
    assigned value by
      A.help('att')

    Display the definitions of all ADPAA class attributes and their currently
    assigned values (Linux only) by
      A.help()

    Display the methods defined for use by an instantiated ADPAA object by
      help(A)

  EXAMPLE:
    Instantiate an ADPAA object
      >>> air_ascent_1Hz = ADPAA()

    Display the defintion of a single ADPAA class attribute and its currently
    assigned value
      >>> A.help('VFREQ')
      VFREQ:
        Time frequency of the data.
        Current value: 1 Hz Data

    Display the definitions of all ADPAA class attributes and their currently
    assigned values
      >>> A.help()
      Opens a scrollable list (Linux only) with all ADPAA class attributes and
      their values.  Instructions to navigate list are included at the top
      of the screen when list opens.

    Display the methods defined for use by an instantiated ADPAA object
      >>> help(A)
      Opens a scrollable list with all the methods within the ADPAA class.


  COPYRIGHT:
    2018 David Delene <delene@aero.und.edu>

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
import numpy as np
import sys
import subprocess

def adpaa_help(self, default_adpaa, att=''):

    # Define text formatting options.
    bold   = '\033[1m'
    normal = '\033[0m'

    # Define the dictionary containing all ADPAA class attributes and their descriptions.
    # All keys are lowercase.  Value is list if attribute needs to be displayed in a
    # different format (capitalized, etc.) or if attribute is natrually capitalized
    # (such as NLHEAD).  Value is string if attribute is in the form to be displayed.
    help_dict = {'nlhead':['NLHEAD','Number of lines (integer) composing the file header. NLHEAD is the\n'+\
                 '  first recorded value on the first line of an exchange file.'],
                 'ffi':['FFI','ASCII file format number. For the UND Citation aircraft data this\n'+\
                 '  will always be 1001.'],
                 'oname':['ONAME','A character string specifying the name(s) of the originator(s) of\n'+\
                 '  the exchange file, last name first. On one line and not exceeding 132\n'+\
                 '  characters.'],
                 'org':['ORG','Character string specifying the organization or affiliation of the\n'+\
                 '  originator of the exchange file. Can include address, phone number, email\n'+\
                 '  address, etc. On one line and not exceeding 132 characters.'],
                 'sname':['SNAME','A character string specifying the source of the measurements or\n'+\
                 '  model results which compose the primary variables, on one line and not\n'+\
                 '  exceeding 132 characters. Can include instrument name, measurement\n'+\
                 '  platform, etc. Required to be defined to write data to file with WriteFile().'],
                 'mname':['MNAME','A character string specifying the name of the field project that\n'+\
                 '  the data were obtained from. Required to be defined to write data to file with WriteFile().'],
                 'ivol':['IVOL','Volume number (integer) of the total number of volumes required to\n'+\
                 '  store a complete dataset, assuming only one file per volume. To be used in\n'+\
                 '  conjunction with VVOL to allow data  exchange of large datasets requiring\n'+\
                 '  more than one volume of the exchange medium (diskette, etc.).'],
                 'vvol':['VVOL','Total number of volumes (integer) required to store the complete\n'+\
                 '  dataset, assuming one file per volume. If VVOL>1 then each volume must\n'+\
                 '  contain a file header with an incremented value for IVOL, and continue the\n'+\
                 '  data records with monotonic independent variable marks.'],
                 'date':['DATE','UT date at which the data within the exchange file begins. For\n'+\
                 '  aircraft data files DATE is the UT date of takeoff. DATE is in the form\n'+\
                 '  YYYY MM DD (year, month, day) with each integer value separated by at least\n'+\
                 '  one space. For example: 2009 09 21.'],
                 'rdate':['RDATE','Date of data reduction or revision, in the same form as DATE.'],
                 'dx':['DX(s)','Interval (real) between values of the s-th independent variable,\n'+\
                 '  X(i,s), i=1,NX(s); in the same units as specified in XNAME(s). DX(s) is\n'+\
                 '  zero for a non-uniform interval. DX(s) is non-zero for a constant interval.\n'+\
                 '  If DX(s) is non-zero then it is required that\n'+\
                 '  NX(s) = (X(NX(s),s)-X(1,s)) / DX(s) + 1. For some file formats the value of\n'+\
                 '  DX also depends on the unbounded independent variable and is expressed as\n'+\
                 '  DX(m,s).'],
                 'xname':['XNAME(s)','A character string giving the name and/or description of the\n'+\
                 '  s-th independent variable, on one line and not exceeding 132 characters.\n'+\
                 '  Include units of measure and order the independent variable names such\n'+\
                 '  that, when reading primary variables from the data records, the most\n'+\
                 '  rapidly varying independent variable is listed first and the most slowly\n'+\
                 '  varying independent variable is listed last. Currently this is\n'+\
                 '  Time [Seconds] from midnight on day aircraft flight started for all UND\n'+\
                 '  exchange files.'],
                 'nv':['NV','Number of primary variables in the exchange file (integer). This number\n'+\
                 '  plus one (for the time value) gives the number of parameters in the data\n'+\
                 '  file.'],
                 'vscal':['VSCAL(n)','Scale factor by which one multiplies recorded values of the n-th\n'+\
                 '  primary variable to convert them to the units specified in VNAME(n).\n'+\
                 '  Currently this is 1 for all UND Citation Aircraft recorded values.'],
                 'vmiss':['VMISS(n)','A quantity indicating missing or erroneous data values for the\n'+\
                 '  n-th primary variable. VMISS(n) must be larger than any "good" data value,\n'+\
                 '  of the n-th primary variable, recorded in the file. The value of VMISS(n)\n'+\
                 '  defined in the file header is the same value that appears in the data\n'+\
                 '  records for missing/bad values of V(X,n). Required to be defined to write data to\n'+\
                 '  file with WriteFile().'],
                 'vname':['VNAME(n)','A character string giving the name and/or description of the\n'+\
                 '  n-th primary variable, on one line and not exceeding 132 characters.\n'+\
                 '  Include units of measure the data will have after multiplying by the n-th\n'+\
                 '  scale factor, VSCAL(n). The order in which the primary variable names are\n'+\
                 '  listed in the file header is the same order in which the primary variables\n'+\
                 '  are read from the data records, and the same order in which scale factors\n'+\
                 '  and missing values for the primary variables are read from the file header\n'+\
                 '  records. Required to be defined to write data to file with WriteFile().'],
                 'nscoml':['NSCOML','Number of special comment lines (integer) within the file header.\n'+\
                 '  Special comments are reserved to note special problems or circumstances\n'+\
                 '  concerning the data within a specific exchange file so they may easily be\n'+\
                 '  found and flagged by those reading the file. If NSCOML=0 then there are\n'+\
                 '  no special comment lines.'],
                 'nncoml':['NNCOML','Number of normal comment lines (integer) within the file header,\n'+\
                 '  including blank lines and data column headers, etc. Normal comments are\n'+\
                 '  those which apply to all of a particular kind of dataset, and can be used\n'+\
                 '  to more completely describe the contents of the file. If NNCOML=0 then\n'+\
                 '  there are no normal comment lines.'],
                 'dtype':['DTYPE',"Version description of the data. Typically either 'Preliminary Data'\n"+\
                 "  or 'Final Data'."],
                 'vfreq':['VFREQ','Time frequency of the data. May need to be manually defined.'],
                 'vdesc':['VDESC','A character string on a single line containing a short description\n'+\
                 '  of each variable in the exchange file. No spaces are allowed in each short\n'+\
                 '  variable description. Required to be defined to write data to file with WriteFile().'],
                 'vunits':['VUNITS','A character string on a single line containing the units of each\n'+\
                 '  variable in the exchange file. No spaces are allowed in each unit\n'+\
                 '  description. Required to be defined to write data to file with WriteFile().'],
                 'name':'The name of the file to read from (for ReadFile methods) or write\n'+\
                 '  to (for WriteFile method).',
                 'data':'A dictionary containing the data within the file. The keys of the\n'+\
                 '  dictionary are the VDESC values and the values of the dictionary are arrays\n'+\
                 '  of floats that correspond to the columns under each VDESC value. Required to be\n'+\
                 '  defined to write data to file with WriteFile().',
                 'mvc':'A dictionary containing the missing value codes for each parameter\n'+\
                 '  within the file.  The keys of the dictionary are the VDESC values and the\n'+\
                 '  values of the dictionary are the corresponding VMISS values.',
                 'format':'An string array containing the output data formats depending on\n'+\
                 '  the format of the corresponding missing value codes.',
                 'start_time':'The start time of the data in Universal Coordinated Time (UTC).',
                 'end_time':'The end time of the data in Universal Coordinated Time (UTC).',
                 'duration':'A string representation of the duration of the file in\n'+\
                 '  hour-minute-second (HH:MM:SS) format.',
                 'jul_date':"The Julian Date based on 'name' if it is in the correct format."}

    # Print information based on user input.
    if att == '': # if no class attribute was given
        out = 'ADPAA Class Attribute Help\n\n'+\
              'A detailed description of each class attribute is given below.  Uppercase\n'+\
              'attributes come directly from the file header and lower case attributes are\n'+\
              'derived from the information within the file.\n\n'+\
              'Use down arrow (or j) to scroll down and up arrow (or k) to scroll up.\n'+\
              'Use / to search and n to show next match or N for previous match.\n'+\
              'Press q to quit. DO NOT USE CTRL+C.\n\n'+\
              "For help with the methods of the ADPAA class, use Python's help() method.\n\n"

        for h in help_dict:
            if type(help_dict[h]) is list:
                # Print attribute name and description.
                out = out + help_dict[h][0]+':\n'
                out = out + '  '+help_dict[h][1]+'\n'

                # Get current value of attribute.
                exec('bar=self.'+h.upper())
                exec('foobar=default_adpaa.'+h.upper())
            else:
                # Print attribute name and description.
                out = out + h+':\n'
                out = out + '  '+help_dict[h]+'\n'

                # Get current value of attribute.
                exec('bar=self.'+h)
                exec('foobar=default_adpaa.'+h)

            # Print current value of attribute.
            if bar == foobar:
                out = out + '  Current value: '+str(bar)+' (default value)\n'
            elif (type(bar) is str) | (type(bar) is float) | (type(bar) is int):
                out = out + '  Current value: '+str(bar)+'\n'
            elif type(bar) is list:
                out = out + '  Current value: list of length '+str(len(bar))+'\n'
            elif type(bar) is dict:
                out = out + '  Current value: dict of length '+str(len(bar))+'\n'
            elif type(bar) is np.ndarray:
                out = out + '  Current value: numpy.ndarray of length '+str(len(bar))+'\n'
            else:
                out = out + '  Current value: unknown\n'

        # Print out help using the native 'less' command--Linux only.
        # This implementation does not allow for quitting the less command with anything other
        # than the 'q' key. If something else is used (i.e. CTRL+C), the terminal will have to
        # be closed in order to properly quit the 'less' command.
        try:
            pager = subprocess.Popen(['less'], stdin=subprocess.PIPE, stdout=sys.stdout)
            pager.communicate(input=out)
            pager.stdin.close()
        except KeyboardInterrupt:
            pager.kill()
            pager.stdin.close()
            print("\nADPAA.help() did not close properly. Please close terminal before continuing.")
            sys.exit(1)

    else: # if a class attribute was given
        if type(help_dict[att.lower()]) is list:
            # Print attribute description.
            print(bold+help_dict[att.lower()][0]+':'+normal)
            print('  '+help_dict[att.lower()][1])

            # Get current value of attribute.
            exec('bar=self.'+att.upper())
            exec('foobar=default_adpaa.'+att.upper())
            if bar == foobar:
                print('  Current value: '+str(bar)+' (default value)')
            else:
                print('  Current value: '+str(bar))
        else:
            # Print attribute description.
            print(bold+att.lower()+':'+normal)
            print('  '+help_dict[att.lower()])

            # Get current value of attribute.
            exec('bar=self.'+att.lower())
            exec('foobar=default_adpaa.'+att.lower())
            if bar == foobar:
                print('  Current value: '+str(bar)+' (default value)')
            else:
                print('  Current value: '+str(bar))
