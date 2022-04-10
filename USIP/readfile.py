#!/usr/bin/env python
"""
  NAME:
    readfile.py
 
  PURPOSE:
    Reads data from a file whose filename is set either by the optional 'name'
    argument or by the ADPAA class attribute 'name' and stores data within the
    ADPAA object the method was called upon.

  CALLS:
    Built-in Python modules numpy and gzip

  MODIFICATIONS:
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/09/23: Written
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2010/28/16: Added comments
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/11/23: 
       Fixed so new-line characters are stripped from all strings.
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/11/29:
       Fixed VNAME so all new-line characters are stripped.
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2016/11/30:
       Added support for files compressed with gzip.
    Nicholas Gapp <nicholas.james.gapp@und.edu> - 2018/07/26:
       Fixed to print filename of file that cannot be read/opened.

  USAGE:
    Create a default APDAA object by
      A = ADPAA()

    Read data from NASA-formatted ASCII file by
      Specifying filename:
        A.ReadFile('input_filename')
      Using filename set by class attribute 'name':
        A.ReadFile()

  EXAMPLE:
    Instantiate an ADPAA object
      >>> air_ascent_1Hz = ADPAA()

    Read in data from the 15_08_08_14_57_48.air.ascent.1Hz file
      >>> air_ascent_1Hz.ReadFile('15_08_08_14_57_48.air.ascent.1Hz')

      --OR--
      >>> air_ascent_1Hz.name = '15_08_08_14_57_48.air.ascent.1Hz'
      >>> air_ascent_1Hz.ReadFile()


  COPYRIGHT:
    2016 David Delene <delene@aero.und.edu>

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
import gzip

def ReadFile(self, filename):

    # Underscored variables are private to the class.
    try:
        # Open the input file for reading.
        # Use the gzip module to read in compressed files.
        if '.gz' in filename:
            _fhand = gzip.open(filename, 'r')
        else:
            _fhand = open(filename, 'r')
        self.name = filename

        # Need skiprows to read the right section of the file.
        #_skiprows = int(_fhand.next().split()[0])
        flines = _fhand.readlines()
        _skiprows = int(flines[0].strip().split()[0])

        # Move file pointer to top of file and grab all header information
        # for parsing.
        _fhand.seek(0)
        _header = [_fhand.readline() for lines in range(_skiprows)]

        # Parse the header information.
        # All new-line characters are stripped from strings.
        self.NLHEAD = _skiprows
        self.FFI    = int(_header[0].split()[1])
        self.ONAME  = _header[1].rstrip()
        self.ORG    = _header[2].rstrip()
        self.SNAME  = _header[3].rstrip()
        self.MNAME  = _header[4].rstrip()
        self.IVOL   = int(_header[5].split()[0])
        self.VVOL   = int(_header[5].split()[1])
        self.DATE   = _header[6][:10]
        self.RDATE  = _header[6][11:].rstrip()
        self.DX     = float(_header[7])
        self.XNAME  = _header[8].rstrip()
        self.NV     = int(_header[9])
        self.VSCAL  = _header[10].split()
        self.VMISS  = _header[11].split()
        self.VNAME  = [_header[12+i].rstrip() for i in range(self.NV)]

        self.NSCOML = int(_header[self.NV+12])
        self.NNCOML = int(_header[self.NV+13])
        self.DTYPE  = _header[self.NV+14].rstrip()
        self.VFREQ  = _header[self.NV+15].rstrip()

        self.VDESC  = _header[self.NV+16].split()
        self.VUNITS = _header[self.NV+17].split()

        # Read data values from file.
        _data = np.loadtxt(filename, dtype='float', \
                           skiprows=self.NLHEAD).T

        # Store data in dictionary.
        # Use "object_name.data['Time']" syntax to access data.
        self.data = dict(zip(self.VDESC, _data))

        # Find starting and ending times of the file and calculate the total
        # file duration.
        _strt = self.data[self.VDESC[0]][0]
        _endt = self.data[self.VDESC[0]][-1]
        self.start_time = self.sfm2hms(_strt)
        self.end_time   = self.sfm2hms(_endt)
        self.duration   = self.sfm2hms(_endt - _strt)

        # Missing Value Code (MVC) access.
        self.mvc = dict(zip(self.VDESC[1:], self.VMISS))

        # Create the string of output formats based on the MVCs.
        self.format = self.create_format_arr(self.VMISS)

        # Close input file.
        _fhand.close()

    except:
        print("Cannot read ASCII file: '"+filename+"'")
        print("Please make sure that the ASCII file is properly formatted.")
