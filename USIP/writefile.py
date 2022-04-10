#!/usr/bin/env python
"""
  NAME:
    writefile.py
 
  PURPOSE:
    Creates a file whose filename is set either by the optional 'name' argument
    or by the ADPAA class attribute 'name' and writes information contained
    within the ADPAA object the method was called upon to that file.

  CALLS:
    Built-in Python modules numpy.ma and re

  MODIFICATIONS:
    Nicholas Gapp - 160923: Written
    Nicholas Gapp - 161028: Added comments
    Nicholas Gapp - 161123: Fixed to print new-line characters where needed and
                            fixed so number of header lines is calculated before
                            data is output to file.

  USAGE:
    Create a default APDAA object by
      A = ADPAA()

    Write data contained within an ADPAA object to a properly formatted file by
      Specifying filename:
        A.WriteFile('output_filename')
      Using filename set by class attribute 'name':
        A.WriteFile()

  EXAMPLE:
    Instantiate an ADPAA object
      >>> air_ascent_1Hz = ADPAA()

    Write data to new file 15_08_08_14_57_48.air.ascent.new.1Hz after making changes
      >>> air_ascent_1Hz.WriteFile('15_08_08_14_57_48.air.ascent.new.1Hz')

      --OR--
      >>> air_ascent_1Hz.name = '15_08_08_14_57_48.air.ascent.new.1Hz'
      >>> air_ascent_1Hz.WriteFile()


  COPYRIGHT:
    2016 David Delene

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

import re
import numpy.ma as ma

def WriteFile(self, filename):

    # Underscored variables are private to the class.
#    try:
        # Update output format array to any changes that may have been made.
        self.format = self.create_format_arr(self.VMISS)

        # Update total number of lines in the header.
        self.NLHEAD = 14 + self.NV + self.NSCOML + self.NNCOML

        # Open output file for writing.
        _fhand = open(filename, 'w')

        # Output header information to output file according to the
        # information within the ADPAA object.
        _fhand.write('%02d %04d\n' % (int(self.NLHEAD), int(self.FFI)))
        _fhand.write(self.ONAME + '\n')
        _fhand.write(self.ORG + '\n')
        _fhand.write(self.SNAME + '\n')
        _fhand.write(self.MNAME + '\n')
        _fhand.write('%01d %01d\n' % (self.IVOL, self.NVOL))
        _fhand.write(self.DATE + ' ' + self.RDATE + '\n')
        _fhand.write(' %11.4f\n' % (self.DX))
        _fhand.write(self.XNAME + '\n')
        _fhand.write(str(self.NV) + '\n')
        for a in range(len(self.VSCAL)):
            _fhand.write(self.format[a+1].format(float(self.VSCAL[a])))
        _fhand.write('\n')
        for a in range(len(self.VMISS)):
            _fhand.write(self.format[a+1].format(float(self.VMISS[a])))
        _fhand.write('\n')
        for field in range(len(self.VNAME)):
            _fhand.write(self.VNAME[field])
            _fhand.write('\n')
        _fhand.write(str(self.NSCOML) + '\n')

        # Add in special comment lines if needed.
        if self.NSCOML > 0:
            for i in range(self.NSCOML):
                _fhand.write(self.SCOMM[i])

        # Continue output of header information to output file.
        _fhand.write(str(self.NNCOML) + '\n')
        _fhand.write(str(self.DTYPE) + '\n')
        _fhand.write(str(self.VFREQ) + '\n')

        # Output variable short names and their respective units.
        # If short names or units are more than 11 characters long,
        # truncate to 11 characters total.
        for b in self.VDESC:
            if len(b) <= 11:
                _fhand.write(' ' + b.ljust(11))
            else:
                _fhand.write(' ' + b[:11])
        _fhand.write('\n')
        for b in self.VUNITS:
            if len(b) <= 11:
                _fhand.write(' ' + b.ljust(11))
            else:
                _fhand.write(' ' + b[:11])
        _fhand.write('\n')

        # Remove missing data masks and substitute in missing value codes.
        for key in self.data:
            if (ma.count_masked(self.data[key]) > 0):
                self.data[key] = self.data[key].filled(fill_value=float(self.VMISS[self.VDESC.index(key)-1]))

        # Output data values to output file.
        for i in xrange(0,len(self.data[self.VDESC[0].strip()])):
            j = 0
            for key in self.VDESC:
                _fhand.write(self.format[j].format(self.data[key.strip()][i]))
                j += 1
            _fhand.write('\n')

        # Close the output file.
        _fhand.close()

#    except:
#        print "Ensure class attribute 'name' is set or ASCII file can be "+\
#              "created."
