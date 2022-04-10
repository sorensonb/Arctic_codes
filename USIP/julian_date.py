#!/usr/bin/env python

"""
NAME:
  julian_date.py

PURPOSE:
  Python package that calculates the Julian day for a given date.

USAGE:
  Importing
    from julian_date import julian_date

  Calculating the Julian day for a given date
    jd = julian_date(yy,mm,dd,hh,mm,ss)

PARAMETERS:
  jd - the calculated Julian day for the given date
  yy - the two-digit year
  mm - the two-digit month
  dd - the two-digit day
  hh - the two-digit hour
  mm - the two-digit minute
  ss - the two-digit second

SUBROUTINES CALLED:
  None

MODIFICATIONS:
  Nicholas Gapp - 160712: Added header, revised comments, and reformatted code.


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

def julian_date(year,month,day,hour,min,sec):

    # Convert the two digit year to a four digit year.
    if (year > 50):
        year = year + 1900
    else:
        year = year + 2000

    # Calculate the terms in the julian date equation.  Note: 'term2' is an
    # intermediary term in the julian day equation.  All terms must be integers
    # except for 'term6' and 'term7'.
    term1 = 367*year
    term2 = (month+9)/12
    term3 = 7*(year+term2)/4
    term4 = 275*month/9
    term5 = day
    term6 = 1721013.5
    term7 = (hour + float(min)/60.0 + float(sec)/3600.0)/24.0

    # Calculate the julian day.
    julianday = float(term1) - float(term3) + float(term4) + float(term5) + \
                float(term6) + float(term7)

    # Return the julian date.
    return julianday
