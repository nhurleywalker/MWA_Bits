#!/usr/bin/env python

# Convert decimal co-ordinates to sexigecimal

from astropy import units as u
from astropy.coordinates import SkyCoord
import sys

c = SkyCoord(ra=float(sys.argv[1])*u.degree, dec=float(sys.argv[2])*u.degree)
print(c.to_string('hmsdms'))
