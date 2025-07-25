#!/usr/bin/env python

#Convert metafits ra, dec to what chgcentre wants

import sys
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits

metafits = fits.open(sys.argv[1])
ra = metafits[0].header["RAPHASE"]
dec = metafits[0].header["DECPHASE"]
coords = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='fk5')
print(coords.to_string('hmsdms'))
