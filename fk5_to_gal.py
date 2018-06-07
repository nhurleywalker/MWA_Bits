#!/usr/bin/env python

#Convert fk5 ra/dec to galactic l/b

import sys
from astropy.coordinates import SkyCoord
from astropy import units as u

a = SkyCoord(sys.argv[1], sys.argv[2], unit=(u.deg, u.deg), frame='fk5')
print a.galactic.l.deg, a.galactic.b.deg
