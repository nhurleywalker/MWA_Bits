#!/usr/bin/env python

# Input is a list of MeerKAT targets as a text file
# and a set of coordinates as some system arguments

import sys
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u

df = pd.read_csv(sys.argv[1], sep='\s+')
print(df)
ras = df["RA(J2000)"]
decs = df["DEC(J2000)"]
coords = SkyCoord(ras, decs, unit=(u.hour, u.deg), frame='fk5')
c0 = SkyCoord(sys.argv[2], sys.argv[3], unit=(u.hour, u.deg), frame='fk5')

for i in range(0, len(ras)):
    print("{0}, {1}, {2:3.0f} arcmin".format(ras[i], decs[i], (c0.separation(coords[i]).deg)*60))
