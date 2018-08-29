#!/usr/bin/env python

import sys
import argparse
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

def measure_pixel_flux(infits, ra, dec):
    if infits is None:
        print "Must supply a filename"
        sys.exit(1)
    else:
        hdu = fits.open(infits)

    w = wcs.WCS(hdu[0].header)
    try:
        location = SkyCoord(ra+" "+dec)#, unit=(u.deg, u.deg))
    except ValueError:
        location = SkyCoord(ra+" "+dec, unit=(u.deg, u.deg))

    cp = np.squeeze(w.wcs_world2pix([[location.ra.value,location.dec.value]],1))
    xp, yp = cp[0], cp[1]

    print hdu[0].data[yp, xp]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Input/output files")
    group1.add_argument("--infits", dest='infits', type=str, default=None,
                        help="The fits image to be cut out from.")
#    group1.add_argument("--psffits", dest='psffits', type=str, default=None,
#                        help="The PSF image (optional)")
    group2 = parser.add_argument_group("Position options")
    group2.add_argument("--ra", dest='ra', type=str, default=None,
                        help="The central RA (use SkyCoord-compatible inputs)")
    group2.add_argument("--dec", dest='dec', type=str, default=None,
                        help="The central Dec (use SkyCoord-compatible inputs)")
    results = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    measure_pixel_flux(infits = results.infits, ra = results.ra, dec = results.dec)

