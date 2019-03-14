#!/usr/bin/env python

import sys
from astropy.io import fits
from numpy import sin, radians
import argparse
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

def wcs_to_pix(ra, dec, fitsimage):
    hdu = fits.open(fitsimage)
    w = wcs.WCS(hdu[0].header)
    x, y = w.wcs_world2pix([[ra,dec]],0).transpose()
#    if hdu[0].header["CTYPE1"] == "GLON":
#        coords = SkyCoord(lon, lat, frame="galactic", unit = (u.deg, u.deg))
#    else:
#        coords = SkyCoord(lon, lat, frame="fk5", unit = (u.deg, u.deg))
    return x, y

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Input/output files")
    group1.add_argument("--ra", dest='ra', type=float, default=None,
                        help="RA / longitude value (deg)")
    group1.add_argument("--dec", dest='outfits', type=float, default=None,
                        help="Dec / latitude value (deg)")
    group1.add_argument("--fitsimage", dest='fitsimage', type=file, default=None,
                        help="Fits image to use")
    results = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    x, y = wcs_to_pix(ra, dec, fitsimage)
    print x, y

