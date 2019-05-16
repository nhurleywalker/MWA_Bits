#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "17/05/2019"

from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy import wcs
import astropy.units as u

from regions import EllipseSkyRegion
from regions import PixCoord
import numpy as np
import sys
import argparse

# Calculate and return the sum or mean or RMS of pixels inside an ellipse in a FITS image

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Files and operations")
    group1.add_argument('--fitsfile', dest='fitsfile', default=None, \
                        help="FITS image to read")
    group1.add_argument('--operation', dest='operation', default=None, \
                        help="Choose sum, mean, or rms")
    group2 = parser.add_argument_group("Ellipse parameters")
    group2.add_argument('--ra', dest='ra', default=None, type=float, \
                        help="RA (decimal deg)")
    group2.add_argument('--dec', dest='dec', default=None, type=float, \
                        help="Dec (decimal deg)")
    group2.add_argument('--major', dest='major', default=None, type=float, \
                        help="Major axis length (deg)")
    group2.add_argument('--minor', dest='minor', default=None, type=float, \
                        help="Minor axis length (deg)")
    group2.add_argument('--pa', dest='pa', default=None, type=float, \
                        help="Position angle (deg)")

    options = parser.parse_args()

    if options.operation != "sum" and options.operation != "mean" and options.operation != "rms":
        print "Need to pick sum, mean or rms"
        sys.exit(1)

    hdu = fits.open(options.fitsfile)
    w = wcs.WCS(hdu[0].header)
    # Get co-ordinate system from image
    data = np.squeeze(hdu[0].data)

    center_sky = SkyCoord(options.ra, options.dec, unit='deg', frame='fk5')
    ellipse_sky = EllipseSkyRegion(center=center_sky,
                                height=2*options.major * u.deg, width=2*options.minor * u.deg,
                                angle=options.pa * u.deg)
    ellipse_pix = ellipse_sky.to_pixel(w)

    xy = np.indices(data.shape, dtype=np.float32)
    x = np.array(xy[1, :])
    y = np.array(xy[0, :])

    ind = ellipse_pix.contains(PixCoord(x, y))

    if options.operation == "sum":
        print np.nansum(data[ind])
    elif options.operation == "mean":
        print np.nanmean(data[ind])
    elif options.operation == "rms":
        print np.nanstd(data[ind])
