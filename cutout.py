#!/usr/bin/env python

import sys
import argparse
import numpy as np

from astropy.io import fits
from astropy import wcs
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord


def cutout(infits, outfits, psffits, ra, dec, size):

    hdu = fits.open(infits)
    w = wcs.WCS(hdu[0].header)

    location = SkyCoord(ra+" "+dec)#, unit=(u.deg, u.deg))
    framesize = u.Quantity(size, u.deg)

    cutout = Cutout2D(hdu[0].data, location, framesize, wcs=w)

    header_new = cutout.wcs.to_header()

    # Read these from the correct PSF image and then put them in the cutout
    if psffits is not None:
        psf = fits.open(psffits)
        wpsf = wcs.WCS(psf[0].header)
        cp = np.squeeze(wpsf.wcs_world2pix([[location.ra.value,location.dec.value,1]],0))
        xp, yp = cp[0], cp[1]
        header_new["BMAJ"] = psf[0].data[0,yp,xp]
        header_new["BMIN"] = psf[0].data[1,yp,xp]
        header_new["BPA"] = psf[0].data[2,yp,xp]

    new = fits.PrimaryHDU(cutout.data,header=header_new) #create new hdu
    newlist = fits.HDUList([new]) #create new hdulist
    newlist.writeto(outfits, overwrite = True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Input/output files")
    group1.add_argument("--infits", dest='infits', type=str, default=None,
                        help="The fits image to be cut out from.")
    group1.add_argument("--outfits", dest='outfits', type=str, default=None,
                        help="The output fits image.")
    group1.add_argument("--psffits", dest='psffits', type=str, default=None,
                        help="The PSF image (optional)")
    group2 = parser.add_argument_group("Cropping options")
    group2.add_argument("--ra", dest='ra', type=str, default=None,
                        help="The central RA (use SkyCoord-compatible inputs)")
    group2.add_argument("--dec", dest='dec', type=str, default=None,
                        help="The central Dec (use SkyCoord-compatible inputs)")
    group2.add_argument("--size", dest='size', type=float, default=1.0,
                        help="The size of the edge of the square cut-out in degrees.")
    results = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    cutout(infits = results.infits, outfits = results.outfits, psffits = results.psffits,
                     ra = results.ra, dec = results.dec, size = results.size)
