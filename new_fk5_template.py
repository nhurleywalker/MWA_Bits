#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "15/08/2019"

import os, sys
from optparse import OptionParser #NB zeus does not have argparse!
import numpy as np
from astropy.io import fits
from astropy import wcs

def new_fk5_template(ra, dec, nx, ny, pixscale, output, overwrite=False):
    if not os.path.exists(output) or overwrite is True:
        w = wcs.WCS(naxis=2)
        # Divisible by 2
        nx = (nx//2)*2
        ny = (ny//2)*2
        w.wcs.crpix = [nx/2, ny/2]
        w.wcs.cdelt = np.array([-pixscale, pixscale])
        w.wcs.crval = [ra, dec]
        w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
        header = w.to_header()
        data = np.zeros([nx, ny], dtype="float32")
        new = fits.PrimaryHDU(data,header=header) #create new hdu
        newlist = fits.HDUList([new]) #create new hdulist
        newlist.writeto(output, overwrite=True)
    return output

if __name__ == '__main__':
    parser = OptionParser(usage = "usage: %prog [options]" +
    """
    Make a new FK5 SIN projected FITS file
    """)
    parser.add_option("--ra", default=180.0, dest="racent", type="float", help="RA centre in decimal deg")
    parser.add_option("--dec", default=0.0, dest="decent", type="float", help="Dec centre in decimal deg")
    parser.add_option("--nx", default=100, dest="nx", type="int", help="Length of RA axis in pixels")
    parser.add_option("--ny", default=100, dest="ny", type="int", help="Length of Dec axis in pixels")
    parser.add_option("--pixscale", default=0.02, dest="pixscale", type="float", help="Pixel scale in degrees")
    parser.add_option("--output", default="template.fits", dest="output", help="Output filename")
    parser.add_option("--overwrite", default=False, action="store_true", dest="overwrite", help="Overwrite existing file (default False)")

    options, args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    else:
        new_fk5_template(options.racent, options.decent, options.nx, options.ny, options.pixscale, options.output, options.overwrite)
