#!/usr/bin/python

# Combine healpix images in an inverse-variance weighted fashion
# Currently assumes some sensible logical filenames are given by the user
# NHW 02/12/2015

from glob import glob
import healpy as hp
import astropy.wcs as wcs
from astropy.table import Table,Column
import os
import numpy as np
import sys
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from optparse import OptionParser

usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('-f','--filename',dest="filename",default=None,
                  help="Input maps to read <FILE>",metavar="FILE")
parser.add_option('-r','--rms',dest="rms",default=None,
                  help="Input rms maps to read <FILE>",metavar="FILE")
parser.add_option('-o','--output',dest="output",default="test.fits",
                  help="Output file <FILE>; default = test.fits",metavar="FILE")
(options, args) = parser.parse_args()

inputs=sorted(glob(options.filename))
noises=sorted(glob(options.rms))

if len(inputs)==0 or len(noises)==0 or (len(inputs)!=len(noises)):
    print "Expecting non-zero, equal number of images and corresponding RMS maps"
    sys.exit(1)

def all_same(items):
    return all(x == items[0] for x in items)

# Check they all have the same nside
nsides=[]
for filename in inputs+noises:
    hdu_in = pyfits.open(filename)
# Healpix is a second hdu table
    header=hdu_in[1].header
    nsides.append(header['NSIDE'])
    hdu_in.close()
if all_same(nsides):
    nside=nsides[0]
else:
    print "All healpix maps must have the same nside dimension for this code to work."
    sys.exit(1)

# Make a blank map
npix=hp.nside2npix(nside)
hpMapFullSky=np.zeros(npix)

# Make a blank weight map
hpMapFullWeight=np.zeros(npix)

print inputs
print noises
for filename,rms in zip(inputs,noises):
    print filename,rms
    hp_map=hp.read_map(filename)
    hp_noise=hp.read_map(rms)
    hp_wmap=hp_map/(hp_noise*hp_noise)
    hp_weight=1/(hp_noise*hp_noise)
# accumulate
    hpMapFullSky+=hp_wmap
    hpMapFullWeight+=hp_weight

# Combine the maps
hpMapFullSky=hpMapFullSky/hpMapFullWeight

#Write to file (assuming Celestial coords; use coord='G' for Galactic)
hp.write_map(options.output,hpMapFullSky,coord='C')

