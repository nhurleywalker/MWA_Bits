import numpy as np
import pickle
import argparse
import os

from astropy.io import fits
from astropy import wcs
import astropy.units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord
from timing import *

wscleanpols = ["XX", "XY", "XYi", "YY"]

def main():

    parser = argparse.ArgumentParser(description="Calculate the primary beam correction and save the correction to the same pickle file (key = 'PB_CORR'). Currently only supports MWA observations. This will overwrite the 'PB_CORR' field in the pickle file.")
    parser.add_argument('ds_file', help="The name of the pickle file to read and write to.")
    parser.add_argument('--metafits', help="The (MWA-style) metafits file associated with this observation. This is used to get the observations''GRIDNUM'. Required for MWA calculations.")
    parser.add_argument('--output_file', help="Write to this other file instead of overwriting the input file.")
    parser.add_argument('--overwrite', action='store_true', help="Write to this other file instead of overwriting the input file")

    args = parser.parse_args()

    output_file = args.output_file or args.ds_file

    # Only proceed if the output file doesn't already have a PB_CORR column in it, OR if overwrite flag is set
    if os.path.exists(output_file):
        try:
            dat = np.load(output_file, allow_pickle=True)
        except:
            response = input(f"Unable to open {output_file} as a pickled dynamic spectrum. Continuing will overwrite this file completely. Do you want to continue? [y/N]: ").strip().lower()
            if response != 'y':
                print("Aborting.")
                exit()
        if 'PB_CORR' in dat.keys():
            if args.overwrite == False:
                print(f"{output_file} already contains a PB_CORR column. Skipping.")
                exit()

    # Read the input file. If the output file *is* the input file, then we've already read it in, above
    if output_file != args.ds_file:
        dat = np.load(args.ds_file, allow_pickle=True)

    coord = SkyCoord("16:27:59.5", "-52:35:04.3", frame='fk5', unit=(u.hour, u.deg))

    if dat['TELESCOPE'] == 'MWA':

        # TODO: read this from the pkl file / check it's correct
        instpols = dat["POLS"]
        if not args.metafits:
            raise ValueError("--metafits option required for MWA PB calculation")

        with fits.open(args.metafits) as hdul:
            gridnum = hdul[0].header["GRIDNUM"]
        freqs = dat['FREQS']
        t = Time(dat['TIMES'][0]/86400, scale='utc', format='mjd')
        finaldynspec = np.empty((len(times),len(freqs),4), dtype='float32')
        for i in range(0, len(times)):
            t = times[i]
            for j in range(0, len(freqs)):
                f = freqs[j]
                for k in range(0, len(instpols)):
                    pol = instpols[k]
                    wpol = wscleanpols[k]
           # Need to create a single pixel FITS file for each instrumental stokes, time, frequency
                    w = wcs.WCS(naxis=2)
                    nx = 1
                    ny = 1
                    nf = 1
                    w.wcs.crpix = [1, 1, 1]
                    w.wcs.cdelt = np.array([-1, 1, 1])
                    w.wcs.crval = [coords.fk5.ra.deg, coords.fk5.dec.deg, freq]
                    w.wcs.ctype = ["RA---SIN", "DEC--SIN", "Hz"]
                    header = w.to_header()
                    new = fits.PrimaryHDU(dat["DS"][i,j,k],header=header) #create new hdu
                    newlist = fits.HDUList([new]) #create new hdulist
                    newlist[0].header["DATE-OBS"] = t.isot
                    # TODO make outputs in line with what wsclean wants (lowercase, xyi instead of yk)
                    output = args.ds_file.replace(".pkl", f"-{i:04d}-{j:04d}-{wpol}-image.fits")
                # TODO fix outputs -- and in the external script I MUST write these to RAM not to /scratch!
                    newlist.writeto(output, overwrite=True)
                os.system(f"beam -2016 -proto {output} -name beam -m {metafits} -ms {ms}")
                os.system(f"pbcorrect {output minus some suffix} image.fits beam stokes")
                if os.path.exists(f"{some smart stokes test}"):
                    for s in stokes:
                        with fits.open(f"{smart naming}{stokes}") as hdu:
                            finaldynspec[i,j,stokes] = hdu[0].data[0,0]
                else:
                   sys.exit(0)


        rX, rY = np.array([beam_lookup_1d(coord.fk5.ra.deg, coord.fk5.dec.deg, gridnum, t, freq) for freq in freqs]).T

        dat['PB_CORR'] = np.array([rX, np.full(rX.shape, np.nan), np.full(rX.shape, np.nan), rY]).T

    else:
        raise NotImplementedError(f"Primary beam correction calculation not implemented for {dat['TELESCOPE']}")

    with open(output_file, 'wb') as pkl:
        pickle.dump(dat, pkl)
