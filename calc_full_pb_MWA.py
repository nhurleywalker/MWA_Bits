import numpy as np
import pickle
import argparse
import os

from astropy.io import fits
from astropy import wcs
import astropy.units as u
from astropy.time import Time
from astropy.constants import c
from astropy.coordinates import SkyCoord

wscleanpols = ["XX", "XY", "XYi", "YY"]
stokes = ["I", "Q", "U", "V"]
beampols = ["xxi" , "xxr", "xyi", "xyr", "yxi", "yxr", "yyi", "yyr"]

parser = argparse.ArgumentParser(description="Calculate the primary beam correction and save the correction to the same pickle file (key = 'PB_CORR'). Currently only supports MWA observations. This will overwrite the 'PB_CORR' field in the pickle file.")
parser.add_argument('ds_file', help="The name of the pickle file to read and write to.")
parser.add_argument('--metafits', help="The (MWA-style) metafits file associated with this observation. This is used to get the observations''GRIDNUM'. Required for MWA calculations.")
parser.add_argument('--ms', help="The measurement set associated with this observation. Unused if metafits file is specified.")
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
    times = Time(dat['TIMES']/86400, scale='utc', format='mjd')
    finaldynspec = np.empty((len(times),len(freqs),4), dtype='float32')
#    print(len(freqs), len(times))
#    print(dat["DS"][40,60,:])
    for i in range(0, len(times)):
        t = times[i]
        for j in range(0, len(freqs)):
            f = freqs[j]
            for k in range(0, len(instpols)):
                pol = instpols[k]
                wpol = wscleanpols[k]
       # Need to create a single pixel FITS file for each instrumental stokes, time, frequency
                w = wcs.WCS(naxis=3)
                nx = 1
                ny = 1
                nf = 1
                w.wcs.crpix = [1, 1, 1]
                w.wcs.cdelt = np.array([-1, 1, 1])
                w.wcs.crval = [coord.fk5.ra.deg, coord.fk5.dec.deg, f]
                w.wcs.ctype = ["RA---SIN", "DEC--SIN", "Hz"]
                header = w.to_header()
# In the absence of any better idea, I will take the real part for now
                print(dat["DS"][i,j,k])
                new = fits.PrimaryHDU(np.real([[dat["DS"][i,j,k]]]),header=header) #create new hdu
                newlist = fits.HDUList([new]) #create new hdulist
                newlist[0].header["DATE-OBS"] = t.isot
                # TODO make outputs in line with what wsclean wants (lowercase, xyi instead of yk)
                outputstem = args.ds_file.replace(".pkl", f"-{i:04d}-{j:04d}")
                outputstemwpol = f"{outputstem}-{wpol}"
                output = f"{outputstemwpol}-image.fits"
                newlist.writeto(output, overwrite=True)
# TODO if I paralellise -- make the beam names specific rather than just 'beam'
            os.system(f"beam -2016 -proto {output} -name beam-{i:04d}-{j:04d} -m {args.metafits} -ms {args.ms}")
            os.system(f"pbcorrect {outputstem} image.fits beam-{i:04d}-{j:04d} {outputstem}")
            if os.path.exists(f"{outputstem}-I.fits"):
                for l in range(0, len(stokes)):
                    s = stokes[l]
                    with fits.open(f"{outputstem}-{s}.fits") as hdu:
                        finaldynspec[i,j,l] = hdu[0].data[0,0]
# Remove fake images of celestial stokes
                    os.remove(f"{outputstem}-{s}.fits")
# Remove beam files
            for b in beampols:
                os.remove(f"beam-{i:04d}-{j:04d}-{b}.fits")
# Remove fake images of instrumental stokes
            for c in wscleanpols:
                os.remove(f"{outputstem}-{c}-image.fits")

#            else:
#               sys.exit(0)


#    rX, rY = np.array([beam_lookup_1d(coord.fk5.ra.deg, coord.fk5.dec.deg, gridnum, t, freq) for freq in freqs]).T

    dat['PB_CORR'] = finaldynspec
    dat['PB_POLS'] = stokes
    with open(output_file, 'wb') as pkl:
        pickle.dump(dat, pkl)

else:
    raise NotImplementedError(f"Primary beam correction calculation not implemented for {dat['TELESCOPE']}")

