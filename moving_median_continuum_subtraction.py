#!/usr/bin/env python

# Perform a moving median continuum subtraction on MWA spectral line data

from astropy.io import fits
import argparse
import sys
import numpy as np
import psutil
# Parallelise the code
import multiprocessing

__author__ = ["Natasha Hurley-Walker"]
__date__ = "2020-04-24"

def _fslice_median(args):
    """
    A shallow wrapper for slice_median

    Parameters
    ----------
    args : list
        A list of arguments for slice_median

    Returns
    -------
    None
    """
    # an easier to debug traceback when multiprocessing
    # thanks to https://stackoverflow.com/a/16618842/1710603
    try:
        return slice_median(*args)
    except:
        import traceback
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

def slice_median(index, chunk):
    median = np.nanmedian(chunk, axis=1)
    return index, median

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Options")
    group1.add_argument("--fitsfile", dest='fitsfile', type=str, default=None,
                        help="Name of the input cube to be continuum-subtracted.")
    group1.add_argument("--binwidth", dest='binwidth', type=float, default=20,
                        help="Number of channels in each bin (default = 20).")
    group1.add_argument('--cores', dest='cores', default=None, type=int,
                        help="Number of cores to use (default = autodetect")
    options = parser.parse_args()

    if len(sys.argv) <= 1 or options.fitsfile is None:
        parser.print_help()
        sys.exit()

    infile = options.fitsfile
    outfile = infile.replace(".fits", "_subtracted.fits")
    contfile = infile.replace(".fits", "_moving_median.fits")
    binwidth = options.binwidth

    if options.cores is None:
        cores = multiprocessing.cpu_count()
    else:
        cores = results.cores

    hdu = fits.open(infile)
    data = hdu[0].data
    nchans = data.shape[1]
    ny = data.shape[2]
    nx = data.shape[3]
# Replace Slice 55 with a copy of Slice 56 to avoid having to use nanmedian
#    data[:,54,:,:] = data[:,55,:,:]

    med = np.zeros(hdu[0].data.shape)

    # We need to order our arguments by an index, since the results will be returned
    # in whatever order the pooled tasks finish
    args = []
# Hopefully this makes views on the data rather than copies!
    for x in range(0, nx):
        for y in range(0, ny):
            n = 0
            for i in range(binwidth/2, nchans - binwidth/2):
                args.append((n, y, x, data[:,i-binwidth/2:i+binwidth/2,y,x]))
                n+=1

#    # start a new process for each median
    pool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    try:
        # chunksize=1 ensures that we only send a single task to each process
        results = pool.map_async(_fslice_median, args, chunksize=1).get(timeout=10000000)
    except KeyboardInterrupt:
        pool.close()
        sys.exit(1)
    pool.close()
    pool.join()

#    for i in range(binwidth/2, nchans - binwidth/2):
#        med[:,i,:,:] = np.median(data[:,i-binwidth/2:i+binwidth/2,:,:], axis=1)

    cindices, yindices, xindices, medians = map(list, zip(*results))
    # Order correctly
    cind = np.argsort(cindices)
    yind = np.argsort(yindices)
    xind = np.argsort(xindices)
    medians = np.array(medians)
    medians = medians[cind,yind,xind]
    print medians.shape
    medians = medians.reshape((medians.shape[1], medians.shape[0], medians.shape[2], medians.shape[3]))
    med[:,binwidth/2:nchans-binwidth/2,:,:] = medians
#    # Flatten list of lists
#    o = [item for sublist in offsets for item in sublist]
#    # Make into array and apply
#    x += np.array(o)

    # Handle the edges with a single median across (half the binwidth) channels
    # TODO: Use the right numpy function (tile?) to make this easier to read
    # TODO: Or come up with a nicer way of handling this, based on feedback from the group
    med_start = np.median(data[:,0:binwidth/2,:,:], axis=1)
    med_end = np.median(data[:,nchans-binwidth/2:nchans,:,:], axis=1)
    for i in range(0, binwidth/2):
        med[:,i,:,:] = med_start
    for i in range(nchans-binwidth/2, nchans):
        med[:,i,:,:] = med_end

    # Re-flag the central channel
#    data[:,54,:,:] = np.nan*np.ones(data[:,54,:,:].shape)
    hdu[0].data = data - med
    hdu.writeto(outfile, overwrite=True)
    hdu[0].data = med 
    hdu.writeto(contfile, overwrite=True)
