#!/usr/bin/env

import sys
import numpy as np
from casacore.tables import table
import argparse

__author__ = "Natasha Hurley-Walker"
__date__ = "2020-05-01"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Input/output files")
    group1.add_argument("--msfile", dest='msfile', type=str, default=None,
                        help="Name of the input visibilities to be continuum-subtracted.")
    group1.add_argument("--datacolumn", dest='datacolumn', type=str, default="CORRECTED_DATA",
                        help="Column to subtract from (default = CORRECTED_DATA).")
    group1.add_argument("--binwidth", dest='binwidth', type=int, default=16,
                        help="Number of channels in each bin (default = 16).")
    options = parser.parse_args()

    if len(sys.argv) <= 1 or options.msfile is None:
        parser.print_help()
        sys.exit()

    binwidth = options.binwidth
    msfile = options.msfile

    mset = table("{0}".format(options.msfile), readonly=True)
    # Have to make a copy first
    mset_sub = mset.copy(msfile.replace(".ms", "_sub.ms"), deep = True)


    pols = [0, 1, 2, 3]
    for pol in pols:
        all_data = mset.getcol(options.datacolumn)
        data = all_data[:,:,pol]
# Free up memory
        del all_data
        all_flag = mset.getcol("FLAG")
        flag = all_flag[:,:,pol]
# Free up memory
        del all_flag

        nchans = data.shape[1] #2 time and baseline, 1 channel, 0 pol

        ## There are also values set to NaN that are not in the flag table
        new_flags = np.isnan(data)
        # If the data is NaN OR the data are flagged then we want to flag it
        total_flags = np.logical_or(flag, new_flags)

        mdata = np.ma.masked_array(data, mask=total_flags)
# Free up memory
        del data
        del flag

        avg = np.zeros(mdata.shape, dtype="complex")

        # Calculate the moving average
        for i in range(binwidth/2, nchans - binwidth/2):
            avg[:,i] = np.ma.mean(mdata[:,i-binwidth/2:i+binwidth/2], axis=1)

        # Handle the edges with a single mean across (half the binwidth) channels
        avg_start = np.ma.mean(mdata[:,0:binwidth/2], axis=1)
        avg_end = np.ma.mean(mdata[:,nchans-binwidth/2:nchans], axis=1)
        for i in range(0, binwidth/2):
            avg[:,i] = avg_start
        for i in range(nchans-binwidth/2, nchans):
            avg[:,i] = avg_end

        mset_sub = table(msfile.replace(".ms", "_sub.ms"), readonly = False)
        new_data = mset.getcol(options.datacolumn)
        new_data[:,:,pol] = mdata - avg
        mset_sub.putcol(options.datacolumn, new_data)
        mset_sub.close()
# Free up memory
        del new_data
        del mdata
        del avg
    mset.close()

    # Old-fashioned binning: too many residuals
    #    for start in np.arange(0, data.shape[1], binwidth):
    #        end = start + binwidth
    #        avg = np.ma.mean(mdata[:,start:end,:], axis=1)
    #        avg_rep = np.repeat(avg[:, np.newaxis, :], binwidth, axis=1)
    #        data[:,start:end,:] = data[:,start:end,:] - avg_rep
