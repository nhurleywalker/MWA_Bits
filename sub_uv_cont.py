#!/usr/bin/env

import sys
import numpy as np
from casacore.tables import table

obsid = sys.argv[1]

mset = table("{0}.ms".format(obsid), readonly=True)

data = mset.getcol("CORRECTED_DATA")
flag = mset.getcol("FLAG")
#real = ms.getdata("real")["real"]
#imag = ms.getdata("imaginary")["imaginary"]
#flag = ms.getdata("flag")["flag"]
#weight = ms.getdata("weight")["weight"]
#mreal = np.ma.masked_array(real, mask=flag)
mdata = np.ma.masked_array(data, mask=flag)

# Assume channel widths of 24
chanwidth = 24

for start in np.arange(0, data.shape[1], chanwidth):
    end = start + chanwidth
    avg = np.ma.mean(mdata[:,start:end,:], axis=1)
    avg_rep = np.repeat(avg[:, np.newaxis, :], chanwidth, axis=1)
    data[:,start:end,:] = data[:,start:end,:] - avg_rep

# Have to make a copy first
mset_sub = mset.copy("{0}_sub.ms".format(obsid), deep = True)
mset_sub = table("{0}_sub.ms".format(obsid), readonly = False)
mset_sub.putcol("CORRECTED_DATA", data)

mset.close()
mset_sub.close()
