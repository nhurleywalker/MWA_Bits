#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import sys
#import os
#from datetime import datetime
from casacore.tables import table #, taql, maketabdesc, makescacoldesc
#from astropy.time import Time
#import astropy.units as u

#Root directory
#root = sys.argv[3]


#Opening ms files
#print(sys.argv)
#sys.stdout.flush()
if sys.argv[1] is not None:
    ms = sys.argv[1]
    mset = table(ms, readonly=True, ack=False)
    time = list(set(mset.getcol("TIME")))
    print(len(time))
    mset.close()
