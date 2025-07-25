#!/usr/bin/env python

import sys
# Take as input a number of channnels, a start and end channel (inclusive), and generate the flagsubbands command

nchan = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

chanstr = f"singularity exec $GXCONTAINER flagsubbands $mst {str(nchan)}"
for chan in range(start, end+1):
     chanstr += f" {str(chan)}"

print(chanstr)
