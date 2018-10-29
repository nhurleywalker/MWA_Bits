#!/usr/bin/env python

import numpy as np
from table_formatting import combine_cols
from astropy.table import Column
from astropy.table import Table

data = 10000000*np.random.rand(20)
errors = np.random.rand(20)*10000#*np.random.rand(20)

print "large numbers,Uniform = True"
col1 = Column(combine_cols(data,errors,uniform=True), name = "large, uniform")
print "large numbers,Uniform = False"
col2 = Column(combine_cols(data,errors,uniform=False), name = "large, mixed")

data = 0.1*np.random.rand(20)
errors = np.random.rand(20)*0.01*np.random.rand(20)

print "small numbers,Uniform = True"
col3 = Column(combine_cols(data,errors,uniform=True), name = "small, uniform")
print "small numbers,Uniform = False"
col4 = Column(combine_cols(data,errors,uniform=False), name = "small, mixed")

t = Table([col1, col2, col3, col4])

t.write("example_table.tex", overwrite=True)

with file('example_table.tex', 'r') as original: data = original.read()
with file('latex_me.tex', 'w') as modified: modified.write("\\documentclass{article} \n\\begin{document}\n" + data + "\n \\end{document}")

print "Now run 'pdflatex latex_me.tex'"
