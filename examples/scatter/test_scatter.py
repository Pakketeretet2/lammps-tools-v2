from lammpstools import block_data
from lammpstools import dump_reader
from lammpstools import data_field
from lammpstools import analysis

import matplotlib.pyplot as plt

import numpy as np
import sys

dname = "dump.melt.bin"
d = dump_reader.dump_reader( dname, file_format = "BIN", dump_format = "LAMMPS" )
d.set_column_headers( ['id', 'type', 'x', 'y', 'z', 'c_pe'] )

q0 = 1e-4
q1 = 4.0
dq = 0.01
qs  = np.arange( q0, q1, dq )

b = d.next_block()
Iavg = np.zeros( len(qs), dtype = float )
N = 0.0
while not d.at_eof():
    #scatter = analysis.scatter.rayleigh_gans( b, qs )
    #print("At t =", b.meta.t, file=sys.stderr)
    #Iavg += scatter
    #N += 1.0
    b = d.next_block()
    break

#Iavg /= N
#for q, I in zip(qs,scatter):
#    print(q,Iavg)
