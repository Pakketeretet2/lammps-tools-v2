from lammpstools import block_data
from lammpstools import dump_reader
from lammpstools import data_field
from lammpstools import analysis

import matplotlib.pyplot as plt

import numpy as np
import sys

class scatter_opts:
    def __init__(self):
        self.q0 = 0.01
        self.q1 = 10.0
        self.dq = 0.01

    

def scatter_dump(d, Natoms, r0, de, scale, opts):
    
    qs  = np.arange( opts.q0, opts.q1, opts.dq )
    Iavg = np.zeros( len(qs), dtype = float )
    N = 0.0
    radius = r0 * np.ones( Natoms )

    for b in d:
        scatter = analysis.scatter.rayleigh_gans( b, qs, radius = radius,
                                                  position_scale = scale,
                                                  d_epsilon = de,
                                                  ids = None )
        print("At t =", b.meta.t, file=sys.stderr)
        Iavg += scatter
        N += 1.0
        b = d.next_block()
    Iavg /= N
    return (qs,Iavg)


de    = 0.01
r0    = 0.5
scale = 1.0
if len(sys.argv) > 1:
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "-o":
            out_file = sys.argv[i+1]
            i += 2
        elif sys.argv[i] == "-de":
            d_epsilon = float(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == "-m":
            if sys.argv[i+1] == "single":
                mode = 'single'
                dname = "../common_dump_files/single_particle.dump"
                d = dump_reader.dump_reader( dname, file_format = "PLAIN",
                                             dump_format = "LAMMPS" )
                d.set_column_headers( ['id', 'type', 'x', 'y', 'z'] )
                d.no_block_data_copy = True
                Natoms = 1

            elif sys.argv[i+1] == "melt":
                mode = 'melt'
                dname = "../common_dump_files/dump.melt.bin"
                d = dump_reader.dump_reader( dname, file_format = "BIN",
                                             dump_format = "LAMMPS" )
                d.set_column_headers( ['id', 'type', 'x', 'y', 'z', 'c_pe'] )
                d.no_block_data_copy = True
                Natoms = 108000
            elif sys.argv[i+1] == "dilute":
                mode = 'melt'
                dname = "../common_dump_files/dump.dilute.bin"
                d = dump_reader.dump_reader( dname, file_format = "BIN",
                                             dump_format = "LAMMPS" )
                d.set_column_headers( ['id', 'type', 'x', 'y', 'z', 'c_pe'] )
                d.no_block_data_copy = True
                Natoms = 108000
                                
            else:
                raise RuntimeError("Unkown mode",sys.argv[i+1])
            i += 2
        elif sys.argv[i] == "-options":
            r0    = float(sys.argv[i+1])
            scale = float(sys.argv[i+2])
            de    = float(sys.argv[i+3])
            i += 4
        else:
            raise RuntimeError("Unkown arg",sys.argv[i])

opts = scatter_opts()
qs, I = scatter_dump(d, Natoms, r0, de, scale, opts)
with open( out_file, "w" ) as fp:
    for q, II in zip(qs,I):
        print(q, II, file = fp)
