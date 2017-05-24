# lammpstools is now a package, you can directly import
# block_data and dump_reader from them.
from lammpstools import block_data
from lammpstools import dump_reader
from lammpstools import data_field   # for data_field_types

# Define a new dump reader for a LAMMPS dump file.
d = dump_reader.dump_reader( "polymer.dump",
                            file_format = "PLAIN", dump_format = "LAMMPS" )

dl = dump_reader.dump_reader( "bondinfo.dump",
                              file_format = "PLAIN", dump_format = "LAMMPS",
                              is_local = True )

# We know all cols in bondinfo.dump are ints, so set the default type:
dl.set_default_column_type( "int" )
# Specify the headers for atom dump file up front: dump_reader
# guesses the type of some cols based on this.
d.set_column_headers( [ 'id', 'mol', 'type', 'x', 'y', 'z' ] )

c = 0

for b, bl in zip(d, dl):
    if c > 0 and c%25 == 0:
        print("At block",c,", t =", b.meta.t, " contains N =",
              b.meta.N, " particles and", bl.meta.N, "bonds")
        i = 1233
        if b.meta.N > i:
            print("Position of particle",i,"is",b.x[i])
        if bl.meta.N > i:
            print("data of bond",i,"is [ ", end = "")
            for j in range(0, len(bl.data)):
                print( bl.data[j][i]," ", end = "")
            print(" ]")
    c += 1
