from lammpstools import transformations
from lammpstools import dump_reader
from lammpstools import block_writers

import sys

dname = "../common_dump_files/assemble.dump"
d = dump_reader.dump_reader( dname, file_format = "PLAIN", dump_format = "LAMMPS" )
d.no_block_data_copy = True
d.set_column_headers( [ "id", "type", "x", "y", "z" ] )

b = d.next_block()

angle  = 3.1415927 / 4.0;
axis   = [ 0, 0, 1 ]
origin = [ 0, 0, 0 ]
block_writers.to_lammps_dump( "test.dump.bin", "wb", b );

for i in range(0,8):
    transformations.rotate_all( b, [0,0,1], [0,0,0], angle )
    block_writers.to_lammps_dump( "test.dump.bin", "ab", b );
