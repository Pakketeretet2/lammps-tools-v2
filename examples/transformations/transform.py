from lammpstools import transformations
from lammpstools import dump_reader
from lammpstools import block_writers

dname = "../common_dump_files/assemble.dump"
d = dump_reader.dump_reader( dname, file_format = "PLAIN", dump_format = "LAMMPS" )
d.no_block_data_copy = True
d.set_column_headers( [ "id", "type", "x", "y", "z" ] )

b = d.next_block()
block_writers.to_lammps_data( "test.data", "w", b );
block_writers.to_lammps_dump( "test.dump.bin", "wb", b );

