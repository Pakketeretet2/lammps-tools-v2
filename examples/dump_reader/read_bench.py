import sys
from lammpstools import dump_reader, block_data

if len(sys.argv) > 1:
    dname = sys.argv[1]
else:
    print("Pass a dump file!", file=sys.stderr)
    sys.exit(-1)

if len(sys.argv) > 2:
    if int(sys.argv[2]) == 0:
        no_copy = True
    else:
        no_copy = False

if dname.endswith(".bin"):
    fformat = "BIN"
else:
    fformat = "PLAIN"

d = dump_reader.dump_reader( dname, file_format = fformat, dump_format = "LAMMPS" )
d.set_column_headers( ['id', 'type', 'x', 'y', 'z'] )
d.no_block_data_copy = no_copy
c = 0
for b in d:
    if c > 0 and (c%50 == 0):
        print( "t =", b.meta.t, "N =",b.meta.N)
        if d.no_block_data_copy:
            print( "Now indexing is normal too, see:", b.x[5][1])
        else:
            print( "Now indexing is normal, see:", b.x[5][1])
    c += 1
