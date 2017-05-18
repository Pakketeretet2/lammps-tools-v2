import dump_reader
import block_data
import sys

d = dump_reader.dump_reader( "A_32.0_R_40_ribbons.dump", 0, 0 )
d.set_column_headers( ["id", "type", "x", "y", "z", "edt" ] )
c = 0
for b in d:
    print(b.meta.t)
    c += 1
    if (c%25 == 0):
        print("At block ", c)
