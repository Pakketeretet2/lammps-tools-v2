#!/usr/bin/env python3
#

from lammpstools import block_writers
from lammpstools import dump_reader
import sys

def print_usage():
    print("Use this like dump-cat -c <column headers> <dump files> "
          " -o <output>,"
          "  where <headers> is a quoted string containing the"
          "  names of the headers in the dump file and <output> is"
          "  an optional output dump file name. If no dump files"
          "  are given, stdcin is read.")


if len(sys.argv) < 2:
    print_usage()
    sys.exit(0)

i = 1
dumps = []
out_fname = None
while i < len(sys.argv):
    arg = sys.argv[i]
    if arg == "-c":
        col_headers = sys.argv[i+1].split(' ')
        i += 2
    elif arg == "-o":
        out_fname = sys.argv[i+1]
        i += 2
    elif arg == "-h":
        print_usage()
        sys.exit(0)
    else:
        dumps.append( arg )
        i += 1

if out_fname is None:
    out_fname = "dump_cat_test.dump"

print("Writing to",out_fname)
print("Writing",len(dumps),"dump files")

for dname in dumps:
    if dname.endswith( ".bin" ):
        fformat = "BIN"
    else:
        fformat = "PLAIN"

    if out_fname.endswith( ".bin" ):
        wmode = "ab"
    else:
        wmode = "a"

    d = dump_reader.dump_reader( dname, file_format = fformat,
            dump_format = "LAMMPS" );
    d.set_column_headers( col_headers )
    d.no_block_data_copy = True
    for b in d:
        block_writers.to_lammps_dump( out_fname, wmode, b )

