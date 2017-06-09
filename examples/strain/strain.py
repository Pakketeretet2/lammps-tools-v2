#from lammpstools import skeletonize, block_data, dump_reader
from lammpstools import *

dname = '../dump_reader/polymer.dump'
d = dump_reader.dump_reader(dname, file_format = "PLAIN", dump_format = "LAMMPS" )
d.set_column_headers( [ 'id', 'mol', 'type', 'x', 'y', 'z' ] )
b = d.next_block()
print("Got N = ", b.meta.N, " particles.")

r0 = 1.0
itype = jtype = 0
method = 1
dims = 3
rc = 1.25
strain = analysis.neighborize.neighbour_strain(b, r0, itype, jtype,
                                               method, dims, rc)

neighs = analysis.neighborize.neigh_list_dist( b, 0, 0, 1, 3, 1.25 )
ids = b.data_by_name("id")

for i in range(0,b.meta.N):
    print("particle",ids[i],"(index",i,") has",len(neighs[i]),"neighbours.")
