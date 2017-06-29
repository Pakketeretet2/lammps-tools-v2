from lammpstools import block_data, dump_reader, data_field
from lammpstools import block_writers, bond_order, neighborize
import numpy as np, sys

dname = 'dump.melt.bin'
d = dump_reader.dump_reader( dname, file_format = "BIN", dump_format = "LAMMPS" )
d.set_column_headers( [ 'id', 'type', 'x', 'y', 'z' ] )
d.no_block_data_copy = True

for b in d:
    print("t = ", b.meta.t )

    psi_n_real = np.zeros( b.meta.N, dtype = np.float64 )
    psi_n_imag = np.zeros( b.meta.N, dtype = np.float64 )

    neighs = neighborize.neigh_list_dist( b, 0, 0, 0, 2, 1.35 )


    psi = bond_order.compute_psi_n( b.get_ref_(), neighs, 6, [0,1,0], psi_n_real, psi_n_imag )
    psi_real = data_field.new_data_field( "psi_r", data_field.DATA_TYPE_DOUBLE, b.meta.N )
    psi_imag = data_field.new_data_field( "psi_i", data_field.DATA_TYPE_DOUBLE, b.meta.N )
    for i in range(0,b.meta.N):
        psi_real[i] = psi_n_real[i]
        psi_imag[i] = psi_n_imag[i]
    b.add_data_field( psi_real )
    b.add_data_field( psi_imag )
    block_writers.to_lammps_dump( "psi.dump.bin", "ab", b )

