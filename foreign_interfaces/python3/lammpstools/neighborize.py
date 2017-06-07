from neighborize_ import *


def neigh_list_dist( block_data_ref, itype, jtype, method, dims, rc,
                     mol_policy = 0, bond_policy = 0, quiet = True ):
    return nearest_neighbours_dist( block_data_ref, itype, jtype,
                                    method, dims, rc, mol_policy,
                                    bond_policy, quiet )

