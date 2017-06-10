import neighborize_


def neigh_list_dist( block_data, itype, jtype, method, dims, rc,
                     mol_policy = 0, bond_policy = 0, quiet = True ):
    return neighborize_.nearest_neighbours_dist( block_data.get_ref(),
                                                 itype, jtype,
                                                 method, dims, rc, mol_policy,
                                                 bond_policy, quiet )

def neighbour_strain( block_data, r0, itype, jtype, method, dims, rc ):
    return neighborize_.neighbour_strain( block_data.get_ref(), r0,
                                          itype, jtype, method, dims, rc )
