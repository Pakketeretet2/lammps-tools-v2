import neighborize_
import block_data_

def neigh_list_dist( block_data, itype, jtype, method, dims, rc,
                     mol_policy = 0, bond_policy = 0, quiet = True ):
    if method == "BIN":
        mmethod = 1
    elif method == "NSQ":
        mmethod = 0
    else:
        raise RuntimeError("Unknown method",method,"!")

    return neighborize_.nearest_neighbours_dist( block_data.get_ref_(),
                                                 itype, jtype,
                                                 mmethod, dims, rc,
                                                 mol_policy,
                                                 bond_policy, quiet )

def neigh_list_dist_indexed( block_data, itype, jtype, method, dims, rc,
                             mol_policy = 0, bond_policy = 0, quiet = True ):
    if method == "BIN":
        mmethod = 1
    elif method == "NSQ":
        mmethod = 0
    else:
        raise RuntimeError("Unknown method",method,"!")

    itype_vec = block_data_.VectorInt(itype)
    jtype_vec = block_data_.VectorInt(jtype)

    return neighborize_.nearest_neighbours_dist_indexed( block_data.get_ref_(),
                                                         itype_vec, jtype_vec,
                                                         mmethod, dims, rc,
                                                         mol_policy,
                                                         bond_policy, quiet )



def neighbour_strain( block_data, r0, itype, jtype, method, dims, rc ):
    return neighborize_.neighbour_strain( block_data.get_ref_(), r0,
                                          itype, jtype, method, dims, rc )

def molecular_connections( b, neighs ):
    return neighborize_.molecular_connections( b.get_ref_(), neighs )

def molecular_network( b, neighs, conns ):
    return neighborize_.molecular_network( b.get_ref_(), neighs, conns )
