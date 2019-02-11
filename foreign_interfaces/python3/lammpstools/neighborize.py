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

    return neighborize_.nearest_neighbors_dist( block_data.get_ref_(),
                                                itype, jtype,
                                                mmethod, dims, rc,
                                                mol_policy,
                                                bond_policy, quiet )

def find_clusters( neigh_list, sort_clusters = False ):
    clusters = neighborize_.find_clusters(neigh_list)
    if sort_clusters:
        clusters.sort( key = len, reverse = True )
    return clusters


def neigh_list_dist_indexed( block_data, i_idxs, j_idxs, method, dims, rc,
                             mol_policy = 0, bond_policy = 0, quiet = True ):
    if method == "BIN":
        mmethod = 1
    elif method == "NSQ":
        mmethod = 0
    else:
        raise RuntimeError("Unknown method",method,"!")

    i_idx_vec = block_data_.VectorInt(i_idxs)
    j_idx_vec = block_data_.VectorInt(j_idxs)

    return neighborize_.nearest_neighbors_dist_indexed( block_data.get_ref_(),
                                                        i_idx_vec, j_idx_vec,
                                                        mmethod, dims, rc,
                                                        mol_policy,
                                                        bond_policy, quiet )



def neighbor_strain( block_data, r0, itype, jtype, method, dims, rc ):
    return neighborize_.neighbor_strain( block_data.get_ref_(), r0,
                                         itype, jtype, method, dims, rc )

def molecular_connections( b, neighs ):
    return neighborize_.molecular_connections( b.get_ref_(), neighs )

def get_empty_neighbor_list():
    return neighborize_.get_empty_neighbor_list()
