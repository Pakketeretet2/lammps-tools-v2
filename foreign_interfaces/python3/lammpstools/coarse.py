import coarsening_, block_data_


def box_to_psi( b, boxes, N_boxes, seed = 123, quiet = False ):
    """ Counts the numbers of atom types 1 and 2 in boxes. """
    return coarsening_.box_to_psi( b.get_ref_(), boxes, N_boxes, seed, quiet )


def data_to_atom( b, boxes, data ):
    """ Maps per-box data onto per-atom data """
    if type(data) == block_data_.VectorInt:
        return coarsening_.map_to_atom_int( b.get_ref_(), boxes, data )

    elif type(data) == block_data_.VectorDouble:
        return coarsening_.map_to_atom_double( b.get_ref_(), boxes, data )
    elif type(data) == block_data_.VectorComplexDouble:
        return coarsening_.map_to_atom_cx_double( b.get_ref_(), boxes, data )
    elif type(data) == block_data_.VectorComplexInt:
        return coarsening_.map_to_atom_cx_int( b.get_ref_(), boxes, data )
    else:
        return None
