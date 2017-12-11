import icosahedra_, block_data_


def dihedral_angles( b, patches1, patches2, conns ):
    """ Returns the dihedral angles. """

    vpatches1 = block_data_.VectorInt( patches1 )
    vpatches2 = block_data_.VectorInt( patches2 )
    return icosahedra_.dihedral_angles( b.get_ref_(), vpatches1,
                                        vpatches2, conns )

def patch_types_to_ids( b, patches1, patches2 ):
    """ Converts patch types to an id list. """
    vpatches1 = block_data_.VectorInt( patches1 )
    vpatches2 = block_data_.VectorInt( patches2 )
    return icosahedra_.patch_types_to_ids( b.get_ref_(), vpatches1, vpatches2 )


def get_mol_connections( b, patches1, patches2, method, dims, rc ):
    """ Determines the molecular connection list. """
    if method == "BIN":
        mmethod = 1
    else:
        mmethod = 0

    vpatches1 = block_data_.VectorInt( patches1 )
    vpatches2 = block_data_.VectorInt( patches2 )
    return icosahedra_.get_mol_connections( b.get_ref_(), vpatches1, vpatches2,
                                            mmethod, dims, rc )

def get_scattering_intensity( b, conns, I0 ):
    """ Calculates an effective, relative scattering intensity. """
    return icosahedra_.get_scattering_intensity( b.get_ref_(), conns, I0 )
