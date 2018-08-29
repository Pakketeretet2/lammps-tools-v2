import density_distribution_
import sys


def box_atoms( b, Nx, dims, cutoff ):
    """ Assigns each atom a box index. """

    dx = 3.0
    boxes = density_distribution_.box_atoms( b.get_ref_(), Nx, dx, dims )

    # You also need to pass back the number of boxes...

    return boxes


def box_count( b, boxes, N_boxes, itype, quiet = True ):
    return density_distribution_.box_count( b.get_ref_(), boxes, N_boxes, itype, quiet )




def density_distribution( b, Nx, itype, dims ):
    """ Determines the distribution of particles of type itype over boxes. """
    return density_distribution_.density_distribution( b.get_ref_(), Nx, itype, dims )
