"""
This module imports the center of mass calculation tools.
"""

import center_of_mass_

def center_of_mass(b):
    """ Returns center of mass of block. """
    com = center_of_mass_.center_of_mass(b.get_ref_())
    return com

def geometric_center(b):
    """ Returns geometric center of block. """
    com = center_of_mass_.geometric_center(b.get_ref_())
    return com
