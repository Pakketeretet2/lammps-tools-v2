import bond_order_

import numpy as np

def psi_n( b, neighs, n, axis ):
    """ Calculates psi n. """

    psi_real = np.zeros( b.meta.N, dtype = np.float64 )
    psi_imag = np.zeros( b.meta.N, dtype = np.float64 )

    psi_avg = bond_order_.compute_psi_n( b.get_ref_(), neighs, n, axis,
                                         psi_real, psi_imag )

    return psi_avg, psi_real, psi_imag
