import scatter_
import conversions_
import numpy as np

def rayleigh_gans( b, qs, ids = None, radius = None, d_epsilon = 1e-2,
                   position_scale = 1.0 ):
    """ Calculates Rayleigh-Gans scattering data for all atoms or given ids. """
    if ids is None:
        ids = b.ids
    else:
        idslist = ids.tolist()
        ids = conversions_.list_to_int_vector(idslist)

    if radius is None:
        radius = np.ones( b.meta.N, dtype = float )

    rradius = conversions_.list_to_double_vector(radius.tolist())
    qqs  = conversions_.list_to_double_vector(qs.tolist())

    return scatter_.rayleigh_gans_( b.get_ref_(), qqs, rradius, position_scale,
                                    d_epsilon, ids )
