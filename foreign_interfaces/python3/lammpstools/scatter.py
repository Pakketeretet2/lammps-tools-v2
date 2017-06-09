import scatter_

def rayleigh_gans( b, ids = None ):
    """ Calculates Rayleigh-Gans scattering data for all atoms or given ids. """
    if ids is None:
        ids = b.ids

    scatter_.rayleigh_gans_( b.get_ref(), ids )
