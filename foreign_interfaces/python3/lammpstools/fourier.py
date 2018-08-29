import fourier_
import block_data

def fft2( Nx, Ny, data ):
    """ Two-dimensional transform of given data on Nx by Ny grid. """
    return fourier_.fft( Nx, Ny, 1, data )


def fft3( Nx, Ny, Nz, data ):
    """ Two-dimensional transform of given data on Nx by Ny by Nz grid. """
    return fourier_.fft( Nx, Ny, Nz, data )


def fft2_int( Nx, Ny, data ):
    """ Two-dimensional transform of given data on Nx by Ny grid. """
    return fourier_.fft_int( Nx, Ny, 1, data )


def fft3_int( Nx, Ny, Nz, data ):
    """ Two-dimensional transform of given data on Nx by Ny by Nz grid. """
    return fourier_.fft_int( Nx, Ny, Nz, data )

def fft3_double( Nx, Ny, Nz, data ):
    """ Two-dimensional transform of given data on Nx by Ny by Nz grid. """
    return fourier_.fft_double( Nx, Ny, Nz, data )


def fft_real( data ):
    """ Grab real part of FFT data. """
    return fourier_.fft_real( data )

def fft_imag( data ):
    """ Grab imaginary part of FFT data. """
    return fourier_.fft_imag( data )

def fft_abs( data ):
    """ Grab absolute value of FFT data. """
    return fourier_.fft_abs( data )

def fft_shift( Nx, Ny, Nz, data ):
    """ Center FFT on zero frequency. """
    return fourier_.fft_shift( Nx, Ny, Nz, data )

def fft_ishift( Nx, Ny, Nz, data ):
    """ Center FFT on zero frequency. """
    return fourier_.fft_ishift( Nx, Ny, Nz, data )

def make_frequency_grid( time_grid ):
    """ Make a frequency grid from time grid. """
    return fourier_.make_frequency_grid( time_grid )

def make_radial_grid( xgrid, ygrid, Nbins ):
    """ Make a radial grid from a 2D grid, linear in r. """
    return fourier_.make_radial_grid( xgrid, ygrid, Nbins )

def make_radial2_grid( xgrid, ygrid, Nbins ):
    """ Make a radial grid from a 2D grid, quadratic in r. """
    return fourier_.make_radial2_grid( xgrid, ygrid, Nbins )


def radial_avg( xgrid, ygrid, data, Nbins ):
    """ Radially average data over Nbins bins, linear in r. """
    # Dispatch right one based on type of data:
    if type(data) == block_data.VectorDouble:
        return fourier_.radial_avg_double( xgrid, ygrid, data, Nbins )
    elif type(data) == block_data.VectorDouble:
        return fourier_.radial_avg_cx_double( xgrid, ygrid, data, Nbins )
    else:
        print( "Type dispatch for ", type(data),
               " failed in radial_avg!", file = sys.stderr )
        return None


def radial2_avg( xgrid, ygrid, data, Nbins ):
    """ Radially average data over Nbins bins, quadratic in r. """
    # Dispatch right one based on type of data:
    if type(data) == block_data.VectorDouble:
        return fourier_.radial2_avg_double( xgrid, ygrid, data, Nbins )
    elif type(data) == block_data.VectorDouble:
        return fourier_.radial2_avg_cx_double( xgrid, ygrid, data, Nbins )
    else:
        print( "Type dispatch for ", type(data),
               " failed in radial_avg!", file = sys.stderr )
        return None
