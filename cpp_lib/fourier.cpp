#include "fourier.hpp"
#include "my_assert.hpp"
#include "constants.hpp"

#ifdef HAVE_ARMADILLO
#include <armadillo>
constexpr const bool use_armadillo = true;
#else
constexpr const bool use_armadillo = false;
#endif // HAVE_ARMADILLO

#include "histogram.hpp"


namespace lammps_tools {

namespace fourier {

using cx_double = std::complex<double>;



template <typename ret_type, typename T>
ret_type data_to_mat( int Nx, int Ny, const std::vector<T> &data );

template <typename ret_type, typename T>
std::vector<T> mat_to_data( int Nx, int Ny, const ret_type &X );

template <typename ret_type, typename T>
ret_type data_to_vec( const std::vector<T> &data );

template <typename ret_type, typename T>
std::vector<T> vec_to_data( const ret_type &X );



std::vector<double> fft_abs( const std::vector<cx_double > &fft )
{
	std::vector<double> a( fft.size() );
	for( std::size_t i = 0; i < fft.size(); ++i ){
		a[i] = std::abs(fft[i]);
	}
	return a;
}


std::vector<double> fft_real( const std::vector<cx_double > &fft )
{
	std::vector<double> re( fft.size() );
	for( std::size_t i = 0; i < fft.size(); ++i ){
		re[i] = std::real(fft[i]);
	}
	return re;
}


std::vector<double> fft_imag( const std::vector<cx_double > &fft )
{
	std::vector<double> im( fft.size() );
	for( std::size_t i = 0; i < fft.size(); ++i ){
		im[i] = std::imag(fft[i]);
	}
	return im;
}


std::vector<cx_double > fft_shift1( const std::vector<cx_double > &fft )
{
	arma::cx_vec X(fft);
	std::size_t Nx = fft.size();
	std::size_t Nlo = Nx/2;

	if( 2*Nlo == Nx ){
		// If Nx is uneven, say 8, we need the following remap:
		// [ 0 1 2 3 4 5 6 7 ] -->
		// [ 4 5 6 7 0 1 2 3 ]
		for( std::size_t i = 0; i < Nlo; ++i ){
			X(i) = fft[i+Nlo];
			X(i+Nlo) = fft[i];
		}
	}else{
		// If Nx is uneven, say 7, we need the following remap:
		// [ 0 1 2 3 4 5 6 ] -->
		// [ 4 5 6 3 0 1 2 ]
		for( std::size_t i = 0; i < Nlo; ++i ){
			X(i) = fft[i+Nlo+1];
			X(i+Nlo+1) = fft[i];
		}
	}
	return vec_to_data<arma::cx_vec, cx_double>( X );
}


std::vector<cx_double > ifft_shift1( const std::vector<cx_double > &fft )
{
	return fft_shift1( fft );
}


std::vector<cx_double > fft_shift2( int  Nx, int Ny,
                                    const std::vector<cx_double > &fft )
{
	arma::cx_mat X = data_to_mat<arma::cx_mat, cx_double>( Nx, Ny, fft );
	arma::cx_mat Xshift(Nx, Ny);

	int Nx_lower = Nx/2;
	int Ny_lower = Ny/2;
	int Nx_upper = Nx_lower;
	int Ny_upper = Nx_lower;

	// If Nx is uneven
	bool x_even = 2*Nx_lower == Nx;
	bool y_even = 2*Ny_lower == Ny;

	if( !x_even ) Nx_upper++;
	if( !y_even ) Ny_upper++;

	for( std::size_t i = 0; i < Nx_lower; ++i ){
		// Quadrant 1 --> Quadrant 4
		for( std::size_t j = 0; j < Ny_lower; ++j ){
			Xshift(i,j) = X(Nx_upper+i, Ny_upper+j);
		}
		// Quadrant 2 --> Quadrant 3
		for( std::size_t j = 0; j < Ny_lower; ++j ){
			Xshift(i+Nx_upper,j) = X(i, Ny_upper+j);
		}
		// Quadrant 3 --> Quadrant 2
		for( std::size_t j = 0; j < Ny_lower; ++j ){
			Xshift(i,j+Ny_upper) = X(Nx_upper+i, j);
		}
		// Quadrant 4 --> Quadrant 1
		for( std::size_t j = 0; j < Ny_lower; ++j ){
			Xshift(Nx_upper + i, Ny_upper + j) = X(i, j);
		}
	}

	return mat_to_data<arma::cx_mat, cx_double>(Nx, Ny, Xshift);
}


std::vector<cx_double > ifft_shift2( int  Nx, int Ny,
                                     const std::vector<cx_double > &fft )
{
	return fft_shift2( Nx, Ny, fft );
}


std::vector<cx_double > fft_shift3( int  Nx, int Ny, int Nz,
                                    const std::vector<cx_double > &fft )
{
	std::vector<cx_double> shift( Nx*Ny*Nz );
	my_runtime_error(__FILE__, __LINE__,
	                 "3D Fourier transform not supported!");
	return shift;
}


std::vector<cx_double > ifft_shift3( int  Nx, int Ny, int Nz,
                                     const std::vector<cx_double > &fft )
{
	std::vector<cx_double> ishift( Nx*Ny*Nz );
	my_runtime_error(__FILE__, __LINE__,
	                 "3D Fourier transform not supported!");
	return ishift;
}


std::vector<cx_double > fft_shift( int  Nx, int Ny, int Nz,
                                   const std::vector<cx_double > &fft )
{
	if( Nz == 1 ){
		if( Ny == 1 ){
			return fft_shift1( fft );
		}else{
			return fft_shift2( Nx, Ny, fft );
		}
	}else{
		return fft_shift3( Nx, Ny, Nz, fft );
	}
}


std::vector<cx_double > ifft_shift( int  Nx, int Ny, int Nz,
                                    const std::vector<cx_double > &fft )
{
	if( Nz == 1 ){
		if( Ny == 1 ){
			return ifft_shift1( fft );
		}else{
			return ifft_shift2( Nx, Ny, fft );
		}
	}else{
		return ifft_shift3( Nx, Ny, Nz, fft );
	}
}



std::vector<cx_double> fft3_impl( int Nx, int Ny, int Nz,
                               const std::vector<double> & )
{
	std::vector<cx_double> fft( Nx*Ny*Nz );
	my_runtime_error(__FILE__, __LINE__,
	                 "3D Fourier transform not supported!");
	return fft;
}



std::vector<cx_double> fft2_impl( int Nx, int Ny,
                                  const std::vector<double> &data )
{
	std::vector<cx_double> fft( Nx*Ny );

#ifdef HAVE_ARMADILLO
	// Convert data to a matrix:
	arma::mat X = data_to_mat<arma::mat, double>( Nx, Ny, data );
	arma::cx_mat Y = fft2( X );
	fft = mat_to_data<arma::cx_mat, cx_double>( Ny, Ny, Y );
#endif
	return fft;
}



std::vector<cx_double> fft1_impl( const std::vector<double> &data )
{
	std::vector<cx_double> fft( data.size() );
#ifdef HAVE_ARMADILLO
	// Convert data to a matrix:
	arma::vec X( data );
	arma::cx_vec fX = arma::fft(X);

	for( std::size_t i = 0; i < data.size(); ++i ){
		fft[i] = fX(i);
	}
#endif

	return fft;
}


std::vector<cx_double> fft3_impl( int Nx, int Ny, int Nz,
                               const std::vector<cx_double> & )
{
	std::vector<cx_double> fft( Nx*Ny*Nz );
	my_runtime_error(__FILE__, __LINE__,
	                 "3D Fourier transform not supported!");
	return fft;
}



std::vector<cx_double> fft2_impl( int Nx, int Ny,
                                  const std::vector<cx_double> &data )
{
	std::vector<cx_double> fft( Nx*Ny );

#ifdef HAVE_ARMADILLO
	// Convert data to a matrix:
	arma::cx_mat X = data_to_mat<arma::cx_mat, cx_double>( Nx, Ny, data );
	arma::cx_mat Y = fft2( X );
	fft = mat_to_data<arma::cx_mat, cx_double>( Ny, Ny, Y );
#endif
	return fft;
}



std::vector<cx_double> fft1_impl( const std::vector<cx_double> &data )
{
	std::vector<cx_double> fft( data.size() );
#ifdef HAVE_ARMADILLO
	// Convert data to a matrix:
	arma::cx_vec X = data;
	arma::cx_vec fX = arma::fft(X);

	for( std::size_t i = 0; i < data.size(); ++i ){
		fft[i] = fX(i);
	}
#endif

	return fft;
}




std::vector<cx_double > fft( int Nx, int Ny, int Nz,
                             const std::vector<cx_double> &data )
{
	return fft_cx_double( Nx, Ny, Nz, data );
}



std::vector<cx_double> fft( int Nx, int Ny, int Nz,
                            const std::vector<double> &data )
{
	return fft_double( Nx, Ny, Nz, data );
}


std::vector<cx_double> fft_double( int Nx, int Ny, int Nz,
                                   const std::vector<double> &data )
{
	if( !use_armadillo ){
		std::vector<cx_double> fft( data.size() );
		my_warning( __FILE__, __LINE__,
		            "Cannot perform FFT without Armadillo1" );
		return fft;
	}else{
		if( Nz == 1 ){
			if( Ny == 1 ){
				return fft1_impl( data );
			}else{
				return fft2_impl( Nx, Ny, data );
			}
		}else{
			return fft3_impl( Nx, Ny, Nz, data );
		}
	}
}


std::vector<cx_double> fft_cx_double( int Nx, int Ny, int Nz,
                                      const std::vector<cx_double> &data )
{
	if( !use_armadillo ){
		std::vector<cx_double> fft( data.size() );
		my_warning( __FILE__, __LINE__,
		            "Cannot perform FFT without Armadillo1" );
		return fft;
	}else{
		if( Nz == 1 ){
			if( Ny == 1 ){
				return fft1_impl( data );
			}else{
				return fft2_impl( Nx, Ny, data );
			}
		}else{
			return fft3_impl( Nx, Ny, Nz, data );
		}
	}
}



std::vector<cx_double> fft_int( int Nx, int Ny, int Nz,
                                const std::vector<int> &data )
{
	std::vector<double> data_d( data.begin(), data.end() );
	return fft_double( Nx, Ny, Nz, data_d );
}


template <typename ret_type, typename T>
ret_type data_to_mat( int Nx, int Ny, const std::vector<T> &data )
{
	ret_type X( Nx, Ny );
	for( std::size_t i = 0; i < data.size(); ++i ){
		std::size_t ix = i % Nx;
		std::size_t iy = i / Nx;
		std::size_t ifull = ix + Nx*iy;
		my_assert( __FILE__, __LINE__, ifull == i,
		           "Incorrect reduction to index!" );
		X( ix, iy ) = data[i];
	}
	return X;

}


template <typename ret_type, typename T>
std::vector<T> mat_to_data( int Nx, int Ny, const ret_type &X )
{
	std::vector<cx_double> data( Nx*Ny );
	for( std::size_t i = 0; i < data.size(); ++i ){
		std::size_t ix = i % Nx;
		std::size_t iy = i / Nx;
		std::size_t ifull = ix + Nx*iy;
		my_assert( __FILE__, __LINE__, ifull == i,
		           "Incorrect reduction to index!" );
		data[i] = X( ix, iy );
	}
	return data;
}


template <typename ret_type, typename T>
ret_type data_to_vec( const std::vector<T> &data )
{
	arma::cx_vec vec( data.size() );
	for( std::size_t i = 0; i < data.size(); ++i ){
		vec(i) = data[i];
	}
	return vec;
}


template <typename ret_type, typename T>
std::vector<T> vec_to_data( const ret_type &X )
{
	std::vector<cx_double> data( X.size() );
	for( std::size_t i = 0; i < X.size(); ++i ){
		data[i] = X(i);
	}
	return data;
}


std::vector<double> make_frequency_grid( const std::vector<double> &x )
{
	// Grid has to be equispaced.
	std::size_t N = x.size();

	double dx = x[1] - x[0];
	for( std::size_t i = 0; i < N - 1; ++i ){
		my_assert( __FILE__, __LINE__,
		           std::fabs(x[i+1]-x[i]-dx) < 1e-12,
		           "Grid needs to be equispaced!" );
	}

	double fmax = 1.0 / dx;
	double fmin = fmax / N;
	double offset = -fmax/2.0;
	std::vector<double> grid(x.size());
	bool N_is_odd = N % 2;

	if( N_is_odd ){
		offset += 0.5*fmin;
		for( std::size_t j = 0; j < N/2; ++j ){
			double fj = offset;
			grid[j] = fj + j*fmin;
		}
		offset -= fmin;
		for( std::size_t j = N/2+1; j < N; ++j ){
			double fj = offset;
			grid[j] = fj + j*fmin;
		}
	}else{
		for( std::size_t j = 0; j < N; ++j ){
			double fj = offset;
			grid[j] = fj + j*fmin;
		}
	}

	return grid;
}




} // namespace fourier


} // lammps_tools
