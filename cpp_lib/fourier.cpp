#include "fourier.hpp"
#include "my_assert.hpp"


#ifdef HAVE_ARMADILLO
#include <armadillo>
constexpr const bool use_armadillo = true;
#else
constexpr const bool use_armadillo = false;
#endif // HAVE_ARMADILLO


namespace lammps_tools {

namespace fourier {

using cx_double = std::complex<double>;


arma::cx_mat data_to_mat( int Nx, int Ny, const std::vector<cx_double> &data );
std::vector<cx_double> mat_to_data( int Nx, int Ny, const arma::cx_mat &X );
arma::cx_vec data_to_vec( const std::vector<cx_double> &data );
std::vector<cx_double> vec_to_data( const arma::cx_vec &X );


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
	return vec_to_data( arma::fft( X ) );
}


std::vector<cx_double > ifft_shift1( const std::vector<cx_double > &fft )
{
	arma::cx_vec X(fft);
	return vec_to_data( arma::ifft( X ) );
}


std::vector<cx_double > fft_shift2( int  Nx, int Ny,
                                    const std::vector<cx_double > &fft )
{
	arma::cx_mat X = data_to_mat( Nx, Ny, fft );
	arma::cx_mat Xshift(Nx, Ny);

	int Nx_lower = Nx/2;
	int Ny_lower = Ny/2;

	// Check. If Nx or Ny are uneven, then it is kinda weird...

	// First quadrant:
	auto XQ1 = X.submat(0, 0, Nx_lower, Ny_lower);
	auto XQ2 = X.submat(Nx_lower, 0, Nx, Ny_lower);
	auto XQ3 = X.submat(0, Ny_lower, Nx_lower, Ny);
	auto XQ4 = X.submat(Nx_lower, Ny_lower, Nx, Ny);

	auto XsQ1 = Xshift.submat(0, 0, Nx_lower, Ny_lower);
	auto XsQ2 = Xshift.submat(Nx_lower, 0, Nx, Ny_lower);
	auto XsQ3 = Xshift.submat(0, Ny_lower, Nx_lower, Ny);
	auto XsQ4 = Xshift.submat(Nx_lower, Ny_lower, Nx, Ny);

	XsQ1 = XQ4;
	XsQ2 = XQ3;
	XsQ3 = XQ2;
	XsQ4 = XQ1;

	return mat_to_data(Nx, Ny, Xshift);
}


std::vector<cx_double > ifft_shift2( int  Nx, int Ny,
                                     const std::vector<cx_double > &fft )
{
	arma::cx_mat X = data_to_mat( Nx, Ny, fft );
	arma::cx_mat Xishift( Nx, Ny ); // , fft );

	int Nx_lower = Nx/2;
	int Ny_lower = Ny/2;

	// Check. If Nx or Ny are uneven, then it is kinda weird...

	// First quadrant:
	auto XQ1 = X.submat(0, 0, Nx_lower, Ny_lower);
	auto XQ2 = X.submat(Nx_lower, 0, Nx, Ny_lower);
	auto XQ3 = X.submat(0, Ny_lower, Nx_lower, Ny);
	auto XQ4 = X.submat(Nx_lower, Ny_lower, Nx, Ny);

	auto XsQ1 = Xishift.submat(0, 0, Nx_lower, Ny_lower);
	auto XsQ2 = Xishift.submat(Nx_lower, 0, Nx, Ny_lower);
	auto XsQ3 = Xishift.submat(0, Ny_lower, Nx_lower, Ny);
	auto XsQ4 = Xishift.submat(Nx_lower, Ny_lower, Nx, Ny);

	XsQ1 = XQ4;
	XsQ2 = XQ3;
	XsQ3 = XQ2;
	XsQ4 = XQ1;

	return mat_to_data(Nx, Ny, Xishift);
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
	arma::mat X( Nx, Ny );

	for( std::size_t i = 0; i < data.size(); ++i ){
		std::size_t ix = i % Nx;
		std::size_t iy = i / Nx;
		std::size_t ifull = ix + Nx*iy;
		my_assert( __FILE__, __LINE__, ifull == i,
		           "Incorrect reduction to index!" );
		X( ix, iy ) = data[i];
	}

	arma::cx_mat Y = fft2( X );



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



std::vector<cx_double> fft_int( int Nx, int Ny, int Nz,
                             const std::vector<int> &data )
{
	std::vector<double> data_d( data.begin(), data.end() );
	return fft_double( Nx, Ny, Nz, data_d );
}



arma::cx_mat data_to_mat( int Nx, int Ny, const std::vector<cx_double> &data )
{
	arma::cx_mat X( Nx, Ny );
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


std::vector<cx_double> mat_to_data( int Nx, int Ny, const arma::cx_mat &X )
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


arma::cx_vec data_to_vec( const std::vector<cx_double> &data )
{
	arma::cx_vec vec( data.size() );
	for( std::size_t i = 0; i < data.size(); ++i ){
		vec(i) = data[i];
	}
	return vec;
}


std::vector<cx_double> vec_to_data( const arma::cx_vec &X )
{
	std::vector<cx_double> data( X.size() );
	for( std::size_t i = 0; i < X.size(); ++i ){
		data[i] = X(i);
	}
	return data;
}





} // namespace fourier

} // lammps_tools
