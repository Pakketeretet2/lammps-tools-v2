#include "fourier.hpp"
#include "random_generator.hpp"
#include "util.hpp"

#include <algorithm>
#include <catch.hpp>
#include <fstream>


TEST_CASE ( "Fourier in one dimension", "[fourier1D]" )
{
	std::size_t N = 55;
	double dx = pi / N;
	double w1 = 0.5;
	double w2 = 0.4;
	std::vector<double> X(N);
	RanGen rng( 1 );

	std::ofstream sig_out( "signal.dat" );
	std::ofstream fft_out( "ffts.dat" );

	for( std::size_t i = 0; i < N; ++i ){
		double x = i*dx;
		X(i) = sin(w1*x) -0.1 + 0.2*rng.uniform() - cos(x);
		sig_out << x << " " << X(i) << "\n";
	}


	std::vector<double> fft_abs, fft_imag, fft_real;
	using lammps_tools::fourier::cx_double;
	std::vector<cx_double> fft = fft( N, 1, 1, X );
	fft_abs = lammps_tools::fourier::fft_abs( fft );
	fft_real = lammps_tools::fourier::fft_real( fft );
	fft_imag = lammps_tools::fourier::fft_imag( fft );

	double idx = 1.0 / dx;
	for( std::size_t i = 0; i < N; ++i ){
		double ix = idx*i;
		fft_out << ix << " " << fft_abs[i] << " " << fft_real[i]
		        << " " << fft_imag[i] << "\n";
	}


}


TEST_CASE ( "Fourier in two dimensions", "[fourier2D]" )
{

}
