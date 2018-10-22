#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../../../cpp_lib/fourier.hpp"

#include "make_vectors_opaque.hpp"

void deprecated_err()
{
	std::cerr << "This function has been moved from fourier to histogram!\n";
	throw -1;
}

PYBIND11_PLUGIN(fourier_) {
	using namespace lammps_tools;
	using fourier::cx_double;

	pybind11::module m("fourier_", "Exposes Fourier transform functions.");

	m.def( "fft_int", &fourier::fft_int,
	       "Perform Fourier transform on discrete data." );

	m.def( "fft_real", &fourier::fft_real,
	       "Grab the real part of fourier transform data." );

	m.def( "fft_imag", &fourier::fft_imag,
	       "Grab the imaginary part of fourier transform data." );

	m.def( "fft_shift", &fourier::fft_shift,
	       "Shift FFT so that zero frequency is in center." );

	m.def( "fft_ishift", &fourier::ifft_shift,
	       "Inverse shift (undo fft_shift" );

	m.def( "make_frequency_grid", &fourier::make_frequency_grid,
	       "Generates a fequency grid for use with fft." );

	m.def( "make_radial_grid", &deprecated_err,
	       "Generates a radial grid from a 2D grid, linear in r." );

	m.def( "make_radial2_grid", &deprecated_err,
	       "Generates a radial grid from a 2D grid, quadratic in r." );

	m.def( "radial_avg_double", &deprecated_err,
	       "Radially average data on a 2D grid, linear in r." );

	m.def( "radial_avg_cx_double", &deprecated_err,
	       "Radially average data on a 2D grid, linear in r." );

	m.def( "radial2_avg_double", &deprecated_err,
	       "Radially average data on a 2D grid, quadratic in r." );

	m.def( "radial2_avg_cx_double", &deprecated_err,
	       "Radially average data on a 2D grid, quadratic in r." );
	return m.ptr();
}
