#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../../../cpp_lib/histogram.hpp"
#include "../../../cpp_lib/fourier.hpp"


#include "make_vectors_opaque.hpp"

PYBIND11_PLUGIN(histogram_) {

	using namespace lammps_tools;

	pybind11::module m("histogram_", "Exposes histogram tools." );

	m.def( "make_histogram_double", &histogram::make_histogram_double,
	       "Makes a histogram of doubles." );

	m.def( "make_histogram_int", &histogram::make_histogram_int,
	       "Makes a histogram of int." );

	using fourier::cx_double;
	m.def( "make_radial_grid", &histogram::make_radial_grid,
	       "Generates a radial grid from a 2D grid, linear in r." );

	m.def( "make_radial2_grid", &histogram::make_radial2_grid,
	       "Generates a radial grid from a 2D grid, quadratic in r." );

	m.def( "radial_avg_double", &histogram::radial_avg<double>,
	       "Radially average data on a 2D grid, linear in r." );

	m.def( "radial_avg_cx_double", &histogram::radial_avg<cx_double>,
	       "Radially average data on a 2D grid, linear in r." );

	m.def( "radial2_avg_double", &histogram::radial2_avg<double>,
	       "Radially average data on a 2D grid, quadratic in r." );

	m.def( "radial2_avg_cx_double", &histogram::radial2_avg<cx_double>,
	       "Radially average data on a 2D grid, quadratic in r." );


	return m.ptr();
}
