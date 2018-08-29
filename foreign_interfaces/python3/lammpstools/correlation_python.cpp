#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../../../cpp_lib/correlation.hpp"

#include "make_vectors_opaque.hpp"

PYBIND11_PLUGIN(correlation_) {
	using namespace lammps_tools;

	pybind11::module m("correlation", "Exposes correlation functions.");

	m.def( "correlate_double", &lammps_tools::correlate::correlate_double,
	       "Correlates data" );

	m.def( "correlate_int", &lammps_tools::correlate::correlate_int,
	       "Correlates data" );


	return m.ptr();
}
