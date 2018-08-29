#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../../../cpp_lib/density_distribution.hpp"

#include "make_vectors_opaque.hpp"


PYBIND11_PLUGIN(density_distribution_) {
	using namespace lammps_tools;

	pybind11::module m("density_distribution_",
	                   "Exposes density distribution functions.");

	m.def("box_atoms", &density::box_atoms,
	      "Boxes atoms in a grid to extract density distribution." );

	m.def("density_distribution", &density::density_distribution,
	      "Determines the distribution of density across boxes.");

	m.def("box_count", &density::box_count,
	      "Count the number of tatoms in each box.");


	return m.ptr();
}
