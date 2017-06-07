#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "../../../cpp_lib/skeletonize.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)

PYBIND11_PLUGIN(skeletonize_) {
	using namespace lammps_tools::skeletonize;

	pybind11::module m("skeletonize", "Exposes skeletonization functions.");

	m.def( "euclidian_distance_transform", &euclidian_distance_transform,
	       "Calculates the Euclidian distance transform for given block." );

	m.def( "neighbour_strain", &neighbour_strain,
	       "Calculates the average inter-neighbour strain" );

	return m.ptr();
}
