#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include <pybind11/numpy.h>

#include "../../../cpp_lib/neighborize.hpp"
#include "../../../cpp_lib/rdf.hpp"
#include "../../../cpp_lib/skeletonize.hpp"


PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)

PYBIND11_PLUGIN(neighborize_) {
	using namespace lammps_tools;


	pybind11::module m("neighborize", "Exposes nearest neighbour functions.");

	m.def( "rdf", &neighborize::rdf,
	       "Calculates the rdf for given block" );
	m.def( "coord", &neighborize::coord,
	       "Calculates the coordination from given rdf" );

	m.def( "nearest_neighbours_dist", &neighborize::nearest_neighs,
	       "Calculates the nearest neighbour list" );




	m.def( "euclidian_distance_transform",
	       &skeletonize::euclidian_distance_transform,
	       "Calculates the Euclidian distance transform for given block." );
	m.def( "neighbour_strain", &skeletonize::neighbour_strain,
	       "Calculates the average inter-neighbour strain" );

	return m.ptr();
}
