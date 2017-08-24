#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../../../cpp_lib/neighborize.hpp"
#include "../../../cpp_lib/rdf.hpp"
#include "../../../cpp_lib/skeletonize.hpp"
#include "../../../cpp_lib/cluster_finder.hpp"

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

	m.def( "nearest_neighbours_dist_indexed",
	       &neighborize::nearest_neighs_indexed,
	       "Calculates nearest neighbour list for specified atoms only." );

	m.def( "euclidian_distance_transform",
	       &skeletonize::euclidian_distance_transform,
	       "Calculates the Euclidian distance transform for given block." );
	m.def( "neighbour_strain", &skeletonize::neighbour_strain,
	       "Calculates the average inter-neighbour strain" );

	m.def( "molecular_connections", &neighborize::get_molecular_connections,
	       "Find the connections between molecules from neigh list." );
	m.def( "molecular_network", &neighborize::make_molecular_networks,
	       "Find the molecular network from the molecule connections." );


	return m.ptr();
}
