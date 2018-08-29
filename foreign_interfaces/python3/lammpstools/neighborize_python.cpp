#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../../../cpp_lib/neighborize.hpp"
#include "../../../cpp_lib/rdf.hpp"
#include "../../../cpp_lib/skeletonize.hpp"
#include "../../../cpp_lib/cluster_finder.hpp"

#include "make_vectors_opaque.hpp"

PYBIND11_PLUGIN(neighborize_) {
	using namespace lammps_tools;

	pybind11::module m("neighborize", "Exposes nearest neighbor functions.");

	m.def( "rdf", &neighborize::rdf,
	       "Calculates the rdf for given block" );
	m.def( "coord", &neighborize::coord,
	       "Calculates the coordination from given rdf" );

	m.def( "nearest_neighbors_dist", &neighborize::nearest_neighs,
	       "Calculates the nearest neighbor list" );

	m.def( "nearest_neighbors_dist_indexed",
	       &neighborize::nearest_neighs_indexed,
	       "Calculates nearest neighbor list for specified atoms only." );

	m.def( "neighbors_to_network",
	       &neighborize::neigh_list_to_network,
	       "Converts the neighbor list to a connectivity network." );

	m.def( "euclidian_distance_transform",
	       &skeletonize::euclidian_distance_transform,
	       "Calculates the Euclidian distance transform for given block." );
	m.def( "neighbor_strain", &skeletonize::neighbor_strain,
	       "Calculates the average inter-neighbor strain" );

	m.def( "molecular_connections", &neighborize::get_molecular_connections,
	       "Find the connections between molecules from neigh list." );

	m.def( "get_empty_neighbor_list",
	       &neighborize::get_empty_neigh_list,
	       "Creates an empty neighbor list to pass to other functions." );

	return m.ptr();
}
