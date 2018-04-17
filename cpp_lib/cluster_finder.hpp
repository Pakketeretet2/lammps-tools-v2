#ifndef CLUSTER_FINDER_HPP
#define CLUSTER_FINDER_HPP

/**
   \file cluster_finder.hpp
*/

#include "neighborize.hpp"

#include <vector>

namespace lammps_tools {

namespace neighborize {


/**
   \brief Converts a neighbor list to a list of clusters.

   The neigh list clusters will contain nc neighbor lists, where nc is the
   number of clusters, and will contain all the clusters, that is, a list of
   all items in the neigh_list conns that are connected to each other.

   \param conns     The neighbor list to process

   \returns a list of clusters.
 */
neigh_list neigh_list_to_clusters( const neigh_list &conns );


/**
   \brief Constructs a list of molecular connections from an atom neighbor list

   \param b            Block data corresponding to the neighbor list
   \param atom_neighs  Atom neighbor list
   \param debug        If true, print some debug info to stderr

 */
neigh_list get_molecular_connections( const block_data &b,
                                      const neigh_list &atom_neighs,
                                      bool debug = false );



} // namespace neighborize

} // namespace lammps_tools

#endif // CLUSTER_FINDER_HPP
