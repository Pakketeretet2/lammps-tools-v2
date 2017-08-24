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
   Finds which molecules are connected to which through given neighbour list.

   The information is stored as follows:
   The vector<int> at index mol_id contains the ids of all the molecules
   mol_id is connected to.

   \param b        Block data to use molecular info from
   \param neighs   Neighbour list to construct network from.
   \param conns    neigh_list that contains the connections between molecules.
   \param networks Contains a list of lists with networks, that is, all
                   molecule ids that are directly or indirectly connected.
*/
void find_molecular_networks( const block_data &b, const neigh_list &neighs,
                              neigh_list &conns,
                              neigh_list &networks );

neigh_list get_molecular_connections( const block_data &b,
                                      const neigh_list &neighs );

neigh_list make_molecular_networks( const block_data &b,
                                    const neigh_list &neighs,
                                    const neigh_list &conns );



} // namespace neighborize

} // namespace lammps_tools

#endif // CLUSTER_FINDER_HPP
