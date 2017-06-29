#ifndef BOND_ORDER_HPP
#define BOND_ORDER_HPP

/**
   \file bond_order.hpp

   Routines for calculating bond order parameters.
*/
#include <vector>

#include "block_data.hpp"
#include "domain.hpp"
#include "geometry.hpp"
#include "neighborize.hpp"
//#include "topology.hpp"

namespace lammps_tools {

namespace order_parameters {

/**
   Calculates bond order parameter psi_n for each particle and the average.

   \param[in]  b        The block data to analyse.
   \param[in]  neighs   The neighbour list that defines all bonded particles.
   \param[in]  n        The order in the bond order.
   \param[in]  axis     The axis relative to which the angles are calculated.
   \param[out] psi_n    Will contain the psi_n for each particle.

   \returns the average psi_n.
*/
double compute_psi_n( const block_data &b,
                      const neighborize::neigh_list &neighs,
                      int n, const point &axis,
                      std::vector<double> &psi_n_real,
	              std::vector<double> &psi_n_imag );

/**
   Calculates all bond angles for all neighbour pairs relative to some axis.

   \param[in]  b           The block data to analyse.
   \param[in]  neighs      The neighbour list that defines all bonded particles.
   \param[in]  axis        The reference axis.
   \param[out] angles      Will contain angles between bonds and axis.
   \param[out] bonds       Will contain the bonds corresponding to angles.
*/
void relative_bond_angles( const block_data &b,
                           const neighborize::neigh_list &neighs,
                           point axis,
                           const std::vector<bond> &bonds,
                           std::vector<double> &angles );

/**
   \brief Calculates angle between two unit vectors.

   \warning This assumes the vectors are unit vectors!
*/
double angle_2pi( const point &v1, const point &v2 );


} // order_parameters

} // lammps_tools



#endif // BOND_ORDER_HPP
