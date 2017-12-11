#ifndef ICOSAHEDRA_HPP
#define ICOSAHEDRA_HPP

/**
   \file icosahedra.hpp

   \brief This file contains some tools for analyzing DNA origami icosahedra.
*/
#include "block_data.hpp"
#include "geometry.hpp"
#include "neighborize.hpp"

#include <vector>

namespace lammps_tools {

namespace icosahedra {


std::vector<double> get_dihedral_angles( const block_data &b,
                                         const std::vector<int> &patches1,
                                         const std::vector<int> &patches2,
                                         const neighborize::neigh_list &conns );


neighborize::neigh_list make_mol_connections( const block_data &b,
                                              const std::vector<int> &patches1,
                                              const std::vector<int> &patches2,
                                              int method, int dims, double rc );

std::vector<point> get_triangle_axes( const block_data &b,
                                      const std::vector<int> &patches1,
                                      const std::vector<int> &patches2 );

std::vector<int> patch_types_to_ids( const block_data &b,
                                     const std::vector<int> &patches1,
                                     const std::vector<int> &patches2 );

double get_scattering_intensity( const block_data &b,
                                 const neighborize::neigh_list &networks,
                                 double I0 );


} // namespace icosahedra

} // namespace lammps_tools


#endif // ICOSAHEDRA_HPP
