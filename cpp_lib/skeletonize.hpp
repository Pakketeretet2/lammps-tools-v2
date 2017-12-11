#ifndef SKELETONIZE_HPP
#define SKELETONIZE_HPP

/**
   \file skeletonize.hpp

   Contains routines for performing skeletonization analysis.
*/

#include <vector>
#include <list>
#include <functional>

#include "neighborize.hpp"

namespace lammps_tools {

namespace skeletonize {

struct ribbon_data
{
	double R;

	std::size_t n_largest;
	double avg_r, var_r, rmin, rmax, max_avg, max_var;

	double area;
};

std::vector<double> get_insideness( const class block_data &b,
                                    const neighborize::neigh_list &neighs );


std::vector<double> euclidian_distance_transform( const class block_data &b,
                                                  const std::vector<double> &insideness,
                                                  double R );

std::vector<double> get_ribbon_widths( const block_data &b,
                                       const neighborize::neigh_list &neighs,
                                       const std::vector<double> &edt,
                                       const std::vector<double> &insideness );

void get_edge( const block_data &b, const std::vector<double> &insideness,
               const std::list<int> &edges );


void skeletonize_edt( const class block_data &b, std::vector<int> &skeleton,
                      const std::vector<double> &insideness, double R );

void skeletonize( const block_data &b, std::vector<int> &skeleton,
                  const std::vector<double> &distance_map, double R );


void get_local_maxima( const std::vector<std::vector<int> > &neighs,
                       const std::vector<double> &field,
                       std::vector<int> &max_indices );

void get_local_maxima( const std::vector<std::vector<int> > & neighs,
                       const std::vector<double> &field,
                       std::vector<int> &max_indices,
                       std::function< double(double) > );

void get_ribbon_data( class block_data &b, ribbon_data &r_data );

std::vector<double> neighbor_strain( class block_data &b, double r0,
                                     int itype, int jtype, int method,
                                     int dims, double rc );


} // namespace skeletonize

} // namespace lammps_tools



#endif // SKELETONIZE_HPP
