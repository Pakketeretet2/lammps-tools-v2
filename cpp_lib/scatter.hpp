#ifndef SCATTER_HPP
#define SCATTER_HPP

/**
   \file scatter.hpp

   \brief Contains routines that generate simulated scattering information.
*/

#include <vector>

namespace lammps_tools {

namespace scatter {

/// Calculates Raygleigh-Gans scattering data from block_data.
void rayleigh_gans( const class block_data &b );
/// Calculates Raygleigh-Gans scattering data from block_data for given ids
void rayleigh_gans( const class block_data &b, const std::vector<int> &ids );

// Ugly crap for pybind11:
void rayleigh_gans_( const class block_data &b, const std::vector<int> &ids )
{ rayleigh_gans(b,ids); }



} // namespace scatter

} // namespace lammps_tools


#endif // SCATTER_HPP
