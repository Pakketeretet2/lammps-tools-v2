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
std::vector<double> rayleigh_gans( const class block_data &b,
                                   const std::vector<double> &qs );

/// Calculates Raygleigh-Gans scattering data from block_data for given ids
std::vector<double> rayleigh_gans( const class block_data &b,
                                   const std::vector<double> &qs,
                                   const std::vector<int> &ids );

// Ugly crap for pybind11:
std::vector<double> rayleigh_gans_( const class block_data &b,
                                    const std::vector<double> &qs,
                                    const std::vector<int> &ids )
{ return rayleigh_gans(b,qs,ids); }



} // namespace scatter

} // namespace lammps_tools


#endif // SCATTER_HPP
