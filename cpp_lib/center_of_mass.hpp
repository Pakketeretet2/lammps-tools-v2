#ifndef CENTER_OF_MASS_HPP
#define CENTER_OF_MASS_HPP

/**
   \file center_of_mass.hpp

   \brief Contains routines that calculate centers of mass and stuff.
*/

#include "types.hpp"


namespace lammps_tools {

class block_data;

point center_of_mass( const block_data &b );
point geometric_center( const block_data &b );

template <typename iter>
point center_of_mass( const block_data &b, iter it, iter end );

template <typename iter>
point geometric_center( const block_data &b, iter it, iter end );


template <typename iter, typename functor>
point weighted_pos_avg( const block_data &b, iter it, iter end, functor w );



} // namespace lammps_tools


#endif // CENTER_OF_MASS_HPP
