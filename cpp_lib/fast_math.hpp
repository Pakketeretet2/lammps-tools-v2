#ifndef FAST_MATH_HPP
#define FAST_MATH_HPP

#include "sincos_lookup.hpp"

namespace lammps_tools {

namespace fast_math {

inline void sincos( double x, double &s, double &c )
{
	static lut_sin_cos table( 512 );
	table.sincos( x, s, c );
}

} // fast_math

} // lammps_tools


#endif // FAST_MATH_HPP
