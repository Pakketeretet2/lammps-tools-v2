#ifndef ZIP_HPP
#define ZIP_HPP

/**
   \file zip.hpp Contains a custom zip iterator struct.

   Copied from
   https://stackoverflow.com/questions/8511035/sequence-zip-function-for-c11
*/

#include <algorithm>
#include <iterator>

namespace lammps_tools {

// Extend the util namespace:
namespace util {

/**
   Applies functor f to zip of containers c1 and c2.
*/
template <typename C1, typename C2, typename F>
void zip_map( const C1 &c1, const C2 &c2, F f )
{
	auto i1 = c1.begin();
	auto i2 = c2.begin();
	for( ; i1 != c1.end() && i2 != c2.end(); ++i1, ++i2 ){
		f( *i1, *i2 );
	}
}


} // namespace lammps_util

} // namespace lammps_tools

#endif // ZIP_HPP
