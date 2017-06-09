#ifndef TYPES_HPP
#define TYPES_HPP

#include <algorithm>

#include <cstdint>

/*!
  \file  types.h
  \brief This file contains type definitions.

  \ingroup cpp_lib
*/

namespace lammps_tools {

typedef int64_t       bigint;    ///< signed long of guaranteed size (64 bits)



/**
   \brief Represents a point in 3D.
*/
struct point
{
	point() : x(0), y(0), z(0) {}
	explicit point(double *v) : x(v[0]), y(v[1]), z(v[2]) {}

	point &operator=( double *v )
	{
		using std::swap;
		point n( v );
		swap( *this, n );
		return *this;
	}

	double x, y, z;
};



}

#endif // TYPES_HPP
