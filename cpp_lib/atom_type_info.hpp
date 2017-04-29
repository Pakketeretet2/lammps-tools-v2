#ifndef ATOM_TYPE_INFO_HPP
#define ATOM_TYPE_INFO_HPP

/**
   \file atom_type_info.hpp

   Definitions of per-atom-type information.
*/

#include <algorithm>

namespace lammps_tools {

/**
   Contains per-atom-type info like mass, pair coeffs, etc.
*/
struct atom_type_info
{
	atom_type_info( int n_types ) : n_types( n_types ){}
	~atom_type_info() {}

	int n_types;

	friend void swap( atom_type_info &f, atom_type_info &s );
};

}// namespace lammps_tools

#endif // ATOM_TYPE_INFO_HPP
