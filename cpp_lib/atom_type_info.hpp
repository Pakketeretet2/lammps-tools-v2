#ifndef ATOM_TYPE_INFO_HPP
#define ATOM_TYPE_INFO_HPP

/**
   \file atom_type_info.hpp

   Definitions of per-atom-type information.
*/

#include <algorithm>
#include <string>
#include <vector>

namespace lammps_tools {

/**
   Contains per-atom-type info like mass, pair coeffs, etc.
*/
struct atom_type_info
{
	std::vector<double> mass;
	std::vector<std::string> type_names;

	atom_type_info();
	explicit atom_type_info( int n_types );

	void set_size( int Ntypes );
	void set_defaults();
};

void swap( atom_type_info &f, atom_type_info &s );


}// namespace lammps_tools

#endif // ATOM_TYPE_INFO_HPP
