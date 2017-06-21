#include "atom_type_info.hpp"
#include <iostream>

namespace lammps_tools {

void swap( atom_type_info &f, atom_type_info &s )
{
	using std::swap;
	swap( s.mass, f.mass );
}


atom_type_info::atom_type_info() : mass()
{
	set_defaults();
}

atom_type_info::atom_type_info( int Ntypes ) : mass(Ntypes + 1)
{
	set_defaults();
}

void atom_type_info::set_size( int Ntypes )
{
	int old_n_types = mass.size();
	mass.resize(Ntypes+1);

	for( int i = old_n_types; i <= Ntypes; ++i ){
		mass[i] = 1.0;
	}
}

void atom_type_info::set_defaults()
{
	for( double &m : mass ){
		m = 1.0;
	}
}


} // namespace lammps_tools
