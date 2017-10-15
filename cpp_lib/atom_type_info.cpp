#include "atom_type_info.hpp"
#include "my_assert.hpp"

#include <iostream>

namespace lammps_tools {

void swap( atom_type_info &f, atom_type_info &s )
{
	using std::swap;
	swap( s.mass, f.mass );
	swap( s.type_names, f.type_names );
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

	if( type_names.empty() ){
		type_names.resize( Ntypes+1 );
	}else{
		std::size_t type_name_size = type_names.size();
		std::size_t n_expect = Ntypes + 1;
		if( type_names.size() != n_expect ){
			std::string warning = "type_names mismatches number of types! ";
			warning += "Resizing invalidates type_names!";
			std::cerr << "type_names is " << type_names.size()
			          << " big but expected " << n_expect << "\n";
			my_warning_if( __FILE__, __LINE__,
			               type_name_size != n_expect, warning );
			type_names.resize( Ntypes+1 );
		}else{
			// Matching size, don't do anything.
			// This is what should happen.
		}
	}

}

void atom_type_info::set_defaults()
{
	for( double &m : mass ){
		m = 1.0;
	}
}




} // namespace lammps_tools
