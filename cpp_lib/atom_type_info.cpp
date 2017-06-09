#include "atom_type_info.hpp"

namespace lammps_tools {

void swap( atom_type_info &f, atom_type_info &s )
{
	using std::swap;
	swap( s, f );
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
	mass.resize(Ntypes);
	set_defaults();
}

void atom_type_info::set_defaults()
{
	for( double &m : mass ){
		m = 1.0;
	}
}


} // namespace lammps_tools
