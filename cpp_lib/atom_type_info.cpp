#include "atom_type_info.hpp"

namespace lammps_tools {

void swap( atom_type_info &f, atom_type_info &s )
{
	using std::swap;
	swap(f.n_types,s.n_types);
}

} // namespace lammps_tools
