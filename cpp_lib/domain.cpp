#include "domain.hpp"

#include <algorithm>

using namespace lammps_tools;

namespace lammps_tools {

void swap( domain &f, domain &s )
{
	using std::swap;

	swap( f.xlo, s.xlo );
	swap( f.xhi, s.xhi );
	
	swap( f.periodic, s.periodic );
}


domain::domain( const domain &o )
	: periodic( o.periodic )
{
	std::copy( o.xlo, o.xlo + 3, xlo );
	std::copy( o.xhi, o.xhi + 3, xhi );
	
}

} // namespace lammps_tools
