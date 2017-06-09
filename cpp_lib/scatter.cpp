#include "block_data.hpp"
#include "block_data_access.hpp"
#include "scatter.hpp"

#include <iostream>
#include <vector>

namespace lammps_tools {

namespace scatter {

void rayleigh_gans( const class block_data &b )
{
	rayleigh_gans( b, get_id( b ) );
}

void rayleigh_gans( const class block_data &b, const std::vector<int> &ids )
{
	std::cerr << "Calculating raygleigh_gans scattering for "
	          << ids.size() << " particles.\n";
	std::vector<double> radius( b.N, 1.0 );

}


} // namespace scatter

} // namespace lammps_tools
