#include "neighborize_nsq.hpp"

namespace lammps_tools {

namespace neighborize {

int neighborizer_nsq::build( neigh_list &neighs,
                             const are_neighbours &criterion )
{
	int total_neighbours = 0;
	for( std::vector<int> &ni : neighs ){
		ni.clear();
	}

	const std::vector<int> &id = data_as<int>(
		b.get_special_field( block_data::ID ) );
	if( !quiet ) std::cerr << "Looping over " << s1.size() << " x "
	                       << s2.size() << " atoms.\n";
	for( int i : s1 ){
		int idi = id[i];
		for( int j : s2 ){

			int idj = id[j];
			if( idi == idj ) continue;

			if( !(criterion( b, i, j ) ) ) continue;
			neighs[i].push_back( j );
			neighs[j].push_back( i );
			total_neighbours += 2;
		}
	}
	return total_neighbours;
}


} // namespace neighborize

} // namespace lammps_tools
