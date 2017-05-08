#include "neighborize_nsq.hpp"

namespace lammps_tools {

namespace neighborize {

double nsq_neighborizer::build( neigh_list &neighs,
                                const are_neighbours &criterion,
                                particle_filter filt,
                                int neigh_est )
{
	int N = b.N;
	
	neighs.resize(N);
	int total_neighbours = 0;
	for( std::vector<int> &ni : neighs ){
		ni.clear();
		if( neigh_est > 0 ) ni.reserve( neigh_est );
	}

	for( int i = 0; i < N; ++i ){
		if( !filt( b, i ) ) continue;
		int idi = id[i];
		for( int j = 0; j < N; ++j ){
			if( !filt( b, j ) ) continue;

			int idj = id[j];
			if( idi >= idj ) continue;

			if( !(type[i] == itype || itype == 0) &&
			     (type[j] == jtype || jtype == 0) ){
				continue;
			}

			if( !(criterion( b, i, j ) ) ) continue;

			neighs[i].push_back(j);
			neighs[j].push_back(i);
			total_neighbours += 2;
		}
	}
	return static_cast<double>( total_neighbours ) / N;
}

	
} // namespace neighborize

} // namespace lammps_tools

