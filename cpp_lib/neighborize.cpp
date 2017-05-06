#include "data_field.hpp"
#include "neighborize.hpp"
#include "neighborize_bin.hpp"
#include "my_assert.hpp"
#include "util.hpp"

#include <algorithm>
#include <cmath>
#include <list>
#include <stdexcept>
#include <vector>

namespace lammps_tools {

namespace neighborize {


double neigh_dist_nsq( neigh_list &neighs,
                       const block_data &b,
                       const std::vector<std::string> &fields,
                       int itype, int jtype, int dims,
                       double rc, particle_filter filt, int neigh_est )
{
	std::vector<int> id, type;
	std::vector<double> x, y, z;
	bool success = grab_common_fields( b, fields, id, type, x, y, z );
	my_assert( __FILE__, __LINE__, success,
	           "Failed to assign required data from block!" );

	double rc2 = rc*rc;
	int N = id.size();
	neighs.resize(N);

	int total_neighbours = 0;
	for( std::vector<int> &ni : neighs ){
		ni.clear();
		if( neigh_est > 0 ) ni.reserve( neigh_est );
	}

	for( int i = 0; i < N; ++i ){
		if( !filt( b, i ) ) continue;
		int idi = id[i];

		double xi[3] = { x[i], y[i], z[i] };

		for( int j = 0; j < N; ++j ){
			if( !filt( b, j ) ) continue;

			int idj = id[j];
			if( idi >= idj ) continue;

			if( !(type[i] == itype || itype == 0) &&
			     (type[j] == jtype || jtype == 0) ){
				continue;
			}

			double xj[3] = { x[j], y[j], z[j] };

			double r[3];
			double r2 = b.dom.dist_2( xi, xj, r );
			if( r2 > rc2 ) continue;

			neighs[i].push_back(j);
			neighs[j].push_back(i);
			total_neighbours += 2;
		}
	}
	return static_cast<double>( total_neighbours ) / N;
}



double neigh_dist_bin( neigh_list &neighs,
                       const block_data &b,
                       const std::vector<std::string> &fields,
                       int itype, int jtype, int dims,
                       double rc, particle_filter filt, int neigh_est )
{
	bin_neighborizer n( b, fields, rc, dims, itype, jtype );
	return n.build( neighs, filt, neigh_est );
}


double make_list_dist( neigh_list &neighs,
                       const block_data &b,
                       const std::vector<std::string> &fields,
                       int itype, int jtype, int method, int dims, double rc,
                       particle_filter filt, int neigh_est )
{
	for( std::vector<int> &ni : neighs ){
		ni.clear();
	}
	if( method == DIST_NSQ ){
		return neigh_dist_nsq( neighs, b, fields, itype, jtype,
		                       dims, rc, filt, neigh_est );
	}else if( method == DIST_BIN ){
		return neigh_dist_bin( neighs, b, fields, itype, jtype,
		                       dims, rc, filt, neigh_est );
	}else{
		my_logic_error( __FILE__, __LINE__,
		                "Unknown dist neighbouring method!" );
	}
}

void verify_unique( int i, int j, const neigh_list &neighs )
{
	if( !util::is_unique(neighs[i], j) ){
		std::cerr << "neighs[i] = \n ";
		for( int k : neighs[i] ){
			std::cerr << " " << k;
		}
		std::cerr << "\nneighs[j] = \n ";
		for( int k : neighs[j] ){
			std::cerr << " " << k;
		}
		std::cerr << "\n";
		my_logic_error( __FILE__, __LINE__,
		                "j not unique neigh of i!" );
	}
}

void verify_contains( int i, int j, const neigh_list &neighs )
{
	if( !util::contains( neighs[j], i ) ){
		my_logic_error( __FILE__, __LINE__,
		                "j was neigh of i but not vice versa!" );
	}
}


void verify_neigh_list( const neigh_list &neighs )
{
	for( std::size_t i = 0; i < neighs.size(); ++i ){
		for( int j : neighs[i] ){
			verify_unique( i, j, neighs );
			verify_contains( i, j, neighs );

		}
	}
}





} // namespace lammps_tools

} // namespace neighborize
