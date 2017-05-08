#include "data_field.hpp"
#include "neighborize.hpp"
#include "neighborize_nsq.hpp"
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



double neigh_dist_bin( neigh_list &neighs,
                       const block_data &b,
                       const std::vector<std::string> &fields,
                       int itype, int jtype, int dims,
                       double rc, particle_filter filt, int neigh_est )
{
	dist_criterion d( rc );
	bin_neighborizer n( b, fields, rc, dims, itype, jtype );
	return n.build( neighs, d, filt, neigh_est );
}


double neigh_dist_nsq( neigh_list &neighs,
                       const block_data &b,
                       const std::vector<std::string> &fields,
                       int itype, int jtype, int dims,
                       double rc, particle_filter filt, int neigh_est )
{
	dist_criterion d( rc );
	nsq_neighborizer n( b, fields, dims, itype, jtype );
	return n.build( neighs, d, filt, neigh_est );
}


double make_list_dist( neigh_list &neighs,
                       const block_data &b,
                       const std::vector<std::string> &fields,
                       int itype, int jtype, int method, int dims,
                       double rc,
                       bool include_mol,
                       bool include_bonds,
                       particle_filter filt,
                       int neigh_est )
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



// Some functors:
bool dist_criterion::operator()( const block_data &b, int i, int j ) const
{
	using dfd = data_field_double;
	const data_field *d_x, *d_y, *d_z;
	
	d_x = static_cast<const dfd*>( b.get_special_field( block_data::X ) );
	d_y = static_cast<const dfd*>( b.get_special_field( block_data::Y ) );
	d_z = static_cast<const dfd*>( b.get_special_field( block_data::Z ) );

	my_assert( __FILE__, __LINE__, d_x && d_y && d_z,
	           "Failed to grab data for x, y and z!" );
	
	const std::vector<double> &x = data_as<double>( d_x );
	const std::vector<double> &y = data_as<double>( d_y );
	const std::vector<double> &z = data_as<double>( d_z );

	if( b.atom_style == block_data::MOLECULAR ){
		const std::vector<int> &mol = data_as<int>(
			b.get_data( b.get_special_field_name(
				            block_data::MOL ) ) );
		if( mol[i] == mol[j] ) return true;
	}

	double xi[3] = { x[i], y[i], z[i] };
	double xj[3] = { x[j], y[j], z[j] };
	double r[3];
	double r2 = b.dom.dist_2( xi, xj, r );

	if( r2 > rc2 ) return false;
	else           return true;
}



} // namespace lammps_tools

} // namespace neighborize
