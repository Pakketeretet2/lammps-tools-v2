#include "data_field.hpp"
#include "neighborize.hpp"
#include "neighborize_nsq.hpp"
#include "neighborize_bin.hpp"
#include "my_assert.hpp"
#include "my_timer.hpp"
#include "util.hpp"

#include <algorithm>
#include <cmath>
#include <list>
#include <stdexcept>
#include <vector>

namespace lammps_tools {

namespace neighborize {

neighborizer::neighborizer( const block_data &b, const std::list<int> &s1,
                            const std::list<int> &s2, int dims )
	: dims(dims), periodic(b.dom.periodic), b(b),
	  xlo{b.dom.xlo[0], b.dom.xlo[1], b.dom.xlo[2]},
	  xhi{b.dom.xhi[0], b.dom.xhi[1], b.dom.xhi[2]},
	  quiet(true), mol_policy(IGNORE), bond_policy(IGNORE), s1(s1), s2(s2)
{ }


double neighborizer::build_list( neigh_list &neighs,
                                 const are_neighbours &criterion )
{
	double avg_neighs = 0.0;
	build( neighs, criterion );


	switch( mol_policy ){
		default:
		case IGNORE:
			break;
		case INCLUDE:
			append_particles_in_mol( neighs );
			break;
		case EXCLUDE:
			remove_particles_in_mol( neighs );
			break;
	}
	switch( bond_policy ){
		default:
		case IGNORE:
			break;
		case INCLUDE:
			my_runtime_error( __FILE__, __LINE__,
			                  "Feature not implemented yet." );
			break;
		case EXCLUDE:
			my_runtime_error( __FILE__, __LINE__,
			                  "Feature not implemented yet." );
			break;
	}

	remove_doubles( neighs );

	for( int i = 0; i < neighs.size(); ++i ){
		avg_neighs += neighs[i].size();
	}
	if( !quiet ){
		std::cerr << "Found " << avg_neighs << " neighs in total for ";
	}
	avg_neighs /= b.N;
	if( !quiet ) std::cerr << avg_neighs << " on average.\n";

	return avg_neighs;
}


int neighborizer::append_particles_in_mol( neigh_list &neighs )
{
	const data_field *mol_ptr = b.get_special_field( block_data::MOL );
	my_assert( __FILE__, __LINE__, mol_ptr,
	           "Appending particles by mol requires special field "
	           "MOL to be set in block_data!" );
	if( !quiet ) std::cerr << "  ....Appending particles in mol...\n";

	// This can be much more efficient...
	// 1. Copy block data to local copy.
	// 2. Sort along molecule.
	// 3. Walk through molecule and done.
	const std::vector<int> &id = data_as<int>(
		b.get_special_field( block_data::ID ) );
	id_map im( id );

	my_timer t(std::cerr);
	if( !quiet ) t.tic();

	block_data local_b( b );
	local_b.sort_along( mol_ptr->name );

	if( !quiet ) t.toc( "Copying and sorting block_data" );

	const std::vector<int> &local_id = data_as<int>(
		local_b.get_special_field( block_data::ID ) );

	int total_neighbours = 0;
	const std::vector<int> &mol = data_as<int>(
		local_b.get_special_field( block_data::MOL ) );

	int i = 0;
	if( !quiet ) t.tic();
	while( i < local_b.N ){
		int current_mol = mol[i];
		int j = i+1;
		my_assert( __FILE__, __LINE__, mol[j] - current_mol <= 1,
		           "Molecule reordering not correct" );
		my_assert( __FILE__, __LINE__, mol[j] - current_mol >= 0,
		           "Molecule reordering not correct" );

		while( mol[j] == current_mol ) ++j;

		// Now the range i .. j-1 marks the current mol.
		for( int k = i; k < j; ++k ){
			for( int l = i; l < k; ++l ){
				// Remember to get the indices the way
				// they were in the original block:
				int idk = local_id[k];
				int idl = local_id[l];
				int idx_k = im[idk];
				int idx_l = im[idl];
				if( !util::contains( neighs[idx_k], idx_l ) ){
					neighs[idx_k].push_back( idx_l );
					neighs[idx_l].push_back( idx_k );
					total_neighbours += 2;
				}
			}
		}
		i = j;
	}
	if( !quiet ) t.toc( "Traversing sorted block_data" );

	return total_neighbours;
}



int neighborizer::remove_particles_in_mol( neigh_list &neighs )
{
	int removed = 0;
	const data_field *mol_ptr = b.get_special_field( block_data::MOL );
	my_assert( __FILE__, __LINE__, mol_ptr,
	           "Appending particles by mol requires special field "
	           "MOL to be set in block_data!" );

	const std::vector<int> &mol = data_as<int>( mol_ptr );

	for( int i = 0; i < neighs.size(); ++i ){
		std::vector<int> &ni = neighs[i];
		int old_size = ni.size();
		auto deleter = [i,mol](int j){ return mol[j] == mol[i]; };
		ni.erase( std::remove_if( ni.begin(), ni.end(), deleter ), ni.end() );
		int new_size = ni.size();
		removed += old_size - new_size;
	}
	return removed;
}





int neighborizer::append_bonded_particles( neigh_list &neighs )
{
	// Do nothing for now.
	return 0;
}




//****************** Stand-alone functions: *******************//


double make_list_dist( neigh_list &neighs,
                       const block_data &b,
                       int itype, int jtype, int method, int dims,
                       double rc,
                       int mol_policy,
                       int bond_policy,
                       bool quiet )
{

	std::list<int> s1, s2;
	const std::vector<int> &type = data_as<int>(
		b.get_special_field( block_data::TYPE ) );

	if( !quiet ) std::cerr << "  ....Filtering out " << b.N
	                       << " atoms looking for types " << itype
	                       << " and " << jtype << ".\n";
	for( int i = 0; i < b.N; ++i ){
		if( !itype || (type[i] == itype)){
			s1.push_back( i );
		}
		if( (type[i] == jtype) || (jtype == 0) ){
			s2.push_back( i );
		}
	}

	if( !quiet ) std::cerr << "  ....Calculating neighbours of " << s1.size()
	                       << " atoms from " << s2.size() << " atoms.\n";

	for( std::vector<int> &ni : neighs ){
		ni.clear();
	}
	neighs.resize( b.N );

	if( method == DIST_NSQ ){
		dist_criterion d( rc );
		neighborizer_nsq n( b, s1, s2, dims );

		n.quiet = quiet;
		n.mol_policy = mol_policy;
		n.bond_policy = bond_policy;

		return n.build_list( neighs, d );
	}else if( method == DIST_BIN ){
		dist_criterion d( rc );
		neighborizer_bin n( b, s1, s2, dims, rc );

		n.quiet = quiet;
		n.mol_policy = mol_policy;
		n.bond_policy = bond_policy;

		return n.build_list( neighs, d );
	}else{
		my_logic_error( __FILE__, __LINE__,
		                "Unknown dist neighbouring method!" );
	}

	return 0.0;
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


void remove_doubles( neigh_list &neighs )
{
	for( int i = 0; i < neighs.size(); ++i ){
		util::remove_doubles( neighs[i] );
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
	           "Failed to grab data for x1, y1 and z1!" );

	const std::vector<double> &x = data_as<double>( d_x );
	const std::vector<double> &y = data_as<double>( d_y );
	const std::vector<double> &z = data_as<double>( d_z );

	double xi[3] = { x[i], y[i], z[i] };
	double xj[3] = { x[j], y[j], z[j] };
	double r[3];
	double r2 = b.dom.dist_2( xi, xj, r );

	if( r2 > rc2 ) return false;
	else           return true;
}




} // namespace lammps_tools

} // namespace neighborize
