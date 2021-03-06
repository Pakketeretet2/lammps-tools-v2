#include "block_data_access.hpp"
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

neighborizer::neighborizer( const block_data &b, const std::vector<int> &s1,
    const std::vector<int> &s2, int dims )
	: dims(dims), periodic(b.dom.periodic), b(b),
	  xlo{b.dom.xlo[0], b.dom.xlo[1], b.dom.xlo[2]},
	  xhi{b.dom.xhi[0], b.dom.xhi[1], b.dom.xhi[2]}, n_atoms(0),
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

	for( std::size_t i = 0; i < neighs.size(); ++i ){
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

		// You relied on there being garbage beyond mol[local_b.N-1]
		// before you idiot. XD
		// Second check makes sure you do not go beyond bounds.
		while( mol[j] == current_mol && j < local_b.N ) ++j;

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

	for( std::size_t i = 0; i < neighs.size(); ++i ){
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

	std::vector<int> s1, s2;
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

	return make_list_dist_indexed( neighs, b, s1, s2, method, dims, rc,
	                               mol_policy, bond_policy, quiet );
}

double make_list_dist_indexed( neigh_list &neighs,
                               const block_data &b,
                               const std::vector<int> &ilist,
                               const std::vector<int> &jlist,
                               int method, int dims,
                               double rc,
                               int mol_policy,
                               int bond_policy,
                               bool quiet )
{
	if( !quiet ) std::cerr << "  ....Calculating neighbours of "
	                       << ilist.size() << " atoms from "
	                       << jlist.size() << " atoms.\n";

	for( std::vector<int> &ni : neighs ){
		ni.clear();
	}
	neighs.resize( b.N );

	if( method == DIST_NSQ ){
		dist_criterion d( rc, dims );
		neighborizer_nsq n( b, ilist, jlist, dims );

		n.quiet = quiet;
		n.mol_policy = mol_policy;
		n.bond_policy = bond_policy;

		return n.build_list( neighs, d );
	}else if( method == DIST_BIN ){
		dist_criterion d( rc, dims );
		neighborizer_bin n( b, ilist, jlist, dims, rc );

		n.quiet = quiet;
		n.mol_policy = mol_policy;
		n.bond_policy = bond_policy;

		return n.build_list( neighs, d );
	}else{
		my_logic_error( __FILE__, __LINE__,
		                "Unknown dist neighbouring method!" );
	}

	double avg_neighs = 0.0;
	double norm = 0.0;
	id_map im( get_id(b) );
	for( int idx : ilist ){
		int i = im[idx];
		avg_neighs += neighs[i].size();

		norm += 1.0;
	}

	return avg_neighs / norm;

}


neigh_list nearest_neighs( const block_data &b,
                           int itype, int jtype, int method, int dims,
                           double rc,
                           int mol_policy, int bond_policy, bool quiet )
{
	neigh_list neighs;
	double avg = make_list_dist( neighs, b, itype, jtype, method, dims, rc,
	                             mol_policy, bond_policy, quiet );
	if( !quiet ) std::cerr << "Average # of neighs = " << avg << "\n";
	return neighs;
}



neigh_list nearest_neighs_indexed( const block_data &b,
                                   const std::vector<int> &ilist,
                                   const std::vector<int> &jlist,
                                   int method, int dims, double rc,
                                   int mol_policy, int bond_policy, bool quiet )
{
	neigh_list neighs;
	double avg = make_list_dist_indexed( neighs, b, ilist, jlist, method, dims, rc,
	                                     mol_policy, bond_policy, quiet );
	if( !quiet ) std::cerr << "Average # of neighs = " << avg << "\n";
	return neighs;
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
	for( std::size_t i = 0; i < neighs.size(); ++i ){
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
	double zi = z[i];
	double zj = z[j];
	if( dims == 2 ){
		zi = zj = 0.0;
	}
	double xi[3] = { x[i], y[i], zi };
	double xj[3] = { x[j], y[j], zj };
	double r[3];
	double r2 = b.dom.dist_2( xi, xj, r );

	if( r2 > rc2 ) return false;
	else           return true;
}


std::vector<bond> neigh_list_to_bonds( const block_data &b,
                                       const neighborize::neigh_list &neighs,
                                       int btype )
{

	std::vector<bond> bonds;
	std::size_t bc = 1;
	const std::vector<int> &id = get_id(b);
	for( std::size_t i = 0; i < neighs.size(); ++i ){
		for( int j : neighs[i] ){
			if( id[i] >= id[j] ) continue;

			bond bb;
			bb.id = bc;
			bb.type = btype;
			bb.particle1 = i;
			bb.particle2 = j;
			bonds.push_back( bb );
			++bc;
		}
	}
	my_assert( __FILE__, __LINE__, bc == bonds.size() + 1,
	           "Incorrect bond count!" );
	return bonds;
}


neigh_list neigh_list_to_network( const neigh_list &neighs, int min_idx )
{
	int N_parts = neighs.size();
	std::vector<bool> part_out( N_parts, false );
	neigh_list networks;
	for( std::size_t i = min_idx; i < part_out.size(); ++i ){
		if( !part_out[i] ){
			std::vector<int> network;
			part_out[i] = true;
			network.push_back( i );

			for( std::size_t it = 0; it < network.size(); ++it ){
				int j = network[it];
				const std::vector<int> &j_neighs = neighs[j];
				for( int o_part : j_neighs ){
					if( part_out[o_part] ) continue;
					part_out[o_part] = true;
					network.push_back(o_part);
				}
			}
			networks.push_back( network );
		}
	}
	return networks;
}


neigh_list get_empty_neigh_list()
{
	neigh_list nl;
	return nl;
}


std::vector<int> all( const lammps_tools::block_data &b )
{
	my_timer timer(std::cerr);
	std::vector<int> a( b.N );
	for( int i = 0; i < b.N; ++i ){
		a[i] = i;
	}
	return a;
}


} // namespace neighborize

} // namespace lammps_tools
