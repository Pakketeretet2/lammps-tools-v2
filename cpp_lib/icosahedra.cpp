#include <algorithm>

#include "block_data.hpp"
#include "block_data_access.hpp"
#include "cluster_finder.hpp"
#include "geometry.hpp"
#include "my_assert.hpp"
#include "my_timer.hpp"
#include "neighborize.hpp"

namespace lammps_tools {

namespace icosahedra {

// Calculates the normalized orientation axis of the triangle.
point triangle_orient_axis( point x11, point x12, point x21, point x22 )
{
	point b01 = x12 - x11;
	point b12 = x22 - x21;

	point ax = util::cross( b01, b12 );
	double n2 = util::dot( ax, ax );
	return ax /= sqrt(n2);
}


// Warning: Assumes the axes are already normalized!
double dihedral_angle( const point &ax1, const point &ax2 )
{
	double cost = util::dot( ax1, ax2 );
	double angle = std::acos( cost );
	return angle;
}

std::vector<int> patch_types_to_ids( const block_data &b,
                                     const std::vector<int> &patches1,
                                     const std::vector<int> &patches2 )
{
	std::vector<int> indices1;
	std::vector<int> indices2;
	const std::vector<int> &types = get_type(b);

	for( std::size_t i = 0; i < b.N; ++i ){
		for( int ti : patches1 ){
			if( types[i] == ti ){
				indices1.push_back(i);
				break;
			}
		}
		for( int tj : patches2 ){
			if( types[i] == tj ){
				indices2.push_back(i);
				break;
			}
		}
	}

	std::vector<int> filter_ids;
	const std::vector<int> &id = get_id(b);
	for( std::size_t idx : indices1 ){
		util::insert_sorted( filter_ids, id[idx] );
	}
	for( std::size_t idx : indices2 ){
		// indices1 are unique, but 2 might have duplicates.
		util::insert_sorted_unique( filter_ids, id[idx] );
	}

	return filter_ids;
}

// Helper function, filters everything but the patch types from neigh lists.
block_data filter_block_data( const block_data &b,
                              const std::vector<int> &patches1,
                              const std::vector<int> &patches2 )
{

	std::vector<int> filter_ids = patch_types_to_ids( b, patches1,
	                                                  patches2 );
	return filter_block( b, filter_ids );
}


neighborize::neigh_list make_mol_connections( const block_data &b,
                                              const std::vector<int> &patches1,
                                              const std::vector<int> &patches2,
                                              int method, int dims, double rc )
{
	const std::vector<int> &filter_types = get_type(b);
	neighborize::neigh_list ns;
	double avg_ns = neighborize::make_list_dist( ns, b, 0, 0,
	                                             method, dims, rc );
	neighborize::neigh_list conns =
		neighborize::get_molecular_connections( b, ns );

	return conns;
}

std::vector<point> get_triangle_axes( const block_data &b,
                                      const std::vector<int> &patches1,
                                      const std::vector<int> &patches2 )
{
	const std::vector<int> &mol = get_mol(b);
	std::size_t N_mols = mol.size() + 1; // mol goes from 0 to N_mols-1
	std::vector<point> mol_axes( N_mols );


	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);
	const std::vector<int> &types = get_type(b);

	// First construct an angle table for each molecule!
	// Grab all points:
	auto idx_to_point = [&x, &y, &z]( std::size_t idx ) -> point {
		return point( x[idx], y[idx], z[idx] ); };

	for( std::size_t mol_i = 1; mol_i < N_mols; ++mol_i ){
		// Find for this molecule the points x11, x12, x21 and x22.
		point x11, x12, x21, x22;
		int pts_out = 0;
		for( std::size_t i = 0; i < b.N; ++i ){
			if( mol[i] == mol_i ){
				if( types[i] == patches1[0] ){
					x11 = idx_to_point( i );
					pts_out += 1;
				}
				if( types[i] == patches1[1] ){
					x21 = idx_to_point( i );
					pts_out += 2;
				}
				if( types[i] == patches2[0] ){
					x12 = idx_to_point( i );
					pts_out += 4;
				}
				if( types[i] == patches2[1] ){
					x22 = idx_to_point( i );
					pts_out += 8;
				}
			}
			if( pts_out == 15 ) break;
		}
		// Now you have x11, x12, x21 and x22 for this molecule.
		mol_axes[mol_i] = triangle_orient_axis( x11, x12, x21, x22 );
	}
	return mol_axes;
}

double get_scattering_intensity( const block_data &b,
                                 const neighborize::neigh_list &conns,
                                 double I0 )
{
	neighborize::neigh_list mol_network =
		neighborize::neigh_list_to_network( conns, 1 );
	double I = 0.0;
	for( const std::vector<int> &nw : mol_network ){
		double s2 = nw.size();
		I += I0 * s2;
	}
	return I;
}


neighborize::neigh_list get_mol_connections( const block_data &b,
                                             const std::vector<int> &patches1,
                                             const std::vector<int> &patches2,
                                             int method, int dims, double rc )
{
	my_timer timer(std::cerr);
	timer.tic();
	neighborize::neigh_list conns = make_mol_connections( b, patches1,
	                                                      patches2, method,
	                                                      dims, rc );
	timer.toc( "    finding molecule connections." );
	return conns;
}



std::vector<double> get_dihedral_angles( const lammps_tools::block_data &b,
                                         const std::vector<int> &patches1,
                                         const std::vector<int> &patches2,
                                         const neighborize::neigh_list &conns )
{

	using namespace neighborize;

	my_assert( __FILE__, __LINE__, patches1.size() == 3,
	           "patches vectors did not contain three types! ");
	my_assert( __FILE__, __LINE__, patches2.size() == 3,
	           "patches vectors did not contain three types! ");

	std::size_t N_mols = conns.size();
	std::vector<point> mol_axes = get_triangle_axes( b, patches1,
	                                                 patches2 );
	std::vector<double> dihedral_angles;
	// Then loop over all molecule connections:
	for( std::size_t mol_i = 1; mol_i < N_mols; ++mol_i ){
		const point &ax1 = mol_axes[mol_i];
		const std::vector<int> &ci = conns[mol_i];
		for( std::size_t mol_j : ci ){
			const point &ax2 = mol_axes[mol_j];
			double angle = dihedral_angle( ax1, ax2 );
			dihedral_angles.push_back( angle );
		}
	}

	return dihedral_angles;

}





} // namespace icosahedra

} // namespace lammps_tools
