#include "block_data_access.hpp"
#include "density_distribution.hpp"
#include "my_timer.hpp"
#include "neighborize_bin.hpp"

#include <cmath>
#include <vector>


namespace lammps_tools {

namespace density {

std::vector<int> box_atoms( const lammps_tools::block_data &b,
                            int Nx, double &dx, int dims )
{
	//my_timer timer( std::cerr );
	double Lx = b.dom.xhi[0] - b.dom.xlo[0];
	double Ly = b.dom.xhi[1] - b.dom.xlo[1];
	double Lz = b.dom.xhi[2] - b.dom.xlo[2];

	// Make sure all lengths are the same.
	my_assert( __FILE__, __LINE__, std::fabs( Lx - Ly ) < 1e-8,
	           "Coarsening analysis only makes sense for square boxes!" );
	if( dims == 3 ){
		my_assert( __FILE__, __LINE__, std::fabs( Lx - Lz ) < 1e-8,
		           "Coarsening analysis only makes sense for square boxes!" );
		my_assert( __FILE__, __LINE__, std::fabs( Ly - Lz ) < 1e-8,
		           "Coarsening analysis only makes sense for square boxes!" );
	}

	dx = Lx / Nx;
	auto all_atoms = neighborize::all(b);
	neighborize::neighborizer_bin nb( b, all_atoms, all_atoms, dims, dx );

	nb.setup_bins();
	nb.bin_atoms();
	int N_bins = nb.n_bins();

	std::vector<int> boxes( b.N );

	for( int i = 0; i < b.N; ++i ){
		int box_idx = nb.get_atom_bin( i );
		if( box_idx < 0 || box_idx >= N_bins ){
			my_warning( __FILE__, __LINE__,
			            "Invalid box idx encountered!" );
		}
		boxes[i] = box_idx;
	}
	//timer.toc( "  Boxing atoms" );

	return boxes;
}



std::vector<int> box_count( const lammps_tools::block_data &b,
                            const std::vector<int> &boxes,
                            int N_boxes, int itype, bool quiet )
{
	//my_timer timer( std::cerr );
	// Boxes is atom-to-box.

	std::vector<int> boxc( N_boxes, 0 );
	const std::vector<int> &types = get_type(b);
	std::size_t icount = 0;
	for( int i = 0; i < b.N; ++i ){
		int box_idx = boxes[i];

		if( types[i] == itype ){
			boxc[box_idx]++;
			icount++;
		}
	}

	// Normalize?
	double box_avg = 0.0;
	for( std::size_t i = 0; i < boxc.size(); ++i ){
		box_avg += boxc[i];
	}
	//timer.toc( "  Box count" );

	return boxc;
}


std::vector<double> density_distribution( const lammps_tools::block_data &b,
                                          int Nx, int itype, int dims )
{
	double dx = 0.0;
	std::vector<int> boxes = box_atoms( b, Nx, dx, dims );
	std::size_t N_boxes = Nx*Nx;
	if( dims == 3 ) N_boxes *= Nx;

	std::vector<double> counts( N_boxes, 0 );

	//my_timer timer(std::cerr);
	//timer.tic();

	// Construct a vector that contains box occupation of itype.
	const std::vector<int> &types = get_type(b);
	int icount = 0;
	for( int i = 0; i < b.N; ++i ){
		if( (itype > 0) && (types[i] != itype) ) continue;
		icount++;
		counts[boxes[i]] += 1.0;
	}

	// You expect each box to occupy N / N_boxes of the total number.
	double n_expect = static_cast<double>(icount) / N_boxes;

	for( std::size_t i = 0; i < counts.size(); ++i ){
		counts[i] /= n_expect;
	}
	//timer.toc( "  Counting atoms per box" );
	return counts;
}


} // namespace density

} // namespace lammps_tools
