#include "block_data.hpp"
#include "block_data_access.hpp"
#include "coarsening.hpp"
#include "density_distribution.hpp"
#include "my_timer.hpp"
#include "neighborize_bin.hpp"

#include "random_generator.hpp"

#include <cmath>
#include <vector>



using namespace lammps_tools;
using namespace coarse;

namespace lammps_tools {

namespace coarse {




std::vector<int> box_to_psi( const lammps_tools::block_data &b,
                             const std::vector<int> &boxes,
                             int N_boxes, int seed, bool quiet )
{
	//my_timer timer( std::cerr );
	// Boxes is atom-to-box.
	std::vector<int> psi( N_boxes, 0 );
	std::vector<int> na = density::box_count( b, boxes, N_boxes, 1, quiet );
	std::vector<int> nb = density::box_count( b, boxes, N_boxes, 2, quiet );

	for( int i = 0; i < N_boxes; ++i ){
		if( na[i] > nb[i] ){
			psi[i] = 1;
		}else if( na[i] < nb[i] ){
			psi[i] = -1;
		}else{
			psi[i] = 0;
		}
	}

	//timer.toc( "  Box to psi" );

	return psi;
}


template <typename T>
std::vector<T> map_to_atom( const lammps_tools::block_data &b,
			    const std::vector<int> &boxes,
			    const std::vector<T> &data )
{
	//my_timer timer( std::cerr );

	std::vector<T> data_per_atom( b.N, 0 );
	for( int i = 0; i < b.N; ++i ){
		int box_i = boxes[i];
		data_per_atom[i] = data[box_i];
	}

	//timer.toc( "  Data to data-per-atom" );

	return data_per_atom;
}



std::vector<double> map_to_atom_double( const lammps_tools::block_data &b,
                                        const std::vector<int> &boxes,
                                        const std::vector<double> &data )
{
	return map_to_atom<double>( b, boxes, data );
}


std::vector<int> map_to_atom_int( const lammps_tools::block_data &b,
                                  const std::vector<int> &boxes,
                                  const std::vector<int> &data )
{
	return map_to_atom<int>( b, boxes, data );
}

std::vector<cx_double> map_to_atom_cx_double( const lammps_tools::block_data &b,
                                              const std::vector<int> &boxes,
                                              const std::vector<cx_double> &data )
{
	return map_to_atom<std::complex<double> >( b, boxes, data );
}


std::vector<cx_int> map_to_atom_cx_int( const lammps_tools::block_data &b,
                                        const std::vector<int> &boxes,
                                        const std::vector<cx_int> &data )
{
	return map_to_atom<std::complex<int> >( b, boxes, data );
}






} // namespace coarse

} // namespace lammps_tools
