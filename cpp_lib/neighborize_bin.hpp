#ifndef NEIGHBORIZE_BIN_HPP
#define NEIGHBORIZE_BIN_HPP

/**
   \file neighborize_bin.hpp
*/

#include <algorithm>
#include <list>
#include <vector>

#include "block_data.hpp"
#include "domain.hpp"
#include "my_assert.hpp"
#include "neighborize.hpp"


namespace lammps_tools {

namespace neighborize {

// For binning it is more convenient to use a class...
class neighborizer_bin : public neighborizer
{
public:

	neighborizer_bin( const block_data &b, const std::vector<int> &s1,
	                  const std::vector<int> &s2, int dims, double rc )
		: neighborizer( b, s1, s2, dims ), atom_to_bin(), bins(),
		  rc(rc), Nx(0), Ny(0), Nz(0), Nbins(0), bin_size(0.0), atoms_binned_( false )
	{}
	virtual ~neighborizer_bin(){}

	int get_atom_bin( int atom ) const;

	void setup_bins();
	void bin_atoms();

	int n_bins() const { return bins.size(); }

	bool atoms_binned() const { return atoms_binned_; }
	bool bins_setup() const { return !atom_to_bin.empty(); }

	// Some helper functions:
	int  xyz_index_to_bin_index( int ix, int iy, int iz ) const;
	int  position_to_bin_index( double x, double y, double z ) const;
	void bin_index_to_xyz_index( int bin, int &ix, int &iy, int &iz ) const;
	void position_to_xyz_index( double x, double y, double z,
	                            int &ix, int &iy, int &iz ) const;


	int  shift_bin_index( int bin, int xinc, int yinc, int zinc ) const;

	template<int dim> std::vector<int> get_nearby_bins( int bin ) const;
	std::vector<int> get_nearby_bins( int bin, int dim ) const;

	// These do the actual work:
	void add_bin_neighs( int i, const std::vector<int> &bin,
	                     neigh_list &neighs,
	                     int &n_neighs,
	                     const are_neighbours &criterion ) const;

	void neigh_bin_atom( int i, neigh_list &neighs,
	                     int &n_neighs,
	                     const are_neighbours &criterion ) const;

	void add_neighs_from_bin_2d( int i, neigh_list &neighs,
	                             int &n_neighs,
	                             const are_neighbours &criterion ) const;

	void add_neighs_from_bin_3d( int i, neigh_list &neighs,
	                             int &n_neighs,
	                             const are_neighbours &criterion ) const;

	const std::vector<int> &get_bin( int i ) const { return bins[i]; }


private:
	virtual int build( neigh_list &neighs,
	                   const are_neighbours &criterion );



	// members:
	std::vector<int> atom_to_bin;
	std::vector<std::vector<int> > bins;


	double rc;

	int Nx, Ny, Nz, Nbins;
	double bin_size;
	bool atoms_binned_;
};



template<int dim> inline
std::vector<int> neighborizer_bin::get_nearby_bins( int bin ) const
{
	if (dim == 2) {
		std::vector<int> loop_idx(9);
		loop_idx[0] = bin;
		loop_idx[ 1] = shift_bin_index(loop_idx[0],  1,  0, 0 );
		loop_idx[ 2] = shift_bin_index(loop_idx[0], -1,  0, 0 );
		loop_idx[ 3] = shift_bin_index(loop_idx[0],  1,  1, 0 );
		loop_idx[ 4] = shift_bin_index(loop_idx[0], -1,  1, 0 );
		loop_idx[ 5] = shift_bin_index(loop_idx[0],  0,  1, 0 );
		loop_idx[ 6] = shift_bin_index(loop_idx[0],  0, -1, 0 );
		loop_idx[ 7] = shift_bin_index(loop_idx[0],  1, -1, 0 );
		loop_idx[ 8] = shift_bin_index(loop_idx[0], -1, -1, 0 );
		return loop_idx;
	} else if (dim == 3) {
		std::vector<int> loop_idx(27);
		loop_idx[0] = bin;

		loop_idx[ 1] = shift_bin_index(loop_idx[0],  1,  0,  0 );
		loop_idx[ 2] = shift_bin_index(loop_idx[0], -1,  0,  0 );
		loop_idx[ 3] = shift_bin_index(loop_idx[0],  1,  1,  0 );
		loop_idx[ 4] = shift_bin_index(loop_idx[0], -1,  1,  0 );
		loop_idx[ 5] = shift_bin_index(loop_idx[0],  0,  1,  0 );
		loop_idx[ 6] = shift_bin_index(loop_idx[0],  0, -1,  0 );
		loop_idx[ 7] = shift_bin_index(loop_idx[0],  1, -1,  0 );
		loop_idx[ 8] = shift_bin_index(loop_idx[0], -1, -1,  0 );

		loop_idx[ 9] = shift_bin_index(loop_idx[0],  0,  0,  1 );
		loop_idx[10] = shift_bin_index(loop_idx[0],  1,  0,  1 );
		loop_idx[11] = shift_bin_index(loop_idx[0], -1,  0,  1 );
		loop_idx[12] = shift_bin_index(loop_idx[0],  1,  1,  1 );
		loop_idx[13] = shift_bin_index(loop_idx[0], -1,  1,  1 );
		loop_idx[14] = shift_bin_index(loop_idx[0],  0,  1,  1 );
		loop_idx[15] = shift_bin_index(loop_idx[0],  0, -1,  1 );
		loop_idx[16] = shift_bin_index(loop_idx[0],  1, -1,  1 );
		loop_idx[17] = shift_bin_index(loop_idx[0], -1, -1,  1 );

		loop_idx[18] = shift_bin_index(loop_idx[0],  0,  0, -1 );
		loop_idx[19] = shift_bin_index(loop_idx[0],  1,  0, -1 );
		loop_idx[20] = shift_bin_index(loop_idx[0], -1,  0, -1 );
		loop_idx[21] = shift_bin_index(loop_idx[0],  1,  1, -1 );
		loop_idx[22] = shift_bin_index(loop_idx[0], -1,  1, -1 );
		loop_idx[23] = shift_bin_index(loop_idx[0],  0,  1, -1 );
		loop_idx[24] = shift_bin_index(loop_idx[0],  0, -1, -1 );
		loop_idx[25] = shift_bin_index(loop_idx[0],  1, -1, -1 );
		loop_idx[26] = shift_bin_index(loop_idx[0], -1, -1, -1 );
		return loop_idx;
	}else{
		my_logic_error_terminate(__FILE__, __LINE__,
		                         "Dimension not implemented");
		return std::vector<int>(0);
	}
}


inline
std::vector<int> neighborizer_bin::get_nearby_bins( int bin, int dim ) const
{
	switch(dim){
		default:
			my_logic_error_terminate(__FILE__, __LINE__,
			                         "Dimension not implemented");
		case 2:
			return get_nearby_bins<2>(bin);

		case 3:
			return get_nearby_bins<3>(bin);
	}
}





} //namespace neighborize

} //namespace lammps_tools

#endif // NEIGHBORIZE_BIN_HPP
