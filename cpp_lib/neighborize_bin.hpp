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
		  rc(rc), n_neighs(0), Nx(0), Ny(0), Nz(0),
		  Nbins(0), bin_size(0.0)
	{}

	~neighborizer_bin(){}


private:
	virtual int build( neigh_list &neighs,
	                   const are_neighbours &criterion );

	// Some helper functions:
	int  xyz_index_to_bin_index( int ix, int iy, int iz );
	int  position_to_bin_index( double x, double y, double z );
	void bin_index_to_xyz_index( int bin, int &ix, int &iy, int &iz );
	void position_to_xyz_index( double x, double y, double z,
	                            int &ix, int &iy, int &iz );


	int  shift_bin_index( int bin, int xinc, int yinc, int zinc );

	// These do the actual work:

	void setup_bins();
	void bin_atoms();

	void add_bin_neighs( int i, const std::vector<int> &bin,
	                     neigh_list &neighs,
	                     const are_neighbours &criterion );

	void neigh_bin_atom( int i, neigh_list &neighs,
	                     const are_neighbours &criterion );

	void add_neighs_from_bin_2d( int i, neigh_list &neighs,
	                             const are_neighbours &criterion );

	void add_neighs_from_bin_3d( int i, neigh_list &neighs,
	                             const are_neighbours &criterion );

	// members:
	std::vector<int> atom_to_bin;
	std::vector<std::vector<int> > bins;


	double rc;

	int n_neighs;
	int Nx, Ny, Nz, Nbins;
	double bin_size;
};

} //namespace neighborize

} //namespace lammps_tools

#endif // NEIGHBORIZE_BIN_HPP
