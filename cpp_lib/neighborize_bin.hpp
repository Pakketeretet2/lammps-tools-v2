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
class bin_neighborizer
{
public:

	bin_neighborizer( const block_data &b,
	                  const std::vector<std::string> &fields,
	                  double rc, int dims, int itype, int jtype )
		: rc(rc), rc2(rc*rc), dims(dims), itype(itype), jtype(jtype),
		  periodic(b.dom.periodic),  b(b),
		  xlo({b.dom.xlo[0], b.dom.xlo[1], b.dom.xlo[2]}),
		  xhi({b.dom.xhi[0], b.dom.xhi[1], b.dom.xhi[2]}),
		  quiet(true), natoms(b.N), n_neighs(0),
		  Nx(0), Ny(0), Nz(0), Nbins(0), bin_size(0.0)

	{
		if( !grab_common_fields( b, fields, id, type, x, y, z ) ){
			my_logic_error( __FILE__, __LINE__,
			                "Failed to grab necessary fields!" );
		}
	}

	double rc, rc2;
	int dims;
	int itype, jtype;
	int periodic;

	const block_data &b;
	double xlo[3];
	double xhi[3];

	std::vector<int> id;
	std::vector<int> type;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

	bool quiet;

	double build( neigh_list &neighs,
	              const are_neighbours &criterion,
	              particle_filter filt = pass_all,
	              int neigh_count_estimate = 100 );

private:
	// Some helper functions:
	int  xyz_index_to_bin_index( int ix, int iy, int iz );
	int  position_to_bin_index( double x, double y, double z );
	void bin_index_to_xyz_index( int bin, int &ix, int &iy, int &iz );
	void position_to_xyz_index( double x, double y, double z,
	                            int &ix, int &iy, int &iz );


	int  shift_bin_index( int bin, int xinc, int yinc, int zinc );

	// These do the actual work:
	void setup_bins();
	void bin_atoms( particle_filter filt );
	void add_bin_neighs( int i, const std::list<int> &bin,
	                     neigh_list &neighs );
	void neigh_bin_atom( int i, neigh_list &neighs );

	void add_neighs_from_bin_2d( int i, neigh_list &neighs );
	void add_neighs_from_bin_3d( int i, neigh_list &neighs );

	// members:
	std::vector<int> atom_to_bin;
	std::vector<std::list<int> > bins;

	int natoms, n_neighs;
	int Nx, Ny, Nz, Nbins;
	double bin_size;
};

} //namespace neighborize

} //namespace lammps_tools

#endif // NEIGHBORIZE_BIN_HPP
