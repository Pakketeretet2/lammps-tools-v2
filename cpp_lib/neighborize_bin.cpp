#include "neighborize_bin.hpp"

#include <cmath>
#include "my_timer.hpp"

namespace lammps_tools {

namespace neighborize {



int bin_neighborizer::xyz_index_to_bin_index( int ix, int iy, int iz )
{
	return ix + Nx*iy + Nx*Ny*iz;
}

void bin_neighborizer::bin_index_to_xyz_index( int bin, int &ix,
                                               int &iy, int &iz )
{
	ix = bin % Nx;
	iy = ( (bin - ix) % (Nx*Ny) ) / Nx;
	iz = (bin - ix - iy*Nx) / (Nx*Ny);
}

void bin_neighborizer::position_to_xyz_index( double x, double y, double z,
                                              int &ix, int &iy, int &iz )
{
	ix = ( x - xlo[0] ) / bin_size;
	iy = ( y - xlo[1] ) / bin_size;
	iz = ( z - xlo[2] ) / bin_size;

	// Adjust for out-of-box coordinates due to how e.g. LAMMPS rewraps.
	if     ( ix >= Nx ) ix -= Nx;
	else if( ix <   0 ) ix += Nx;
	if     ( iy >= Ny ) iy -= Ny;
	else if( iy <   0 ) iy += Ny;
	if     ( iz >= Nz ) iz -= Nz;
	else if( iz <   0 ) iz += Nz;

	if( ix >= Nx || iy >= Ny || iz >= Nz ||
	    ix < 0 || iy < 0 || iz < 0 ){
		std::cerr << "  Bin index for particle at ( " << x << ", " << y
		          << ", " << z << " ) invalid!\n";
		std::cerr << "  Bin index " << xyz_index_to_bin_index(ix,iy,iz)
		          << "; ( " << ix << ", " << iy << ", " << iz
		          << " ) not valid!\n";

		my_logic_error(__FILE__, __LINE__, "Incorrect bin index!" );
	}
}

int bin_neighborizer::position_to_bin_index( double x, double y, double z )
{
	int ix, iy, iz;
	position_to_xyz_index( x, y, z, ix, iy, iz );
	return ix + iy*Nx + iz*Nx*Ny;
}


int bin_neighborizer::shift_bin_index( int bin, int xinc, int yinc, int zinc )
{
	int ix = bin % Nx;
	int iy = ( (bin - ix) % (Nx*Ny) ) / Nx;
	int iz = (  bin - ix - iy*Nx) / (Nx*Ny);

	ix += xinc;
	iy += yinc;
	iz += zinc;

	if( periodic ){
		if     ( ix < 0 )    ix += Nx;
		else if( ix >=  Nx ) ix -= Nx;

		if     ( iy < 0 )    iy += Ny;
		else if( iy >= Ny )  iy -= Ny;

		if     ( iz < 0 )    iz += Nz;
		else if( iz >= Nz )  iz -= Nz;

		my_assert( __FILE__, __LINE__,
		           (ix >= 0) && (iy >= 0) && (iz >= 0) &&
		           (ix < Nx) && (iy < Ny) && (iz < Nz),
		           "Illegal bin index after correction!" );
	}

	return xyz_index_to_bin_index( ix, iy, iz );
}

void bin_neighborizer::setup_bins()
{
	if( atom_to_bin.size() > 0 ){
		atom_to_bin.clear();
	}
	atom_to_bin.resize( natoms );
	double pad = 0.1*rc;
	bin_size = rc + pad;

	Nx = std::floor( (xhi[0] - xlo[0]) / bin_size );
	Ny = std::floor( (xhi[1] - xlo[1]) / bin_size );
	Nz = std::floor( (xhi[2] - xlo[2]) / bin_size );

	if( Nx < 3 ){
		Nx = 3;
	}
	if( Ny < 3 ){
		Ny = 3;
	}
	if( (Nz < 3) && (dims == 3) ){
		Nx = 3;
	}
	Nbins = Nx*Ny*Nz;

	// Reset the bin size so that it matches exactly:
	double bin_size_x = (xhi[0] - xlo[0]) / Nx;
	double bin_size_y = (xhi[1] - xlo[1]) / Ny;
	double bin_size_z = (xhi[2] - xlo[2]) / Nz;
	double tol = 1e-8;
	my_assert( __FILE__, __LINE__, (bin_size_x - bin_size_y) < tol,
	           "Incorrect bin dimensions for x!" );
	my_assert( __FILE__, __LINE__, (bin_size_x - bin_size_z) < tol,
	           "Incorrect bin dimensions for y!" );
	my_assert( __FILE__, __LINE__, (bin_size_y - bin_size_z) < tol,
	           "Incorrect bin dimensions for z!" );
	bin_size = bin_size_x;

	if( !quiet ) std::cerr << "  ....Making " << Nbins << " bins...\n";
	try {
		bins.resize( Nbins );
	}catch( std::bad_alloc ){
		std::cerr << "  Failed to allocate " << Nbins
		          << " bins! Use NSQ binning instead!\n";
		throw;
	}
}

void bin_neighborizer::bin_atoms( particle_filter filt )
{
	if( !quiet ) std::cerr << "  ....Binning atoms...\n";
	for( std::size_t i = 0; i < id.size(); ++i ){
		if( !filt( b, i ) ) continue;

		int bin = position_to_bin_index( x[i], y[i], z[i] );
		bins[bin].push_back(i);
		atom_to_bin[i] = bin;
		if( !quiet ) std::cerr << "    ....Atom " << i
		                       << " in bin " << bin << "\n";
	}
}

// Adds to the neighbour list of i the atoms in bin that are in range.
void bin_neighborizer::add_bin_neighs( int i, const std::list<int> &bin,
                                       neigh_list &neighs )
{
	double xi[3];
	xi[0] = x[i];
	xi[1] = y[i];
	xi[2] = z[i];

	for( int j : bin ){
		if( id[i] >= id[j] ) continue;

		if( !(type[i] == itype || itype == 0) &&
		     (type[j] == jtype || jtype == 0) ){
			continue;
		}

		// Don't know if you have to check if i is already in j's
		// bin or not...

		double xj[3];
		double r[3];
		xj[0] = x[j];
		xj[1] = y[j];
		xj[2] = z[j];

		double r2 = b.dom.dist_2( xi, xj, r );
		if( r2 > rc2 ) continue;

		if( !quiet ) std::cerr << "    ....atoms " << i << " and "
		                       << j << " in bins " << atom_to_bin[i]
		                       << " and " << atom_to_bin[j]
		                       << " are neighbours!\n";

		neighs[i].push_back(j);
		neighs[j].push_back(i);
		n_neighs += 2;
	}
}

void bin_neighborizer::add_neighs_from_bin_2d( int i, neigh_list &neighs )
{
	int loop_idx[9];
	loop_idx[0] = atom_to_bin[i];
	loop_idx[ 1] = shift_bin_index(loop_idx[0],  1,  0, 0 );
	loop_idx[ 2] = shift_bin_index(loop_idx[0], -1,  0, 0 );
	loop_idx[ 3] = shift_bin_index(loop_idx[0],  1,  1, 0 );
	loop_idx[ 4] = shift_bin_index(loop_idx[0], -1,  1, 0 );
	loop_idx[ 5] = shift_bin_index(loop_idx[0],  0,  1, 0 );
	loop_idx[ 6] = shift_bin_index(loop_idx[0],  0, -1, 0 );
	loop_idx[ 7] = shift_bin_index(loop_idx[0],  1, -1, 0 );
	loop_idx[ 8] = shift_bin_index(loop_idx[0], -1, -1, 0 );

	for( int j = 0; j < 9; ++j ){
		int bin_index = loop_idx[j];
		if( bin_index < 0 || bin_index >= Nbins ){
			if( !periodic ){
				continue;
			}else{
				my_logic_error( __FILE__, __LINE__,
				                "Illegal bin index!" );
			}
		}
		const std::list<int> &bin = bins[bin_index];
		add_bin_neighs( i, bin, neighs );
	}
}

void bin_neighborizer::add_neighs_from_bin_3d( int i, neigh_list &neighs )
{
	int loop_idx[27];
	loop_idx[0] = atom_to_bin[i];
	loop_idx[ 1] = shift_bin_index(loop_idx[0],  1,  0, 0 );
	loop_idx[ 2] = shift_bin_index(loop_idx[0], -1,  0, 0 );
	loop_idx[ 3] = shift_bin_index(loop_idx[0],  1,  1, 0 );
	loop_idx[ 4] = shift_bin_index(loop_idx[0], -1,  1, 0 );
	loop_idx[ 5] = shift_bin_index(loop_idx[0],  0,  1, 0 );
	loop_idx[ 6] = shift_bin_index(loop_idx[0],  0, -1, 0 );
	loop_idx[ 7] = shift_bin_index(loop_idx[0],  1, -1, 0 );
	loop_idx[ 8] = shift_bin_index(loop_idx[0], -1, -1, 0 );

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

	for( int j = 0; j < 27; ++j ){
		int bin_index = loop_idx[j];
		if( bin_index < 0 || bin_index >= Nbins ){
			if( !periodic ){
				continue;
			}else{
				my_logic_error( __FILE__, __LINE__,
				                "Illegal bin index!" );
			}
		}
		const std::list<int> &bin = bins[bin_index];
		add_bin_neighs( i, bin, neighs );
	}
}


void bin_neighborizer::neigh_bin_atom( int i, neigh_list &neighs )
{
	if( dims == 2 ){
		add_neighs_from_bin_2d( i, neighs );
	}else{
		add_neighs_from_bin_3d( i, neighs );
	}
}




double bin_neighborizer::build( neigh_list &neighs,
                                const are_neighbours &criterion,
                                particle_filter filt,
                                int neigh_count_estimate )
{
	my_timer m;
	m.tic();
	if( !quiet ) std::cerr << "Building neigh list for "
	                       << id.size() << " particles.\n";
	neighs.resize( id.size() );
	for( std::vector<int> &ni : neighs ){
		ni.clear();
		if( neigh_count_estimate > 0 )
			ni.reserve( neigh_count_estimate );
	}
	n_neighs = 0;
	m.toc("Clearing neigh list");

	// 0. Allocates bins and atom_to_bin containers,
	//    sets number of bins and bin size
	m.tic();
	setup_bins();
	m.toc("Setting up bins");

	// 1. Put each atom that is not filtered in a bin.
	m.tic();
	bin_atoms( filt );
	m.toc("Binning atoms");

	// 2. Loop over all atoms ids, and add all atoms in their bin
	//    or adjacent bins to their neighbour list, if within rc.
	m.tic();
	for( std::size_t i = 0; i < id.size(); ++i ){
		// Skip i here too since it's not indexed in atom_to_bin:
		if( !filt( b, i ) ) continue;

		neigh_bin_atom( i, neighs );
	}
	m.toc("Neighborizing");


	if( !quiet ){
		std::cerr << "Done! Counted " << n_neighs << " neighbours.\n";
	}

	return static_cast<double>( n_neighs ) / natoms;
}


} //namespace neighborize

} //namespace lammps_tools
