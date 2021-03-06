#include "neighborize_bin.hpp"

#include <cmath>
#include "my_timer.hpp"

#define USE_OPENMP 1

#ifdef USE_OPENMP
#include <omp.h>
constexpr const bool use_openmp = true;
#else
constexpr const bool use_openmp = false;
// Just to make stuff compile:
int omp_get_thread_num(){ return 0; }
#endif


namespace lammps_tools {

namespace neighborize {


int neighborizer_bin::xyz_index_to_bin_index( int ix, int iy, int iz ) const
{
	return ix + Nx*iy + Nx*Ny*iz;
}

void neighborizer_bin::bin_index_to_xyz_index( int bin, int &ix,
                                               int &iy, int &iz ) const
{
	ix = bin % Nx;
	iy = ( (bin - ix) % (Nx*Ny) ) / Nx;
	if( dims == 3 ){
		iz = (bin - ix - iy*Nx) / (Nx*Ny);
	}else{
		iz = 0;
	}
}

void neighborizer_bin::position_to_xyz_index( double x, double y, double z,
                                              int &ix, int &iy, int &iz ) const
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

	if( dims == 2 ) iz = 0;

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

int neighborizer_bin::position_to_bin_index( double x, double y, double z ) const
{
	int ix, iy, iz;
	position_to_xyz_index( x, y, z, ix, iy, iz );
	if( dims == 3 ){
		return ix + iy*Nx + iz*Nx*Ny;
	}else{
		return ix + iy*Nx;
	}
}


int neighborizer_bin::shift_bin_index( int bin, int xinc, int yinc, int zinc ) const
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
	if( dims == 2 ) iz = 0;

	return xyz_index_to_bin_index( ix, iy, iz );
}

void neighborizer_bin::setup_bins()
{
	if( atom_to_bin.size() > 0 ){
		atom_to_bin.clear();
	}
	atom_to_bin.resize( b.N );
	bin_size = rc;

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
		Nz = 3;
	}else if( dims == 2 ){
		Nz = 1;
	}

	Nbins = Nx*Ny*Nz;

	// Reset the bin size so that it matches exactly:
	double bin_size_x = (xhi[0] - xlo[0]) / Nx;
	double bin_size_y = (xhi[1] - xlo[1]) / Ny;
	double bin_size_z = (xhi[2] - xlo[2]) / Nz;

	bin_size = ( bin_size_x > bin_size_y ) ? bin_size_x : bin_size_y;
	bin_size = ( bin_size_z > bin_size )   ? bin_size_z : bin_size;


	if( !quiet ) std::cerr << "  ....Bin size = " << bin_size << ".\n";
	if( !quiet ) std::cerr << "  ....Making " << Nbins << " bins...\n";
	try {
		bins.resize( Nbins );
	}catch( std::bad_alloc & ){
		std::cerr << "  Failed to allocate " << Nbins
		          << " bins! Use NSQ binning instead!\n";
		throw;
	}

}

void neighborizer_bin::bin_atoms()
{
	if( atoms_binned_ ){
		// This is a no-op.
		if( !quiet ) std::cerr << "  ....Atoms already binned!\n";
		return;
	}

	if( !quiet ) std::cerr << "  ....Binning "
	                       << s2.size() << " atoms...\n";

	const std::vector<double> &x = data_as<double>(
		b.get_special_field( block_data::X ) );
	const std::vector<double> &y = data_as<double>(
		b.get_special_field( block_data::Y ) );
	const std::vector<double> &z = data_as<double>(
		b.get_special_field( block_data::Z ) );
	// std::cerr << "x[3] = " << x[3] << ", " << y[3] << ", " << z[3] << "\n";

	// Bin only the atoms in the second group. This way, you can
	//    later loop over the first group, get their bin on-the-fly,
	//    and automatically you only compare against atoms from the
	//    second group. If they are both all, it reduces to the
	//    standard approach anyway.
	if( use_openmp ){
		std::vector<int> s2_vec( s2.size() );
		int c = 0;
		for( int j : s2 ){
			s2_vec[c++] = j;
		}
		for( std::size_t i = 0; i < s2_vec.size(); ++i ){
			int j = s2_vec[i];
			if( !quiet ){
				std::cerr << "Binning atom " << j << " at ( "
				          << x[j] << ", " << y[j] << ", "
				          << z[j] << " ), ";
			}
			int bin = position_to_bin_index( x[j], y[j], z[j] );
			if( !quiet ){
				int ix, iy, iz;
				bin_index_to_xyz_index( bin, ix, iy, iz );
				std::cerr << "bin = " << ix << ", " << iy
				          << ", " << iz << "\n";
			}

			bins[bin].push_back( j );
			atom_to_bin[j] = bin;
		}
	}else{
		for( int j : s2 ){

			if( !quiet ){
				std::cerr << "Binning atom " << j << " at ( "
				          << x[j] << ", " << y[j] << ", "
				          << z[j] << " ), ";
			}
			int bin = position_to_bin_index( x[j], y[j], z[j] );
			if( !quiet ){
				int ix, iy, iz;
				bin_index_to_xyz_index( bin, ix, iy, iz );
				std::cerr << "bin = " << ix << ", " << iy
				          << ", " << iz << "\n";
			}

			bins[bin].push_back( j );
			atom_to_bin[j] = bin;
		}
	}

	atoms_binned_ = true;
}

// Adds to the neighbour list of i the atoms in bin that are in range.
void neighborizer_bin::add_bin_neighs( int i, const std::vector<int> &bin,
                                       neigh_list &neighs,
                                       int &n_neighs,
                                       const are_neighbours &criterion ) const
{
	const std::vector<int> &id = data_as<int>(
		b.get_special_field( block_data::ID ) );

	for( int j : bin ){
		if( id[i] == id[j] ) continue;

		// Don't know if you have to check if i is already in j's
		// bin or not...
		if( !criterion( b, i, j ) ) continue;

		if( !quiet ) std::cerr << "    ....atoms " << i << " and "
		                       << j << " in bins " << atom_to_bin[i]
		                       << " are neighbours!\n";

		neighs[i].push_back( j );
		neighs[j].push_back( i );
		n_neighs += 2;
	}
}


void neighborizer_bin::add_neighs_from_bin_2d( int i, neigh_list &neighs,
                                               int &n_neighs,
                                               const are_neighbours &criterion ) const
{
	const std::vector<double> &x = data_as<double>(
		b.get_special_field( block_data::X ) );
	const std::vector<double> &y = data_as<double>(
		b.get_special_field( block_data::Y ) );
	const std::vector<double> &z = data_as<double>(
		b.get_special_field( block_data::Z ) );
	const std::vector<int> &id = data_as<int>(
		b.get_special_field( block_data::ID ) );
	double xi = x[i];
	double yi = y[i];
	double zi = z[i];

	int my_bin = position_to_bin_index( xi, yi, zi );
	std::vector<int> loop_idx = get_nearby_bins( my_bin, 2 );

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
		const std::vector<int> &bin = bins[bin_index];
		if( !quiet ){
			std::cerr << "  ....Checking bins of atom " << i
			          << " ( id = " << id[i] << ", bin "
			          << bin_index << " )\n";
		}
		add_bin_neighs( i, bin, neighs, n_neighs, criterion );
	}
}


void neighborizer_bin::add_neighs_from_bin_3d( int i, neigh_list &neighs,
                                               int &n_neighs,
                                               const are_neighbours &criterion ) const
{
	const std::vector<double> &x = data_as<double>(
		b.get_special_field( block_data::X ) );
	const std::vector<double> &y = data_as<double>(
		b.get_special_field( block_data::Y ) );
	const std::vector<double> &z = data_as<double>(
		b.get_special_field( block_data::Z ) );
	const std::vector<int> &id = data_as<int>(
		b.get_special_field( block_data::ID ) );
	double xi = x[i];
	double yi = y[i];
	double zi = z[i];
	int my_bin = position_to_bin_index( xi, yi, zi );
	std::vector<int> loop_idx = get_nearby_bins<3>(my_bin);

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
		const std::vector<int> &bin = bins[bin_index];
		if( !quiet ){
			int ix, iy, iz;
			bin_index_to_xyz_index( bin_index, ix, iy, iz );
			std::cerr << "  ....Checking bins of atom " << i
			          << " ( id = " << id[i] << ", bin "
			          << ix << ", " << iy << ", " << iz << " )\n";
		}
		add_bin_neighs( i, bin, neighs, n_neighs, criterion );
	}
}


void neighborizer_bin::neigh_bin_atom( int i, neigh_list &neighs,
                                       int &n_neighs,
                                       const are_neighbours &criterion ) const
{
	if( dims == 2 ){
		add_neighs_from_bin_2d( i, neighs, n_neighs, criterion );
	}else{
		add_neighs_from_bin_3d( i, neighs, n_neighs, criterion );
	}
}




int neighborizer_bin::build( neigh_list &neighs,
                             const are_neighbours &criterion )
{
	my_timer m(std::cerr);
	if( !quiet ) m.tic();

	for( std::vector<int> &ni : neighs ){
		ni.clear();
	}
	int n_neighs = 0;
	if( !quiet ) m.toc("  Clearing neigh list");

	// 0. Allocates bins and atom_to_bin containers,
	//    sets number of bins and bin size
	if( !quiet ) m.tic();
	setup_bins();
	if( !quiet ) m.toc("  Setting up bins");

	// 1. Put each atom that is not filtered in a bin.
	if( !quiet ) m.tic();
	bin_atoms();
	if( !quiet ) m.toc("  Binning atoms");

	// 2. Loop over atoms in first group, get their bin on-the-fly,
	//    and neighbourize the atoms in that and adjacent bins, which
	//    are only those of group 2. :D
	if( !quiet ) std::cerr << "  ....Looping over "
	                       << s1.size() << " atoms...\n";
	if( !quiet ) m.tic();
	if( use_openmp ){
		std::vector<int> s1_vec( s1.size() );
		int c = 0;
		for( int j : s1 ){
			s1_vec[c++] = j;
		}
		for( std::size_t j = 0; j < s1_vec.size(); ++j ){
			int i = s1_vec[j];
			neigh_bin_atom( i, neighs, n_neighs, criterion );
		}
	}else{
		for( int i : s1 ){
			neigh_bin_atom( i, neighs, n_neighs, criterion );
		}
	}
	if( !quiet ) m.toc("  Neighborizing");


	return n_neighs;
}


int neighborizer_bin::get_atom_bin( int atom ) const
{
	if( !bins_setup() ){
		my_warning( __FILE__, __LINE__,
		            "Bins not built yet! Call setup_bins() first!" );
		return -1;
	}

	return( atom_to_bin[atom] );
}





} //namespace neighborize

} //namespace lammps_tools
