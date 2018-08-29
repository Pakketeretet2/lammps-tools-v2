#include "block_data_access.hpp"
#include "constants.hpp"
#include "correlation.hpp"
#include "my_timer.hpp"
#include "neighborize.hpp"

#include <cmath>

using namespace lammps_tools;
using namespace correlate;

namespace lammps_tools {

namespace correlate {


template <typename T>
std::vector<double> correlate_impl( const lammps_tools::block_data &b,
                                    const std::vector<T> &data,
                                    double x0, double x1, double dx, int dims )
{
	// x1 is the max, so you can bin the atoms
	std::cerr << "Correlating data between " << x0 << " and "
	          << x1 << "...\n";
	my_timer timer( std::cerr );
	std::vector<int> all_vec = neighborize::all(b);
	/*
	neighborize::neigh_list neighs = nearest_neighs( b, 0, 0,
	                                                 neighborize::DIST_BIN,
	                                                 dims, x1 );
	*/
	int method = neighborize::DIST_BIN;
	neighborize::neigh_list neighs = neighborize::nearest_neighs(b, 0, 0,
	                                                             method,
	                                                             dims, x1);
	if( method == neighborize::DIST_NSQ ){
		timer.toc( "Building neighbor list (NSQ)" );
	}else{
		timer.toc( "Building neighbor list (BIN)" );
	}


	double Lgrid = x1 - x0;
	int n_bins = 1 + Lgrid / dx;

	std::vector<double> corr( n_bins );
	std::vector<double> counts( n_bins );

	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);

	const double pi = constants::pi;

	timer.tic();
	// For each atom, loop over its neighbors and
	// correlate the data between them.
	for( int i = 0; i < b.N; ++i ){
		double xi[3] = { x[i], y[i], z[i] };
		T di = data[i];
		// You also need to take into account the
		// self-correlation in this case.

		for( int j : neighs[i] ){
			if( i > j ) continue; // No double counting.
			T dj = data[j];

			double xj[3] = { x[j], y[j], z[j] };
			double rij[3];
			double r2 = b.dom.dist_2( xi, xj, rij );
			double r  = std::sqrt(r2);

			int bin = ( r - x0 ) / dx;
			if( bin >= n_bins || bin < 0 ) continue;

			double didj = di*dj;
			corr[bin] += didj;

			counts[bin] += 1.0;
		}
		// Count self-correlation too:
		corr[0]   += di*di;
		counts[0] += 1.0;
	}
	timer.toc( "Correlating data" );

	// Normalize:
	for( int bin = 0; bin < n_bins; ++bin ){
		double r = x0 + dx * bin;
		double bin_r = r + 0.5*dx;
		double bin_V = 2*pi*bin_r*bin_r;
		if( dims == 3 ) bin_V *= 2.0*bin_r/3.0;
		if( counts[bin] > 0 ){
			// double factor = 1.0 / (counts[bin]*bin_V);
			double factor = 1.0 / counts[bin];
			corr[bin] *= factor;
		}
	}

	return corr;
}



std::vector<double> correlate_int( const lammps_tools::block_data &b,
                                   const std::vector<int> &data,
                                   double x0, double x1, double dx, int dims )
{
	return correlate_impl<int>( b, data, x0, x1, dx, dims );
}

std::vector<double> correlate_double( const lammps_tools::block_data &b,
                                      const std::vector<double> &data,
                                      double x0, double x1, double dx, int dims )
{
	return correlate_impl<double>( b, data, x0, x1, dx, dims );
}




} // namespace correlate

} // namespace lammps_tools
