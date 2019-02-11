#include "block_data_access.hpp"
#include "constants.hpp"
#include "correlation.hpp"
#include "my_timer.hpp"
#include "neighborize.hpp"
#include "neighborize_bin.hpp"


#include <cmath>

using namespace lammps_tools;
using namespace correlate;

namespace lammps_tools {

namespace correlate {


template <typename T>
void correlate_impl( const lammps_tools::block_data &b,
                     const std::vector<T> &data,
                     std::vector<double> &Cr,
                     double x0, double x1, double dx, int dims )
{
	// x1 is the max, so you can bin the atoms
	std::cerr << "Correlating data between " << x0 << " and "
	          << x1 << "...\n";
	double Lgrid = x1 - x0;
	int n_bins = Lgrid / dx;

	Cr.resize( n_bins, 0.0 );
	std::vector<double> counts( n_bins );

	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);

	const double pi = constants::pi;

	my_timer timer( std::cerr );
	std::vector<int> a = neighborize::all(b);
	neighborize::neighborizer_bin nb( b, a, a, 2, x1 );
	nb.setup_bins();
	nb.bin_atoms();
	timer.toc("Binning atoms");


	timer.tic();
	std::cerr << "Correlating " << nb.n_bins() << " bins...\n";

	// Loop over the bins rather than the atoms.
	for( int bin_i = 0; bin_i < nb.n_bins(); ++bin_i ){
		std::vector<int> loop_idx = nb.get_nearby_bins<2>(bin_i);
		for( int bin_j : loop_idx ){
			if (bin_j < bin_i) continue;

			// Loop over all particles in bin_i:
			for( int i : nb.get_bin(bin_i) ){
				double xi[3] = { x[i], y[i], z[i] };
				T di = data[i];
				// You also need to take into account the
				// self-correlation in this case.

				for( int j : nb.get_bin(bin_j) ){
					T dj = data[j];
					double xj[3] = { x[j], y[j], z[j] };
					double rij[3];
					double r2 = b.dom.dist_2( xi, xj, rij );
					double r  = std::sqrt(r2);

					int bin = ( r - x0 ) / dx;
					if( bin >= n_bins || bin < 0 ) continue;

					// You will encounter i == j twice.
					double didj = di*dj * ( ( i == j ) ? 0.5 : 1.0 );
					Cr[bin] += didj;

					counts[bin] += 1.0;
				}
			}
		}
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
			Cr[bin] *= factor;
		}
	}

}



void correlate_int( const lammps_tools::block_data &b,
                    const std::vector<int> &data,
                    std::vector<double> &Cr,
                    double x0, double x1, double dx, int dims )
{
	correlate_impl<int>( b, data, Cr, x0, x1, dx, dims );
}

void correlate_double( const lammps_tools::block_data &b,
                       const std::vector<double> &data,
                       std::vector<double> &Cr,
                       double x0, double x1, double dx, int dims )
{
	correlate_impl<double>( b, data, Cr, x0, x1, dx, dims );
}




} // namespace correlate

} // namespace lammps_tools
