#include "histogram.hpp"
#include "fourier.hpp"
#include <cmath>

namespace lammps_tools {

namespace histogram {

template <typename T>
std::vector<double> make_histogram( const std::vector<T> &data,
                                    const T &y0, const T &y1, int N_bins )
{
	std::vector<double> distr( N_bins + 1 );
	double dx = (y1 - y0) / N_bins;

	double normalize = 1.0 / data.size();
	for( const T &t : data ){
		int bin = ( t - y0 ) / dx;

		if( bin < 0 || bin > N_bins ) continue;

		distr[bin]++;
	}

	for( double &v : distr ){
		v *= normalize;
	}
	return distr;
}

// Declare specializations:
template std::vector<double> make_histogram<int>( const std::vector<int> &data,
                                                  const int &y0, const int &y1, int N_bins );
template std::vector<double> make_histogram<double>( const std::vector<double> &data,
                                                     const double &y0, const double &y1, int N_bins );



std::vector<double> make_histogram_double( const std::vector<double> &data,
                                           double y0, double y1, int N_bins )
{
	return make_histogram<double>( data, y0, y1, N_bins );
}

std::vector<double> make_histogram_int( const std::vector<int> &data,
                                        int y0, int y1, int N_bins )
{
	return make_histogram<int>( data, y0, y1, N_bins );
}


std::vector<double> make_radial_grid( const std::vector<double> &xgrid,
                                      const std::vector<double> &ygrid,
                                      int Nbins )
{
	std::size_t Nx = xgrid.size();
	std::size_t Ny = ygrid.size();

	double rmax = std::sqrt( xgrid[0]*xgrid[0] + ygrid[0]*ygrid[0] );
	double dx = rmax / Nbins;
	std::vector<double> grid(Nbins);
	for( int bin = 0; bin < Nbins; ++bin ){
		double r = bin*dx;
		grid[bin] = r;
	}
	return grid;
}


template <typename T>
std::vector<T> radial_avg( const std::vector<double> &xgrid,
                           const std::vector<double> &ygrid,
                           const std::vector<T> &y, int Nbins )
{
	std::size_t Nx = xgrid.size();
	std::size_t Ny = ygrid.size();

	double rmax = std::sqrt( xgrid[0]*xgrid[0] + ygrid[0]*ygrid[0] );
	double dx = rmax / Nbins;
	std::vector<T> avg( Nbins, 0 );
	std::vector<double> counts( Nbins, 0 );

	for( std::size_t ix = 0; ix < Nx; ++ix ){
		double xi = xgrid[ix];
		for( std::size_t iy = 0; iy < Ny; ++iy ){
			double yi = ygrid[iy];
			double rij = std::sqrt( xi*xi + yi*yi );
			int bin = rij / dx;
			if( bin >= Nbins || bin < 0 ) continue;

			avg[bin] += y[ix*Nx + iy];
			counts[bin] += 1.0;
		}
	}

	for( int bin = 0; bin < Nbins; ++bin ){
		if( counts[bin] > 0 ){
			avg[bin] /= counts[bin];
		}
	}
	return avg;
}



std::vector<double> make_radial2_grid( const std::vector<double> &xgrid,
                                       const std::vector<double> &ygrid,
                                       int Nbins )
{
	std::size_t Nx = xgrid.size();
	std::size_t Ny = ygrid.size();

	double rmax2 = xgrid[0]*xgrid[0] + ygrid[0]*ygrid[0];
	double dx2 = rmax2 / Nbins;
	std::vector<double> grid2(Nbins);
	for( int bin = 0; bin < Nbins; ++bin ){
		double r2 = bin*dx2;
		grid2[bin] = r2;
	}
	return grid2;
}



template <typename T>
std::vector<T> radial2_avg( const std::vector<double> &xgrid,
                            const std::vector<double> &ygrid,
                            const std::vector<T> &y, int Nbins )
{
	std::size_t Nx = xgrid.size();
	std::size_t Ny = ygrid.size();

	double rmax2 = xgrid[0]*xgrid[0] + ygrid[0]*ygrid[0];
	double dx2 = rmax2 / Nbins;
	std::vector<T> avg( Nbins, 0 );
	std::vector<double> counts( Nbins, 0 );

	for( std::size_t ix = 0; ix < Nx; ++ix ){
		double xi = xgrid[ix];
		for( std::size_t iy = 0; iy < Ny; ++iy ){
			double yi = ygrid[iy];
			double rij2 = xi*xi + yi*yi;
			int bin = rij2 / dx2;
			if( bin >= Nbins || bin < 0 ) continue;

			avg[bin] += y[ix*Nx + iy];
			counts[bin] += 1.0;
		}
	}

	for( int bin = 0; bin < Nbins; ++bin ){
		if( counts[bin] > 0 ){
			avg[bin] /= counts[bin];
		}
	}
	return avg;
}



template std::vector<double> radial_avg<double>(
	const std::vector<double> &xgrid,
	const std::vector<double> &ygrid,
	const std::vector<double> &y, int Nbins );



template std::vector<double> radial2_avg<double>(
	const std::vector<double> &xgrid,
	const std::vector<double> &ygrid,
	const std::vector<double> &y, int Nbins );


// Instantiate for cx_double:
template std::vector<fourier::cx_double> radial_avg<fourier::cx_double>(
	const std::vector<double> &xgrid,
	const std::vector<double> &ygrid,
	const std::vector<fourier::cx_double> &y, int Nbins );


template std::vector<fourier::cx_double> radial2_avg<fourier::cx_double>(
	const std::vector<double> &xgrid,
	const std::vector<double> &ygrid,
	const std::vector<fourier::cx_double> &y, int Nbins );



std::vector<double> radial_avg( const std::vector<double> &xgrid,
                                const std::vector<double> &ygrid,
                                const std::vector<int> &y, int Nbins )
{
	std::vector<double> y_as_double( y.begin(), y.end() );
	return radial_avg<double>( xgrid, ygrid, y_as_double, Nbins );
}


std::vector<double> radial2_avg( const std::vector<double> &xgrid,
                                 const std::vector<double> &ygrid,
                                 const std::vector<int> &y, int Nbins )
{
	std::vector<double> y_as_double( y.begin(), y.end() );
	return radial2_avg<double>( xgrid, ygrid, y_as_double, Nbins );
}




} // namespace histogram

} // namespace lammps_tools
