#include <cmath>

#include "block_data_access.hpp"
#include "constants.hpp"
#include "rdf.hpp"


namespace lammps_tools {

namespace neighborize {

using lammps_tools::constants::pi;
void compute_rdf( const block_data &b, int Nbins, double r0, double r1,
                  int dims, int itype, int jtype,
                  std::vector<double> &rdf, std::vector<double> &coord )
{
	neigh_list neighs;
	make_list_dist( neighs, b, itype, jtype, DIST_BIN, dims, r1 );
	compute_rdf_with_neighs( b, Nbins, r0, r1, dims, neighs, rdf, coord );
}

void compute_rdf_with_neighs( const block_data &b, int Nbins,
                              double r0, double r1, int dims,
                              const neigh_list &neighs,
                              std::vector<double> &rdf,
                              std::vector<double> &coord )
{
	rdf.resize( Nbins );
	coord.resize( Nbins );

	double dr = ( r1 - r0 ) / ( Nbins - 1 );
	for( int i = 0; i < Nbins; ++i ){
		rdf[i] = coord[i] = 0.0;
	}

	// TODO: Think of normalisation.
	double total_neighs = 0.0;
	for( std::size_t i = 0; i < neighs.size(); ++i ){
		total_neighs += neighs[i].size();
	}
	double avg_neighs = total_neighs / b.N;


	const std::vector<int> &ids = get_id( b );
	const std::vector<double> &x = get_x( b );
	const std::vector<double> &y = get_y( b );
	const std::vector<double> &z = get_z( b );

	double rc2 = r1*r1;
	int adds = 0;

	for( int i = 0; i < neighs.size(); ++i ){
		const std::vector<int> &ni = neighs[i];
		int idi = ids[i];
		for( int j : ni ){
			int idj = ids[j];
			double xi[3] = { x[i], y[i], z[i] };
			double xj[3] = { x[j], y[j], z[j] };

			double r[3];
			double r2 = b.dom.dist_2( xi, xj, r );
			if( r2 > rc2 ) continue;

			double rr = std::sqrt( r2 );
			int bini = ( rr - r0 ) / dr;
			if( bini < 0 || bini >= Nbins ) continue;

			rdf[bini] += 1.0;
			++adds;
		}
	}

	double Lx = b.dom.xhi[0] - b.dom.xlo[0];
	double Ly = b.dom.xhi[1] - b.dom.xlo[1];
	double Lz = b.dom.xhi[2] - b.dom.xlo[2];

	double V = Lx*Ly;
	if( dims != 2 ) V *= Lz;

	double volume_scale;
	if( dims != 2 ){
		volume_scale = 4.0*pi / (3.0 * V);
	}else{
		volume_scale = pi / V;
	}
	volume_scale *= b.N;
	coord[0] = 0.0;
	// V is already set, it is total volume for 3D or area for 2D.
	for( int bin = 0; bin < Nbins; ++bin ){
		double bin_r   = r0 + dr*bin;
		double bin_r_p = bin_r + dr;
		double bin_r3 = bin_r*bin_r*bin_r;
		double bin_r_p3 = bin_r_p*bin_r_p*bin_r_p;
		double n_ideal = volume_scale*( bin_r_p3 - bin_r3 );

		rdf[bin] /= n_ideal;

		if( bin == 0 ){
			coord[bin] = 0.0;
		}else{
			coord[bin] = coord[bin-1] + rdf[bin]*n_ideal;
		}
	}
	std::ofstream out("rdf_test.dat");
	for( int bin = 0; bin < Nbins; ++bin ){
		out << bin << " " << r0 + dr*bin << " " << rdf[bin]
		    << " " << coord[bin] << "\n";
	}
}


std::vector<double> rdf( const block_data &b, int Nbins, double r0, double r1,
                         int dims, int itype, int jtype )
{
	std::vector<double> rrdf( Nbins );

	double dr = ( r1 - r0 ) / ( Nbins - 1 );
	for( int i = 0; i < Nbins; ++i ){
		rrdf[i] = 0.0;
	}

	neigh_list neighs;
	make_list_dist( neighs, b, itype, jtype, DIST_BIN, dims, r1 );

	// TODO: Think of normalisation.
	double total_neighs = 0.0;
	for( std::size_t i = 0; i < neighs.size(); ++i ){
		total_neighs += neighs[i].size();
	}
	double avg_neighs = total_neighs / b.N;


	const std::vector<int> &ids = get_id( b );
	const std::vector<double> &x = get_x( b );
	const std::vector<double> &y = get_y( b );
	const std::vector<double> &z = get_z( b );

	double rc2 = r1*r1;
	int adds = 0;

	for( int i = 0; i < neighs.size(); ++i ){
		const std::vector<int> &ni = neighs[i];
		int idi = ids[i];
		for( int j : ni ){
			int idj = ids[j];
			double xi[3] = { x[i], y[i], z[i] };
			double xj[3] = { x[j], y[j], z[j] };

			double r[3];
			double r2 = b.dom.dist_2( xi, xj, r );
			if( r2 > rc2 ) continue;

			double rr = std::sqrt( r2 );
			int bini = ( rr - r0 ) / dr;
			if( bini < 0 || bini >= Nbins ) continue;

			rrdf[bini] += 1.0;
			++adds;
		}
	}

	double Lx = b.dom.xhi[0] - b.dom.xlo[0];
	double Ly = b.dom.xhi[1] - b.dom.xlo[1];
	double Lz = b.dom.xhi[2] - b.dom.xlo[2];

	double V = Lx*Ly;
	if( dims != 2 ) V *= Lz;

	double volume_scale;
	if( dims != 2 ){
		volume_scale = 4.0*pi / (3.0 * V);
	}else{
		volume_scale = pi / V;
	}
	volume_scale *= b.N;

	// V is already set, it is total volume for 3D or area for 2D.
	for( int bin = 0; bin < Nbins; ++bin ){
		double bin_r   = r0 + dr*bin;
		double bin_r_p = bin_r + dr;
		double bin_r3 = bin_r*bin_r*bin_r;
		double bin_r_p3 = bin_r_p*bin_r_p*bin_r_p;
		double n_ideal = volume_scale*( bin_r_p3 - bin_r3 );

		rrdf[bin] /= n_ideal;
	}
	return rrdf;
}

std::vector<double> coord( const block_data &b, double r0, double r1, int dims,
                           const std::vector<double> &rrdf )
{
	int Nbins = rrdf.size();
	std::vector<double> ccoord( rrdf.size() );

	double Lx = b.dom.xhi[0] - b.dom.xlo[0];
	double Ly = b.dom.xhi[1] - b.dom.xlo[1];
	double Lz = b.dom.xhi[2] - b.dom.xlo[2];

	double dr = ( r1 - r0 ) / ( Nbins - 1 );

	double V = Lx*Ly;
	if( dims != 2 ) V *= Lz;

	double volume_scale;
	if( dims != 2 ){
		volume_scale = 4.0*pi / (3.0 * V);
	}else{
		volume_scale = pi / V;
	}
	volume_scale *= b.N;

	ccoord[0] = 0.0;
	for( std::size_t bin = 1; bin < ccoord.size(); ++bin ){
		double bin_r   = r0 + dr*bin;
		double bin_r_p = bin_r + dr;
		double bin_r3 = bin_r*bin_r*bin_r;
		double bin_r_p3 = bin_r_p*bin_r_p*bin_r_p;
		double bin_vol = bin_r_p3 - bin_r3;
		ccoord[bin] = ccoord[bin-1] + rrdf[bin]*volume_scale*bin_vol;
	}

	return ccoord;
}


} // namespace lammps_tools

} // namespace neighborize
