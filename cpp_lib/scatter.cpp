#include "block_data.hpp"
#include "block_data_access.hpp"
#include "id_map.hpp"
#include "scatter.hpp"
#include "fast_math.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace lammps_tools {

namespace scatter {

constexpr const bool use_lut = true;
//constexpr const bool use_lut = false;


inline double bessel_j1( double x )
{
	double inv_x  = 1.0/x;
	if( use_lut ){
		double s, c;
		// sincos.sincos( x, s, c );
		fast_math::sincos( x, s, c );
		return inv_x * ( s * inv_x - c );
	}else{
		double s = std::sin(x);
		double c = std::cos(x);
		return inv_x * ( s * inv_x - c );
	}
}


void get_position( const block_data &b, int idx, double ri[3], double scale )
{
	double xx = get_x(b)[idx];
	double yy = get_y(b)[idx];
	double zz = get_z(b)[idx];

	ri[0] = xx*scale;
	ri[1] = yy*scale;
	ri[2] = zz*scale;
}


std::vector<double> rayleigh_gans( const class block_data &b,
                                   const std::vector<double> &qs,
                                   const std::vector<double> &radius,
                                   double position_scale, double d_epsilon0 )
{
	return rayleigh_gans( b, qs, radius, position_scale,
	                      d_epsilon0, get_id( b ) );
}

std::vector<double> rayleigh_gans( const class block_data &b,
                                   const std::vector<double> &qs,
                                   const std::vector<double> &radius,
                                   double position_scale, double d_epsilon0,
                                   const std::vector<int> &ids )
{
	// In units of micrometers:
	id_map im( get_id( b ) );

	double d_epsilon_02 = d_epsilon0 * d_epsilon0;

	std::size_t Nqs = qs.size();
	std::vector<double> FXre( Nqs, 0.0 );
	std::vector<double> FXim( Nqs, 0.0 );
	std::vector<double> Iq  ( Nqs, 0.0 );

	// Assume detector is far away, and at 90 degrees from box.
	double n [3] = { 0.0, 1.0, 0.0 };
	double n0[3] = { 1.0, 0.0, 0.0 };
	// Wave vector for 500nm in 1 / micrometers
	double qdir[3] = { n[0] - n0[0], n[1] - n0[1], n[2] - n0[2] };
	double inv_n = 1.0 / std::sqrt( util::norm2<3>( qdir ) );
	qdir[0] *= inv_n;
	qdir[1] *= inv_n;
	qdir[2] *= inv_n;


	for( int idi : ids ){
		int idx = im[ idi ];
		double Ri = radius[idx];
		double Ri2 = Ri*Ri;
		double Ri3 = Ri2*Ri;
		// You need to convert the positions somewhow!
		double ri[3];
		get_position( b, idx, ri, position_scale );

		std::size_t bin = 0;
		for( double q : qs ){
			double qq[3] = { q*qdir[0], q*qdir[1], q*qdir[2] };
			double qRi = q * Ri;
			double FFq = Ri3 * bessel_j1( qRi ) / qRi;
			double qdotr = util::dot<3>( qq, ri );
			double exp_part_real, exp_part_imag;
			if( use_lut ){
				fast_math::sincos( -qdotr,
				                   exp_part_imag,
				                   exp_part_real );
			}else{
				exp_part_real = std::cos( -qdotr );
				exp_part_imag = std::sin( -qdotr );
			}

			FXre[bin] += FFq * exp_part_real;
			FXim[bin] += FFq * exp_part_imag;
			++bin;
		}
	}
	// Determine intensity:
	for( std::size_t i = 0; i < Iq.size(); ++i ){
		double abs_val_square = FXre[i]*FXre[i] + FXim[i]*FXim[i];
		Iq[i] = d_epsilon_02 * abs_val_square;
	}
	return Iq;
}


} // namespace scatter

} // namespace lammps_tools
