#ifndef SINCOS_LOOKUP_HPP
#define SINCOS_LOOKUP_HPP

#include <cmath>
#include <vector>

namespace lammps_tools {

namespace fast_math {

static constexpr const double pi =
	3.1415926535897932384626433832795028841971693993751058209749446;
static constexpr const double pi2 = 2*pi;



struct lut_sin_cos
{
	static constexpr const bool interpolate = false; // true;
	explicit lut_sin_cos( int N ) : dt( pi2 / static_cast<double>(N) ),
	                                hdt( 0.5*dt ),
	                                sin_tab(N+1), cos_tab(N+1)
	{
		build();
	}
	void build()
	{
		std::size_t N = sin_tab.size();
		double t = 0.0;
		for( std::size_t i = 0; i < N; ++i ){
			sin_tab[i] = std::sin(t);
			cos_tab[i] = std::cos(t);
			t += dt;
		}
	}

	void sincos( double x, double &s, double &c )
	{
		x = std::fmod( x, pi2 );
		//if( x <= 0 ) x += pi2; // Branch, but might be jump on sign bit.
		//x += std::signbit(x)*pi2;
		if( std::signbit(x) ) x += pi2;

		if( interpolate ){
			double rem = fmod( x, dt );
			int    bin = x / dt;
			double xs  = rem / dt;
			double xm  = 1.0 - xs;
			s = sin_tab[bin]*xm + xs*sin_tab[bin+1];
			c = cos_tab[bin]*xm + xs*cos_tab[bin+1];
		}else{
			int    bin = (x + hdt) / dt;
			s = sin_tab[bin];
			c = cos_tab[bin];
		}
	}

	double dt, hdt;
	std::vector<double> sin_tab, cos_tab;
};

} // namespace fast_math

} // namespace lammps_tools

#endif // SINCOS_LOOKUP_HPP
