#ifndef FAST_MATH_HPP
#define FAST_MATH_HPP

#include "sincos_lookup.hpp"

namespace lammps_tools {

namespace fast_math {

void sincos( double x, double &s, double &c )
{
	static lut_sin_cos table( 512 );
	table.sincos( x, s, c );
}


double sqrt( double x )
{
	// initial guess:
	if( x < 1e-4 ) return 0.0;
	/*
	int factors = 0;
	while( x > 4.0 ){
		x /= 4.0;
		++factors;
	}
	int corr = 1 << factors;
	*/

	// Smart initial guess: take x, divide exponent by 2.

	// Nice union for this, see
	// https://stackoverflow.com/questions/15685181/how-to-get-the-sign-
	//     mantissa-and-exponent-of-a-floating-point-number
	union double_int {
		double d;
		struct {
			unsigned long int mantissa : 52;
			unsigned long int exponent : 11;
			unsigned long int sign     : 1;
		} parts;
	};

	double_int xx;
	xx.d = x;
	//std::cerr << "Before division, x as double is " << xx.d << "\n";
	//std::cerr << "                 x has exponent " << xx.parts.exponent << "\n";
	xx.parts.exponent /= 2;
	xx.parts.exponent += 512;
	double y = xx.d;
	//std::cerr << "After  division, x as double is " << y << "\n";
	//std::cerr << "                 x has exponent " << xx.parts.exponent << "\n";
	//xx.parts.sign = 1;
	//std::cerr << "Setting the sign bit leads to " << xx.d << "\n";
	//std::cerr << "\n\n";


	y -= 0.5*( y - x / y );
	y -= 0.5*( y - x / y );
	y -= 0.5*( y - x / y );
	y -= 0.5*( y - x / y );

	return y;
}

// 1 / sqrt(x)
double inv_sqrt( double x )
{
	if( x > 1e8 ) return 0; // 1e-4 precision.

	int factors = 0;
	while( x < 0.25 ){
		x *= 4;
		++factors;
	}
	double corr = 1.0 / static_cast<double>( 1 << factors );
	double y = 0.618;
	double yy = y*y;
	y -= y*( 1.5 - x * yy );
	yy = y*y;
	y -= y*( 1.5 - x * yy );
	yy = y*y;
	y -= y*( 1.5 - x * yy );
	yy = y*y;
	y -= y*( 1.5 - x * yy );
	yy = y*y;

	return y * corr;

}

} // fast_math

} // lammps_tools


#endif // FAST_MATH_HPP
