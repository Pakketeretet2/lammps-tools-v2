#ifndef TRIANGULATE_HPP
#define TRIANGULATE_HPP


#include <vector>
#include <cmath>

#include "block_data.hpp"
#include "my_assert.hpp"

namespace lammps_tools {

namespace triangulate {


inline void copy_arr( double *dest, const double *src, int N )
{
	for( int i = 0; i < N; ++i ){
		dest[i] = src[i];
	}
}



struct triangle
{
	int i1, i2, i3;
	double x1[3], x2[3], x3[3];

	triangle( int j1, int j2, int j3,
	          const double y1[3], const double y2[3], const double y3[3] )
		: i1(-1), i2(-1), i3(-1), x1{0,0,0}, x2{0,0,0}, x3{0,0,0}
	{
		// Find the right mapping, which is indices sorted
		// from low to high.
		if ( j1 < j2 ){
			if( j1 < j3 ){
				// Order is j1 < (j2,j3)
				i1 = j1;
				copy_arr( x1, y1, 3 );

				if( j2 < j3 ){
					i2 = j2;
					i3 = j3;
					copy_arr( x2, y2, 3 );
					copy_arr( x3, y3, 3 );
				}else{
					i2 = j3;
					i3 = j2;
					copy_arr( x2, y3, 3 );
					copy_arr( x3, y2, 3 );
				}
			}else{
				// Order is j3 < (j1,j2):
				i1 = j3;
				copy_arr( x1, y3, 3 );
				if( j1 < j2 ){
					i2 = j1;
					i3 = j2;
					copy_arr( x2, y1, 3 );
					copy_arr( x3, y2, 3 );
				}else{
					i2 = j2;
					i3 = j1;
					copy_arr( x2, y2, 3 );
					copy_arr( x3, y1, 3 );
				}
			}
		}else{
			if( j2 < j3 ){
				// Order is j2 < (j1,j3)
				i1 = j2;
				copy_arr( x1, y2, 3 );
				if( j1 < j3 ){
					i2 = j1;
					i3 = j3;
					copy_arr( x2, y1, 3 );
					copy_arr( x3, y3, 3 );
				}else{
					i2 = j3;
					i3 = j1;
					copy_arr( x2, y3, 3 );
					copy_arr( x3, y1, 3 );
				}
			}else{
				// Order is j3 < (j1,j2)
				i1 = j3;
				copy_arr( x1, y3, 3 );
				if( j1 < j2 ){
					i2 = j1;
					i3 = j2;
					copy_arr( x2, y1, 3 );
					copy_arr( x3, y2, 3 );
				}else{
					i2 = j2;
					i3 = j1;
					copy_arr( x2, y2, 3 );
					copy_arr( x3, y1, 3 );
				}
			}
		}
		my_assert( __FILE__, __LINE__, i1 <= i2,
		           "Id order error!" );
		my_assert( __FILE__, __LINE__, i1 <= i3,
		           "Id order error!" );

	}


	bool operator==( const triangle &o ) const
	{
		return (i1==o.i1) && (i2==o.i2) && (i3==o.i3);
	}


	double area() const
	{
		// area using Heron's formula:
		double dx1[3];
		double dx2[3];
		double dx3[3];
		dx1[0] = x1[0] - x2[0];
		dx1[1] = x1[1] - x2[1];
		dx1[2] = x1[2] - x2[2];

		dx2[0] = x1[0] - x3[0];
		dx2[1] = x1[1] - x3[1];
		dx2[2] = x1[2] - x3[2];

		dx3[0] = x3[0] - x2[0];
		dx3[1] = x3[1] - x2[1];
		dx3[2] = x3[2] - x2[2];

		double a = std::sqrt(dx1[0]*dx1[0] + dx1[1]*dx1[1] + dx1[2]*dx1[2]);
		double b = std::sqrt(dx2[0]*dx2[0] + dx2[1]*dx2[1] + dx2[2]*dx2[2]);
		double c = std::sqrt(dx3[0]*dx3[0] + dx3[1]*dx3[1] + dx3[2]*dx3[2]);

		double s = 0.5*( a + b + c );
		double sa = s - a;
		double sb = s - b;
		double sc = s - c;

		return std::sqrt( s*sa*sb*sc );
	}
};


void triangulate_block( class block_data &b, double rc, int periodic,
                        int dims, int method, std::vector<triangle> &triangles );

double triangulation_area( class block_data &b, std::vector<triangle> &triangles );


} // namespace triangulate

} // namespace lammps_tools


#endif // TRIANGULATE_HPP
