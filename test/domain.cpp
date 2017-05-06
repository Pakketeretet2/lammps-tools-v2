#include <catch.hpp>

#include "domain.hpp"

TEST_CASE( "Tests if periodic boundaries work as expected.", "[domain_pbc]" ){
	using lammps_tools::domain;
	domain dom;
	dom.xlo[0] = dom.xlo[1] = dom.xlo[2] = 0.0;
	dom.xhi[0] = dom.xhi[1] = dom.xhi[2] = 4.0;

	double x1[3] = {    2,   1, 0 };
	double x2[3] = {    1,   1, 0 };
	double x3[3] = {  3.9, 3.9, 0 };
	double x4[3] = {    0,   1, 1 };

	/* 
	   Non-periodic distances:
	   x1 --> x2  = {    1,    0,  0 }
	   x1 --> x3  = { -1.1, -2.1,  0 }
	   x1 --> x4  = {    2,    0, -1 }
	   x2 --> x3  = { -2.1, -2.1,  0 }
	   x2 --> x4  = {    1,    0, -1 }
	   x3 --> x4  = { 3.9,   2.9, -1 }
	*/
	dom.periodic = 0;
	double r[3];
	dom.dist_2( x1, x2, r );
	REQUIRE( r[0] == Approx( 1.0 ) );
	REQUIRE( r[1] == Approx( 0.0 ) );
	REQUIRE( r[2] == Approx( 0.0 ) );
	
	dom.dist_2( x1, x3, r );
	REQUIRE( r[0] == Approx( -1.9 ) );
	REQUIRE( r[1] == Approx( -2.9 ) );
	REQUIRE( r[2] == Approx(  0.0 ) );
	
	dom.dist_2( x1, x4, r );
	REQUIRE( r[0] == Approx( 2.0 ) );
	REQUIRE( r[1] == Approx( 0.0 ) );
	REQUIRE( r[2] == Approx( -1.0 ) );
	
	dom.dist_2( x2, x3, r );
	REQUIRE( r[0] == Approx( -2.9 ) );
	REQUIRE( r[1] == Approx( -2.9 ) );
	REQUIRE( r[2] == Approx(  0.0 ) );
	
	dom.dist_2( x2, x4, r );
	REQUIRE( r[0] == Approx(  1.0 ) );
	REQUIRE( r[1] == Approx(  0.0 ) );
	REQUIRE( r[2] == Approx( -1.0 ) );
	
	dom.dist_2( x3, x4, r );
	REQUIRE( r[0] == Approx(  3.9 ) );
	REQUIRE( r[1] == Approx(  2.9 ) );
	REQUIRE( r[2] == Approx( -1.0 ) );

	/*
	   Periodic distances:
	   x1 --> x2  = {    1,    0,  0 }
	   x1 --> x3  = { -1.9,  1.1,  0 }
	   x1 --> x4  = {    2,    0,  0 }
	   x2 --> x3  = {  1.9,  1.9,  0 }
	   x2 --> x4  = {    1,    0, -1 }
	   x3 --> x4  = { -0.1,  -1.1, -1 }
	*/
	dom.periodic = domain::BIT_X + domain::BIT_Y + domain::BIT_Z;
	dom.dist_2( x1, x2, r );
	REQUIRE( r[0] == Approx( 1.0 ) );
	REQUIRE( r[1] == Approx( 0.0 ) );
	REQUIRE( r[2] == Approx( 0.0 ) );
	
	dom.dist_2( x1, x3, r );
	REQUIRE( r[0] == Approx( -1.9 ) );
	REQUIRE( r[1] == Approx(  1.1 ) );
	REQUIRE( r[2] == Approx(  0.0 ) );
	
	dom.dist_2( x1, x4, r );
	REQUIRE( r[0] == Approx(  2.0 ) );
	REQUIRE( r[1] == Approx(  0.0 ) );
	REQUIRE( r[2] == Approx( -1.0 ) );
	
	dom.dist_2( x2, x3, r );
	REQUIRE( r[0] == Approx(  1.1 ) );
	REQUIRE( r[1] == Approx(  1.1 ) );
	REQUIRE( r[2] == Approx(  0.0 ) );
	
	dom.dist_2( x2, x4, r );
	REQUIRE( r[0] == Approx(  1.0 ) );
	REQUIRE( r[1] == Approx(  0.0 ) );
	REQUIRE( r[2] == Approx( -1.0 ) );
	
	dom.dist_2( x3, x4, r );
	REQUIRE( r[0] == Approx( -0.1 ) );
	REQUIRE( r[1] == Approx( -1.1 ) );
	REQUIRE( r[2] == Approx( -1.0 ) );
	
}
