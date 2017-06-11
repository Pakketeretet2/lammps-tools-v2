#include <catch.hpp>

#include "geometry.hpp"

TEST_CASE ( "Quaternions", "[quat]" )
{
	using lammps_tools::quat;
	using lammps_tools::point;

	quat one( 1.0, 2.0, 3.0, 4.0 );
	quat two( 1.0, 1.0, 1.0, 0.0 );

	quat oneptwo = one + two;
	quat onemtwo = one - two;
	quat twomone = two - one;

	quat one_con = one.conj();
	quat two_con = two.conj();
	quat one_inv = one.inv();
	quat two_inv = two.inv();

	quat oneone   = one*one;
	quat oneone_i = one*one_inv;
	quat twotwo   = two*two;
	quat twotwo_i = two*two_inv;
	quat onetwo   = one*two;
	quat twoone   = two*one;

	REQUIRE( one_con[0] == Approx(  1.0 ) );
	REQUIRE( one_con[1] == Approx( -2.0 ) );
	REQUIRE( one_con[2] == Approx( -3.0 ) );
	REQUIRE( one_con[3] == Approx( -4.0 ) );

	REQUIRE( two_con[0] == Approx(  1.0 ) );
	REQUIRE( two_con[1] == Approx( -1.0 ) );
	REQUIRE( two_con[2] == Approx( -1.0 ) );
	REQUIRE( two_con[3] == Approx(  0.0 ) );

	REQUIRE( oneone_i[0] == Approx( 1.0 ) );
	REQUIRE( oneone_i[1] == Approx( 0.0 ) );
	REQUIRE( oneone_i[2] == Approx( 0.0 ) );
	REQUIRE( oneone_i[3] == Approx( 0.0 ) );

	REQUIRE( oneone_i[0] == Approx( 1.0 ) );
	REQUIRE( oneone_i[1] == Approx( 0.0 ) );
	REQUIRE( oneone_i[2] == Approx( 0.0 ) );
	REQUIRE( oneone_i[3] == Approx( 0.0 ) );

	REQUIRE( oneptwo[0] == Approx( 2 ) );
	REQUIRE( oneptwo[1] == Approx( 3 ) );
	REQUIRE( oneptwo[2] == Approx( 4 ) );
	REQUIRE( oneptwo[3] == Approx( 4 ) );

	REQUIRE( onemtwo[0] == Approx( 0 ) );
	REQUIRE( onemtwo[1] == Approx( 1 ) );
	REQUIRE( onemtwo[2] == Approx( 2 ) );
	REQUIRE( onemtwo[3] == Approx( 4 ) );

	REQUIRE( twomone[0] == Approx( -0 ) );
	REQUIRE( twomone[1] == Approx( -1 ) );
	REQUIRE( twomone[2] == Approx( -2 ) );
	REQUIRE( twomone[3] == Approx( -4 ) );

	REQUIRE( onetwo[0] == Approx( -4 ) );
	REQUIRE( onetwo[1] == Approx( -1 ) );
	REQUIRE( onetwo[2] == Approx(  8 ) );
	REQUIRE( onetwo[3] == Approx(  3 ) );

        REQUIRE( twoone[0] == Approx( -4 ) );
	REQUIRE( twoone[1] == Approx(  7 ) );
	REQUIRE( twoone[2] == Approx(  0 ) );
	REQUIRE( twoone[3] == Approx(  5 ) );



}
