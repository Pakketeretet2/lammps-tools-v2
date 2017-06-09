#include <catch.hpp>

#include "block_data.hpp"
#include "block_data_access.hpp"
#include "center_of_mass.hpp"
#include "types.hpp"


TEST_CASE( "Calculating centers of mass works", "[center_of_mass]" )
{
	using lammps_tools::block_data;
	using lammps_tools::data_field_int;
	using lammps_tools::data_field_double;

	lammps_tools::block_data b( 4 );
	std::vector<double> x = { -1, -1,  1,  1 };
	std::vector<double> y = { -1,  1, -1,  1 };
	std::vector<double> z = {  0,  1,  2,  3 };
	std::vector<int> id   = { 1,2,3,4 };
	std::vector<int> type = { 1, 1, 2, 2 };

	b.set_ntypes( 2 );
	b.ati.mass[1] = 1.0;
	b.ati.mass[2] = 2.0;

	b.add_field( data_field_int(   "id", id ), block_data::ID );
	b.add_field( data_field_int( "type", type ), block_data::TYPE );
	b.add_field( data_field_double( "x", x ), block_data::X );
	b.add_field( data_field_double( "y", y ), block_data::Y );
	b.add_field( data_field_double( "z", z ), block_data::Z );

	lammps_tools::point p1 = geometric_center( b );
	lammps_tools::point p2 = center_of_mass( b );


	REQUIRE( p1.x == Approx(0.0) );
	REQUIRE( p1.y == Approx(0.0) );
	REQUIRE( p1.z == Approx(1.5) );

	REQUIRE( p2.x == Approx(2.0/6.0) );
	REQUIRE( p2.y == Approx(0.0) );
	REQUIRE( p2.z == Approx(11.0/6.0) );

}
