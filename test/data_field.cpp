#include "data_field.hpp"
#include <catch.hpp>

using namespace lammps_tools;

TEST_CASE ( "Data indicators are properly set up.", "[data_field_data_types]" ) {
	data_field_int    ids( "id" );
	data_field_int    types( "types" );
	data_field_double x( "x" );

	REQUIRE( x.type() == data_field::DOUBLE );
	REQUIRE( ids.type() == data_field::INT );
	REQUIRE( types.type() == data_field::INT );

	REQUIRE( x.type() != data_field::INT );
	REQUIRE( ids.type() != data_field::DOUBLE );
	REQUIRE( types.type() != data_field::DOUBLE );
}


TEST_CASE ( "Data field creation works as expected.", "[data_field_creation]" ) {
	data_field_int a( "a" );
	data_field_int b( "b", 5 );
	data_field_double c( "c", 5 );

	INFO( "Type of a is " << a.type() );
	INFO( "Type of b is " << b.type() );
	INFO( "Type of c is " << c.type() );

	b[0] = b[1] = b[2] = b[3] = 2;
	b[4] = 3;

	c[0] = 1.2;
	c[1] = 0.3 + c[0];
	c[2] = 0.3 + c[1];
	c[3] = 0.3 + c[2];
	c[4] = 0.3 + c[3];
	INFO( "a is named " << a.name );
	INFO( "b is named " << b.name );
	INFO( "c is named " << c.name );

	REQUIRE( c[4] == Approx(2.4) );
}

TEST_CASE ( "Data field copying works as expected.", "[data_field_copy]" ) {
	data_field_double c( "c", 5 );

	c[0] = 1.2;
	c[1] = 0.3 + c[0];
	c[2] = 0.3 + c[1];
	c[3] = 0.3 + c[2];
	c[4] = 0.3 + c[3];

	data_field *copy_of_c = copy( &c );
	REQUIRE( copy_of_c != nullptr );
	REQUIRE( copy_of_c->type() == c.type() );
	REQUIRE( copy_of_c->size() == c.size() );

	delete copy_of_c;
}

TEST_CASE ( "Data field copy constructor-by-ref.", "[data_field_constructor]" ) {
	data_field_double c( "c", 5 );

	c[0] = 1.2;
	c[1] = 0.3 + c[0];
	c[2] = 0.3 + c[1];
	c[3] = 0.3 + c[2];
	c[4] = 0.3 + c[3];

	data_field_double copy_of_c( c );
	REQUIRE( copy_of_c.type() == c.type() );
	REQUIRE( copy_of_c.size() == c.size() );

	for( int i = 0; i < 5; ++i ){
		REQUIRE( copy_of_c[i] == c[i] );
	}
}


TEST_CASE ( "Data field's data_as works as expected.", "[data_field_data_as]" ) {
	data_field_double c( "c", 5 );
	data_field_int d( "d", 5 );

	c[0] = 1.2;
	c[1] = 0.3 * c[0];
	c[2] = 0.3 * c[1];
	c[3] = 0.3 * c[2];
	c[4] = 0.3 * c[3];

	d[0] = 1;
	d[1] = 2 * d[0];
	d[2] = 2 * d[1];
	d[3] = 2 * d[2];
	d[4] = 2 * d[3];


	const std::vector<double> &c_data = data_as<double>( &c );
	const std::vector<int> &d_data    = data_as<int>( &d );

	for( int i = 0; i < 5; ++i ){
		REQUIRE( c_data[i] == c[i] );
		REQUIRE( d_data[i] == d[i] );
	}
}
