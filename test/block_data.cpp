#include "block_data.hpp"
#include <catch.hpp>

TEST_CASE ( "block_data assignment works correctly.", "[block_data_assignment]" ) {
	block_data b1;
	data_field_double d( "data", 3 );
	b1.set_natoms( d.size() );
	d[0] = 1.337;
	d[1] = d[0]*2;
	d[2] = d[1]*2;
	SECTION( "data field can be added" ){
		REQUIRE( b1.n_data_fields() == 0);
		b1.add_field(d);
		REQUIRE( b1.n_data_fields() == 1);
		REQUIRE( b1.get_data( "data" ) != nullptr );
		REQUIRE( b1.get_data( "data" )->size() == d.size() );
	}
}

TEST_CASE ( "block_data copy constructor works correctly.", "[block_data_copy_constructor]" ) {
	data_field_double d( "data", 3 );
	d[0] = 1.337;
	d[1] = d[0]*2;
	d[2] = d[1]*2;
	block_data b1( d.size() );
	b1.add_field(d);
	block_data b2(b1);
	REQUIRE( b1.n_data_fields() == 1 );
	REQUIRE( b2.n_data_fields() == b1.n_data_fields() );
	const std::vector<std::string> &names = b1.get_data_names();
	for( const std::string &n : names ){
		REQUIRE( b1.get_data( n ) != nullptr );
		REQUIRE( b2.get_data( n ) != nullptr );
		// Should _NOT_ point to same data:
		REQUIRE( b1.get_data(n) != b2.get_data(n) );
		// Should be same data:
		
		const data_field_double &df1 =
			dynamic_cast<const data_field_double&>( *b1.get_data( n ) );
		const data_field_double &df2 =
			dynamic_cast<const data_field_double&>( *b2.get_data( n ) );
		REQUIRE( df1.name == df2.name );
		
		for( std::size_t i = 0; i < df1.size(); ++i ){
			REQUIRE( df1[i] == df2[i] );
		}
	}
}


