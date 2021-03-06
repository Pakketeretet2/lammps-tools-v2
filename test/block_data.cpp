#include "block_data.hpp"
#include "block_data_access.hpp"
#include "data_field.hpp"
#include "dump_reader_lammps.hpp"

#include <catch.hpp>

#include <memory>

TEST_CASE ( "block_data constructor works correctly.", "[block_data_constructor]" ) {

	using namespace lammps_tools;

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

	using namespace lammps_tools;

	data_field_double d( "data", 3 );
	data_field_int id( "id", 3 );
	id[0] = 1;
	id[1] = 3;
	id[2] = 2;
	d[0] = 1.337;
	d[1] = d[0]*2;
	d[2] = d[1]*2;
	block_data b1( d.size() );
	b1.add_field( d );
	b1.add_field( id, block_data::ID );
	block_data b2(b1);
	REQUIRE( b1.n_data_fields() == 2 );
	REQUIRE( b2.n_data_fields() == b1.n_data_fields() );
	REQUIRE( b1.get_special_field( block_data::ID) );
	REQUIRE( b2.get_special_field( block_data::ID) );

	for( std::size_t i = 0; i < b1.n_data_fields(); ++i ){
		const std::string &n = b1[i].name;
		REQUIRE( b1.get_data( n ) != nullptr );
		REQUIRE( b2.get_data( n ) != nullptr );
		// Should _NOT_ point to same data:
		REQUIRE( b1.get_data(n) != b2.get_data(n) );
		// Should be same data:
		if( n == "data" ){
			const data_field_double &df1 =
				dynamic_cast<const data_field_double&>( *b1.get_data( n ) );
			const data_field_double &df2 =
				dynamic_cast<const data_field_double&>( *b2.get_data( n ) );
			REQUIRE( df1.name == df2.name );
			for( std::size_t i = 0; i < df1.size(); ++i ){
				REQUIRE( df1[i] == df2[i] );
			}

		}else{
			const data_field_int &df1 =
				dynamic_cast<const data_field_int&>( *b1.get_data( n ) );
			const data_field_int &df2 =
				dynamic_cast<const data_field_int&>( *b2.get_data( n ) );
			REQUIRE( df1.name == df2.name );
			for( std::size_t i = 0; i < df1.size(); ++i ){
				REQUIRE( df1[i] == df2[i] );
			}

		}

	}
}



TEST_CASE ( "block_data assignment works correctly.", "[block_data_assignment]" ) {

	using namespace lammps_tools;

	block_data b1;
	data_field_double d( "data", 3 );
	b1.set_natoms( d.size() );
	d[0] = 1.337;
	d[1] = d[0]*2;
	d[2] = d[1]*2;
	b1.add_field(d);
	block_data b2 = b1;
	const data_field_double *d1 = static_cast<const data_field_double*>( b1.get_data("data") );
	const data_field_double *d2 = static_cast<const data_field_double*>( b2.get_data("data") );
	REQUIRE( d1 != nullptr );
	REQUIRE( d2 != nullptr );
	REQUIRE( d1 != d2 );

	const data_field_double &dd1 = *d1;
	const data_field_double &dd2 = *d2;


	REQUIRE( dd1[0] == Approx(1.337) );
	REQUIRE( dd1[1] == Approx(2.674) );
	REQUIRE( dd1[2] == Approx(5.348) );
	REQUIRE( dd1[0] == dd2[0] );
	REQUIRE( dd1[1] == dd2[1] );
	REQUIRE( dd1[2] == dd2[2] );
}



TEST_CASE ( "block_data remove_field works correctly.", "[block_data_add_remove]" ) {

	using namespace lammps_tools;
	using dfd = data_field_double;

	block_data b;
	dfd d1( "data", 3 );
	dfd d2( "data_2", 3 );
	dfd  x( "x", 3 );
	dfd  y( "y", 3 );

	b.set_natoms( d1.size() );
	d1[0] = 1.337;
	d1[1] = d1[0]*2;
	d1[2] = d1[1]*2;
	d2[0] = -12;
	d2[1] = -18;
	d2[2] =   4;
	x[0] = 0.1;
	x[1] = 1.1;
	x[2] = 2.1;
	y[0] = y[1] = y[2] = 3;

	b.add_field(d1);
	b.add_field(d2);
	b.add_field( x, block_data::X );
	b.add_field( y, block_data::Y );


	const dfd *dd1 = static_cast<const dfd*>( b.get_data("data") );
	const dfd *dd2 = static_cast<const dfd*>( b.get_data("data_2") );
	const dfd *dx  = static_cast<const dfd*>( b.get_data("x") );
	const dfd *dy  = static_cast<const dfd*>( b.get_data("y") );

	REQUIRE( dd1 );
	REQUIRE( dd2 );
	REQUIRE( dx );
	REQUIRE( dy );

	REQUIRE( b.get_special_field_name( block_data::X ) == "x" );
	REQUIRE( b.get_special_field_name( block_data::Y ) == "y" );
	REQUIRE( b.get_special_field_name( block_data::Z ) == "" );
	REQUIRE( b.n_data_fields() == 4 );

	REQUIRE( b[0].name == "data" );
	REQUIRE( b[1].name == "data_2" );
	REQUIRE( b[2].name == "x" );
	REQUIRE( b[3].name == "y" );


	int special_field;
	std::shared_ptr<data_field> tmp( b.remove_field( "x", special_field ) );

	REQUIRE( tmp.get() == dx );
	REQUIRE( special_field == block_data::X );

	REQUIRE( tmp->name == "x" );
	REQUIRE( b.n_data_fields() == 3 );
	REQUIRE( b.get_special_field_name( block_data::X ) == "" );
	REQUIRE( b.get_special_field_name( block_data::Y ) == "y" );

	REQUIRE( b[0].name == "data" );
	REQUIRE( b[1].name == "data_2" );
	REQUIRE( b[2].name == "y" );

	data_field *tmp_y = b.remove_field( "y", special_field );

	REQUIRE( tmp_y == dy );
	REQUIRE( special_field == block_data::Y );

	REQUIRE( tmp_y->name == "y" );
	REQUIRE( b.n_data_fields() == 2 );
	REQUIRE( b.get_special_field_name( block_data::X ) == "" );
	REQUIRE( b.get_special_field_name( block_data::Y ) == "" );

	REQUIRE( b[0].name == "data" );
	REQUIRE( b[1].name == "data_2" );

	delete tmp_y;

	b.add_field( *tmp, block_data::X );


	REQUIRE( b.n_data_fields() == 3 );
	REQUIRE( b.get_special_field_name( block_data::X ) == "x" );
	REQUIRE( b.get_special_field_name( block_data::Y ) == "" );

	REQUIRE( b[0].name == "data" );
	REQUIRE( b[1].name == "data_2" );
	REQUIRE( b[2].name == "x" );
}




TEST_CASE ( "Convenience accessors work.", "[block_data_convenience]" ) {
	using namespace lammps_tools;
	using dfd = data_field_double;

	block_data b;
	dfd d1( "data", 3 );
	dfd d2( "data_2", 3 );
	dfd  x( "x", 3 );
	dfd  y( "y", 3 );

	b.set_natoms( d1.size() );
	d1[0] = 1.337;
	d1[1] = d1[0]*2;
	d1[2] = d1[1]*2;
	d2[0] = -12;
	d2[1] = -18;
	d2[2] =   4;
	x[0] = 0.1;
	x[1] = 1.1;
	x[2] = 2.1;
	y[0] = y[1] = y[2] = 3;

	b.add_field(d1);
	b.add_field(d2);
	b.add_field( x, block_data::X );
	b.add_field( y, block_data::Y );

	std::vector<double>       &x_rw = get_x_rw( b );
	const std::vector<double> &y_ro = get_y( b );

	REQUIRE( x_rw[0] == Approx(x[0]) );
	REQUIRE( x_rw[1] == Approx(x[1]) );
	REQUIRE( x_rw[2] == Approx(x[2]) );
	REQUIRE( x_rw[3] == Approx(x[3]) );
	REQUIRE( y_ro[0] == Approx(y[0]) );
	REQUIRE( y_ro[1] == Approx(y[1]) );
	REQUIRE( y_ro[2] == Approx(y[2]) );
	REQUIRE( y_ro[3] == Approx(y[3]) );

}
