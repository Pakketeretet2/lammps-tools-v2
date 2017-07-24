// Tests the C interface separately.
#include <catch.hpp>

#include "data_field.hpp"
#include "enums.hpp"
#include "lammps_tools.h"



TEST_CASE ( "Tests the creation of a dump_reader and getting a block.", "[interface_dump_reader]" )
{
	using dfd = lammps_tools::data_field_double;

	const char *dname = "lammps_dump_file_test.dump.bin";
	int fformat = lammps_tools::FILE_FORMAT_BIN;
	int dformat = lammps_tools::DUMP_FORMAT_LAMMPS;

	lt_dump_reader_handle reader = lt_new_dump_reader( dname, fformat, dformat );
	int n_col = 6;
	char **col_headers = new char*[n_col];
	for( int i = 0; i < n_col; ++i ){
		col_headers[i] = new char[256];
	}
	sprintf(col_headers[0], "id" );
	sprintf(col_headers[1], "type" );
	sprintf(col_headers[2], "x" );
	sprintf(col_headers[3], "y" );
	sprintf(col_headers[4], "z" );
	sprintf(col_headers[5], "c_pe" );

	for( int i = 0; i < n_col; ++i ){
		lt_set_col_header( reader, i, col_headers[i] );
	}


	for( int i = 0; i < n_col; ++i ){
		delete [] col_headers[i];
	}
	delete [] col_headers;


	REQUIRE( lt_dump_reader_status(reader) == IS_GOOD );

	lt_block_data_handle block;
	std::cerr << "Block data is located at " << block.bd << ".\n";
	for( int i = 0; i < 3; ++i ){
		int status = lt_get_next_block( reader, &block );

		status = lt_get_next_block( reader, &block );
		REQUIRE( status == 0 );
		REQUIRE( block.bd != nullptr );
		REQUIRE( block.bd->N == 4000 );

		const lammps_tools::data_field *f_x = block.bd->get_data( "x" );
		const dfd *df_x = static_cast<const dfd*>( f_x );
		REQUIRE( df_x != nullptr );
		const dfd &x = *df_x;
		const std::vector<double> &xx =
			lammps_tools::data_as<double>( block.bd->get_data("x") );

		// No need for Approx, they should be _exactly_ the same!
		REQUIRE( x[0] == xx[0] );
	}
	lt_delete_dump_reader( reader );
}



TEST_CASE ( "Tests getting, adding and changing data fields.", "[interface_data_fields]" )
{
	using dfd = lammps_tools::data_field_double;
	const char *dname = "lammps_dump_file_test.dump.bin";
	int fformat = lammps_tools::FILE_FORMAT_BIN;
	int dformat = lammps_tools::DUMP_FORMAT_LAMMPS;

	lt_dump_reader_handle reader = lt_new_dump_reader( dname, fformat, dformat );
	int n_col = 6;
	char **col_headers = new char*[n_col];
	for( int i = 0; i < n_col; ++i ){
		col_headers[i] = new char[256];
	}
	sprintf(col_headers[0], "id" );
	sprintf(col_headers[1], "type" );
	sprintf(col_headers[2], "x" );
	sprintf(col_headers[3], "y" );
	sprintf(col_headers[4], "z" );
	sprintf(col_headers[5], "c_pe" );

	for( int i = 0; i < n_col; ++i ){
		lt_set_col_header( reader, i, col_headers[i] );
	}


	for( int i = 0; i < n_col; ++i ){
		delete [] col_headers[i];
	}
	delete [] col_headers;


	REQUIRE( lt_dump_reader_status(reader) == IS_GOOD );

	lt_block_data_handle block;
	std::cerr << "Block data is located at " << block.bd << ".\n";
	int status = lt_get_next_block( reader, &block );
	status = lt_get_next_block( reader, &block );
	REQUIRE( status == 0 );

	const lammps_tools::data_field *f_x = block.bd->get_data( "x" );
	lammps_tools::data_field *f_z = block.bd->get_data_rw( "z" );


	const dfd *df_x = static_cast<const dfd*>( f_x );
	dfd *df_z = static_cast<dfd*>( f_z );

	REQUIRE( df_x != nullptr );
	REQUIRE( df_z != nullptr );

	// Test setting data:
	std::cerr << "df_x[2] = " << (*df_x)[2] << "\n";
	std::cerr << "df_z[2] = " << (*df_z)[2] << "\n";

	// Change stuff:
	lt_data_field_handle dfh_x( df_x );
	lt_data_field_handle dfh_z( df_z );

	double df_x2 = lt_data_field_get_indexed_double_data( &dfh_x, 2 );
	double df_z2 = lt_data_field_get_indexed_double_data( &dfh_z, 2 );

	double df_x_old = df_x2;

	REQUIRE( df_x2 == (*df_x)[2] );
	REQUIRE( df_z2 == (*df_z)[2] );

	int status_x = lt_data_field_set_indexed_double_data( &dfh_x, 2, 1337 );
	int status_z = lt_data_field_set_indexed_double_data( &dfh_z, 2, -666 );

	REQUIRE( status_x != 0 );
	REQUIRE( status_z == 0 );

	df_x2 = lt_data_field_get_indexed_double_data( &dfh_x, 2 );
	df_z2 = lt_data_field_get_indexed_double_data( &dfh_z, 2 );

	REQUIRE( df_x2 == (*df_x)[2] );
	REQUIRE( df_z2 == (*df_z)[2] );
	REQUIRE( df_x2 == df_x_old );
	REQUIRE( df_z2 == -666 );

	lt_delete_dump_reader( reader );
}
