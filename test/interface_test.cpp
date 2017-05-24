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
	std::cerr << "Block data is located at " << block.bd.get() << ".\n";
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
