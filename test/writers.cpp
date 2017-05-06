#include "readers.hpp"
#include "writers.hpp"

#include "dump_reader_lammps.hpp"
#include "util.hpp"

#include <algorithm>
#include <catch.hpp>
#include <fstream>
#include <memory>


TEST_CASE ( "LAMMPS data file gets written correctly.", "[write_lammps_data]" ) {
	using namespace lammps_tools;

	std::ifstream in( "lammps_data_file_test.data" );
	REQUIRE( in );
	int status = -1;

	block_data b = readers::block_data_from_lammps_data( in, status );
	writers::block_to_lammps_data( "lammps_data_file_test_out.data", b );

}

TEST_CASE ( "LAMMPS dump file gets written correctly.", "[write_lammps_dump]" ) {
	using lammps_tools::writers::block_to_lammps_dump;
	using lammps_tools::readers::block_data_from_lammps_data;
	using lammps_tools::readers::dump_reader_lammps;
	using lammps_tools::readers::make_dump_reader_lammps;

	std::ifstream in( "lammps_data_file_test.data" );
	REQUIRE( in );
	int status = -1;

	lammps_tools::block_data b = block_data_from_lammps_data( in, status );
	block_to_lammps_dump( "lammps_dump_file_test_out.dump", b, 0 );
	block_to_lammps_dump( "lammps_dump_file_test_out.dump.gz", b, 1 );

	// For the binary, read a binary file and write it too to check:
	std::vector<std::string> headers =
		{ "id", "type", "x", "y", "z", "c_pe" };
	std::string fname = "lammps_dump_file_test.dump.bin";
	std::unique_ptr<dump_reader_lammps> d(
		make_dump_reader_lammps( fname, 2, headers ) );
	lammps_tools::block_data b2;

	std::ofstream out( "lammps_dump_file_test_out.dump.bin",
	                   std::ios::binary );
	status = d->next_block( b2 );
	while( status == 0 ){
		block_to_lammps_dump( out, b2, 2 );
		status = d->next_block( b2 );
	}
}
