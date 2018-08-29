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


TEST_CASE ( "HOOMD dump file file gets written correctly.", "[write_hoomd_gsd]" )
{
	std::string fname = "triangles.gsd";

	using namespace lammps_tools;
	using namespace readers;
	using namespace writers;

	std::shared_ptr<dump_reader> r(
		make_dump_reader( fname, FILE_FORMAT_BIN, DUMP_FORMAT_HOOMD ) );
	block_data b;
	std::cerr << "Reading block...\n";
	int status = r->next_block(b);
	REQUIRE( status == 0 );
	REQUIRE( b.tstep == 0 );
	REQUIRE( b.N == 98560 );
	REQUIRE( b.N_types == 8 );

	// Try to write, make sure no info is lost?
	std::cerr << "At block " << b.tstep << ", N = " << b.N << "\n";

	block_to_hoomd_gsd( "triangles_copy.gsd", b, "wb" );
	block_to_hoomd_gsd( "triangles_copy2.gsd", b, "ab" );

	int nblock = 0;
	while( r->next_block(b) == 0 && nblock < 250 ){
		std::cerr << "At block " << b.tstep << ", N = " << b.N << "\n";
		block_to_hoomd_gsd( "triangles_copy.gsd", b, "ab" );
		++nblock;
	}
}
