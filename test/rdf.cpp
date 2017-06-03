#include <catch.hpp>

#include "block_data.hpp"
#include "dump_reader_lammps.hpp"
#include "my_timer.hpp"
#include "neighborize.hpp"
#include "rdf.hpp"
#include "util.hpp"
#include "writers.hpp"

#include <string>

TEST_CASE( "Calculating RDF works", "[neigh_rdf]" )
{
	using namespace lammps_tools;
	using namespace lammps_tools::readers;
	using namespace lammps_tools::neighborize;

	std::string dname = "lammps_dump_file_test.dump.bin";
	std::vector<std::string> headers = { "id", "type",
	                                     "x", "y", "z", "c_pe" };
	std::vector<std::string> atom_headers = { "id", "type",
	                                          "x", "y", "z" };
	std::vector<std::vector<int> > neighs;
	std::unique_ptr<dump_reader_lammps> d(
		make_dump_reader_lammps( dname, 2 ) );


	d->set_column_headers( headers );
	d->set_column_header_as_special(   "id", block_data::ID );
	d->set_column_header_as_special( "type", block_data::TYPE );
	d->set_column_header_as_special(    "x", block_data::X );
	d->set_column_header_as_special(    "y", block_data::Y );
	d->set_column_header_as_special(    "z", block_data::Z );

	block_data b;
	for( int i = 0; i < 4; ++i ){
		int status = d->next_block(b);
		REQUIRE( status == 0 );
	}
	std::vector<double> rdf, coords;
	lammps_tools::neighborize::compute_rdf( b, 251, 0.0, 4.0, 3, 0, 0, rdf, coords );

}
