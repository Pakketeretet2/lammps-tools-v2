#include "dump_reader_lammps.hpp"
#include "enums.hpp"
#include "id_map.hpp"
#include "readers.hpp"
#include "util.hpp"

#include <algorithm>
#include <catch.hpp>
#include <fstream>


TEST_CASE ( "LAMMPS data file gets read correctly.", "[read_lammps_data]" )
{
	using namespace lammps_tools;

	std::ifstream in( "lammps_data_file_test.data" );
	REQUIRE( in );
	int status = -1;

	block_data b = readers::block_data_from_lammps_data( in, status, false );
	REQUIRE( b.n_data_fields() == 11 );
	// These entries should _not_ be found:
	REQUIRE( b.get_data( "xlo" ) == nullptr );
	REQUIRE( b.get_data( "xhi" ) == nullptr );
	REQUIRE( b.get_data( "atom_style" ) == nullptr );
	REQUIRE( b.get_data( "timestep" ) == nullptr );


	// These entries should all be found:
	REQUIRE( b.get_data( "id" )   != nullptr );
	REQUIRE( b.get_data( "type" ) != nullptr );
	REQUIRE( b.get_data( "x" )    != nullptr );
	REQUIRE( b.get_data( "y" )    != nullptr );
	REQUIRE( b.get_data( "z" )    != nullptr );
	REQUIRE( b.get_data( "vx" )   != nullptr );
	REQUIRE( b.get_data( "vy" )   != nullptr );
	REQUIRE( b.get_data( "vz" )   != nullptr );
	REQUIRE( b.get_data( "ix" )   != nullptr );
	REQUIRE( b.get_data( "iy" )   != nullptr );
	REQUIRE( b.get_data( "iz" )   != nullptr );

	// And these should be set:
	REQUIRE( b.get_special_field_name(block_data::ID) != "" );
	REQUIRE( b.get_special_field_name(block_data::TYPE) != "" );
	REQUIRE( b.get_special_field_name(block_data::X) != "" );
	REQUIRE( b.get_special_field_name(block_data::Y) != "" );
	REQUIRE( b.get_special_field_name(block_data::Z) != "" );
	REQUIRE( b.get_special_field_name(block_data::VX) != "" );
	REQUIRE( b.get_special_field_name(block_data::VY) != "" );
	REQUIRE( b.get_special_field_name(block_data::VZ) != "" );
	REQUIRE( b.get_special_field_name(block_data::IX) != "" );
	REQUIRE( b.get_special_field_name(block_data::IY) != "" );
	REQUIRE( b.get_special_field_name(block_data::IZ) != "" );
}



TEST_CASE ( "LAMMPS plain text dump file gets read correctly.", "[read_lammps_dump_plain]" )
{
	using namespace lammps_tools;
	using namespace readers;

	std::string fname = "lammps_dump_file_test.dump";
	std::vector<std::string> headers = { "id", "type", "x", "y", "z", "c_pe" };
	std::shared_ptr<dump_reader> r( make_dump_reader( fname, FILE_FORMAT_PLAIN, DUMP_FORMAT_LAMMPS ) );

	REQUIRE( r );
	REQUIRE( r->good() );

	double lo = 0.0;
	double hi = 1.6795961913825074e+01;

	block_data b;

	for( int i = 0; i < 2; ++i ){
		int status = r->next_block( b );
		break;
		REQUIRE( status == 0 );
		REQUIRE( b.N == 4000 );

		REQUIRE( b.dom.xlo[0] == Approx(lo) );
		REQUIRE( b.dom.xlo[1] == Approx(lo) );
		REQUIRE( b.dom.xlo[2] == Approx(lo) );

		REQUIRE( b.dom.xhi[0] == Approx(hi) );
		REQUIRE( b.dom.xhi[1] == Approx(hi) );
		REQUIRE( b.dom.xhi[2] == Approx(hi) );

		REQUIRE( b.get_data( "x" ) );
		REQUIRE( b.get_data( "y" ) );
		REQUIRE( b.get_data( "z" ) );
		REQUIRE( b.get_data( "id" ) );
		REQUIRE( b.get_data( "type" ) );

		const std::vector<double> &x = data_as<double>( b.get_data("x") );
		const std::vector<double> &y = data_as<double>( b.get_data("y") );
		const std::vector<double> &z = data_as<double>( b.get_data("z") );
		const std::vector<int> &id   = data_as<int>( b.get_data("id") );
		const std::vector<int> &type = data_as<int>( b.get_data("type") );

		const std::vector<double> &pe = data_as<double>( b.get_data("c_pe") );

		if( i == 0 ){
			REQUIRE( b.tstep == 0 );

			REQUIRE( id[6] == 7 );
			REQUIRE( type[6] == 1 );
			REQUIRE( x[6] == Approx(  2.51939  ) );
			REQUIRE( y[6] == Approx(  0 ) );
			REQUIRE( z[6] == Approx( 0.839798  ) );
			REQUIRE( pe[6] == Approx( -6.77377 ).epsilon(1e-4) );
		} else if( i == 1 ){
			REQUIRE( b.tstep == 50 );

			REQUIRE( id[2888] == 3066 );
			REQUIRE( type[2888] == 1 );
			REQUIRE( x[2888] == Approx( 10.7546 ) );
			REQUIRE( y[2888] == Approx( 11.0767 ) );
			REQUIRE( z[2888] == Approx( 11.6878 ) );

			REQUIRE( pe[2888] == Approx( -5.25933 ).epsilon(1e-4) );
		}
	}
}


TEST_CASE ( "LAMMPS binary dump file gets read correctly.", "[read_lammps_dump_bin]" )
{

	using namespace lammps_tools;
	using namespace readers;

	std::string fname = "lammps_dump_file_test.dump.bin";
	std::vector<std::string> headers = { "id", "type", "x", "y", "z", "c_pe" };
	std::shared_ptr<dump_reader> r( make_dump_reader_lammps( fname, FILE_FORMAT_BIN, headers ) );
	dump_reader_lammps *rl = static_cast<dump_reader_lammps*>(r.get());
	rl->set_column_header_as_special( "id", block_data::ID );
	rl->set_column_header_as_special( "type", block_data::TYPE );
	rl->set_column_header_as_special( "x", block_data::X );
	rl->set_column_header_as_special( "y", block_data::Y );
	rl->set_column_header_as_special( "z", block_data::Z );

	REQUIRE( r );
	REQUIRE( r->good() );

	double lo = 0.0;
	double hi = 1.6795961913825074e+01;

	block_data b;

	for( int i = 0; i < 2; ++i ){
		int status = r->next_block( b );
		REQUIRE( status == 0 );
		REQUIRE( b.N == 4000 );

		REQUIRE( b.dom.xlo[0] == Approx(lo) );
		REQUIRE( b.dom.xlo[1] == Approx(lo) );
		REQUIRE( b.dom.xlo[2] == Approx(lo) );

		REQUIRE( b.dom.xhi[0] == Approx(hi) );
		REQUIRE( b.dom.xhi[1] == Approx(hi) );
		REQUIRE( b.dom.xhi[2] == Approx(hi) );

		REQUIRE( b.get_data( "x" ) );
		REQUIRE( b.get_data( "y" ) );
		REQUIRE( b.get_data( "z" ) );
		REQUIRE( b.get_data( "id" ) );
		REQUIRE( b.get_data( "type" ) );

		const std::vector<double> &x = data_as<double>( b.get_data("x") );
		const std::vector<double> &y = data_as<double>( b.get_data("y") );
		const std::vector<double> &z = data_as<double>( b.get_data("z") );
		const std::vector<int> &id   = data_as<int>( b.get_data("id") );
		const std::vector<int> &type = data_as<int>( b.get_data("type") );

		const std::vector<double> &pe = data_as<double>( b.get_data("c_pe") );

		if( i == 0 ){
			REQUIRE( b.tstep == 0 );

			REQUIRE( id[6] == 7 );
			REQUIRE( type[6] == 1 );
			REQUIRE( x[6] == Approx(  2.51939  ) );
			REQUIRE( y[6] == Approx(  0 ) );
			REQUIRE( z[6] == Approx( 0.839798  ) );
			REQUIRE( pe[6] == Approx( -6.77377 ).epsilon(1e-4) );
		} else if( i == 1 ){
			REQUIRE( b.tstep == 50 );

			REQUIRE( id[2888] == 3066 );
			REQUIRE( type[2888] == 1 );
			REQUIRE( x[2888] == Approx( 10.7546 ) );
			REQUIRE( y[2888] == Approx( 11.0767 ) );
			REQUIRE( z[2888] == Approx( 11.6878 ) );

			REQUIRE( pe[2888] == Approx( -5.25933 ).epsilon(1e-4) );
		}
	}
}



TEST_CASE ( "LAMMPS gzipped text dump file gets read correctly.", "[read_lammps_dump_gzip]" )
{
	using namespace lammps_tools;
	using namespace readers;

	std::string fname = "lammps_dump_file_test.dump.gz";
	// std::vector<std::string> headers = { "id", "type", "x", "y", "z", "c_pe" };
	std::shared_ptr<dump_reader> r( make_dump_reader( fname, FILE_FORMAT_GZIP, DUMP_FORMAT_LAMMPS ) );

	REQUIRE( r );
	REQUIRE( r->good() );

	double lo = 0.0;
	double hi = 1.6795961913825074e+01;

	block_data b;

	for( int i = 0; i < 2; ++i ){
		int status = r->next_block( b );
		REQUIRE( status == 0 );
		REQUIRE( b.N == 4000 );

		REQUIRE( b.dom.xlo[0] == Approx(lo) );
		REQUIRE( b.dom.xlo[1] == Approx(lo) );
		REQUIRE( b.dom.xlo[2] == Approx(lo) );

		REQUIRE( b.dom.xhi[0] == Approx(hi) );
		REQUIRE( b.dom.xhi[1] == Approx(hi) );
		REQUIRE( b.dom.xhi[2] == Approx(hi) );

		REQUIRE( b.get_data( "x" ) );
		REQUIRE( b.get_data( "y" ) );
		REQUIRE( b.get_data( "z" ) );
		REQUIRE( b.get_data( "id" ) );
		REQUIRE( b.get_data( "type" ) );

		const std::vector<double> &x = data_as<double>( b.get_data("x") );
		const std::vector<double> &y = data_as<double>( b.get_data("y") );
		const std::vector<double> &z = data_as<double>( b.get_data("z") );
		const std::vector<int> &id   = data_as<int>( b.get_data("id") );
		const std::vector<int> &type = data_as<int>( b.get_data("type") );

		const std::vector<double> &pe = data_as<double>( b.get_data("c_pe") );

		if( i == 0 ){
			REQUIRE( b.tstep == 0 );

			REQUIRE( id[6] == 7 );
			REQUIRE( type[6] == 1 );
			REQUIRE( x[6] == Approx(  2.51939  ) );
			REQUIRE( y[6] == Approx(  0 ) );
			REQUIRE( z[6] == Approx( 0.839798  ) );
			REQUIRE( pe[6] == Approx( -6.77377 ).epsilon(1e-4) );
		} else if( i == 1 ){
			REQUIRE( b.tstep == 50 );

			REQUIRE( id[2888] == 3066 );
			REQUIRE( type[2888] == 1 );
			REQUIRE( x[2888] == Approx( 10.7546 ) );
			REQUIRE( y[2888] == Approx( 11.0767 ) );
			REQUIRE( z[2888] == Approx( 11.6878 ) );

			REQUIRE( pe[2888] == Approx( -5.25933 ).epsilon(1e-4) );
		}
	}
}


TEST_CASE ( "Dump readers count correct number of blocks.", "[dump_reader_number_of_blocks]" )
{
	using namespace lammps_tools;
	using namespace readers;

	std::vector<std::string> dumps = { "lammps_dump_file_test.dump",
	                                   "lammps_dump_file_test.dump.bin",
	                                   "lammps_dump_file_test.dump.gz" };
	std::vector<int> dformats   = { DUMP_FORMAT_LAMMPS, DUMP_FORMAT_LAMMPS, DUMP_FORMAT_LAMMPS };
	std::vector<int> fformats   = { FILE_FORMAT_PLAIN, FILE_FORMAT_BIN, FILE_FORMAT_GZIP };
	std::vector<int> dump_sizes = { 51, 51, 51 };
	std::vector<std::string> cols = { "id", "type", "x", "y", "z", "c_pe" };

	// Verify that any extensions to this test are done properly:
	REQUIRE( dump_sizes.size() == dumps.size() );
	REQUIRE( dformats.size() == dumps.size() );
	REQUIRE( fformats.size() == dumps.size() );

	for( std::size_t i = 0; i < dumps.size(); ++i ){
		std::string dname = dumps[i];
		int size = dump_sizes[i];
		int dformat = dformats[i];
		int fformat = fformats[i];

		dump_reader *d = make_dump_reader( dname, fformat, dformat );
		std::cerr << "Reading file " << dname << ".\n";
		if( util::str_contains( dname, "lammps_dump_" ) ){
			dump_reader_lammps* ptr;
			ptr = static_cast<dump_reader_lammps*>( d );
			REQUIRE( ptr != nullptr );
			ptr->set_column_headers( cols );
		}

		REQUIRE( number_of_blocks( *d ) == size );

		delete d;
	}
}



TEST_CASE ( "Setting data type explicitly works.", "[read_lammps_dump_explicit_data_type]" )
{
	using namespace lammps_tools;
	using namespace readers;

	std::string fname = "lammps_dump_file_test.dump.gz";

	std::shared_ptr<dump_reader_lammps> r( make_dump_reader_lammps( fname, FILE_FORMAT_GZIP ) );

	std::vector<std::string> headers = { "id", "type", "x", "y", "z", "c_pe" };
	std::vector<int> special_fields  = { block_data::ID, block_data::TYPE, block_data::X,
	                                     block_data::Y, block_data::Z, block_data::UNKNOWN };
	r->set_column_headers( headers );
	auto set_header_field = [&r]( const std::string &h, int t ) {
		if( t == block_data::UNKNOWN ){
			std::cerr << "Warning: Ignoring special field UNKNOWN!\n";
			return;
		}
		std::cerr << "Setting " << h << " to " << t << "\n";
		if( !r->set_column_header_as_special( h, t ) ){
			my_runtime_error(__FILE__,__LINE__,
			                 "Failed to find header!" );
		}
	};
	util::zip_map( headers, special_fields, set_header_field );
	std::cerr << "Set all fields.\n";

	REQUIRE( r );
	REQUIRE( r->good() );

	double lo = 0.0;
	double hi = 1.6795961913825074e+01;

	block_data b;

	int status = r->next_block( b );
	REQUIRE( status == 0 );
	REQUIRE( b.N == 4000 );

	REQUIRE( b.dom.xlo[0] == Approx(lo) );
	REQUIRE( b.dom.xlo[1] == Approx(lo) );
	REQUIRE( b.dom.xlo[2] == Approx(lo) );

	REQUIRE( b.dom.xhi[0] == Approx(hi) );
	REQUIRE( b.dom.xhi[1] == Approx(hi) );
	REQUIRE( b.dom.xhi[2] == Approx(hi) );

	REQUIRE( b.get_data( "x" ) );
	REQUIRE( b.get_data( "y" ) );
	REQUIRE( b.get_data( "z" ) );
	REQUIRE( b.get_data( "id" ) );
	REQUIRE( b.get_data( "type" ) );

	const std::vector<double> &x = data_as<double>( b.get_special_field( block_data::X ) );
	const std::vector<int> &id = data_as<int>( b.get_special_field( block_data::ID ) );
	id_map im1( id );
	REQUIRE( x[im1[127]] == Approx(2.51939) );

	r->set_column_type( "x", data_field::INT );
	status = r->next_block( b );

	const std::vector<int> &id2 = data_as<int>( b.get_special_field( block_data::ID ) );
	id_map im2( id2 );

	const std::vector<int> &x2 = data_as<int>( b.get_special_field( block_data::X ) );
	REQUIRE( x2[im2[127]] == 2 );
}
