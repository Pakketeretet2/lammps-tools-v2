#include "readers.hpp"

#include <algorithm>
#include <catch.hpp>
#include <fstream>

TEST_CASE ( "LAMMPS data file gets read correctly.", "[read_lammps_data]" ) {
	std::ifstream in( "lammps_data_file_test.data" );
	REQUIRE( in );
	int status = -1;

	block_data b = data_readers::block_data_from_lammps_data( in, status );
	REQUIRE( b.get_data_names().size() == 8 );
	// These entries should _not_ be found:
	REQUIRE( b.get_data( "xlo" ) == nullptr );
	REQUIRE( b.get_data( "xhi" ) == nullptr );
	REQUIRE( b.get_data( "atom_style" ) == nullptr );
	REQUIRE( b.get_data( "timestep" ) == nullptr );
	
	
	// These entries should all be found:
	REQUIRE( b.get_data( "id" ) != nullptr );
	REQUIRE( b.get_data( "type" ) != nullptr );
	REQUIRE( b.get_data( "x" ) != nullptr );
	REQUIRE( b.get_data( "y" ) != nullptr );
	REQUIRE( b.get_data( "z" ) != nullptr );
	REQUIRE( b.get_data( "vx" ) != nullptr );
	REQUIRE( b.get_data( "vy" ) != nullptr );
	REQUIRE( b.get_data( "vz" ) != nullptr );
}



TEST_CASE ( "LAMMPS plain text dump file gets read correctly.", "[read_lammps_dump_plain]" ) {

	using namespace dump_readers;
	using dfd = data_field_double;
	using dfi = data_field_int;

	std::string fname = "lammps_dump_file_test.dump";
	// std::vector<std::string> headers = { "id", "type", "x", "y", "z", "c_pe" };
	std::shared_ptr<dump_reader> r( make_dump_reader( fname, LAMMPS, PLAIN ) );
	
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


TEST_CASE ( "LAMMPS binary dump file gets read correctly.", "[read_lammps_dump_bin]" ) {

	using namespace dump_readers;
	using dfd = data_field_double;
	using dfi = data_field_int;
	
	std::string fname = "lammps_dump_file_test.dump.bin";
	std::vector<std::string> headers = { "id", "type", "x", "y", "z", "c_pe" };
	std::shared_ptr<dump_reader> r( make_dump_reader_lammps( fname, BIN, headers ) );
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



TEST_CASE ( "LAMMPS gzipped text dump file gets read correctly.", "[read_lammps_dump_gzip]" ) {

	using namespace dump_readers;
	using dfd = data_field_double;
	using dfi = data_field_int;

	std::string fname = "lammps_dump_file_test.dump.gz";
	// std::vector<std::string> headers = { "id", "type", "x", "y", "z", "c_pe" };
	std::shared_ptr<dump_reader> r( make_dump_reader( fname, LAMMPS, GZIP ) );
	
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
