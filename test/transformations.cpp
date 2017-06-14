#include <catch.hpp>
#include <string>
#include <vector>

#include "block_data.hpp"
#include "block_data_access.hpp"
#include "data_field.hpp"
#include "dump_reader_lammps_bin.hpp"
#include "geometry.hpp"
#include "transformations.hpp"
#include "writers.hpp"


TEST_CASE( "Tests if transformations work.", "[transformations]" ){
	using namespace lammps_tools;

	using transformations::shift;
	using transformations::shift_all;

	using readers::dump_reader_lammps_bin;

	dump_reader_lammps_bin d( "lammps_dump_file_test.dump.bin" );
	std::vector<std::string> headers    = {"id", "type", "x",
	                                       "y", "z", "c_pe"};
	std::vector<int> special_types = { block_data::ID,
	                                   block_data::TYPE,
	                                   block_data::X,
	                                   block_data::Y,
	                                   block_data::Z };
	d.set_column_headers( headers );
	for( std::size_t i = 0; i < special_types.size(); ++i ){
		d.set_column_header_as_special( headers[i], special_types[i] );
	}


	block_data b;
	int status = d.next_block(b);
	REQUIRE( status == 0 );
	std::vector<int> image_flags( b.N, 0 );
	b.add_field( data_field_int( "ix", image_flags ), block_data::IX );
	b.add_field( data_field_int( "iy", image_flags ), block_data::IY );
	b.add_field( data_field_int( "iz", image_flags ), block_data::IZ );


	const double *xlo = b.dom.xlo;
	const double *xhi = b.dom.xhi;
	double L[3];
	L[0] = xhi[0] - xlo[0];
	L[1] = xhi[1] - xlo[1];
	L[2] = xhi[2] - xlo[2];

	SECTION( "Shift one" ){
		int target = 1234;
		REQUIRE( target < b.N - 1 );

		const std::vector<double> &x = get_x(b);
		const std::vector<double> &y = get_y(b);
		const std::vector<double> &z = get_z(b);

		const std::vector<int> &ix  = data_as<int>(
			b.get_special_field( block_data::IX ) );
		const std::vector<int> &iy  = data_as<int>(
			b.get_special_field( block_data::IY ) );
		const std::vector<int> &iz  = data_as<int>(
			b.get_special_field( block_data::IZ ) );


		double dx = L[0];
		double dy = -L[1];
		double dz = L[2];

		std::vector<int> ids(1);
		ids[0] = get_id(b)[target];
		lammps_tools::point delta( dx, dy, dz );
		std::cerr << "Shifting points...\n";
		block_data b_shift = shift( b, delta, ids );


		const std::vector<int> &ixn = data_as<int>(
			b_shift.get_special_field( block_data::IX ) );
		const std::vector<int> &iyn = data_as<int>(
			b_shift.get_special_field( block_data::IY ) );
		const std::vector<int> &izn = data_as<int>(
			b_shift.get_special_field( block_data::IZ ) );



		// Particle should still be in box:
		double xn = get_x(b_shift)[target];
		double yn = get_y(b_shift)[target];
		double zn = get_z(b_shift)[target];

		double tmp[3];
		double xt[3];
		double xtn[3];

		xt[0] = x[target];
		xt[1] = y[target];
		xt[2] = z[target];

		xtn[0] = xn;
		xtn[1] = yn;
		xtn[2] = zn;

		// Shift back to same position.
		REQUIRE( b.dom.dist_2( xtn, xt, tmp ) == Approx(0.0) );
		// Image flags should be different.
		REQUIRE( ixn[target] == ix[target] + 1 );
		REQUIRE( iyn[target] == ix[target] - 1 );
		REQUIRE( izn[target] == ix[target] + 1 );

		REQUIRE( ixn[target+1] == ix[target+1] );
		REQUIRE( iyn[target+1] == ix[target+1] );
		REQUIRE( izn[target+1] == ix[target+1] );

	}


	SECTION( "Shift many" ){
		point delta( 3.4, 5.2, -4.1 );
		std::cerr << "Shifting points...\n";
		block_data b_shift = shift_all( b, delta );


		const std::vector<double> &x = get_x(b);
		const std::vector<double> &y = get_y(b);
		const std::vector<double> &z = get_z(b);


		const std::vector<double> &xn = get_x(b_shift);
		const std::vector<double> &yn = get_y(b_shift);
		const std::vector<double> &zn = get_z(b_shift);
		std::cerr << "Calculating " << b.N
		          << "^2 distances, this can take a while...\n";
		int done = 0;
		int total = b.N*(b.N-1) / 2.0;

		double *xhi = b.dom.xhi;
		double *xlo = b.dom.xlo;

		for( int i = 0; i < b.N; ++i ){
			double xi [3] = {  x[i],  y[i],  z[i] };
			double xni[3] = { xn[i], yn[i], zn[i] };
			REQUIRE( xi[0] != xni[0] );
			REQUIRE( xi[1] != xni[1] );
			REQUIRE( xi[2] != xni[2] );
			// Check the rewrap:
			REQUIRE( xni[0] <  xhi[0] );
			REQUIRE( xni[0] >= xlo[0] );
			REQUIRE( xni[1] <  xhi[1] );
			REQUIRE( xni[1] >= xlo[1] );
			REQUIRE( xni[2] <  xhi[2] );
			REQUIRE( xni[2] >= xlo[2] );



			for( int j = 0; j < i; ++j ){
				double xj [3] = {  x[j],  y[j],  z[j] };
				double xnj[3] = { xn[j], yn[j], zn[j] };
				double tmp[3];
				double r2_old = b.dom.dist_2( xi, xj, tmp );
				double r2_new = b.dom.dist_2( xni, xnj, tmp );
				REQUIRE( r2_old == Approx( r2_new ) );
				++done;
				if( done % 500000 == 0 ){
					std::cerr << "  Done with " << done
					          << " of of " << total
					          << " distances...\n";
				}
			}
		}
		std::cerr << "Done!\n";

		std::ofstream shift_out( "lammps_dump_file_test_shift_out.dump.bin",
		                         std::ios::binary );
		writers::block_to_lammps_dump( shift_out, b_shift,
		                               FILE_FORMAT_BIN );
	}

	SECTION( "Tests if rotating works." ){
		using transformations::rotate_all;

		point origin( 0.5*( b.dom.xhi[0] + b.dom.xlo[0] ),
		              0.5*( b.dom.xhi[1] + b.dom.xlo[1] ),
		              0.5*( b.dom.xhi[2] + b.dom.xlo[2] ) );
		double quart_pi = 3.1415927 / 4.0;
		point axis( 0, 0, 1.0 );
		block_data b_center = transformations::rotate_all( b, axis, origin, quart_pi );

		std::ofstream center_out( "lammps_dump_file_test_rotate_out.dump.bin",
		                          std::ios::binary );
		writers::block_to_lammps_dump( center_out, b_center, FILE_FORMAT_BIN );
	}

	SECTION( "Tests if centering works." ){
		using namespace lammps_tools;

		using transformations::shift;
		using transformations::shift_all;


		point origin( 0.5*( b.dom.xhi[0] + b.dom.xlo[0] ),
		              0.5*( b.dom.xhi[1] + b.dom.xlo[1] ),
		              0.5*( b.dom.xhi[2] + b.dom.xlo[2] ) );
		block_data b_center = transformations::center_box_on( b, origin );

		std::ofstream center_out( "lammps_dump_file_test_center_out.dump.bin",
		                          std::ios::binary );
		writers::block_to_lammps_dump( center_out, b_center, FILE_FORMAT_BIN );
	}
}
