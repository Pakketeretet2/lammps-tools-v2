#include <catch.hpp>

#include "block_data.hpp"
#include "block_data_access.hpp"
#include "dump_reader_lammps_bin.hpp"
#include "geometry.hpp"
#include "transformations.hpp"


TEST_CASE( "Tests if atom shifting works.", "[transformation_shift]" ){
	using lammps_tools::domain;
	using lammps_tools::block_data;
	using lammps_tools::readers::dump_reader_lammps_bin;

	dump_reader_lammps_bin d( "lammps_dump_file_test.dump.bin" );
	d.set_column_headers( {"id", "type", "x", "y", "z", "c_pe"} );
	block_data b;
	int status = d.next_block(b);
	REQUIRE( status == 0 );

	const double *xlo = b.dom.xlo;
	const double *xhi = b.dom.xhi;
	double L[3];
	L[0] = xhi[0] - xlo[0];
	L[1] = xhi[1] - xlo[1];
	L[2] = xhi[2] - xlo[2];

	int target = 1234;
	REQUIRE( target < b.N );
	double dx = L[0];
	double dy = -L[1];
	double dz = L[2];

	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);

	std::vector<int> ids(1);
	ids[0] = get_id(b)[target];
	lammps_tools::point delta( dx, dy, dz );
	std::cerr << "Shifting points...\n";
	lammps_tools::transformations::shift( b, delta, ids );

}
