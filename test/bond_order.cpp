#include <catch.hpp>

#include "block_data.hpp"
#include "block_data_access.hpp"
#include "bond_order.hpp"
#include "constants.hpp"
#include "data_field.hpp"
#include "fast_math.hpp"
#include "geometry.hpp"
#include "topology.hpp"

TEST_CASE( "Bond order angle calculation makes sense", "[bond_order_angle]" ) {

	std::vector<double> e_x = { 1.0, 0.0, 0.0 };
	std::vector<double> e_y = { 0.0, 1.0, 0.0 };

	int Nangles = 720;
	double dt = lammps_tools::constants::pi2 / Nangles;
	int Nparticles = Nangles + 1;
	std::vector<double> x(Nparticles), y(Nparticles), z( Nparticles, 0.0 );
	double r0 = 1.0;
	double theta = 0.0;
	std::vector<int> id(Nparticles), type(Nparticles, 1);
	x[0] = 0.0;
	y[0] = 0.0;
	z[0] = 0.0;
	id[0] = 1;
	type[0] = 2;
	for( int i = 1; i < Nparticles; ++i ){
		double c, s;
		s = std::sin(theta);
		c = std::cos(theta);

		x[i] = r0 * c;
		y[i] = r0 * s;
		id[i] = i+1;

		lammps_tools::point v1( x[i], y[i], 0.0 );
		lammps_tools::point v2( e_x[0], e_x[1], e_x[2] );
		double a_test = lammps_tools::order_parameters::angle_2pi( v2, v1 );
		REQUIRE( a_test == Approx( theta ) );

		theta += dt;
	}

	lammps_tools::block_data b( Nparticles );
	lammps_tools::data_field_int d_id( "id", id );
	lammps_tools::data_field_int d_type( "type", type );
	lammps_tools::data_field_double d_x( "x", x );
	lammps_tools::data_field_double d_y( "y", y );
	lammps_tools::data_field_double d_z( "z", z );

	b.add_field(   d_id, lammps_tools::block_data::ID );
	b.add_field( d_type, lammps_tools::block_data::TYPE );
	b.add_field(    d_x, lammps_tools::block_data::X );
	b.add_field(    d_y, lammps_tools::block_data::Y );
	b.add_field(    d_z, lammps_tools::block_data::Z );

	lammps_tools::neighborize::neigh_list neighs;
	double avg = make_list_dist( neighs, b, 2, 1,
	                             lammps_tools::neighborize::DIST_NSQ, 2,
	                             1.35, 0, 0, false );

	std::vector<lammps_tools::bond> bonds =
		lammps_tools::neighborize::neigh_list_to_bonds( b, neighs, 1 );

	std::vector<double> angles;
	lammps_tools::order_parameters::relative_bond_angles( b, neighs, e_x,
	                                                      bonds, angles );
	for( int i = 0; i < Nangles; ++i ){
		std::cerr << "Angle " << i << " = " << angles[i] << "\n";
	}


}
