#include "block_data.hpp"
#include "block_data_access.hpp"
#include "domain.hpp"
#include "writers.hpp"

#include <catch2/catch.hpp>

#include <memory>

TEST_CASE ( "Test if image flag remapping works..", "[reconstruct_image_flag]" )
{
	using namespace lammps_tools;
	
	int n_atoms = 6;
	block_data b(n_atoms);
	b.dom.xlo[0] = b.dom.xlo[1] = b.dom.xlo[2] = -4;
	b.dom.xhi[0] = b.dom.xhi[1] = b.dom.xhi[2] =  4;

	b.dom.periodic = 7;

	{
		// Local block so I can reuse the same names...
		std::vector<double> x =    { 3.5, 4.0,  4.5, -3.0,  -2.5, -2.0 };
		std::vector<double> y =    { 2.0, 2.25, 2.5, 2.75, 3.0,  3.25 };
		std::vector<double> z =    { 1.0, 1.0,  1.0, 1.0,  1.0,  1.0 };
		std::vector<int>    id =   { 1,  2,    3,   4,    5,    6 };
		std::vector<int>    type = { 1,  1,    2,   2,    1,    1 };
		std::vector<int>    mol =  { 1,  1,    1,   1,    1,    1 };

		data_field_double dx("x", x);
		data_field_double dy("y", y);
		data_field_double dz("z", z);

		data_field_int did("id", id);
		data_field_int dmol("mol", mol);
		data_field_int dtype("type", type);
		
		b.add_field(did, block_data::ID);
		b.add_field(dmol, block_data::MOL);
		b.add_field(dtype, block_data::TYPE);
		b.add_field(dx, block_data::X);
		b.add_field(dy, block_data::Y);
		b.add_field(dz, block_data::Z);
	}

	std::vector<int> image_x(b.N);
	std::vector<int> image_y(b.N);
	std::vector<int> image_z(b.N);

	b.dom.reconstruct_image_flags(b, image_x, image_y, image_z);

	REQUIRE(image_x[0] == 0);
	REQUIRE(image_x[1] == 0);
	REQUIRE(image_x[2] == 0);
	REQUIRE(image_x[3] == 1);
	REQUIRE(image_x[4] == 1);
	REQUIRE(image_x[5] == 1);

	REQUIRE(image_y[0] == 0);
	REQUIRE(image_y[1] == 0);
	REQUIRE(image_y[2] == 0);
	REQUIRE(image_y[3] == 0);
	REQUIRE(image_y[4] == 0);
	REQUIRE(image_y[5] == 0);

	REQUIRE(image_z[0] == 0);
	REQUIRE(image_z[1] == 0);
	REQUIRE(image_z[2] == 0);
	REQUIRE(image_z[3] == 0);
	REQUIRE(image_z[4] == 0);
	REQUIRE(image_z[5] == 0);

	b.set_ntypes(2);
	b.add_field(data_field_int("ix", image_x), block_data::IX);
	b.add_field(data_field_int("iy", image_y), block_data::IY);
	b.add_field(data_field_int("iz", image_z), block_data::IZ);

	writers::block_to_lammps_data("recreate_img_flags.data", b);
	
}
