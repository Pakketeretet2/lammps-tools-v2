#include <catch.hpp>

#include "block_data.hpp"
#include "dump_reader_lammps.hpp"
#include "my_timer.hpp"
#include "neighborize.hpp"
#include "util.hpp"
#include "writers.hpp"

#include <string>

TEST_CASE( "Neighbour list works as expected", "[neigh_list_dist]" ) {
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

	double rc[5] = { 1.35, 2.0, 3.0, 3.5, 5.0 };
	double count[3][5] = { { 12.000, 18.000, 86.000, 140.000, 428.000 },
	                       {  9.2385, 26.706, 94.452, 150.528, 441.148 },
	                       {  9.234, 26.738, 94.478, 150.534, 441.133 } };

	int neigh_est = 0;
	std::vector<block_data> blocks;
	block_data b;
	int frame = 0;
	while( d->next_block(b) == 0 ){
		if( frame % 25 == 0 ){
			blocks.push_back(b);
		}
		++frame;
	}
	int n_rcs = 5;
	REQUIRE( blocks.size() == 3 );

	double eps = 1e-3;
	bool nsq = true;
	bool bin = true;

	if( nsq ){
		my_timer timer_inner( std::cerr );
		timer_inner.tic();
		int j = 0;
		for( const block_data &b : blocks ){
			for( int i = 0; i < n_rcs; ++i ){
				double rrc = rc[i];
				std::cerr << "At NSQ, rc = " << rrc << "\n";

				double avg = make_list_dist( neighs, b,
				                             0, 0, DIST_NSQ,
				                             3, rrc, 0, 0, true );
				REQUIRE( avg ==
			        	 Approx(count[j][i]).epsilon(eps) );

			}
			++j;
		}
		timer_inner.toc("  NSQ neighbouring");
	}
	if( bin ){
	        my_timer timer_inner( std::cerr );
	        timer_inner.tic();
	        int j = 0;
	        for( const block_data &b : blocks ){
			for( int i = 0; i < n_rcs; ++i ){
				double rrc = rc[i];
				std::cerr << "At BIN, rc = " << rrc << "\n";
				timer_inner.tic();
				double avg = make_list_dist( neighs, b,
				                             0, 0, DIST_BIN,
				                             3, rrc, 0, 0, true );
				REQUIRE( avg ==
				         Approx(count[j][i]).epsilon(eps) );
			}
			++j;
		}
		timer_inner.toc("  BIN neighbouring");
	}
}

TEST_CASE( "Binning neighbour list works as expected", "[neigh_list_dist_bin]" ) {
	using namespace lammps_tools;
	using namespace lammps_tools::readers;
	using namespace lammps_tools::neighborize;

	int N = 8;
	data_field_int id( "id", N );
	data_field_int type( "type", N );
	data_field_double x( "x", N );
	data_field_double y( "y", N );
	data_field_double z( "z", N );
	block_data b(N);
	b.tstep = 0;
	b.N_types = 1;
	b.atom_style = block_data::ATOMIC;
	domain dom;
	dom.xlo[0] = dom.xlo[1] = dom.xlo[2] = 0.0;
	dom.xhi[0] = dom.xhi[1] = dom.xhi[2] = 9.0;
	b.dom = dom;

	double rc = 3.0;

	x[0] = 5.0;
	x[1] = 4.1;
	x[2] = 3.9;
	x[3] = 6.5;
	x[4] = 7.0;
	x[5] = 7.5;
	x[6] = 8.4;
	x[7] = 0.6;

	y[0] = 5.0;
	y[1] = 5.5;
	y[2] = 3.9;
	y[3] = 4.0;
	y[4] = 6.5;
	y[5] = 3.5;
	y[6] = 1.0;
	y[7] = 1.0;

	for( int i = 0; i < N; ++i ){
		z[i] = 0.0;
		id[i] = i+1;
		type[i] = 1;
	}

	b.add_field( id, block_data::ID );
	b.add_field( type, block_data::TYPE );
	b.add_field( x, block_data::X );
	b.add_field( y, block_data::Y );
	b.add_field( z, block_data::Z );

	std::vector<std::vector<int> > neighs;

	double avg = make_list_dist( neighs, b, 0, 0, DIST_BIN, 3, rc );
}



TEST_CASE( "Type discriminator works", "[neigh_list_types]" ) {
using namespace lammps_tools;
	using namespace lammps_tools::readers;
	using namespace lammps_tools::neighborize;

	int N = 8;
	data_field_int id( "id", N );
	data_field_int type( "type", N );
	data_field_double x( "x", N );
	data_field_double y( "y", N );
	data_field_double z( "z", N );
	block_data b(N);
	b.tstep = 0;
	b.N_types = 1;
	b.atom_style = block_data::ATOMIC;
	domain dom;
	dom.xlo[0] = dom.xlo[1] = dom.xlo[2] = 0.0;
	dom.xhi[0] = dom.xhi[1] = dom.xhi[2] = 9.0;
	dom.periodic = domain::BIT_X + domain::BIT_Y + domain::BIT_Z;
	b.dom = dom;

	double rc = 3.0;

	x[0] = 5.0;
	x[1] = 4.1;
	x[2] = 3.9;
	x[3] = 6.5;
	x[4] = 7.0;
	x[5] = 7.5;
	x[6] = 8.4;
	x[7] = 0.6;

	y[0] = 5.0;
	y[1] = 5.5;
	y[2] = 3.9;
	y[3] = 4.0;
	y[4] = 6.5;
	y[5] = 3.5;
	y[6] = 1.0;
	y[7] = 1.0;

	type[0] = 1;
	type[1] = 2;
	type[2] = 1;
	type[3] = 2;
	type[4] = 1;
	type[5] = 2;
	type[6] = 1;
	type[7] = 1;

	for( int i = 0; i < N; ++i ){
		z[i] = 0.0;
		id[i] = i+1;
	}

	b.add_field( id, block_data::ID );
	b.add_field( type, block_data::TYPE );
	b.add_field( x, block_data::X );
	b.add_field( y, block_data::Y );
	b.add_field( z, block_data::Z );

	lammps_tools::writers::block_to_lammps_dump( "neighs_input.dump", b, 0 );

	std::vector<std::vector<int> > neighs;

	std::vector<data_field_int> ncs;
	data_field_int nc1( "nc", N );
	data_field_int nc2( "nc_12", N );
	data_field_int nc3( "nc_11", N );
	data_field_int nc4( "nc_21", N );
	data_field_int nc5( "nc_22", N );
	ncs.push_back(nc1);
	ncs.push_back(nc2);
	ncs.push_back(nc3);
	ncs.push_back(nc4);
	ncs.push_back(nc5);

	std::vector<std::vector<int> > types = {{0,0},
	                                        {1,2},
	                                        {1,1},
	                                        {2,1},
	                                        {2,2}};

	std::vector<std::vector<std::vector<bool> > > are_neighs(5);
 	for( int i = 0; i < 5; ++i ){
		are_neighs[i].resize(8);
		for( int j = 0; j < 8; ++j ){
			are_neighs[i][j].resize(8);
		}
	}
	int i = 0;
	are_neighs[i][0][1] = true;
	are_neighs[i][0][2] = true;
	are_neighs[i][0][3] = true;
	are_neighs[i][0][4] = true;
	are_neighs[i][0][5] = true;

	are_neighs[i][1][0] = true;
	are_neighs[i][1][2] = true;
	are_neighs[i][1][3] = true;

	are_neighs[i][2][0] = true;
	are_neighs[i][2][1] = true;
	are_neighs[i][2][3] = true;

	are_neighs[i][3][0] = true;
	are_neighs[i][3][1] = true;
	are_neighs[i][3][2] = true;
	are_neighs[i][3][4] = true;
	are_neighs[i][3][5] = true;

	are_neighs[i][4][0] = true;
	are_neighs[i][4][3] = true;

	are_neighs[i][5][0] = true;
	are_neighs[i][5][3] = true;
	are_neighs[i][5][6] = true;

	are_neighs[i][6][5] = true;
	are_neighs[i][6][7] = true;

	are_neighs[i][7][6] = true;

	i = 1;
	are_neighs[i][0][1] = true;
	are_neighs[i][0][3] = true;
	are_neighs[i][0][5] = true;

	are_neighs[i][1][0] = true;
	are_neighs[i][1][2] = true;

	are_neighs[i][2][1] = true;
	are_neighs[i][2][3] = true;

	are_neighs[i][3][0] = true;
	are_neighs[i][3][2] = true;
	are_neighs[i][3][4] = true;

	are_neighs[i][4][3] = true;

	are_neighs[i][5][0] = true;
	are_neighs[i][5][6] = true;

	are_neighs[i][6][5] = true;

	i = 2;
	are_neighs[i][0][2] = true;
	are_neighs[i][0][4] = true;

	are_neighs[i][2][0] = true;

	are_neighs[i][4][0] = true;

	are_neighs[i][6][7] = true;

	are_neighs[i][7][6] = true;

	are_neighs[3] = are_neighs[1];

	i = 4;
	are_neighs[i][1][3] = true;

	are_neighs[i][3][1] = true;
	are_neighs[i][3][5] = true;

	are_neighs[i][5][3] = true;

	int method = DIST_BIN;
	// int method = DIST_NSQ;
	double avg = 0.0;
	for( int k = 0; k < 5; ++k ){
		int itype = types[k][0];
		int jtype = types[k][1];

		avg = make_list_dist( neighs, b, itype, jtype, method, 3, rc );
		for( int i = 0; i < b.N; ++i ){
			for( int j = 0; j < b.N; ++j ){
				if( util::contains( neighs[i], j ) ){
					REQUIRE(  are_neighs[k][i][j] );
				}else{
					REQUIRE( !are_neighs[k][i][j] );
				}
			}
			ncs[k][i] = neighs[i].size();
		}
		b.add_field( ncs[k] );
	}

	lammps_tools::writers::block_to_lammps_dump( "neighs_output.dump", b, 0 );

}



TEST_CASE( "Neighbour list can work with in-molecule connections", "[neigh_list_mol]" ) {
	using namespace lammps_tools;
	using namespace lammps_tools::readers;
	using namespace lammps_tools::neighborize;

	std::string dname = "lammps_dump_file_test_icosahedron.dump";

	std::vector<std::string> all_headers = { "id", "mol", "type",
	                                         "x", "y", "z" };

	std::unique_ptr<dump_reader_lammps> d(
		make_dump_reader_lammps( dname, 0 ) );

	d->set_column_headers( all_headers );

	d->set_column_header_as_special(   "id", block_data::ID );
	d->set_column_header_as_special(  "mol", block_data::MOL );
	d->set_column_header_as_special( "type", block_data::TYPE );
	d->set_column_header_as_special(    "x", block_data::X );
	d->set_column_header_as_special(    "y", block_data::Y );
	d->set_column_header_as_special(    "z", block_data::Z );

	std::ofstream out( "neigh_counts.dump.bin", std::ios::binary );
	block_data b;
	double rc = 2.0;
	int dims = 3;
	int mol_policy = neighborizer::INCLUDE;
	int bond_policy = neighborizer::IGNORE;

	while( d->next_block(b) == 0 ){
		const std::vector<int> &id =
			data_as<int>( b.get_special_field( block_data::ID ) );
		neigh_list neighs;
		double avg_neighs = make_list_dist( neighs, b, 2, 3, DIST_BIN,
		                                    dims, rc, mol_policy,
		                                    bond_policy, true );
		data_field_int ncs( "ncs", b.N );
		for( std::size_t i = 0; i < neighs.size(); ++i ){
			ncs[i] = neighs[i].size();
			if( neighs[i].size() == 0 ) continue;

		}
		b.add_field( ncs );
		lammps_tools::writers::block_to_lammps_dump( out, b, 2 );
	}

}
