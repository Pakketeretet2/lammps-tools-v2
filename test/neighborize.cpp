#include <catch.hpp>

#include "block_data.hpp"
#include "dump_reader_lammps.hpp"
#include "my_timer.hpp"
#include "neighborize.hpp"

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
				                             atom_headers,
				                             0, 0, DIST_NSQ,
				                             3, rrc );
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
				                             atom_headers,
				                             0, 0, DIST_BIN,
				                             3, rrc );
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

	std::vector<std::string> atom_headers = { "id", "type", "x", "y", "z" };
	std::vector<std::vector<int> > neighs;

	double avg = make_list_dist( neighs, b, atom_headers, 0, 0, DIST_BIN, 3, rc );
}
