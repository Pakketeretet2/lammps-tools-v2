#include <catch.hpp>

#include "block_data.hpp"
#include "cluster_finder.hpp"
#include "dump_reader_lammps.hpp"
#include "my_timer.hpp"
#include "neighborize.hpp"
#include "util.hpp"
#include "writers.hpp"

#include <list>
#include <fstream>
#include <vector>


TEST_CASE( "Cluster analysis on triangles", "[cluster_triangles]" ) {
	using namespace lammps_tools;
	using namespace lammps_tools::readers;
	using namespace lammps_tools::neighborize;

	std::string dname = "icosahedron_0.dump.bin";

	std::vector<std::string> common_headers = { "id", "type",
	                                            "x", "y", "z" };

	std::vector<std::string> all_headers = { "id", "mol", "type",
	                                         "x", "y", "z" };

	std::unique_ptr<dump_reader_lammps> d(
		make_dump_reader_lammps( dname, 2 ) );

	d->set_column_headers( all_headers );

	d->set_column_header_as_special(   "id", block_data::ID );
	d->set_column_header_as_special(  "mol", block_data::MOL );
	d->set_column_header_as_special( "type", block_data::TYPE );
	d->set_column_header_as_special(    "x", block_data::X );
	d->set_column_header_as_special(    "y", block_data::Y );
	d->set_column_header_as_special(    "z", block_data::Z );

	std::vector<block_data> blocks;
	block_data tmp;
	int c = 0;
	my_timer t(std::cerr);
	t.tic();
	// int n_blocks = 20;
	bigint last_time = 0;
	while( d->next_block(tmp) == 0 ){
		blocks.push_back(tmp);
		++c;
		int n_skip = 25;
		if( c % n_skip == 0 ){
			std::string msg = "Read ";
			msg += std::to_string(n_skip);
			msg += " blocks of ";
			msg += std::to_string( tmp.N );
			msg += " atoms each";
			t.toc(msg);
			t.tic();
			std::cerr << "At " << c << " blocks...\n";
		}
		last_time = tmp.tstep;
	}

	double rc = 2.0;
	int dims = 3;
	int mol_policy = neighborizer::INCLUDE;
	int bond_policy = neighborizer::IGNORE;
	int n_block = 0;

	std::ofstream net_out( "networks.dat" );
	std::ofstream out( "icosahedron_0_colours.dump.bin", std::ios::binary );
	for( block_data &b : blocks ){
		t.tic();
		++n_block;

		const std::vector<int> &type = data_as<int>(
			b.get_special_field( block_data::TYPE ) );
		const std::vector<int> &mol = data_as<int>(
			b.get_special_field( block_data::MOL ) );

		neigh_list neighs;
		double avg_neighs = make_list_dist( neighs, b,
		                                    2, 3, DIST_BIN, dims, rc,
		                                    mol_policy, bond_policy );
		t.toc("Neighbourising");

		t.tic();
		neigh_list conns;
		std::list<std::list<int> > networks;

		find_molecular_networks ( b, neighs, conns, networks );
		t.toc("Finding networks");

		for( int i = 0; i < b.N; ++i ){

			if( b.tstep == 0 ){
				// 153 in a molecule so 152 neighbours.
				if( type[i] == 1 ){
					REQUIRE( neighs[i].size() == 152 );
				}
				if( type[i] == 2 ){
					REQUIRE( neighs[i].size() == 152 );
				}
				if( type[i] == 3 ){
					REQUIRE( neighs[i].size() == 152 );
				}
			}
		}


		t.tic();
		std::vector<int> id_to_color( b.N );
		data_field_int colour( "cluster_size", b.N );

		for( int i = 0; i < b.N; ++i ){
			colour[i] = 0;
			for( const std::list<int> &n : networks ){
				if( util::contains( n, mol[i] ) ){
					colour[i] = n.size();
					break;
				}
			}
		}
		b.add_field( colour );
		t.toc("Colouring");

		t.tic();
		std::vector<int> histogram( 30 );
		for( const std::list<int> &n : networks ){
			int s = n.size();
			if( s > 29 ) continue;

			histogram[s]++;
		}
		int s = 0;
		for( int c : histogram ){
			net_out << b.tstep << " " << s << " " << c
			        << " " << c*s << "\n";
			++s;
		}
		net_out << "\n";
		t.toc("Dumping to histogram");

		t.tic();
		writers::block_to_lammps_dump_bin( out, b );
		t.toc("Write to file");

	}

}
