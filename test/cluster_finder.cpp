#include <catch.hpp>

#include "block_data.hpp"
#include "cluster_finder.hpp"
#include "dump_reader_lammps.hpp"
#include "my_timer.hpp"
#include "neighborize.hpp"
#include "readers.hpp"
#include "util.hpp"
#include "writers.hpp"

#include <list>
#include <fstream>
#include <vector>



TEST_CASE( "Cluster analysis on triangles (smaller set)", "[cluster_triangles_smaller]" ) {
	using namespace lammps_tools;
	using namespace lammps_tools::readers;
	using namespace lammps_tools::neighborize;

	std::string dname = "triangle_neighs_test.data";
	int status;
	block_data b = block_data_from_lammps_data( dname, status );

	REQUIRE( status == 0 );

	neigh_list neighs;

	const std::vector<int> &id = data_as<int>(
		b.get_special_field( block_data::ID ) );
	const std::vector<int> &type = data_as<int>(
		b.get_special_field( block_data::TYPE ) );
	const std::vector<int> &mol = data_as<int>(
		b.get_special_field( block_data::MOL ) );

	const std::vector<double> &x = data_as<double>(
		b.get_special_field( block_data::X ) );
	const std::vector<double> &y = data_as<double>(
		b.get_special_field( block_data::Y ) );
	const std::vector<double> &z = data_as<double>(
		b.get_special_field( block_data::Z ) );


	double rc = 2.0;
	int dims = 3;
	int mol_policy = neighborizer::IGNORE;
	int bond_policy = neighborizer::IGNORE;

	double avg_neighs = make_list_dist( neighs, b,
	                                    2, 3, DIST_BIN, dims, rc,
	                                    mol_policy, bond_policy );
	id_map im( id );
	std::vector<int> patches = { 52094, 52083, 52077,
	                             52117, 52111, 52100 };
	for( int idi : patches ){
		int i = im[idi];
		REQUIRE( (type[i] == 2 || type[i] == 3) );
		REQUIRE( neighs[i].size() == 1 );
	}
}



TEST_CASE( "Cluster analysis on triangles", "[cluster_triangles]" ) {
	using namespace lammps_tools;
	using namespace lammps_tools::readers;
	using namespace lammps_tools::neighborize;

	std::string dname = "icosahedron_0_first_last.dump.bin";

	std::vector<std::string> common_headers = { "id", "type",
	                                            "x", "y", "z" };

	std::vector<std::string> all_headers = { "id", "mol", "type",
	                                         "x", "y", "z" };

	std::unique_ptr<dump_reader_lammps> d(
		make_dump_reader_lammps( dname, FILE_FORMAT_BIN ) );

	d->set_column_headers( all_headers );

	d->set_column_header_as_special(   "id", block_data::ID );
	d->set_column_header_as_special(  "mol", block_data::MOL );
	d->set_column_header_as_special( "type", block_data::TYPE );
	d->set_column_header_as_special(    "x", block_data::X );
	d->set_column_header_as_special(    "y", block_data::Y );
	d->set_column_header_as_special(    "z", block_data::Z );

	std::vector<block_data> blocks;
	block_data tmp;
	// int n_blocks = 20;
	while( d->next_block(tmp) == 0 ){
		blocks.push_back(tmp);
	}

	double rc = 2.0;
	int dims = 3;
	int mol_policy = neighborizer::INCLUDE;
	int bond_policy = neighborizer::IGNORE;
	int n_block = 0;

	for( block_data &b : blocks ){
		++n_block;

		const std::vector<int> &id = data_as<int>(
			b.get_special_field( block_data::ID ) );
		const std::vector<int> &type = data_as<int>(
			b.get_special_field( block_data::TYPE ) );
		const std::vector<int> &mol = data_as<int>(
			b.get_special_field( block_data::MOL ) );

		neigh_list neighs;
		double avg_neighs = make_list_dist( neighs, b,
		                                    2, 3, DIST_BIN, dims, rc,
		                                    mol_policy, bond_policy );
		neigh_list conns;
		std::list<std::list<int> > networks;

		find_molecular_networks ( b, neighs, conns, networks );
		if( b.tstep == 0 ){
			for( int i = 0; i < b.N; ++i ){

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
		}else{
			// Only two blocks, so.
			int target_mol = 341;
			const auto &mol_conns = conns[target_mol];
			std::cerr << target_mol << " connects to ";
			for( int mol_id : mol_conns ){
				std::cerr << " " << mol_id;
			}
			std::cerr << "\n";
			id_map im( id );
			std::vector<int> patches = { 52094, 52083, 52077,
			                             52117, 52111, 52100 };
			for( int i : patches ){
				int idx = im[ i ];
				const std::vector<int> &ni = neighs[idx];
				std::cerr << "Particle " << i << " (index = "
				          << idx << ") has neighs:";
				for( int j : ni ){
					std::cerr << " " << id[j];
				}
				std::cerr << "\n";
			}


			REQUIRE( util::contains( mol_conns, 15 ) );
			REQUIRE( util::contains( mol_conns, 260 ) );
			REQUIRE( util::contains( mol_conns, 325 ) );
		}
	}

}
