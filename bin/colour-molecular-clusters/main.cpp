#include <iomanip>
#include <iostream>
#include <list>
#include <memory>
#include <string>

#include "cluster_finder.hpp"
#include "dump_reader_lammps.hpp"
#include "my_timer.hpp"
#include "neighborize.hpp"
#include "readers.hpp"
#include "util.hpp"
#include "writers.hpp"

using namespace lammps_tools;
using neighborize::make_list_dist_indexed;

std::string out_file;
std::list<std::string> dumps;
int n_patches;

int parse_args( int argc, char **argv )
{
	int i = 1;
	while( i < argc ){
		std::string arg = argv[i];
		if( arg == "-o" || arg == "--output" ){
			out_file = argv[i+1];
			i += 2;
		}else if( arg == "-np" || arg == "--n-patches" ){
			n_patches  = std::atoi( argv[i+1] );
			i += 2;
		}else{
			dumps.push_back( arg );
			++i;
		}
	}
	return 0;
}


void colour_blocks( readers::dump_reader_lammps *d,
                    std::ostream &hist, std::ostream &out,
                    std::ostream &cluster_info,
                    std::ostream &subunit_info,
                    int fformat, int max_cluster_size )
{
	block_data b;
	int block_count = 0;
	my_timer t(std::cerr);
	t.tic();

	while( d->next_block(b) == 0 ){

		double rc = 2.0;

		const std::vector<int> &id = data_as<int>(
			b.get_special_field( block_data::ID ) );
		const std::vector<int> &type = data_as<int>(
			b.get_special_field( block_data::TYPE ) );
		const std::vector<int> &mol = data_as<int>(
			b.get_special_field( block_data::MOL ) );

		neighborize::neigh_list neighs;
		neighborize::neigh_list neighs_patches_only;

		int method = neighborize::DIST_BIN;
		int dims = 3;

		int policy_include = neighborize::neighborizer::INCLUDE;
		int policy_ignore  = neighborize::neighborizer::IGNORE;

		std::vector<int> first, second;
		for( int i = 0; i < b.N; ++i ){
			if( type[i] > 1 ){
				first.push_back(i);
				second.push_back(i);
			}
		}

		double avg_neighs = make_list_dist_indexed( neighs, b,
		                                            first,
		                                            second,
		                                            method,
		                                            dims, rc,
		                                            policy_include,
		                                            policy_ignore );
		make_list_dist_indexed( neighs_patches_only, b, first, second,
		                        method, dims, rc,
		                        policy_ignore, policy_ignore );

		neighborize::neigh_list conns;
		std::list<std::list<int> > networks;

		neighborize::find_molecular_networks( b, neighs, conns,
		                                      networks );
		int max_mol = *std::max_element( mol.begin(),
		                                 mol.end() );

		// Store mol_id --> size_of_cluster:
		std::vector<int> mol_to_cluster_size( max_mol + 1 );
		std::vector<int> size_distro( max_cluster_size + 1 );
		for( const std::list<int> &li : networks ){
			for( int mol_id : li ){
				mol_to_cluster_size[mol_id] = li.size();
			}
			int ss = li.size();
			if( ss < max_cluster_size ) size_distro[ss]++;
		}

		// Output for the histogram:
		for( int cs = 1; cs < max_cluster_size; ++cs ){
			hist << std::setw(8) << b.tstep << "\t"
			     << std::setw(8) << cs << "\t"
			     << std::setw(8) << size_distro[cs] << "\n";
		}
		hist << "\n";


		// Colour each molecule according to the size of
		// the cluster it's in.
		data_field_int cluster_size( "cluster_size", b.N );
		for( int i = 0; i < b.N; ++i ){
			int mol_id = mol[i];
			cluster_size[i] = mol_to_cluster_size[mol_id];
		}
		b.add_field( cluster_size );


		// Colour each molecule according to the number of
		// connections it has to other molecules.
		std::vector<std::vector<int> > mol_to_patches( max_mol + 1 );
		for( int i = 1; i <= max_mol; ++i ){
			mol_to_patches[i].reserve( 6 );
		}

		for( int i = 0; i < b.N; ++i ){
			if( type[i] > 1 ){
				int mol_id = mol[i];
				mol_to_patches[mol_id].push_back(i);
			}
		}

		std::vector<int> mol_to_connection_count( max_mol + 1 );
		for( int mol_id = 1; mol_id <= max_mol; ++mol_id ){
			int connections = 0;
			my_assert( __FILE__, __LINE__,
			           mol_to_patches[mol_id].size() == n_patches,
			           "Failed to find all patches in mol!" );

			for( int patch_idx : mol_to_patches[mol_id] ){
				const std::vector<int> &neighs_i
					= neighs_patches_only[patch_idx];
				connections += neighs_i.size();
			}
			mol_to_connection_count[mol_id] = connections;
		}
		// Colour this entire molecule this colour.
		data_field_int connections( "connections", b.N );
		for( int i = 0; i < b.N; ++i ){
			int mol_id = mol[i];
			connections[i] = mol_to_connection_count[mol_id];
		}
		b.add_field( connections );


		// Output the per-assembly info:
		int cluster_idx = 0;
		for( std::vector<int> &cluster : conns ){
			if( cluster.empty() ) continue;
			int mol_id = cluster.front();
			int cs = mol_to_cluster_size[mol_id];
			cluster_info << std::setw(8) << b.tstep << "\t"
			             << std::setw(8) << cluster_idx << "\t"
			             << std::setw(8) << cs << "\n";
			++cluster_idx;
		}
		cluster_info << "\n";

		// Output the per-subunit info:
		for( int mol_id = 1; mol_id <= max_mol; ++mol_id ){
			int cs = mol_to_cluster_size[mol_id];
			int cc = mol_to_connection_count[mol_id];
			subunit_info << std::setw(8) << b.tstep << "\t"
			             << std::setw(8) << mol_id << "\t"
			             << std::setw(8) << cs << "\t"
			             << std::setw(8) << cc << "\n";
		}
		subunit_info << "\n";
		// Write output
		writers::block_to_lammps_dump( out, b, fformat );
		if( (block_count > 0) && (block_count % 25 == 0) ){
			double d_time = t.toc( "Colouring 25 blocks" );
			std::cerr << "At t = " << b.tstep << "...\n";
			std::cerr << "Average neighs = " << avg_neighs
			          << "\n";
			double rate = d_time / 25.0;
			std::cerr << "Rate is " << rate << " ms/block"
			          << " ( " << 1000.0 / rate
			          << " blocks / s ).\n";
			t.tic();

		}
		++block_count;
	}
}



int main( int argc, char **argv )
{
	std::ostream *out;
	std::unique_ptr<std::ofstream> of;
	n_patches = 6;

	out_file = "";
	if( parse_args( argc, argv ) ){
		std::cerr << "Error parsing args!\n";
		return -1;
	}
	if( out_file.empty() ){
		out = &std::cout;
	}else{
		if( util::ends_with( out_file, ".bin" ) ){
			auto bin = std::ios::binary;
			std::ofstream *tmp = new std::ofstream( out_file, bin );
			of = std::unique_ptr<std::ofstream>( tmp );
		}else{
			std::ofstream *tmp = new std::ofstream( out_file );
			of = std::unique_ptr<std::ofstream>( tmp );
		}
		if( !of ){
			std::cerr << "Failed to open file!\n";
			return -2;
		}else{
			out = of.get();
		}
	}

	std::ofstream hist( "histogram.dat" );
	int max_cluster_size = 30;

	std::ofstream cluster_info( "clusters.dat" );
	std::ofstream subunit_info( "subunits.dat" );

	std::cerr << "Colouring " << dumps.size() << " dump file(s)...\n";

	std::vector<std::string> common_headers = { "id", "type",
	                                            "x", "y", "z" };

	std::vector<std::string> all_headers = { "id", "mol", "type",
	                                         "x", "y", "z" };

	for( const std::string &dname : dumps ){
		int fformat = FILE_FORMAT_PLAIN;

		if( util::ends_with( dname, ".bin" ) ){
			fformat = FILE_FORMAT_BIN;
		}else if( util::ends_with( dname, ".gz" ) ){
			fformat = FILE_FORMAT_GZIP;
		}

		std::unique_ptr<readers::dump_reader_lammps> d(
			readers::make_dump_reader_lammps( dname, fformat ) );

		d->set_column_headers( all_headers );
		d->set_column_header_as_special(   "id", block_data::ID );
		d->set_column_header_as_special(  "mol", block_data::MOL );
		d->set_column_header_as_special( "type", block_data::TYPE );
		d->set_column_header_as_special(    "x", block_data::X );
		d->set_column_header_as_special(    "y", block_data::Y );
		d->set_column_header_as_special(    "z", block_data::Z );

		colour_blocks( d.get(), hist, *out, cluster_info,
		               subunit_info, fformat, max_cluster_size );
	}
}
