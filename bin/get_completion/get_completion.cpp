/**
   \file get_completion.cpp

   This program extracts the concentrations of different sized aggregates.

   Link against lammpstools.
   Compile with -llammpstools -I(Clara_dir)
*/

#include "clara.hpp"

#include "block_data.hpp"
#include "block_data_access.hpp"
#include "cluster_finder.hpp"
#include "dump_reader.hpp"
#include "neighborize.hpp"

#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <vector>

#ifdef DO_FILTER_BLOCK
constexpr const bool do_filter = true;
#else
constexpr const bool do_filter = false;
#endif // DO_FILTER_BLOCK

using namespace lammps_tools;
using namespace lammps_tools::neighborize;

enum dump_formats {
	GNUPLOT,
	OLD
};

void dump_header_gnuplot( std::ostream &out )
{
	out << "# time  size  count\n";
}


void dump_population_gnuplot( std::ostream &out, double time,
                              const std::vector<int> &sizes,
                              const std::vector<int> &counts )
{
	out << std::fixed << std::setprecision(0);
	for( std::size_t i = 0; i < counts.size(); ++i ){
		out << std::setw(10) << time << " "
		    << std::setw(4) << sizes[i] << " "
		    << std::setw(6) << counts[i] << "\n";
	}
	out << "\n";
}


void dump_header( std::ostream &out, const std::vector<int> &sizes )
{
	out << "# time ";
	for( std::size_t i = 0; i < sizes.size(); ++i ){
		out << " n" << sizes[i];
	}
	out << " max_size\n";
}


void dump_population( std::ostream &out, double time,
                      const std::vector<int> &counts,
                      std::size_t max_size )
{
	out << std::fixed << std::setprecision(0);
	out << time;
	for( std::size_t i = 0; i < counts.size(); ++i ){
		out << " " << std::setw(8) << counts[i];
	}
	out << " " << std::setw(8) << max_size << "\n";

}



block_data filter_block( const block_data &b )
{
	const std::vector<int> &types = get_type( b );
	const std::vector<int> &ids = get_id(b);

	std::vector<int> id;
	id.reserve(b.N/10);

	const std::vector<int> type_i = { 3 , 5, 7 };
	const std::vector<int> type_j = { 4 , 6, 8 };


	for( std::size_t idx = 0; idx < b.N; ++idx ){
		int ttype = types[idx];

		if( ttype == type_i[0] ||
		    ttype == type_i[1] ||
		    ttype == type_i[2] ||
		    ttype == type_j[0] ||
		    ttype == type_j[1] ||
		    ttype == type_j[2] ){
			id.push_back(ids[idx]);
		}
	}
	return filter_block(b,id);
}



void analyze_dump( readers::dump_reader *d, std::vector<int> sizes,
                   double dt, double rc, std::ostream &out,
                   uint64_t start_frame, int output_format )
{
	block_data b, b2;

	int bc = 0;

	std::vector<int> type_i = { 3 , 5, 7 };
	std::vector<int> type_j = { 4 , 6, 8 };

	std::map< int, int > size2count;
	int idx = 0;
	for( int s : sizes ){
		size2count[ s ] = idx;
		++idx;
	}

	double old_t = 0.0;

	switch( output_format ){
		case OLD:
		default:
			dump_header( out, sizes );
			break;
		case GNUPLOT:
			dump_header_gnuplot( out );
	}


	while( d->next_block( b2 ) == 0 ){
		if (do_filter) {
			b = filter_block(b2);
		} else {
			b = b2;
		}
		b.dom.periodic = 7;
		double time = b.tstep * dt;
		// std::cerr << "time = " << time << ", old_time = " << old_t << "\n";
		if( (bc != 0) && (time <= old_t) ){
			// This means there was two trajectories
			// in the dump file.
			std::cerr << "Got block with mismatched time! "
			          << "Ignoring frame " << bc << " (t = "
			          << b.tstep << ")!\n";
			++bc;
			continue;
		}

		const std::vector<int> &mols = get_mol( b );
		const std::vector<int> &types = get_type( b );

		int nmols = *std::max_element( mols.begin(), mols.end() );
		int atoms_per_mol = b.N / nmols;
		int target_atom = 25180;
		int target_mol = target_atom / atoms_per_mol;

		/*
		  Analysis scheme:

		  1. Construct molecular network.
		  2. Count clusters.
		  3. Make distribution of clusters sizes.
		*/
		std::vector<int> ilist, jlist;

		for( std::size_t idx = 0; idx < b.N; ++idx ){
			int ttype = types[idx];

			if( ttype == type_i[0] ||
			    ttype == type_i[1] ||
			    ttype == type_i[2] ){
				ilist.push_back(idx);
			}
			if( ttype == type_j[0] ||
			    ttype == type_j[1] ||
			    ttype == type_j[2] ){
				jlist.push_back(idx);
			}
		}

		double scale_factor = static_cast<double>(b.N) / nmols;
		neigh_list nl;
		double avg_neighs = make_list_dist_indexed( nl, b, ilist, jlist,
		                                            DIST_BIN, 3, rc );

		// Step 1:  Construct a list of molecular clusters
		auto mol_nlist = get_molecular_connections( b, nl, false );

		// Step 2:  Now we need to convert the molecular connections
		//          to actual clusters.
		auto mol_clusters = neigh_list_to_clusters( mol_nlist );

		// Step 3: Make distribution:
		std::vector<int> counts( sizes.size(), 0 );
		int max_size = 0;
		for( const auto &mneighs : mol_clusters ){
			int size = mneighs.size();
			if( size > max_size ) max_size = size;

			if( size2count.count( size ) ){
				++counts[ size2count[ size ] ];
			}
		}

		switch( output_format ){
			default:
			case GNUPLOT:
				dump_population_gnuplot( out, time,
				                         sizes, counts );
				break;
			case OLD:
				dump_population( out, time, counts, max_size );
				break;
		}


		++bc;
		old_t = time;
		if( bc % 50 == 0 ){
			std::cerr << "  At block " << bc << "...\n";
		}

	}
	std::cerr << "Done with analysis on " << bc << " blocks!\n";
}








int main( int argc, char **argv )
{
	std::string dump_name = "";
	std::string file_format = "bin";
	std::string dump_format = "hoomd";
	std::string out_fname = "-";
	std::string dump_format_str = "OLD";

	double rc = 1.25;
	bool print_help = false;
	uint64_t start_frame = 0;
	bool dump_all_sizes = false;
	int max_capsid_size = 80;
	int output_format = OLD;


	if( argc < 2 ){
		std::cerr << "Pass a dump file to process!\n";
		return -1;
	}


	auto cli = clara::Opt( dump_format, "dump-format" )
		["-d"]["--dump-format"]( "Dump format of dump file" )
		| clara::Opt( file_format, "file-format" )
		["-f"]["--file-format"]( "File format of dump file" )
		| clara::Opt( rc, "rc" )
		["-r"]["--cutoff"]( "Cut-off for the neighbor list construction" )
		| clara::Arg( dump_name, "dump-name" )
		| clara::Opt( dump_all_sizes, "dump_all_sizes" )
		["-a"]["--all-sizes"]( "Dump all sizes instead of selection?" )
		| clara::Help( print_help )
		( "Display usage information" )
		| clara::Opt( start_frame, "start-frame" )
		["-s"]["--start"]( "Frame to start analysis at." )
		| clara::Opt( dump_format_str, "dump-format" )
		["--output-format"]( "Output format." )
		| clara::Opt( max_capsid_size, "max_capsid_size" )
		["-m"]["--max-capsid-size"]( "Largest capsid size to dump" )
		| clara::Opt( out_fname, "output" )
		["-o"]["--output"]( "Output file name (\"-\" for stdout)" );

	auto result = cli.parse( clara::Args( argc, argv ) );
	if( !result ){
		std::cerr << "Error in command line args: "
		          << result.errorMessage() << "\n";
		return -1;
	}

	if( print_help ){
		std::cerr << cli << "\n";
		return 0;
	}

	if( dump_format_str == "OLD" ){
		output_format = OLD;
	}else if( dump_format_str == "GNUPLOT" ){
		output_format = GNUPLOT;
	}else{
		std::cerr << "Unknown output format " << dump_format_str
		          << "! Aborting!\n";
		return -1;
	}

	std::cerr << "**  Options:  **\n";
	std::cerr << "    Dump file:     " << dump_name << "\n"
	          << "    File format:   " << file_format << "\n"
		  << "    Dump format:   " << dump_format << "\n"
	          << "    Cut-off:       " << rc << "\n"
	          << "    Output file:   " << out_fname << "\n"
	          << "    Output format: " << dump_format_str << "\n"
	          << "    Start-frame:   " << start_frame << "\n\n";

	int fformat = FILE_FORMAT_PLAIN;
	int dformat = DUMP_FORMAT_LAMMPS;

	if( dump_name.empty() ){
		std::cerr << "Pass a dump file!\n";
		return -1;
	}

	if( file_format == "plain" ){
		fformat = FILE_FORMAT_PLAIN;
	}else if( file_format == "bin" ){
		fformat = FILE_FORMAT_BIN;
	}else{
		std::cerr << "Unrecognized file format \"" << file_format
		          << "\"!\n";
	}

	if( dump_format == "lammps" ){
		dformat = DUMP_FORMAT_LAMMPS;
	}else if( dump_format == "hoomd" ){
		dformat = DUMP_FORMAT_HOOMD;
	}else{
		std::cerr << "Unrecognized dump format \"" << dump_format
		          << "\"!\n";
	}

	std::ostream *out = &std::cout;
	std::ofstream out_file;
	if( out_fname != "-" && !out_fname.empty() ){
		out_file.open( out_fname );
		out = &out_file;
	}

	double dt = 0.0025;

	std::unique_ptr<readers::dump_reader> d(
		readers::make_dump_reader( dump_name, fformat, dformat ) );

	// Interesting sizes:
	std::vector<int> sizes = { 1, 2, 3, 5, 20, 60, 80, 8, 10, 12, 15, 18, 19 };
	if( dump_all_sizes ){
		sizes.resize(max_capsid_size);
		for( int i = 0; i < max_capsid_size; ++i ){
			sizes[i] = i+1;
		}
	}

	analyze_dump( d.get(), sizes, dt, rc, *out, start_frame, output_format );



	return 0;
}
