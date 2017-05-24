#include <cmath>
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
#include "skeletonize.hpp"
#include "util.hpp"
#include "writers.hpp"

using namespace lammps_tools;
using neighborize::make_list_dist;

std::string out_file;
std::list<std::string> dumps;
double R;

int parse_args( int argc, char **argv )
{
	int i = 1;
	while( i < argc ){
		std::string arg = argv[i];
		if( arg == "-o" || arg == "--output" ){
			out_file = argv[i+1];
			i += 2;
		}else if( arg == "-R" || arg == "--radius" ){
			R = std::stof( argv[i+1] );
			std::cerr << "R is now " << R << "\n";
			i += 2;
		}else{
			dumps.push_back( arg );
			++i;
		}
	}
	return 0;
}



void skeletonize_blocks( readers::dump_reader_lammps *d,
                         std::ostream &out, int fformat,
                         std::ostream &widths_out,
                         double width_thresh = 1.0 )
{
	block_data b;
	int block_count = 0;
	my_timer t(std::cerr);

	std::vector<std::string> remove_these = { "c_stress[1]", "c_stress[2]",
	                                          "c_stress[3]", "c_stress[4]",
	                                          "c_stress[5]", "c_stress[6]" };

	t.tic();
	std::list<double> widths_outer;
	std::list<double> diams;

	while( d->next_block(b) == 0 ){
		double rc = 1.3;
		int dims = 3;

		// Remove some fields:
		data_field *df;
		for( const std::string &n : remove_these ){
			int special_field = block_data::UNKNOWN;
			df = b.remove_field( n, special_field );
			delete df;
		}

		// Construct the Euclidian distance transform:
		neighborize::neigh_list neighs;
		neighborize::make_list_dist( neighs, b, 0, 0,
		                             neighborize::DIST_BIN, dims, rc );
		data_field_double edt( "edt" );
		std::vector<double> inside;
		inside = skeletonize::get_insideness( b, neighs );
		edt.get_data_rw() =
			skeletonize::euclidian_distance_transform( b, inside, R );

		// Reduce the EDT to its local maxima and find
		// the distances to the edge.
		std::vector<double> widths =
			skeletonize::get_ribbon_widths( b, neighs,
			                                edt.get_data(),
			                                inside );
		auto filter_width = [width_thresh]( double w )
			{ return false; };


		widths.erase( std::remove_if( widths.begin(), widths.end(),
		                              filter_width ), widths.end() );

		double width_mean = util::mean( widths );
		double width_var  = util::var( widths );
		double diameter = *std::max_element( widths.begin(),
		                                     widths.end() );
		widths_out << std::setw(8) << b.tstep << "\t"
		           << std::setw(8) << diameter << "\t"
		           << std::setw(8) << width_mean << "\t"
		           << std::setw(8) << std::sqrt(width_var) << "\n";

		widths_outer.push_back( width_mean );
		diams.push_back( diameter );

		// Write to dump file:
		b.add_field( edt );
		writers::block_to_lammps_dump( out, b, fformat );

		if( block_count > 0 && (block_count % 25 == 0) ){
			double t_toc = t.toc( "Getting EDT of 25 blocks" );
			double dt = t_toc;
			double rate = 25.0 * 1000.0 / dt;
			std::cerr << "Rate: " << rate << " blocks/s.\n";
			std::cerr << "At block " << block_count << ".\n\n";

			t.tic();
		}
		++block_count;
	}

	double avg_width = util::mean( widths_outer );
	double var_width = util::var( widths_outer );
	std::cerr << "\n<w> = " << avg_width << "\n";
	double avg_diam  = util::mean( diams );
	double var_diam = util::var( diams );

	std::ofstream avg_width_out( "WIDTH_AVG" );
	std::ofstream avg_diam_out( "DIAM_AVG" );
	avg_width_out << avg_width << " " << std::sqrt(var_width) << "\n";
	avg_diam_out << avg_diam << " " << std::sqrt(var_diam) << "\n";
}



int main( int argc, char **argv )
{
	std::ostream *out;
	std::unique_ptr<std::ofstream> of;

	R = -1;
	out_file = "";
	if( parse_args( argc, argv ) ){
		std::cerr << "Error parsing args!\n";
		return -1;
	}
	int out_fformat = FILE_FORMAT_PLAIN;
	if( out_file.empty() ){
		out = &std::cout;
	}else{
		if( util::ends_with( out_file, ".bin" ) ){
			auto bin = std::ios::binary;
			std::ofstream *tmp = new std::ofstream( out_file, bin );
			of = std::unique_ptr<std::ofstream>( tmp );
			out_fformat = FILE_FORMAT_BIN;
			std::cerr << "Writing to binary file "
			          << out_file << "\n";
		}else{
			std::ofstream *tmp = new std::ofstream( out_file );
			of = std::unique_ptr<std::ofstream>( tmp );
			std::cerr << "Writing to text file "
			          << out_file << "\n";
		}
		if( !of ){
			std::cerr << "Failed to open file!\n";
			return -2;
		}else{
			out = of.get();
		}
	}

	if( R <= 0 ){
		std::cerr << "Set radius of template from command line!\n";
		return -2;
	}

	std::vector<std::string> common_headers = { "id", "type",
	                                            "x", "y", "z" };

	std::vector<std::string> all_headers = { "id", "type",
	                                         "x", "y", "z", "c_stress[1]",
	                                         "c_stress[2]", "c_stress[3]",
	                                         "c_stress[4]", "c_stress[5]",
	                                         "c_stress[6]" };
	std::ofstream widths_out( "widths.dat" );

	for( const std::string &dname : dumps ){
		int in_fformat = FILE_FORMAT_PLAIN;


		if( util::ends_with( dname, ".bin" ) ){
			in_fformat = FILE_FORMAT_BIN;
		}else if( util::ends_with( dname, ".gz" ) ){
			in_fformat = FILE_FORMAT_GZIP;
		}

		std::unique_ptr<readers::dump_reader_lammps> d(
			readers::make_dump_reader_lammps( dname, in_fformat ) );

		d->set_column_headers( all_headers );
		d->set_column_header_as_special(   "id", block_data::ID );
		// d->set_column_header_as_special(  "mol", block_data::MOL );
		d->set_column_header_as_special( "type", block_data::TYPE );
		d->set_column_header_as_special(    "x", block_data::X );
		d->set_column_header_as_special(    "y", block_data::Y );
		d->set_column_header_as_special(    "z", block_data::Z );

		skeletonize_blocks( d.get(), *out,  out_fformat,
		                    widths_out, 0.2 );
	}
}
