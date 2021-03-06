#include <iostream>
#include <memory>
#include <vector>

#include "writers_lammps.hpp"
#include "dump_reader_lammps.hpp"
#include "util.hpp"

void print_usage()
{
	std::cerr << "Use this like dump-cat -c <column headers> <dump files> "
	          << " -o <output>,\n"
	          << "  where <headers> is a quoted string containing the\n"
	          << "  names of the headers in the dump file and <output> is\n"
	          << "  an optional output dump file name. If no dump files\n"
	          << "  are given, stdcin is read.\n";
}


int main( int argc, char **argv )
{
	using namespace lammps_tools;

	std::vector<std::string> dumps;
	std::string headers;
	std::string out_file = "-";
	bool to_local = false;
	bool silent = false;
	bool ignore_first = false;

	int i = 1;
	while( i < argc ){
		std::string a( argv[i] );
		if( a[0] == '-' ){
			if( a == "-c" || a == "--column-headers" ){
				headers = argv[i+1];
				i += 2;
			}else if( a == "-o" || a == "--output" ){
				out_file = argv[i+1];
				i += 2;
			}else if( a == "-l" || a == "--local" ){
				to_local = true;
				i += 1;
			}else if( a == "-h" || a == "--help" ){
				print_usage();
				return 1;
			}else if( a == "-i" || a == "--ignore-first" ){
				ignore_first = true;
				i += 1;
			}else if( a == "-s" || a == "--silent" ){
				silent = true;
				i += 1;
			}
		}else{
			dumps.push_back( a );
			i += 1;
		}
	}

	bool read_stdin = dumps.empty();
	if( read_stdin ){
		std::cerr << "Reading from stdin...\n";
	}else{
		std::cerr << "Catting ";
		for( const std::string &s : dumps ){
			std::cerr << " " << s;
		}
		std::cerr << "\n";
	}
	std::cerr << "Headers are " << headers << ".\n";
	if( ignore_first ){
		std::cerr << "Dropping first frame of subsequent dump file!\n";
	}

	std::ofstream *out_fstream = nullptr;
	std::ostream *out = nullptr;
	int out_fformat = FILE_FORMAT_PLAIN;

	if( out_file == "-" ){
		out = &std::cout;
	}else{
		if( util::ends_with( out_file, ".bin") ){
			out_fformat = FILE_FORMAT_BIN;
			out_fstream = new std::ofstream( out_file, std::ios::binary );
			out = out_fstream;
		}else if( util::ends_with( out_file, ".gz" ) ){
			std::cerr << "Writing to gzip is not supported by "
			          << "dump-cat. Instead, stream text output "
			          << "through gzip!\n";
			return -2;
		}else{
			std::cerr << "Assuming file " << out_file
			          << " is plain text dump file.\n";
			out_fformat = FILE_FORMAT_PLAIN;
			out_fstream = new std::ofstream( out_file );
			out = out_fstream;
		}
	}

	int n_writes = 0;
	auto status_print = [&n_writes, silent]	{
		if( !silent && n_writes > 0  && n_writes % 50 == 0 ){
			std::cerr << "At block " << n_writes << "...\n";
		}
		n_writes++; };

	if( read_stdin ){
		std::istream &in( std::cin );
		int fformat = FILE_FORMAT_PLAIN;
		block_data b;
		if( to_local ){
			std::unique_ptr<readers::dump_reader_lammps> d(
				readers::make_dump_reader_lammps( in, readers::dump_reader_lammps::LOCAL ) );
			d->set_column_headers( util::split( headers ) );
			while( d->next_block(b) == 0 ){
				writers::block_to_lammps_dump( *out, b, out_fformat, to_local );
				status_print();
			}
		}else{
			std::unique_ptr<readers::dump_reader_lammps> d(
				readers::make_dump_reader_lammps( in ) );
			d->set_column_headers( util::split( headers ) );
			while( d->next_block(b) == 0 ){
				writers::block_to_lammps_dump( *out, b, out_fformat, to_local );
				status_print();
			}
		}

	}else{
		int dump_idx = 0;
		for( const std::string &dname : dumps ){
			int fformat = FILE_FORMAT_PLAIN;
			if( util::ends_with( dname, ".bin" ) ){
				fformat = FILE_FORMAT_BIN;
			}else if( util::ends_with( dname, ".gz" ) ){
				fformat = FILE_FORMAT_GZIP;
			}
			std::cerr << "Catting " << dname << ".\n";
			block_data b;
			if( to_local ){
				std::unique_ptr<readers::dump_reader_lammps> d(
					readers::make_dump_reader_lammps( dname, fformat, readers::dump_reader_lammps::LOCAL ) );
				d->set_column_headers( util::split( headers ) );
				if( ignore_first && dump_idx > 0 ){
					d->next_block(b);
				}

				while( d->next_block(b) == 0 ){
					writers::block_to_lammps_dump( *out, b,
					                               out_fformat,
					                               to_local );
					status_print();
				}
			}else{
				std::unique_ptr<readers::dump_reader_lammps> d(
					readers::make_dump_reader_lammps( dname, fformat ) );
				d->set_column_headers( util::split( headers ) );

				if( ignore_first && dump_idx > 0 ){
					d->next_block(b);
				}

				while( d->next_block(b) == 0 ){
					writers::block_to_lammps_dump( *out, b,
					                               out_fformat,
					                               to_local );
					status_print();
				}
			}
			++dump_idx;
		}
	}
	if( out_fstream ) delete out_fstream;
}
