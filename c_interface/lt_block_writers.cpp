#include "lt_block_writers.h"

#include <iostream>
#include <string>
#include <cstring>


int lt_block_writers_lammps_data( const char *fname, const char *w_mode,
                                  lt_block_data_handle bdh )
{
	std::ofstream out;
	std::ios_base::openmode wmode = std::ios_base::out;
	if( c_like_open_out_file( fname, w_mode, out, wmode ) ){
		std::cerr << "An error occured opening file " << fname
		          << " with write mode " << w_mode << "! Aborting...\n";
		return -1;
	}

	// Writing data to binary is an error:
	if( wmode & std::ios_base::binary ){
		std::cerr << "Cannot write binary to LAMMPS data files!\n";
		return -1;
	}

	std::cerr << "Writing " << bdh.bd->N_types << " types to " << fname << "\n";
	lammps_tools::writers::block_to_lammps_data( out, *bdh.bd );
	return 0;
}

int lt_block_writers_lammps_dump( const char *fname, const char *w_mode,
                                  lt_block_data_handle bdh )
{
	std::ofstream out;
	std::ios_base::openmode wmode = std::ios_base::out;
	if( c_like_open_out_file( fname, w_mode, out, wmode ) ){
		std::cerr << "An error occured opening file " << fname
		          << " with write mode " << w_mode << "! Aborting...\n";
		return -1;
	}

	// Determine the output format:
	int fformat = lammps_tools::FILE_FORMAT_UNSET;
	if( wmode & std::ios_base::binary ){
		fformat = lammps_tools::FILE_FORMAT_BIN;
	}else if( wmode & std::ios_base::out ){
		fformat = lammps_tools::FILE_FORMAT_PLAIN;
	}else{
		std::cerr << "No clue which write mode is given in format "
		          << w_mode << "! Aborting!\n";
		return -1;
	}

	lammps_tools::writers::block_to_lammps_dump( out, *bdh.bd, fformat );
	return 0;
}

int c_like_open_out_file( const char *fname, const char *w_mode,
                          std::ofstream &out, std::ios_base::openmode &wmode )
{
	// Go over the write mode and set appropriate bits.
	// We already know you want out, not in:
	bool write_type_set = false;
	int len = std::strlen(w_mode);
	for( int i = 0; i < len; ++i ){
		char c = w_mode[i];
		switch( c ){
			case 'r':
				std::cerr << "Don't use these functions "
				          << "to read files!\n";
				return -2;
				break;
			case 'a':
				if( write_type_set ){
					std::cerr << "You cannot use multiple "
					          << "write modes!\n";
					return -3;
				}
				write_type_set = true;
				wmode |= std::ios_base::app;
				break;
			case 'w':
				if( write_type_set ){
					std::cerr << "You cannot use multiple "
					          << "write modes!\n";
					return -3;
				}
				write_type_set = true;
				wmode |= std::ios_base::out;
				break;
			case '+':
				std::cerr << "Updating writes not supported!\n";
				return -4;
				break;
			case 'b':
				wmode |= std::ios_base::binary;
				break;
			default:
				std::cerr << "Unknown write option " << c
				          << "!\n";
				return -1;
				break;
		}
	}

	out = std::ofstream( fname, wmode );
	return 0;
}
