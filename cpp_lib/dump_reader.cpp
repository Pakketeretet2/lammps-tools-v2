#include "dump_reader.hpp"

#include "dump_reader_lammps_bin.hpp"
#include "dump_reader_lammps_gzip.hpp"
#include "dump_reader_lammps_plain.hpp"

#include <memory> // Smart pointers.

namespace lammps_tools {

namespace readers {

const char *fformat_to_str( int file_format )
{
	switch( file_format ){
		default:
			return "UNKNOWN!";
		case FILE_FORMAT_PLAIN:
			return "PLAIN TEXT";
		case FILE_FORMAT_GZIP:
			return "GZIPPED TEXT";
		case FILE_FORMAT_BIN:
			return "BINARY";
	}
}


const char *dformat_to_str( int dformat )
{
	switch( dformat ){
		default:
			return "UNKOWN!";
		case DUMP_FORMAT_LAMMPS:
			return "LAMMPS";
		case DUMP_FORMAT_HOOMD:
			return "HOOMD";
		case DUMP_FORMAT_NAMD:
			return "NAMD";
	}
}


dump_reader *make_dump_reader( const std::string &fname,
                               int fformat, int dformat )
{
	dump_reader *reader = nullptr;

	if( dformat == DUMP_FORMAT_LAMMPS ){
		reader = make_dump_reader_lammps( fname, fformat );
	}else if( dformat == DUMP_FORMAT_HOOMD ){
		if( fformat == FILE_FORMAT_BIN ){
			// reader = dump_reader_hoomd_gsd( fname );
		}
	}else if( dformat == DUMP_FORMAT_NAMD ){
		if( fformat == FILE_FORMAT_BIN ){
			// reader = dump_reader_namd_dcd( fname );
		}
	}
	if( !reader ){
		std::string msg = "Failed to construct reader; file format = ";
		msg += fformat_to_str(fformat);
		msg += ", dump format = ";
		msg += dformat_to_str(dformat);
		my_runtime_error( __FILE__, __LINE__, msg );
	}

	return reader;
}

dump_reader *make_dump_reader( std::istream &input,
                               int fformat, int dformat )
{
	dump_reader *reader = nullptr;

	if( dformat == DUMP_FORMAT_LAMMPS ){
		if( fformat == FILE_FORMAT_PLAIN ){
			reader = make_dump_reader_lammps( input );
		}
	}

	if( !reader ){
		std::string msg = "Failed to construct reader; file format = ";
		msg += fformat_to_str(fformat);
		msg += ", dump format = ";
		msg += dformat_to_str(dformat);
		my_runtime_error( __FILE__, __LINE__, msg );
	}

	return reader;
}





std::size_t number_of_blocks( dump_reader &dr )
{
	std::size_t c = 0;
	block_data b;
	while( dr.next_block( b ) == 0 ) ++c;
	return c;
}

} // namespace readers

} // namespace lammps_tools
