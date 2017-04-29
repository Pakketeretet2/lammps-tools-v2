#include "dump_reader.hpp"

#include "dump_reader_lammps_bin.hpp"
#include "dump_reader_lammps_gzip.hpp"
#include "dump_reader_lammps_plain.hpp"

#include <memory> // Smart pointers.


using namespace dump_readers;


const char *dump_readers::fformat_to_str( int file_format )
{
	switch( file_format ){
		default:
			return "UNKNOWN!";
		case PLAIN:
			return "PLAIN TEXT";
		case GZIP:
			return "GZIPPED TEXT";
		case BIN:
			return "BINARY";
	}
}


const char *dump_readers::dformat_to_str( int dformat )
{
	switch( dformat ){
		default:
			return "UNKOWN!";
		case LAMMPS:
			return "LAMMPS";
		case HOOMD:
			return "HOOMD";
		case NAMD:
			return "NAMD";
	}
}


dump_reader *dump_readers::make_dump_reader( const std::string &fname,
                                             int dformat, int fformat )
{
	dump_reader *reader = nullptr;

	if( dformat == LAMMPS ){
		reader = make_dump_reader_lammps( fname, fformat );
	}else if( dformat == HOOMD ){
		if( fformat == BIN ){
			// reader = dump_reader_hoomd_gsd( fname );			
		}
	}else if( dformat == NAMD ){
		if( fformat == BIN ){
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

dump_reader *dump_readers::make_dump_reader( std::istream &input,
                                             int dformat, int fformat )
{
	dump_reader *reader = nullptr;
	
	if( dformat == LAMMPS ){
		if( fformat == PLAIN ){
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




dump_reader *dump_readers::make_dump_reader_lammps( const std::string &fname,
                                                    int fformat,
                                                    std::vector<std::string> headers )
{
	dump_reader_lammps *reader = nullptr;
	if( fformat == PLAIN ){
		reader = new dump_reader_lammps_plain( fname );
	}else if( fformat == BIN ){
		reader = new dump_reader_lammps_bin( fname );
		if( headers.empty() ){
			my_logic_error( __FILE__, __LINE__, "Column headers are "
			                "required for binary LAMMPS dump files!" );
		}
	}else if( fformat == GZIP ){
		reader = new dump_reader_lammps_gzip( fname );
	}
	if( reader ) reader->set_column_headers( headers );
	return reader;
}

dump_reader *dump_readers::make_dump_reader_lammps( std::istream &input,
                                                    std::vector<std::string> headers )
{
	dump_reader_lammps *reader = nullptr;
	reader = new dump_reader_lammps_plain( input );
	if( reader ) reader->set_column_headers( headers );
	return reader;
	
}

dump_reader *dump_readers::make_dump_reader_lammps( std::istream &input )
{
	std::vector<std::string> empty;
	return make_dump_reader_lammps( input, empty );
}

dump_reader *dump_readers::make_dump_reader_lammps( const std::string &fname,
                                                    int fformat )
{
	std::vector<std::string> empty;
	return make_dump_reader_lammps( fname, fformat, empty );	
}




std::size_t dump_readers::dump_reader::block_count()
{
	std::size_t c = 0;
	int status = 0;
	while( !status ){
		block_data b;
		status = next_block(b);
		++c;
	}
	return c;
}
