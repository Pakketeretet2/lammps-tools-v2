#include "dump_reader.hpp"

#include "dump_reader_lammps.hpp"
#include "dump_reader_lammps_bin.hpp"
#include "dump_reader_lammps_gzip.hpp"
#include "dump_reader_lammps_plain.hpp"
#include "dump_reader_xyz.hpp"

#include <memory> // Smart pointers.

namespace lammps_tools {

namespace readers {


int dump_reader::next_block( block_data &block, bool warn_if_no_special )
{
	int status = get_next_block( block );
	if( status ) return status;

	if( warn_if_no_special && (block.n_special_fields() == 0) ){
		dump_reader_lammps *drl;
		drl = dynamic_cast<dump_reader_lammps*>( this );
		if( !drl || (drl && drl->dump_style != drl->LOCAL) ){
			std::cerr << "Dump style is " << drl->dump_style << "\n";
			my_warning( __FILE__, __LINE__,
			            "Block has no special fields!" );
		}

	}


	// You might have read in atom types. In this case, make sure
	// the block data knows that.
	if( auto df = block.get_special_field( block_data::TYPE ) ){
		const std::vector<int> &types = data_as<int>( df );
		int ntypes = *std::max_element( types.begin(), types.end() );
		block.set_ntypes( ntypes );
	}
	return status;
}


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
	}else if( dformat == DUMP_FORMAT_LAMMPS_LOCAL ){
		reader = make_dump_reader_lammps( fname, fformat, dump_reader_lammps::LOCAL );
	}else if( dformat == DUMP_FORMAT_HOOMD ){
		if( fformat == FILE_FORMAT_BIN ){
			// reader = dump_reader_hoomd_gsd( fname );
		}
	}else if( dformat == DUMP_FORMAT_NAMD ){
		if( fformat == FILE_FORMAT_BIN ){
			// reader = dump_reader_namd_dcd( fname );
		}
	}else if( dformat == DUMP_FORMAT_XYZ ){
		if( fformat == FILE_FORMAT_PLAIN ){
			reader = new dump_reader_xyz( fname );
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
	}else if( dformat == DUMP_FORMAT_LAMMPS_LOCAL ){
		if( fformat == FILE_FORMAT_PLAIN ){
			reader = make_dump_reader_lammps( input, dump_reader_lammps::LOCAL );
		}
	}else if( dformat == DUMP_FORMAT_XYZ ){
		if( fformat == FILE_FORMAT_PLAIN ){
			reader = new dump_reader_xyz( input );
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
