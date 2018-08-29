#include "dump_reader.hpp"

#include "dump_reader_hoomd_gsd.hpp"
#include "dump_reader_lammps.hpp"
#include "dump_reader_lammps_bin.hpp"
#include "dump_reader_lammps_gzip.hpp"
#include "dump_reader_lammps_plain.hpp"
#include "dump_reader_xyz.hpp"

#ifdef THREADED_READ_BLOCKS
#include "../dependencies/readerwriterqueue/readerwriterqueue.h"
constexpr const bool threaded_read_blocks = true;
#else
constexpr const bool threaded_read_blocks = false;

// Dummy class:
namespace moodycamel {

template <typename T, size_t MAX_BLOCK_SIZE>
class ReaderWriterQueue
{
public:
	ReaderWriterQueue(){};

	T t;
};

}// namespace moodycamel
#endif // THREADED_READ_BLOCKS

#include <memory> // Smart pointers.

namespace lammps_tools {

namespace readers {

dump_reader::dump_reader()
	: quiet(true), read_blocks(nullptr), read_started(false)
{
	if( threaded_read_blocks ){
		read_blocks = new moodycamel::ReaderWriterQueue<block_data>;
	}
}

dump_reader::~dump_reader()
{
	if( read_blocks ){
		delete read_blocks;
	}
}


/**
   \brief adds type names that reflect the integer types of the particles.

   \param block The block to modify.
*/
void set_type_names_to_type( block_data &block )
{
	block.ati.type_names[0] = "__UNUSED__";
	for( std::size_t i = 1; i < block.ati.type_names.size(); ++i ){
		block.ati.type_names[i] = std::to_string( i );
	}
}



int dump_reader::next_block( block_data &block, bool warn_if_no_special )
{
	if( threaded_read_blocks ){
		return next_block_thr_impl( block, warn_if_no_special );
	}else{
		return next_block_impl( block, warn_if_no_special );
	}
}


int dump_reader::next_block_thr_impl( block_data &block,
				      bool warn_if_no_special )
{
	return next_block_impl( block, warn_if_no_special );
}


int dump_reader::next_block_impl( block_data &block, bool warn_if_no_special )
{
	int status = get_next_block( block );

	if( status ){
		if( check_eof() ){
			std::cerr << "Reached EOF.\n";
		}else{
			std::cerr << "Encountered status = " << status
			          << " in next_block!\n";
			my_warning( __FILE__, __LINE__,
			            "Status != 0 while getting block" );
		}
		return status;
	}

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
		bool fix_names = false;
		if( block.ati.type_names.empty() ) fix_names = true;

		block.set_ntypes( ntypes );
		if( fix_names ) set_type_names_to_type( block );
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
			reader = new dump_reader_hoomd_gsd( fname );
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
