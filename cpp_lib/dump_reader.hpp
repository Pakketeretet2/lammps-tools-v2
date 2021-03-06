#ifndef DUMP_READER_HPP
#define DUMP_READER_HPP

/**
   \file dump_reader.hpp

   Contains declaration of general dump_reader interface and helper functions.
*/

#include "block_data.hpp"
#include "enums.hpp"

#include <iosfwd>
#include <memory>

#ifdef THREADED_READ_BLOCKS
namespace moodycamel {
template <typename T, size_t MAX_BLOCK_SIZE> class ReaderWriterQueue;
}
#else
namespace moodycamel {
template <typename T, size_t MAX_BLOCK_SIZE = 512> class ReaderWriterQueue;
}
#endif


namespace lammps_tools {

/// Contains functions and classes that are related to reading dump files.
namespace readers {


/// A generic class for reading in dump files.
class dump_reader
{
public:
	/// Allocates the queue.
	dump_reader();

	/// Empty destructor:
	virtual ~dump_reader();

	/**
	   Attempts to read in the next block from file.

	   \param block   Will contain the new block_data on success, and shall
	                  be  unchanged upon failure.

           \returns 0 on success, positive if EOF reached, negative on failure.
	*/
	int next_block( block_data &block, bool warn_if_no_special = false );

	/// Checks if the internal file is at eof:
	bool eof()  const { return check_eof(); }

	/// Checks if the internal file is good:
	bool good() const { return check_good(); }

	/// If true, do not print output to stderr
	bool quiet;

	/// Some dump readers need to specify headers for the columns.
	virtual void set_column_headers(const std::vector<std::string> &headers)
	{}

	/// Skips n blocks, might not actually be fast.
	virtual int skip_n_blocks( uint n );


	/// Skips to specific block from current block.
	virtual int skip_to_block( uint n, uint current );
	

private:
	virtual int  get_next_block( block_data &block ) = 0;
	virtual bool check_eof()  const = 0;
	virtual bool check_good() const = 0;

	/// Implementation of serial next_block
	int next_block_impl( block_data &block, bool warn_if_no_special );

	/// Implementation of threaded next_block
	int next_block_thr_impl( block_data &block, bool warn_if_no_special );


	/// Contains a buffer to store the blocks in.
	moodycamel::ReaderWriterQueue<block_data,512> *read_blocks;
	/// If true, the first read was called and a thread will start
	/// filling read_blocks.
	bool read_started;
};




/**
   \brief Pretty-prints given file format.

   \param   file_format The file format to pretty-print.

   \returns A pretty-printed string of the file format.
*/
const char *fformat_to_str( int file_format );

/**
   \brief Pretty-prints given dump format.

   \param   dump_format The file format to pretty-print.

   \returns A pretty-printed string of the dump format.
*/
const char *dformat_to_str( int dformat );


/**
   \brief Constructs correct dump_reader for file of given file and dump format.

   \param fname    Name of the dump file.
   \param fformat  File format (see \p FILE_FORMATS).
   \param dformat  Dump format (see \p DUMP_FORMATS).

   \returns a pointer to the open dump_reader object, or nullptr on failure.

   \warning The dump_reader object is heap-allocated. Be sure to delete it
            or wrap it in a smart pointer.

*/
dump_reader* make_dump_reader( const std::string &fname,
                               int fformat, int dformat );


/**
   \brief Constructs correct dump_reader for input stream of given dump format.

   \param fname    Name of the dump file.
   \param fformat  File format (see \p FILE_FORMATS).
   \param dformat  Dump format (see \p DUMP_FORMATS).

   \returns a pointer to the open dump_reader object, or nullptr on failure.

   \warning The dump_reader object is heap-allocated. Be sure to delete it
            or wrap it in a smart pointer.
*/
dump_reader* make_dump_reader( std::istream &input,
                               int fformat, int dformat );



/**
   \brief Counts the number of blocks left in given dump reader.

   \warning This function is _destructive_ in the sense that the
            state of dump_reader is changed. If you want to know
            the number of blocks of a dump file you are going to
            read, you should make a dump_reader, call this function,
            and open a new dump_reader to re-read the dump file.
*/
std::size_t number_of_blocks( dump_reader &dr );


} // namespace readers

} // namespace lammps_tools

#endif // DUMP_READER_HPP
