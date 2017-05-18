#ifndef DUMP_READER_HPP
#define DUMP_READER_HPP

/**
   \file dump_reader.hpp

   Contains declaration of general dump_reader interface and helper functions.
*/

#include "block_data.hpp"

#include <iosfwd>
#include <memory>


namespace lammps_tools {

/// Contains functions and classes that are related to reading dump files.
namespace readers {

class dump_reader_lammps;

/// Specifies various file formats
enum FILE_FORMATS {
	PLAIN   = 0,  ///< Plain text
	GZIP    = 1,  ///< Gzipped plain text
	BIN     = 2   ///< Binary
};

/// Specifies various dump formats
enum DUMP_FORMATS {
	LAMMPS  = 0,   ///< LAMMPS format
	HOOMD   = 1,   ///< HOOMD GSD format
	NAMD    = 2    ///< NAMD DCD format
};


/// A generic class for reading in dump files.
class dump_reader
{
public:
	/// Empty constructor; this just defines the interface.
	dump_reader() : quiet(true), current_block(-1) {}

	/// Empty destructor:
	virtual ~dump_reader(){}

	/**
	   Attempts to read in the next block from file.

	   \param block   Will contain the new block_data on success, and shall
	                  be  unchanged upon failure.

           \returns 0 on success, positive if EOF reached, negative on failure.
	*/
	int next_block( block_data &block ) { return get_next_block( block ); }

	/// Checks if the internal file is at eof:
	bool eof()  const { return check_eof(); }

	/// Checks if the internal file is good:
	bool good() const { return check_good(); }

	bool quiet;

private:
	int current_block;

	virtual int  get_next_block( block_data &block ) = 0;
	virtual bool check_eof()  const = 0;
	virtual bool check_good() const = 0;

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
   \param dformat  Dump format (see \p DUMP_FORMATS).
   \param fformat  File format (see \p FILE_FORMATS).

   \returns a pointer to the open dump_reader object, or nullptr on failure.

   \warning The dump_reader object is heap-allocated. Be sure to delete it
            or wrap it in a smart pointer.

*/
dump_reader* make_dump_reader( const std::string &fname,
                               int dformat, int fformat );


/**
   \brief Constructs correct dump_reader for input stream of given dump format.

   \param fname    Name of the dump file.
   \param dformat  Dump format (see \p DUMP_FORMATS).
   \param fformat  File format (see \p FILE_FORMATS).

   \returns a pointer to the open dump_reader object, or nullptr on failure.

   \warning The dump_reader object is heap-allocated. Be sure to delete it
            or wrap it in a smart pointer.
*/
dump_reader* make_dump_reader( std::istream &input,
                               int dformat, int fformat );


/**
   \brief Constructs a LAMMPS dump reader from file name.

   \param fname    Name of the dump file.
   \param fformat  File format (see \p FILE_FORMATS).
   \param header   Optionally you can pass the expected headers.
                   If the file is binary these are assumed to be
                   the column headers. For other file formats these
                   _have_ to match the column headers in the file
                   and are used for verification.

  \returns a pointer to the open dump_reader object, or nullptr on failure.

  \warning The dump_reader object is heap-allocated. Be sure to delete it
           or wrap it in a smart pointer.
*/
dump_reader_lammps *make_dump_reader_lammps( const std::string &fname, int fformat,
                                             std::vector<std::string> headers );

/**
   \brief Constructs a LAMMPS dump reader from input file stream.

   \param fname    Name of the dump file.
   \param header   Optionally you can pass the expected headers.
                   If the file is binary these are assumed to be
                   the column headers. For other file formats these
                   _have_ to match the column headers in the file
                   and are used for verification.

  \returns a pointer to the open dump_reader object, or nullptr on failure.

  \warning The dump_reader object is heap-allocated. Be sure to delete it
           or wrap it in a smart pointer.
*/
dump_reader_lammps *make_dump_reader_lammps( std::istream &input,
                                             std::vector<std::string> headers );


/**
   \brief Constructs a LAMMPS dump reader from file name.
   \overloads make_dump_reader_lammps
*/
dump_reader_lammps *make_dump_reader_lammps( const std::string &fname, int fformat );


/**
   \brief Constructs a LAMMPS dump reader from input file stream.
   \overloads make_dump_reader_lammps
*/
dump_reader_lammps *make_dump_reader_lammps( std::istream &input );


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
