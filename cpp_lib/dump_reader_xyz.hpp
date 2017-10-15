#ifndef DUMP_READER_XYZ_HPP
#define DUMP_READER_XYZ_HPP

/**
   \file dump_reader_xyz.hpp

   Contains declaration for a dump reader for XYZ files.
*/


#include "dump_reader.hpp"
#include "block_data.hpp"

#include <iosfwd>
#include <memory>


namespace lammps_tools {

class block_data;

/// Contains functions and classes that are related to reading dump files.
namespace readers {


/// A generic class for reading in dump files.
class dump_reader_xyz : public dump_reader
{
public:
	/// Constructor that takes file name
	dump_reader_xyz( const std::string &fname );
	/// Constructor that takes input stream
	dump_reader_xyz( std::istream &istream );
	/// Empty destructor
	virtual ~dump_reader_xyz();

private:
	virtual int  get_next_block( block_data &block );
	virtual bool check_eof()  const;
	virtual bool check_good() const;

	std::ifstream in_file;
	std::istream &in;
	bool no_lammps_warned;
};

} // namespace readers

} // namespace lammps_tools


#endif //DUMP_READER_XYZ_HPP
