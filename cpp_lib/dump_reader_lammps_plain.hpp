#ifndef DUMP_READER_LAMMPS_PLAIN_HPP
#define DUMP_READER_LAMMPS_PLAIN_HPP

/**
   \file dump_reader_lammps_gzip.hpp
   
   Declaration of dump reader for lammps plain text dump files.
*/

#include "dump_reader_lammps.hpp"

#include <string>
#include <iosfwd>

namespace lammps_tools {

namespace dump_readers {

class dump_reader_lammps_plain : public dump_reader_lammps
{
public:
	/// Initialises dump reader from file. 
	dump_reader_lammps_plain( const std::string &fname );
	/// Initialises dump reader from input stream.
	dump_reader_lammps_plain( std::istream &istream );
	/// Cleanup:
	virtual ~dump_reader_lammps_plain();

private:
	virtual int  get_next_block( block_data &block );
	virtual bool check_eof()  const;
	virtual bool check_good() const;

	virtual bool get_line( std::string &line );

	int next_block_meta( block_data &block, std::string &last_line );
	int next_block_body( block_data &block, const std::string &last_line );

	std::istream  *in;
	std::ifstream *in_file;
};

} // namespace dump_readers

} // namespace lammps_tools

#endif // DUMP_READER_LAMMPS_PLAIN_HPP
