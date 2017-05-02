#ifndef DUMP_READER_LAMMPS_HPP
#define DUMP_READER_LAMMPS_HPP

// LAMMPS dump readers come in three flavours: Plain, gzip and bin.
// They all share this general interface.

#include <iosfwd>

#include "dump_reader.hpp"

namespace lammps_tools {

class block_data;

namespace dump_readers {

class dump_reader_lammps : public dump_reader
{
public:
	enum dump_styles { ATOMIC,
	                   CUSTOM };
	
	/// Empty constructor
	dump_reader_lammps(){}
	

	/// Empty destructor
	virtual ~dump_reader_lammps() {}

	/// Sets up a vector containing the expected columns.
	void set_column_headers( const std::vector<std::string> &headers );

	/// Returns a vector containing the expected column headers.
	const std::vector<std::string> &get_column_headers() const;
	
	
private:
	virtual int  get_next_block( lammps_tools::block_data &block ) = 0;
	virtual bool check_eof()  const = 0;
	virtual bool check_good() const = 0;

	std::vector<std::string> column_headers;

};


/// Checks if given column header corresponds to an integer quantity.
bool is_int_data_field( const std::string header );

} // namespace dump_readers

} // namespace lammps_tools

#endif // DUMP_READER_LAMMPS_HPP
