#ifndef DUMP_READER_LAMMPS_PLAIN_HPP
#define DUMP_READER_LAMMPS_PLAIN_HPP

/**
   \file dump_reader_lammps_gzip.hpp

   Declaration of dump reader for lammps plain text dump files.
*/

#include "dump_reader_lammps.hpp"

#include <memory>
#include <string>
#include <iosfwd>

namespace lammps_tools {

namespace readers {

class dump_reader_lammps_plain : public dump_reader_lammps
{
public:

	/// Initialises dump reader from file.
	dump_reader_lammps_plain( const std::string &fname,
	                          int dump_style = dump_reader_lammps::CUSTOM );
	/// Initialises dump reader from input stream.
	dump_reader_lammps_plain( std::istream &istream,
	                          int dump_style = dump_reader_lammps::CUSTOM );
	/// Cleanup:
	virtual ~dump_reader_lammps_plain();

private:
	virtual int  get_next_block( block_data &block );
	virtual bool check_eof()  const;
	virtual bool check_good() const;

	virtual bool get_line( std::string &line );

	int next_block_meta( block_data &block, std::string &last_line );
	int next_block_body( block_data &block, const std::string &last_line );

	// Extracts the custom names for each data field from given line.
	void set_custom_data_fields( block_data &block,
	                             const std::string &line,
	                             std::vector<std::string> &headers,
	                             std::vector<data_field*> &data_fields );

	// Set custom fields for dump local.
	void set_local_data_fields( block_data &block,
	                            const std::string &line,
	                            std::vector<std::string> &headers,
	                            std::vector<data_field*> &data_fields );



	// Reads the info of the current block into block_data from file.
	void append_data_to_fields( block_data &block,
	                            std::vector<data_field*> &data_fields );

	dump_reader_lammps_plain &operator=(dump_reader_lammps_plain&) = delete;
	dump_reader_lammps_plain(dump_reader_lammps_plain&) = delete;

	std::ifstream *in_file;
	std::istream *in;
};

} // namespace readers

} // namespace lammps_tools

#endif // DUMP_READER_LAMMPS_PLAIN_HPP
