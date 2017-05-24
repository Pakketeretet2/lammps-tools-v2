#ifndef DUMP_READER_LAMMPS_BIN_HPP
#define DUMP_READER_LAMMPS_BIN_HPP

/**
   \file dump_reader_lammps_bin.hpp

   Declarations for binary lammps dump file reader, which is mostly
   based on the binary2txt tool that comes with lammps.
*/

#include "dump_reader_lammps.hpp"

#include <cstdio>
#include <string>
#include <memory>

namespace lammps_tools {

namespace readers {

/**
   A dump reader for binary LAMMPS dump files.
*/
class dump_reader_lammps_bin : public dump_reader_lammps
{
public:
	/**
	   Initialises dump reader from file.

	   \param fname Name of dump file.

	   If no headers are provided, it assumes a standard atomic format.
	*/
	dump_reader_lammps_bin( const std::string &fname,
	                        int dump_style = dump_reader_lammps::CUSTOM );

	/**
	   Initialises dump reader from file.

	   \param fname Name of dump file.
	   \param h     Vector containing the columns contained in the file.
	*/
	dump_reader_lammps_bin( const std::string &fname,
	                        std::vector<std::string> h,
	                        int dump_style = dump_reader_lammps::CUSTOM );


	/// Cleanup:
	virtual ~dump_reader_lammps_bin();

private:
	virtual int  get_next_block( block_data &block );
	virtual bool check_eof()  const;
	virtual bool check_good() const;

	int next_block_meta( block_data &block, int &size_one, int &nchunk );
	int next_block_body( block_data &block, int size_one, int nchunk );

	dump_reader_lammps_bin( dump_reader_lammps_bin& ) = delete;
	dump_reader_lammps_bin &operator=(dump_reader_lammps_bin&) = delete;

	std::FILE *in;
};

} // namespace dump_readers

} // namespace lammps_tools


#endif // DUMP_READER_LAMMPS_BIN_HPP
