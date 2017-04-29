#ifndef DUMP_READER_LAMMPS_BIN_HPP
#define DUMP_READER_LAMMPS_BIN_HPP

#include "dump_reader_lammps.hpp"

#include <string>
#include <cstdio>

// This code is mostly based on the binary2txt tool provided with LAMMPS.

namespace dump_readers {

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
	dump_reader_lammps_bin( const std::string &fname );

	/**
	   Initialises dump reader from file. 

	   \param fname Name of dump file.
	   \param h     Vector containing the columns contained in the file.
	*/
	dump_reader_lammps_bin( const std::string &fname,
	                        std::vector<std::string> h );
	

	/// Cleanup:
	virtual ~dump_reader_lammps_bin();
	
private:
	virtual int  get_next_block( block_data &block );
	virtual bool check_eof()  const;
	virtual bool check_good() const;

	int next_block_meta( block_data &block, int &size_one, int &nchunk );
	int next_block_body( block_data &block, int size_one, int nchunk );

	std::FILE* in;
};

}



#endif // DUMP_READER_LAMMPS_BIN_HPP
