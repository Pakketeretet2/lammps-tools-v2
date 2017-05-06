#ifndef DUMP_READER_LAMMPS_GZIP_HPP
#define DUMP_READER_LAMMPS_GZIP_HPP

/**
   \file dump_reader_lammps_gzip.hpp

   Declaration of dump reader for lammps gzip dump files.
*/

#include "dump_reader_lammps_plain.hpp"

#ifdef HAVE_BOOST_GZIP
#  include <boost/iostreams/filter/gzip.hpp>
#  include <boost/iostreams/filtering_stream.hpp>
#endif

#include <fstream>

namespace lammps_tools {

namespace readers {

class dump_reader_lammps_gzip : public dump_reader_lammps_plain
{
public:
	/// Initialises dump reader from file.
	dump_reader_lammps_gzip( const std::string &fname );

	/// Cleanup:
	virtual ~dump_reader_lammps_gzip();

private:
	virtual bool get_line( std::string &line );
	std::ifstream infile;
#ifdef HAVE_BOOST_GZIP
	boost::iostreams::filtering_istream in;
#else
	// Just to trick the compiler. Using this class _will_ lead to a
	// runtime error if GZIP is not compiled in.
	std::ifstream in;
#endif

};

} // namespace readers

} // namespace lammps_tools


#endif // DUMP_READER_LAMMPS_GZIP_HPP
