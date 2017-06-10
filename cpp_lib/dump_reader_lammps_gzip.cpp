#include "dump_reader_lammps_gzip.hpp"

#ifdef HAVE_BOOST_GZIP
#  include <boost/iostreams/filter/gzip.hpp>
#  include <boost/iostreams/filtering_stream.hpp>
#endif

using namespace lammps_tools;
using namespace readers;

namespace lammps_tools {

namespace readers {

#ifdef HAVE_BOOST_GZIP
dump_reader_lammps_gzip::dump_reader_lammps_gzip( const std::string &fname, int dump_style )
	: dump_reader_lammps_plain( fname, dump_style ),
	  infile( fname, std::ios_base::in |std::ios_base::binary ), in()
{
	my_assert( __FILE__, __LINE__, util::file_exists( fname ),
	           "Dump file does not exist!" );

	in.push( boost::iostreams::gzip_decompressor() );
	in.push( infile );

}
#else
dump_reader_lammps_gzip::dump_reader_lammps_gzip( const std::string &fname, int dump_style )
	: dump_reader_lammps_plain( fname, dump_style ),
	  infile( fname, std::ios_base::in |std::ios_base::binary ), in(fname)
{
	my_logic_error( __FILE__, __LINE__, "Gzipped files are not supported "
	                "without boost support! Recompile with HAVE_BOOST_GZIP "
	                "defined and boost installed or gunzip file!" );
}
#endif


dump_reader_lammps_gzip::~dump_reader_lammps_gzip()
{}


bool dump_reader_lammps_gzip::get_line( std::string &line )
{
	if( std::getline( in, line ) ){
		return true;
	}else{
		return false;
	}
}

} // namespace readers

} // namespace lammps_tools
