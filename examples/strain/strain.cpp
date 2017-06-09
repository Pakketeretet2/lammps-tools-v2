#include "dump_reader_lammps_bin.hpp"
#include "skeletonize.hpp"
#include "writers.hpp"

#include <memory>
#include <string>


int main( int argc, char **argv )
{
	using namespace lammps_tools;

	std::string dname = "edt.dump.bin";
	readers::dump_reader_lammps_bin d( dname );
	d.set_column_headers( { "id", "type", "x", "y", "z", "edt" } );

	d.set_column_header_as_special( "id", block_data::ID );
	d.set_column_header_as_special( "type", block_data::TYPE );
	d.set_column_header_as_special( "x", block_data::X );
	d.set_column_header_as_special( "y", block_data::Y );
	d.set_column_header_as_special( "z", block_data::Z );
	block_data b;
	int status = d.next_block(b);
	double rc = 1.25;
	double r0 = 1.0;
	std::ofstream dump_out( "strain.dump.bin", std::ios::binary );
	while( status == 0 ){
		std::vector<double> strain =
			skeletonize::neighbour_strain( b, r0, 1, 1,
			                               neighborize::DIST_BIN,
			                               3, rc );
		b.add_field( data_field_double( "strain", strain ) );
		writers::block_to_lammps_dump( dump_out, b, FILE_FORMAT_BIN );
		status = d.next_block(b);
       }
}
