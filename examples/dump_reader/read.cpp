/*
  This C++ program is functionally identical to the Python
  script in terms of functionality. Use these to see how
  much slower (or faster) Python is.

  Compile as
  $CC -O3 -std=c++11 read.cpp -I../../cpp_lib -llammpstools -o read_cpp
*/

#include "block_data.hpp"
#include "block_data_access.hpp"
#include "dump_reader_lammps_plain.hpp"
#include "data_field.hpp"

int main( int argc, char **argv )
{
	using namespace lammps_tools;

	using readers::dump_reader_lammps_plain;
	using readers::dump_reader_lammps;

	dump_reader_lammps_plain dl( "bondinfo.dump", dump_reader_lammps::LOCAL );
	dump_reader_lammps_plain d ( "polymer.dump",  dump_reader_lammps::CUSTOM );

	block_data b, bl;
	std::vector<std::string> headers = {"id", "mol", "type", "x", "y", "z"};
	std::vector<int> types = { block_data::ID, block_data::MOL,
	                           block_data::TYPE, block_data::X,
	                           block_data::Y, block_data::Z };
	for( std::size_t i = 0; i < headers.size(); ++i ){
		d.set_column_header( i, headers[i], types[i] );
	}

	dl.default_col_type = data_field::INT;

	int c = 0;
	int status = d.next_block(b);
	status |= dl.next_block(bl);
	while( status == 0 ){
		if( c > 0 && c % 25 == 0 ){

			std::cout << "At block " << c << ", t = " << b.tstep
			          << " contains N = " << b.N << " particles and "
			          << bl.N << " bonds\n";
			int i = 1233;
			if( b.N > i ){
				const std::vector<double> &x = get_x( b );
				const std::vector<double> &y = get_y( b );
				const std::vector<double> &z = get_z( b );

				std::cout << "Position of particle " << i << " is [ "
				          << x[i] << " " << y[i] << " " << z[i]
				          << " ]\n";
			}
			if( bl.N > i ){
				std::cout << "data of bond " << i << " is [ ";
				for( std::size_t j = 0; j < bl.n_data_fields(); ++j ){
					const std::vector<int> &data =
						data_as<int>( &bl[j] );
					std::cout << data[i] << " ";
				}
				std::cout << "]\n";
			}
		}

		status = d.next_block(b);
		status |= dl.next_block(bl);
		++c;
	}
	return 0;
}
