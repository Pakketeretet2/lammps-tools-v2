#include "data_field.hpp"
#include "dump_reader_xyz.hpp"

#include <iostream>
#include <fstream>
#include <string>

#include "enums.hpp"
#include "types.hpp"


using namespace lammps_tools;
using namespace readers;

namespace lammps_tools {

namespace readers {


dump_reader_xyz::dump_reader_xyz( const std::string &fname )
	: in_file(fname), in(in_file), no_lammps_warned(false)
{ }

dump_reader_xyz::dump_reader_xyz( std::istream &istream )
	: in(istream), no_lammps_warned(false)
{ }


dump_reader_xyz::~dump_reader_xyz()
{ }


int dump_reader_xyz::get_next_block( block_data &block )
{
	std::string line;
	if( !std::getline(in,line) ){
		// Probably at eof?
		if( eof() ) return 1;
		else return -1;
	}
	block_data tmp;

	std::stringstream ss(line);
	bigint natoms;
	ss >> natoms;
	tmp.set_natoms(natoms);

	ss.str("");
	ss.clear();
	if( !std::getline(in,line) ) return -1;


	ss.str(line);
	std::string first, second;
	ss >> first >> second >> tmp.tstep;
	if( first != "Atoms." && second != "Timestep:" ){
		if( !no_lammps_warned ){
			my_warning( __FILE__, __LINE__,
			            "XYZ file not written by LAMMPS!" );
		}
	}
	ss.clear();

	data_field_double x("x", natoms);
	data_field_double y("y", natoms);
	data_field_double z("z", natoms);
	data_field_int id("id", natoms);
	data_field_int type("type", natoms);

	int max_type = 1;
	bigint i = 0;
	while( std::getline(in,line) ){
		ss.str(line);

		ss >> type[i] >> x[i] >> y[i] >> z[i];
		id[i] = i + 1;
		max_type = type[i] > max_type ? type[i] : max_type;

		++i;
		ss.clear();
		if( i == natoms ) break;
	}

	my_assert( __FILE__, __LINE__, i == natoms,
	           "Incorrect number of atoms read!" );


	if( !good() ) return -1;
	if( eof() ) return 1;

	// Copy all data to the block:
	tmp.set_ntypes( max_type );

	tmp.add_field( id, block_data::special_fields::ID );
	tmp.add_field( type, block_data::special_fields::TYPE );
	tmp.add_field( x, block_data::special_fields::X );
	tmp.add_field( y, block_data::special_fields::Y );
	tmp.add_field( z, block_data::special_fields::Z );

	block = tmp;


	return 0;
}

bool dump_reader_xyz::check_eof() const
{
	return in.eof();
}

bool dump_reader_xyz::check_good() const
{
	return in.good();
}



} // namespace readers

} // namespace lammps_tools
