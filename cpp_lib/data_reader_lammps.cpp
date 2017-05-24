#include "enums.hpp"
#include "id_map.hpp"
#include "readers.hpp"
#include "util.hpp"

#include <string>
#include <sstream>


namespace lammps_tools {


/// Contains stuff that has to do with reading data.
namespace readers {

/**
   Reads in header info from data file in given stream.
   It assumes the file is correctly formatted!

   \param in Input stream to read from.
   \param b  block_data to read into.

   \returns 0 on success, non-zero integer otherwise.
*/
int get_header_info( std::istream &in, block_data &b,
                     std::string &last_line, bool quiet )
{
	std::string line;
	std::vector<std::string> body_headers =
		{ "Atoms", "Velocities", "Masses", "Bonds",
		  "Angles", "Impropers", "Dihedrals", "Pair", "PairIJ" };

	while( in ){
		std::getline(in,line);

		if( line.empty() ) continue;

		std::vector<std::string> words = util::split(line);
		// Search for keywords:
		if( words.size() >= 2 ){
			if( words[1] == "atoms" ){
				b.set_natoms( std::stoul( words[0] ) );
			}else if( words[1] == "atom" && words[2] == "types" ){
				b.N_types = std::stoi( words[0] );
			}
		}
		if( words.size() >= 4 ){
			if( words[2] == "xlo" && words[3] == "xhi" ){
				b.dom.xlo[0] = std::stof( words[0] );
				b.dom.xhi[0] = std::stof( words[1] );
			}else if( words[2] == "ylo" && words[3] == "yhi" ){
				b.dom.xlo[1] = std::stof( words[0] );
				b.dom.xhi[1] = std::stof( words[1] );
			}else if( words[2] == "zlo" && words[3] == "zhi" ){
				b.dom.xlo[2] = std::stof( words[0] );
				b.dom.xhi[2] = std::stof( words[1] );
			}
		}

		if( std::find( body_headers.begin(), body_headers.end(),
		               words[0] ) != body_headers.end() ){
			// Encountered a body section.
			if( !quiet ) std::cerr << "  ....Body keyword encountered.\n";
			last_line = line;
			return 0;
		}
	}

	return -1;
}


int read_data_atoms( std::istream &in, block_data &b, bool quiet )
{
	std::string line;
	std::getline(in,line);
	std::vector<std::string> words = util::split(line);
	// Check size.
	std::size_t n_cols = words.size();
	std::size_t n_col_target = 5;

	if( b.atom_style == ATOM_STYLE_MOLECULAR ){
		n_col_target = 6;
	}
	my_assert( __FILE__, __LINE__,
	           n_cols == n_col_target || n_cols == n_col_target + 3,
	           "Incorrect column count for atomic!" );

	if( !quiet ) std::cerr << "    ....Allocating storage for "
	                       << b.N << " worth of atoms.\n";

	data_field_int id  ( "id",   b.N );
	data_field_int type( "type", b.N );
	data_field_double x( "x",    b.N );
	data_field_double y( "y",    b.N );
	data_field_double z( "z",    b.N );

	data_field_int mol( "mol" );
	data_field_int ix( "ix" );
	data_field_int iy( "iy" );
	data_field_int iz( "iz" );

	bool has_image_flags = (n_cols == n_col_target + 3);

	if( has_image_flags ){
		ix.resize(b.N);
		iy.resize(b.N);
		iz.resize(b.N);
	}
	if( b.atom_style == ATOM_STYLE_MOLECULAR ){
		mol.resize(b.N);
	}

	if( b.atom_style == ATOM_STYLE_ATOMIC ){
		for( bigint i = 0; i < b.N; ++i ){
			std::stringstream ss(line);
			ss >> id[i] >> type[i] >> x[i] >> y[i] >> z[i];

			if( has_image_flags ){
				ss >> ix[i] >> iy[i] >> iz[i];
			}

			std::getline(in,line);
		}
	}else if( b.atom_style == ATOM_STYLE_MOLECULAR ){
		for( bigint i = 0; i < b.N; ++i ){
			std::stringstream ss(line);
			ss >> id[i] >> mol[i] >> type[i] >> x[i] >> y[i] >> z[i];

			if( has_image_flags ){
				ss >> ix[i] >> iy[i] >> iz[i];
			}

			std::getline(in,line);
		}
	}
	if( !quiet ) std::cerr << "    ....Adding fields.\n";

	b.add_field(   id, block_data::ID );
	if( b.atom_style == ATOM_STYLE_MOLECULAR ){
		b.add_field( mol, block_data::MOL );
	}
	b.add_field( type, block_data::TYPE );
	b.add_field(    x, block_data::X );
	b.add_field(    y, block_data::Y );
	b.add_field(    z, block_data::Z );

	if( has_image_flags ){
		b.add_field( ix, block_data::IX );
		b.add_field( iy, block_data::IY );
		b.add_field( iz, block_data::IZ );
	}

	return 0;
}


int read_data_atoms_velocities( std::istream &in, block_data &b, bool quiet )
{
	std::string line;
	std::getline(in,line);
	std::vector<std::string> words = util::split(line);
	// Check size.
	std::size_t n_cols = words.size();
	my_assert( __FILE__, __LINE__, n_cols == 4,
	           "Incorrect column count for velocities!" );

	const data_field *id_base = b.get_data( "id" );
	my_assert( __FILE__, __LINE__, id_base != nullptr,
	           "Block data did not contain atom ids!" );
	const data_field_int *ids = static_cast< const data_field_int* >( id_base );
	id_map im( ids->get_data() );

	data_field_double vx( "vx", b.N );
	data_field_double vy( "vy", b.N );
	data_field_double vz( "vz", b.N );

	for( bigint i = 0; i < b.N; ++i ){
		std::stringstream ss(line);
		int id;
		ss >> id;
		int idx = im [ id ];
		ss >> vx[idx] >> vy[idx] >> vz[idx];
		std::getline(in,line);
	}

	b.add_field( vx, block_data::VX );
	b.add_field( vy, block_data::VY );
	b.add_field( vz, block_data::VZ );

	return 0;
}



/**
   Reads in header info from data file in given stream.
   It assumes the file is correctly formatted!

   \param in Input stream to read from.
   \param b  block_data to read into.

   \returns 0 on success, non-zero integer otherwise.
*/
int get_data_body( std::istream &in, block_data &b,
                   const std::string &last_line, bool quiet )
{
	std::string line = last_line;
	int status = 0;
	while( in ){
		std::vector<std::string> words = util::split(line);
		std::string keyword = words[0];
		if( keyword == "Masses" ){
			if( !quiet ) std::cerr << "    ....Reading Masses...\n";
			std::getline(in,line);
			std::getline(in,line);
			for( int i = 0; i < b.N_types; ++i ){
				std::stringstream ss(line);
				int type;
				double mass;
				ss >> type >> mass;
				std::getline(in,line);
			}
			std::getline(in,line);

		}else if( keyword == "Atoms" ){
			if( !quiet ) std::cerr << "    ....Reading atoms...\n";
			std::getline(in,line);
			if( words[2] == "atomic" ){
				b.atom_style = ATOM_STYLE_ATOMIC;
			}else if( words[2] == "molecular" ){
				b.atom_style = ATOM_STYLE_MOLECULAR;
			}
			status = read_data_atoms( in, b, quiet );
			my_assert( __FILE__, __LINE__, status == 0,
			           "Error reading in positions!" );
			if( !quiet ) std::cerr << "    ....Read in " << b.N
			                       << " atoms.\n";
			std::getline(in,line);

		}else if( keyword == "Velocities" ){
			if( !quiet ) std::cerr << "    ....Reading velocities...\n";
			std::getline(in,line);
			status = read_data_atoms_velocities( in, b, quiet );
			my_assert( __FILE__, __LINE__, status == 0,
			           "Error reading in velocities!" );
			std::getline(in,line);
		}else if( keyword == "Bonds" ){

		}else if( keyword == "Angles" ){

		}else if( keyword == "Impropers" ){

		}else if( keyword == "Dihedrals" ){

		}else if( keyword == "Pair" ){
			if( !quiet ) std::cerr << "    ....Reading pair coeffs...\n";
			std::getline(in,line);
			std::getline(in,line);
			for( int i = 0; i < b.N_types; ++i ){
				std::stringstream ss(line);
				// Determine the number of coeffs.
				int ncoeffs = util::split(line).size() - 1;
				std::vector<double> coeffs(ncoeffs);
				for( int j = 0; j < ncoeffs; ++j ){
					ss >> coeffs[j];
				}
				std::getline(in,line);
			}
			std::getline(in,line);
		}else if( keyword == "PairIJ" ){

		}else{
			std::cerr << "I don't know how to handle \""
			          << line << "\"!\n";
			my_runtime_error( __FILE__, __LINE__,
			                  "Unrecognized line encountered!" );
		}


		if( !quiet ){
			std::cerr << "    ....Block_data now has "
			          << b.n_data_fields() << " data fields:";
			for( std::size_t i = 0; i < b.n_data_fields(); ++i ){
				const std::string &header = b[i].name;
				std::cerr << " " << header;
			}
			std::cerr << "\n";
		}
	}
	my_assert( __FILE__, __LINE__, b.n_data_fields() > 0,
	           "Data file had no data fields?" );
	return 0;
}

block_data block_data_from_lammps_data( const std::string &fname, int &status,
                                        bool quiet )
{
	std::ifstream in( fname );
	return block_data_from_lammps_data( in, status, quiet );
}


block_data block_data_from_lammps_data( std::istream &in, int &status,
                                        bool quiet )
{
	if( !quiet ) std::cerr << "Reading LAMMPS data....\n";
	block_data b;
	std::string line;
	// Ignore first line, get_header_info skips blanks so no worries.
	std::getline(in,line);

	while( in ){
		std::string last_line = "";
		status = get_header_info( in, b, last_line, quiet );

		my_assert( __FILE__, __LINE__, status == 0,
		           "Failure in reading data file!" );
		if( !quiet ) std::cerr << "  ....Attempting to read body, starting from "
		                       << last_line << ".\n";
		status = get_data_body( in, b, last_line, quiet );
	}

	if( !quiet ) std::cerr << "Succesfully read LAMMPS data....\n\n";
	return b;
}

} // namespace readers

} // namespace lammps_tools
