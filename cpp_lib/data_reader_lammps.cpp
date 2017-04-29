#include "id_map.hpp"
#include "readers.hpp"
#include "util.hpp"

#include <string>
#include <sstream>

namespace data_readers {

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
		std::stringstream ss(line);
		bool line_handled = false;
		
		if( line.empty() ) continue;
		
		std::vector<std::string> words = split(line);
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


int read_data_atoms_atomic( std::istream &in, block_data &b, bool quiet )
{
	std::string line;
	std::getline(in,line);
	std::vector<std::string> words = split(line);
	// Check size.
	std::size_t n_cols = words.size();
	my_assert( __FILE__, __LINE__, n_cols == 5 || n_cols == 8,
	           "Incorrect column count for atomic!" );
	if( !quiet ) std::cerr << "    ....Allocating storage for " << b.N << " worth of atoms.\n";

	data_field_int id  ( "id",   b.N );
	data_field_int type( "type", b.N );
	data_field_double x( "x",    b.N );
	data_field_double y( "y",    b.N );
	data_field_double z( "z",    b.N );

	data_field_int ix( "iz" );
	data_field_int iy( "iy" );
	data_field_int iz( "iz" );

	if( n_cols == 8 ){
		ix.resize(b.N);
		iy.resize(b.N);
		iz.resize(b.N);
	}

	for( std::size_t i = 0; i < b.N; ++i ){
		std::stringstream ss(line);
		ss >> id[i] >> type[i] >> x[i] >> y[i] >> z[i];

		if( n_cols == 8 ){
			ss >> ix[i] >> iy[i] >> iz[i];
		}
		
		std::getline(in,line);
	}

	
	b.add_field( id );
	b.add_field( type );
	b.add_field( x );
	b.add_field( y );
	b.add_field( z );

	return 0;
}

int read_data_atoms_molecular( std::istream &in, block_data &b, bool quiet )
{
	return 0;
}


int read_data_atoms_velocities( std::istream &in, block_data &b, bool quiet )
{
	std::string line;
	std::getline(in,line);
	std::vector<std::string> words = split(line);
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
	
	for( std::size_t i = 0; i < b.N; ++i ){
		std::stringstream ss(line);
		int id;
		ss >> id;
		int idx = im [ id ];
		ss >> vx[idx] >> vy[idx] >> vz[idx];
		std::getline(in,line);
	}

	b.add_field( vx );
	b.add_field( vy );
	b.add_field( vz );

	
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
		std::vector<std::string> words = split(line);
		std::string keyword = words[0];
		if( keyword == "Masses" ){
			if( !quiet ) std::cerr << "    ....Reading Masses...\n";
			std::getline(in,line);
			std::getline(in,line);
			for( std::size_t i = 0; i < b.N_types; ++i ){
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
				b.atom_style = block_data::ATOMIC;
				status = read_data_atoms_atomic( in, b, quiet );
			}else if( words[2] == "molecular" ){
				b.atom_style = block_data::MOLECULAR;
				status = read_data_atoms_molecular( in, b, quiet );
			}
			my_assert( __FILE__, __LINE__, status == 0,
			           "Error reading in positions!" );
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
			for( std::size_t i = 0; i < b.N_types; ++i ){
				std::stringstream ss(line);
				int type;
				// Determine the number of coeffs.
				int ncoeffs = split(line).size() - 1;
				std::vector<double> coeffs(ncoeffs);
				for( int j = 0; j < ncoeffs; ++j ){
					ss >> coeffs[j];
				}
				std::getline(in,line);
			}
			std::getline(in,line);
		}else if( keyword == "PairIJ" ){
			
		}else{
			
		}

		
		if( !quiet ){
			std::cerr << "    ....Block_data now has "
			          << b.n_data_fields() << " data fields:";
			for( const std::string &header : b.get_data_names() ){
				std::cerr << " " << header;
			}
			std::cerr << "\n";
		}
		
	}
	return 0;
}


/**
   Reads in a block_data from a LAMMPS data file.
   It assumes the file is correctly formatted!

   \param in      Input stream to read from.
   \param status  Will contain 0 on success, non-zero integer otherwise.

   \returns A new block_data containing the info from the given data file.
*/
block_data block_data_from_lammps_data( std::istream &in, int &status,
                                        bool quiet )
{
	if( !quiet ) std::cerr << "Reading LAMMPS data....\n";
	block_data b;
	std::string line;
	// Ignore first two lines.
	std::getline(in,line);
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

} // namespace data_readers
