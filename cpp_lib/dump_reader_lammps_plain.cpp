#include "dump_reader_lammps_plain.hpp"
#include "types.hpp"
#include "util.hpp"

#include <fstream>



using namespace dump_readers;




dump_reader_lammps_plain::dump_reader_lammps_plain( const std::string &fname )
	: in_file( nullptr ), in( nullptr )
{
	in_file = new std::ifstream( fname );
	my_assert( __FILE__, __LINE__, in_file,
	           "Failed to open input stream!" );
	in = in_file;
}

dump_reader_lammps_plain::dump_reader_lammps_plain( std::istream &istream )
	: in_file( nullptr ), in( nullptr )
{
	in = &istream;
	my_assert( __FILE__, __LINE__, in,
	           "Failed to open input stream!" );
}

dump_reader_lammps_plain::~dump_reader_lammps_plain()
{
	// Works for either case because in always points to the open stream.
	if( in_file ) delete in_file;
}

int dump_reader_lammps_plain::get_next_block( block_data &block )
{
	if( !quiet ) std::cerr << "Reading block from LAMMPS dump file....\n";
	
	std::string last_line = "";
	block_data tmp_block;
	if( !quiet ) std::cerr << "  ....Reading block meta....\n";
	int status = next_block_meta( tmp_block, last_line );
	if( status ){
		std::cerr << "!! An error occured while reading meta !!\n";
		return status;
	}
	if( !quiet ) std::cerr << "  ....Reading block body....\n";
	status = next_block_body( tmp_block, last_line );
	if( status ) return status;

	// At this point, tmp_block should be in a 100% correct state.
	
	
	block = tmp_block;
	
	return 0;
}

bool dump_reader_lammps_plain::check_eof() const
{
	return in->eof();	
}
bool dump_reader_lammps_plain::check_good() const
{
	return in->good();
}


int dump_reader_lammps_plain::next_block_meta( block_data &block,
                                               std::string &last_line )
{
	std::string line;
	bigint tstep, N;
	std::string boxline;
	double xlo[3], xhi[3];

	while( get_line( line ) ){
		if( line == "ITEM: TIMESTEP" ){
			get_line( line );
			tstep = std::stoul( line );
			if( !quiet ) std::cerr << "    ....t = " << tstep << "\n";
		}else if( line == "ITEM: NUMBER OF ATOMS" ){
			get_line( line );
			N = std::stoul( line );
			if( !quiet ) std::cerr << "    ....N = " << N << "\n";
		}else if( starts_with( line, "ITEM: BOX BOUNDS " ) ){
			boxline = line;
			get_line( line );
			std::stringstream dims( line );
			dims >> xlo[0] >> xhi[0];
			
			get_line( line );
			dims.str(""); dims.clear();
			dims.str( line );
			dims >> xlo[1] >> xhi[1];

			get_line( line );
			dims.str(""); dims.clear();
			dims.str( line );
			dims >> xlo[2] >> xhi[2];
			if( !quiet )
				std::cerr << "    ....box = [ " << xlo[0]
				          << ", " << xhi[0] << " ] x [ "
				          << xlo[1] << ", " << xhi[1]
				          << " ] x [ " << xlo[2] << ", "
				          << xhi[2] << " ].\n";
		}else if( starts_with( line, "ITEM: ATOMS" ) ){
			// Stop there.
			last_line = line;
			if( check_good() ){

				block.tstep = tstep;
				block.N = N;
				block.dom.xlo[0] = xlo[0];
				block.dom.xlo[1] = xlo[1];
				block.dom.xlo[2] = xlo[2];

				block.dom.xhi[0] = xhi[0];
				block.dom.xhi[1] = xhi[1];
				block.dom.xhi[2] = xhi[2];

				block.dom.periodic = 0;
				std::stringstream bl( boxline );
				std::string tmp;
				bl >> tmp >> tmp >> tmp; // ITEM: BOX BOUNDS
				bl >> tmp;
				if( tmp == "pp" ){
					block.dom.periodic += domain::BIT_X;
				}
				bl >> tmp;
				if( tmp == "pp" ){
					block.dom.periodic += domain::BIT_Y;
				}
				bl >> tmp;
				if( tmp == "pp" ){
					block.dom.periodic += domain::BIT_Z;
				}
				return 0;
			}else{
				std::cerr << "!! File no longer good after "
				          << "encountering ATOMS !!";
				return -1;
			}
		}else{
			std::cerr << "!! Encountered line '" << line
			          << "' and have no clue what to do !!\n";
			return -1;
		}
	}
	if( check_eof() ){
		return 1;
	}else{
		return -1;
	}
}

int dump_reader_lammps_plain::next_block_body( block_data &block,
                                               const std::string &last_line )
{
	using dfi = data_field_int;
	using dfd = data_field_double;
		
	std::string line = last_line;

	block.set_natoms( block.N );
	int dump_style = ATOMIC;

	const std::vector<std::string> &col_headers = get_column_headers();


	if( starts_with( line, "ITEM: ATOMS" ) ){
		if( !quiet ) std::cerr << "    ....Reading atoms....\n";
		// Figure out which column maps which.
		// Read out the next block.N lines.
		std::vector<std::string> headers;
		std::vector<data_field*> data_fields;

		
		if( line == "ITEM: ATOMS" ){
			dump_style = ATOMIC;
			headers.push_back("id");
			headers.push_back("type");
			headers.push_back("x");
			headers.push_back("y");
			headers.push_back("z");

			dfi *id = new dfi( "id",   block.N );
			dfi *type = new dfi( "type", block.N );
			dfd *x = new dfd( "x", block.N );
			dfd *y = new dfd( "y", block.N );
			dfd *z = new dfd( "z", block.N );

			data_fields.push_back( id );
			data_fields.push_back( type );
			data_fields.push_back( x );
			data_fields.push_back( y );
			data_fields.push_back( z );
			
		}else{
			dump_style = CUSTOM;
			std::vector<std::string> words = split(line);
			for( int i = 2; i < words.size(); ++i ){
				std::string w = words[i];
				if( !col_headers.empty() ){
					my_assert( __FILE__, __LINE__,
					           contains( col_headers, w ),
					           "Found column not in given "
					           "column headers!" );
				}
				
				headers.push_back(w);
				if( !quiet ) std::cerr << "      ....Current header is \"" << w << "\"....\n";
				// Depending on the keyword, you want to take
				// either an int or a double.
				if( is_int_data_field( w ) ){
					dfi *new_field = new dfi( w, block.N );
					data_fields.push_back( new_field );
				}else{
					// Assume it's a double.
					dfd *new_field = new dfd( w, block.N );
					data_fields.push_back( new_field );					
				}

				if( w == "mol" ){
					// Assume atom_style is molecular.
					block.atom_style = block_data::MOLECULAR;
				}
			}
		}
		std::size_t n_cols = headers.size();
		my_assert( __FILE__, __LINE__, n_cols == data_fields.size(),
		           "# of data fields does not match # of columns!" );

		for( int i = 0; i < block.N; ++i ){
			get_line( line );
			
			std::stringstream ss( line );
			for( int j = 0; j < n_cols; ++j ){
				int type = data_fields[j]->type();
				if( type == data_field::INT ){
					std::vector<int> &vec =
						data_as_rw<int>( data_fields[j] );
					ss >> vec[i];
				}else if( type == data_field::DOUBLE ){
					std::vector<double> &vec =
						data_as_rw<double>( data_fields[j] );
					ss >> vec[i];
				}
			}
		}

		// std::cerr << "Block has atom_style " << block.atom_style << ".\n";

		for( data_field *df : data_fields ){
			block.add_field( *df );
			delete df;
		}
		
		return 0;
		
	}else{
		// std::cerr << "Encountered unknown header!\n";
		// std::cerr << line << "\n";
		return -1;
	}

	
}

bool dump_reader_lammps_plain::get_line( std::string &line )
{
	if( std::getline( *in, line ) ){
		return true;
	}else{
		return false;
	}
}
