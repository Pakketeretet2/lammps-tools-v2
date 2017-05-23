#include "dump_reader_lammps_plain.hpp"
#include "types.hpp"
#include "util.hpp"
#include "my_assert.hpp"

#include <fstream>


using namespace lammps_tools;
using namespace readers;


using lammps_tools::util::starts_with;
using lammps_tools::util::ends_with;
using lammps_tools::util::split;
using lammps_tools::util::contains;

namespace lammps_tools {

namespace readers {

dump_reader_lammps_plain::dump_reader_lammps_plain( const std::string &fname,
                                                    int dump_style )
	: dump_reader_lammps( dump_style ), in_file( nullptr ), in( nullptr )
{
	in_file = new std::ifstream( fname );
	my_assert( __FILE__, __LINE__, in_file,
	           "Failed to open input stream!" );
	in = in_file;
}

dump_reader_lammps_plain::dump_reader_lammps_plain( std::istream &istream,
                                                    int dump_style )
	: dump_reader_lammps( dump_style ), in_file( nullptr ), in( nullptr )
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
	if( status > 0 ){
		if( !quiet ) std::cerr << "....EOF reached....\n";
		return status;
	}else if( status < 0 ){
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
	double xlo[3] = {0,0,0};
	double xhi[3] = {0,0,0};
	tstep = N = 0;

	bool setup_box_and_get_body = false;

	while( get_line( line ) ){
		if( line == "ITEM: TIMESTEP" ){
			get_line( line );
			tstep = std::stoul( line );
			if( !quiet ) std::cerr << "    ....t = " << tstep << "\n";
		}else if( line == "ITEM: NUMBER OF ENTRIES" ){
			my_assert( __FILE__, __LINE__, dump_style == LOCAL,
			           "Got NUMBER OF ENTRIES in dump_style atom!" );
			get_line( line );
			N = std::stoul( line );
			if( !quiet ) std::cerr << "    ....N = " << N << "\n";
		}else if( line == "ITEM: NUMBER OF ATOMS" ){
			my_assert( __FILE__, __LINE__, dump_style != LOCAL,
			           "Got NUMBER OF ATOMS in dump_style local!" );
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
			my_assert( __FILE__, __LINE__, dump_style != LOCAL,
			           "Got ITEM: ATOMS in dump_style local!" );

			// Stop there.
			last_line = line;
			setup_box_and_get_body = true;
		}else if( starts_with( line, "ITEM: ENTRIES" ) ){
			my_assert( __FILE__, __LINE__, dump_style == LOCAL,
			           "Got NUMBER OF ENTRIES in dump_style atom!" );
			if( !quiet ) std::cerr << "Got number of entries!\n";

			// Stop there.
			last_line = line;
			setup_box_and_get_body = true;
		}else{
			std::cerr << "!! Encountered line '" << line
			          << "' and have no clue what to do !!\n";
			return -1;
		}

		if( setup_box_and_get_body ){
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
				          << "Atom or Entry section! !!\n";
				return -1;
			}
		}

	}
	if( check_eof() ){
		return 1;
	}else{
		return -1;
	}
}

void dump_reader_lammps_plain::set_custom_data_fields(
	block_data &block, const std::string &line,
	std::vector<std::string> &headers,
	std::vector<data_field*> &data_fields )
{
	using dfi = data_field_int;
	using dfd = data_field_double;

	if( dump_style == CUSTOM ){
		my_assert( __FILE__, __LINE__,
		           util::starts_with( line, "ITEM: ATOMS" ),
		           "Wrong line passed to set_custom_data_fields!" );
	}else if( dump_style == LOCAL ){
		my_assert( __FILE__, __LINE__,
		           util::starts_with( line, "ITEM: ENTRIES" ),
		           "Wrong line passed to set_custom_data_fields!" );

	}
	std::vector<std::string> words = split(line);
	std::vector<std::string> found_headers( words.begin()+2, words.end() );

	const std::vector<std::string> &col_headers = get_column_headers();

	if( col_headers.empty() ){
		// Add word to column headers.
		set_column_headers( found_headers );
		for( const std::string &w : found_headers ){
			// Assume they are doubles, unless is_int_data_field(w)
			// is true
			if( is_int_data_field(w) ){
				set_column_type( w, data_field::INT );
			}else{
				set_column_type( w, data_field::DOUBLE );
			}
		}
	}

	for( const std::string &w : found_headers ){
		if( !contains( col_headers, w ) ){
			std::cerr << "Column header " << w << " encountered "
			          << " but was not in column headers!\n";
			my_runtime_error(__FILE__, __LINE__,
			                 "Column header mismatch");
		}

		headers.push_back(w);
		if( !quiet ) std::cerr << "      ....Current header is \""
		                       << w << "\"....\n";
		// Depending on the keyword, you want to take
		// either an int or a double.


		int type = get_column_type( w );

		if( type == data_field::INT ){
			dfi *new_field = new dfi( w, block.N );
			data_fields.push_back( new_field );
		}else{
			// Assume it's a double.
			dfd *new_field = new dfd( w, block.N );
			data_fields.push_back( new_field );
		}

		if( header_to_special_field.count(w) ){
			if( header_to_special_field[w] == block_data::MOL ){
				block.atom_style = ATOM_STYLE_MOLECULAR;
			}
		}
	}
}


void dump_reader_lammps_plain::append_data_to_fields(
	block_data &block, std::vector<data_field*> &data_fields )
{
	std::size_t n_cols = data_fields.size();
	std::string line;
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
}



int dump_reader_lammps_plain::next_block_body(
	block_data &block, const std::string &last_line )
{
	using dfi = data_field_int;
	using dfd = data_field_double;

	std::string line = last_line;

	block.set_natoms( block.N );

	if( starts_with( line, "ITEM: ATOMS" ) ||
	    starts_with( line, "ITEM: ENTRIES" ) ){


		// Figure out which column maps which.
		// Read out the next block.N lines.
		std::vector<std::string> headers;
		std::vector<data_field*> data_fields;


		if( line == "ITEM: ATOMS" ){
			if( !quiet ) std::cerr << "    ....Reading atoms....\n";
			my_assert( __FILE__, __LINE__,
			           dump_style == ATOMIC,
			           "Inconsistent dump style for atomic dump!" );

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
		}else if( util::starts_with( line, "ITEM: ATOMS " ) ){
			my_assert( __FILE__, __LINE__,
			           dump_style == CUSTOM,
			           "Inconsistent dump style for custom dump!" );
			if( !quiet ) std::cerr << "    ....Reading atoms....\n";
			set_custom_data_fields( block, line, headers,
			                        data_fields );
		}else if( util::starts_with( line, "ITEM: ENTRIES" ) ){

			my_assert( __FILE__, __LINE__,
			           dump_style == LOCAL,
			           "Inconsistent dump style for local dump!" );

			if( !quiet ) std::cerr << "    ....Reading entries....\n";
			set_custom_data_fields( block, line, headers,
			                        data_fields );
		}


		std::size_t n_cols = headers.size();
		my_assert( __FILE__, __LINE__, n_cols == data_fields.size(),
		           "# of data fields does not match # of columns!" );
		append_data_to_fields( block, data_fields );


		if( dump_style == ATOMIC ){
			block.add_field( *data_fields[0], block_data::ID );
			block.add_field( *data_fields[1], block_data::TYPE );
			block.add_field( *data_fields[2], block_data::X );
			block.add_field( *data_fields[3], block_data::Y );
			block.add_field( *data_fields[4], block_data::Z );

			for( data_field *df : data_fields ) delete df;
		}else{
			add_custom_data_fields( data_fields, block );
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

} // namespace readers

} // namespace lammps_tools
