#include "dump_reader_hoomd_gsd.hpp"

#ifdef HAVE_GSD
#include <gsd.h>
#endif // HAVE_GSD

#include <cstdlib>
#include <cstring>
#include <unistd.h> // For fdopen

using namespace lammps_tools;
using namespace readers;

namespace lammps_tools {

namespace readers {

#ifdef HAVE_GSD
// Some helper functions. Maybe might be useful to put in GSD itself?
// Not portable though probably.
bool gsd_at_eof( gsd_handle *gh )
{
	FILE *fp = fdopen( gh->fd, "r" );
	return std::feof( fp );
}

bool gsd_good( gsd_handle *gh )
{
	FILE *fp = fdopen( gh->fd, "r" );
	return !std::ferror( fp );
}



#endif // HAVE_GSD



dump_reader_hoomd_gsd::dump_reader_hoomd_gsd( const std::string &fname )
	: status(0), gh( nullptr ), current_frame(-1),
	  eof_(true), good_(false), optional_data_found(0)
{

#ifndef HAVE_GSD

	std::string msg = "Cannot read HOOMD GSD files without having ";
	msg += "GSD support compiled in! Enable USE_GSD and recompile!";
	my_runtime_error( __FILE__, __LINE__, msg );

#else

	std::cerr << "Attempting to open gsd file " << fname << ".\n";
	gh = new gsd_handle;

	status = gsd_open( gh, fname.c_str(), GSD_OPEN_READONLY );
	my_assert( __FILE__, __LINE__, status == 0,
	           "Error opening gsd file!" );

	eof_  = false;
	good_ = true;

	max_frame = gsd_get_nframes( gh );
	double size_in_MB = gh->file_size / (1024.0*1024.0);
	std::cerr << "Opened gsd file with " << max_frame
	          << " frames, file size is " << size_in_MB << " MB.\n";

#endif
}

dump_reader_hoomd_gsd::~dump_reader_hoomd_gsd()
{
	// prevents compiler complaint about deleting forward-declared struct:
#ifdef HAVE_GSD
	if( gh ){
		std::cerr << "Closing file @ " << gh << ".\n";
		gsd_close( gh );
		delete gh;
	}
#endif
}

bool dump_reader_hoomd_gsd::check_eof()  const
{
	return eof_;
}

bool dump_reader_hoomd_gsd::check_good() const
{
	return good_;
}

#ifndef HAVE_GSD

// Just to make it compile:
template <typename T>
int dump_reader_hoomd_gsd::get_chunk_data( const char *name, T *dest,
                                           uint size, T *store )
{ return -1; }

int dump_reader_hoomd_gsd::get_next_block( block_data &block )
{ return -1; }

#else



template <typename T>
int dump_reader_hoomd_gsd::get_chunk_data( const char *name, T *dest,
                                           uint size, T *store )
{
	const gsd_index_entry *idx = gsd_find_chunk( gh, current_frame, name );

	if( !idx ){
		// This is not strictly an error, because GSD may store some
		// simulation information in the first frame only. Copy store
		// to dest:
		std::copy( store, store + size, dest );
		return 2;
	}

	// The chunk was definitely present, so try to read it:
	int status = gsd_read_chunk( gh, dest, idx );
	if( status == -1 )
		my_warning( __FILE__, __LINE__, "I/O failure!" );
	else if( status == -2 )
		my_warning( __FILE__, __LINE__, "Invalid input!" );
	else if( status == -3 )
		my_warning( __FILE__, __LINE__, "Invalid data file!" );
	else if( status != 0 ){
		my_warning( __FILE__, __LINE__, "Generic non-zero!" );
	}

	if( status == 0 ){
		// In this case, you succesfully read the data. Also store it:
		std::copy( dest, dest + size, store );
		return 0;
	}

	// In this case, there was some sort of error. Check for EOF:
	if( gsd_at_eof( gh ) ){
		eof_ = true;
		return 1;
	}else{
		good_ = false;
		return status;
	}
}

template <typename T>
int dump_reader_hoomd_gsd::get_chunk_data( const char *name, T &dest, T &store )
{
	const gsd_index_entry *idx = gsd_find_chunk( gh, current_frame, name );

	if( !idx ){
		// This is not strictly an error, because GSD may store some
		// simulation information in the first frame only. Copy store
		// to dest:
		dest = store;
		return 2;
	}

	// The chunk was definitely present, so try to read it:
	my_assert( __FILE__, __LINE__, dest.size() > 0,
	           "Container appears to have zero storage room!" );

	int status = gsd_read_chunk( gh, dest.data(), idx );

	if( status == -1 )
		my_warning( __FILE__, __LINE__, "I/O failure!" );
	else if( status == -2 )
		my_warning( __FILE__, __LINE__, "Invalid input!" );
	else if( status == -3 )
		my_warning( __FILE__, __LINE__, "Invalid data file!" );
	else if( status != 0 ){
		my_warning( __FILE__, __LINE__, "Generic non-zero!" );
	}

	if( status == 0 ){
		// In this case, you succesfully read the data. Also store it:
		store = dest;
		return 0;
	}

	// In this case, there was some sort of error. Check for EOF:
	if( gsd_at_eof( gh ) ){
		eof_ = true;
		return 1;
	}else{
		good_ = false;
		return status;
	}
}




int dump_reader_hoomd_gsd::get_type_names( std::vector<std::string> &type_names )
{
	std::size_t n_types = type_names.size() - 1;


	int buff_size = n_types * gsd::TYPE_BUFFER_SIZE;
	char *type_names_buffer = new char[ buff_size ];
	char *type_names_stored = new char[ buff_size ](); // Make sure to 0!

	if( ! (type_names_buffer && type_names_stored) ){
		my_runtime_error( __FILE__, __LINE__,
		                  "Failed to allocate storage for type names" );
	}

	if( !store_typenames.empty() ){
		uint longest_word_length = 0;
		// Makes sure you skip the first, unused type, because that
		// one is also not in the buffer in the file.
		for( std::size_t i = 1; i < store_typenames.size(); ++i ){
			const std::string &s = store_typenames[i];
			std::size_t l = s.length();
			if( l > longest_word_length )
				longest_word_length = l;
		}

		// Write the words in the right formatting:
		int stride = longest_word_length + 1; // +1 for trailing 0.
		char *start = type_names_stored;
		for( std::size_t i = 1; i < store_typenames.size(); ++i ){
			std::strcpy( start, store_typenames[i].c_str() );
			start += stride;
		}
	}

	int status = get_chunk_data( "particles/types", type_names_buffer,
	                             buff_size, type_names_stored );
	if( status != 0 && status != 2 ) return status;

	// If you get here, status was good.
	uint current_name = 0;
	char *name = type_names_buffer;
	int size_left = buff_size;

	type_names[0] = "__UNUSED__";
	std::string next_name;
	bool at_end_of_word = false;

	while( current_name < n_types ){
		utf16_char next_char;
		next_char.d = 0;
		int chars = util::next_utf16_char( name, next_char, size_left );
		my_assert( __FILE__, __LINE__, chars > 0,
		           "Failed to read next utf16 char" );
		if( chars == 1 ){
			// Check for '\0':
			if( next_char.c[0] == '\0' ){
				if( at_end_of_word ){
					// Keep scrolling until you hit
					// first non-0 char.
					size_left -= chars;
					name += chars;
					continue;
				}

				type_names[current_name+1] = next_name;
				next_name = "";
				++current_name;
				at_end_of_word = true;
			}else{
				at_end_of_word = false;
				next_name += next_char.c[0];
			}
		} else if( chars == 2 ){
			at_end_of_word = false;
			utf8_char utf8;
			util::down_cast_utf16_char( next_char, utf8 );
			next_name += utf8.c[0];
			next_name += utf8.c[1];
		} else {
			my_warning( __FILE__, __LINE__,
			            "Got > UTF-8 character, dropping it!" );
		}


		size_left -= chars;
		name += chars;
	}

	store_typenames = type_names;

	delete [] type_names_buffer;
	delete [] type_names_stored;
	return 0;
}


template <int data_type, typename T>
int dump_reader_hoomd_gsd::add_optional_data( block_data &b, T &data, uint n_fields,
                                              const std::vector<std::string> &names )
{
	// If the data could not be found this block, we defaulted
	// to store_data. If the data was _not_ in the first run encountered,
	// store_data is empty, and hence, so will the data. Therefore, we can
	// check data.empty() to make sure this optional data is not present at all.
	if( data.empty() ){
		data.resize( n_fields * b.N );
	}

	std::size_t stride = names.size();

	if( n_fields != stride ){
		std::cerr << "Stride " << stride << " != # of fields "
		          << n_fields << "\n";
		my_runtime_error( __FILE__, __LINE__,
		                  "Stride != Number of fields" );
	}

	my_assert( __FILE__, __LINE__, stride * b.N == data.size(),
	           "Data sizes mismatch!" );

	if( data_type == data_field::DOUBLE ){
		// You need this constructor and copy per-element because
		// data might need to be reinterpreted. Also, the data
		// might contain more than one per-atom data, and those
		// need to be split to different data fields.
		for( std::size_t nf = 0; nf < n_fields; ++nf ){
			data_field_double d( names[nf], b.N );
			for( bigint i = 0; i < b.N; ++i ){
				std::size_t index = stride*i + nf;
				d[i] = data[index];
			}

			b.add_field( d );
		}
	}else if( data_type == data_field::INT ){
		for( std::size_t nf = 0; nf < n_fields; ++nf ){
			data_field_int d( names[nf], b.N );
			for( bigint i = 0; i < b.N; ++i ){
				std::size_t index = stride*i + nf;
				d[i] = data[index];
			}
			b.add_field( d );
		}
	}else{
		my_logic_error( __FILE__, __LINE__,
		                "Unknown data type encountered!" );
		return -2;
	}

	return 0;
}



int dump_reader_hoomd_gsd::get_next_block( block_data &block )
{
	++current_frame;
	if( current_frame == max_frame ){
		// All frames exhausted, this is as good as EOF:
		eof_ = true;
		return 1;
	}
	optional_data_found = 0;

	// const gsd_index_entry *entry;

	// Read out all the chunks you want to add to your dumpfile.
	block_data tmp;

	uint64_t tstep;
	uint8_t  dims;
	float    box[6];
	uint32_t N;

	bool default_type = false;


	int status = get_chunk_data( "configuration/step", &tstep,
	                             1, &store_tstep );
	if( status != 0 && status != 2 ){
		if( status == 1 ){
			std::cerr << "Encountered EOF\n";
		}
		return status;
	}


	status = get_chunk_data("configuration/dimensions", &dims,
	                        1, &store_dims );
	if( status != 0 && status != 2 ) return status;

	status = get_chunk_data( "configuration/box", box, 6, store_box );
	if( status != 0 && status != 2 ) return status;

	status = get_chunk_data( "particles/N", &N, 1, &store_N );
	if( status != 0 && status != 2 ) return status;


	tmp.tstep = tstep;
	tmp.dom.xlo[0] = -0.5*box[0];
	tmp.dom.xlo[1] = -0.5*box[1];
	tmp.dom.xlo[2] = -0.5*box[2];
	tmp.dom.xhi[0] =  0.5*box[0];
	tmp.dom.xhi[1] =  0.5*box[1];
	tmp.dom.xhi[2] =  0.5*box[2];

	tmp.set_natoms(N);

	// Read the data chunks:
	std::vector<float> xx(3*N);
	std::vector<uint32_t> type_ids(N);
	std::vector<int32_t> body(N);

	status = get_chunk_data( "particles/position", xx, store_x );
	if( status != 0 && status != 2 ) return status;

	status = get_chunk_data( "particles/body", body, store_body );
	if( status != 0 && status != 2 ) return status;

	status = get_chunk_data( "particles/typeid", type_ids, store_type_ids );
	if( status != 0 && status != 2 ) return status;


	// Set the types:

	if( default_type ){
		// This means each particle has HOOMD type 0,
		// which is our internal type 1.
		for( std::size_t i = 0; i < N; ++i ){
			type_ids[i] = 0;
		}
	} // Else we read everything correctly already.

	// Determine the number of types (+1 because Hoomd is 0-indexed).
	unsigned int n_types = *std::max_element( type_ids.begin(),
	                                          type_ids.end() ) + 1;

	// Extract the type names:
	std::vector<std::string> type_names( n_types + 1 );
	status = get_type_names( type_names );
	if( status ){
		good_ = false;
		return -1;
	}

	// Copy everything to the block_data format:

	optional_data_found[BODY] = !body.empty();

	// Convert everything here to data_fields:
	using dfd = data_field_double;
	using dfi = data_field_int;

	dfd x( "x", N );
	dfd y( "y", N );
	dfd z( "z", N );

	dfi id( "id", N );
	dfi mol( "mol", N );
	dfi type( "type", N );

	for( std::size_t i = 0; i < N; ++i ){
		x[i] = xx[3*i];
		y[i] = xx[3*i+1];
		z[i] = xx[3*i+2];
		// Offset with one because HOOMD uses 0-indexing:
		type[i] = type_ids[i] + 1;
		id[i] = i+1;
		if( optional_data_found[BODY] ){
			// Same for mol type:
			mol[i] = body[i] + 1;
		}
	}

	// Set the type names:
	tmp.ati.type_names = type_names;
	tmp.atom_style = lammps_tools::ATOM_STYLE_ATOMIC;

	tmp.add_field( id, block_data::special_fields::ID );
	if( optional_data_found[BODY] ){
		tmp.atom_style = lammps_tools::ATOM_STYLE_MOLECULAR;
		tmp.add_field( mol, block_data::special_fields::MOL );
	}
	tmp.add_field( type, block_data::special_fields::TYPE );
	tmp.add_field( x, block_data::special_fields::X );
	tmp.add_field( y, block_data::special_fields::Y );
	tmp.add_field( z, block_data::special_fields::Z );


	// ****    Check for additional fields that might be present:    ****
	std::vector<float> v( 3*N );
	status = get_chunk_data( "particles/velocity", v, store_v );
	if( status != 0 && status != 2 ) return status;
	if( status == 0 ) optional_data_found[VELOCITY] = 1;

	std::vector<float> o( 4*N );
	status = get_chunk_data( "particles/orientation", o, store_orient );
	if( status != 0 && status != 2 ) return status;
	if( status == 0 ) optional_data_found[ORIENTATION] = 1;

	std::vector<float> moment_inertia(3*N);
	status = get_chunk_data( "particles/moment_inertia", moment_inertia,
	                         store_moment_inertia );
	if( status != 0 && status != 2 ) return status;
	if( status == 0 ) optional_data_found[MOMENT_INERTIA] = 1;

	constexpr const int double_type = data_field::DOUBLE;

	if( optional_data_found[VELOCITY] ){
		add_optional_data<double_type>( tmp, v, 3, { "v.x", "v.y", "vz" } );
	}

	if( optional_data_found[ORIENTATION] ){
		add_optional_data<double_type>( tmp, o, 4,
		                                { "orientation.x", "orientation.y",
	                                          "orientation.z", "orientation.w" } );
	}
	if( optional_data_found[MOMENT_INERTIA] ){
		add_optional_data<double_type>( tmp, moment_inertia, 3,
		                                { "mom_inertia.x", "mom_inertia.y",
		                                  "mom_inertia.z" } );
	}

	block = tmp;

	return 0;

}


int64_t dump_reader_hoomd_gsd::file_size() const
{
	return gh->file_size;
}



#endif // HAVE_GSD



// ***********    Declare all instantiations for templated members:    *********
#define MAKE_FORWARD_TEMPLATE_DECL_PTR( TYPE ) \
	template \
	int dump_reader_hoomd_gsd::get_chunk_data<TYPE> \
	( const char *, TYPE *, uint, TYPE * );

#define MAKE_FORWARD_TEMPLATE_DECL_REF( TYPE ) \
	template  \
	int dump_reader_hoomd_gsd::get_chunk_data<TYPE> \
	( const char *, TYPE &, TYPE & );


MAKE_FORWARD_TEMPLATE_DECL_PTR(uint64_t)
MAKE_FORWARD_TEMPLATE_DECL_PTR(uint32_t)
MAKE_FORWARD_TEMPLATE_DECL_PTR(uint8_t)
MAKE_FORWARD_TEMPLATE_DECL_PTR(float)

MAKE_FORWARD_TEMPLATE_DECL_REF(std::vector<float>)
MAKE_FORWARD_TEMPLATE_DECL_REF(std::vector<uint8_t>)
MAKE_FORWARD_TEMPLATE_DECL_REF(std::vector<uint32_t>)



#undef MAKE_FORWARD_TEMPLATE_DECL_REF
#undef MAKE_FORWARD_TEMPLATE_DECL_PTR


// ***********    End of declarations of templated members.     ****************



} // namespace readers

} // lammps_tools
