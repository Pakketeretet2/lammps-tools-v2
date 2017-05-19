#include "util.hpp"
#include "writers_lammps.hpp"

#include <fstream>

#ifdef HAVE_BOOST_GZIP
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#endif // HAVE_BOOST_GZIP

#include <cstdio>
#include <cstring>

#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif


#if defined(LAMMPS_SMALLBIG)
typedef int tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#elif defined(LAMMPS_SMALLSMALL)
typedef int tagint;
typedef int bigint;
#define BIGINT_FORMAT "%d"
#else /* LAMMPS_BIGBIG */
typedef int64_t tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#endif


namespace lammps_tools {

namespace writers {


void block_to_lammps_data( const std::string &fname, const block_data &b )
{
	std::ofstream out( fname );
	block_to_lammps_data( out, b );
}


void block_to_lammps_data( std::ostream &out, const block_data &b )
{
	std::vector<int> id, type, mol, ix, iy, iz;
	std::vector<double> x, y, z;

	// Set the important fields.
	std::vector<std::string> names(6);
	names[0] = b.get_special_field_name(block_data::ID);
	names[1] = b.get_special_field_name(block_data::TYPE);
	names[2] = b.get_special_field_name(block_data::X);
	names[3] = b.get_special_field_name(block_data::Y);
	names[4] = b.get_special_field_name(block_data::Z);
	names[5] = b.get_special_field_name(block_data::MOL);

	std::vector<std::string> image_flag_names(3);
	image_flag_names[0] = b.get_special_field_name(block_data::IX);
	image_flag_names[1] = b.get_special_field_name(block_data::IY);
	image_flag_names[2] = b.get_special_field_name(block_data::IZ);

	bool has_image_flags = !( image_flag_names[0].empty() ||
	                          image_flag_names[1].empty() ||
	                          image_flag_names[2].empty() );

	bool success = grab_field_as<int>( b, names[0],id );
	success &= grab_field_as<int>   ( b, names[1], type );
	success &= grab_field_as<double>( b, names[2], x );
	success &= grab_field_as<double>( b, names[3], y );
	success &= grab_field_as<double>( b, names[4], z );

	if( b.atom_style == block_data::MOLECULAR ){
		success &= grab_field_as<int>( b, names[5], mol );
	}

	if( has_image_flags ){
		success &= grab_field_as<int>( b, image_flag_names[0], ix );
		success &= grab_field_as<int>( b, image_flag_names[1], iy );
		success &= grab_field_as<int>( b, image_flag_names[2], iz );
	}

	my_assert( __FILE__, __LINE__, success,
	           "Failure to grab essential fields for write to data!" );


	out << "LAMMPS data file via lammps-tools\n\n";
	out << b.N << " atoms\n";
	int n_types = b.N_types;
	out << n_types << " atom types\n\n";
	const char *words[3][2] = { { "xlo", "xhi"}, {"ylo", "yhi"}, {"zlo", "zhi"} };
	for( int dim = 0; dim < 3; ++dim ){
		out << b.dom.xlo[dim] << " " << b.dom.xhi[dim]
		  << " " << words[dim][0] << " " << words[dim][1] << "\n";
	}
	out << "\n";

	std::string atom_style = "atomic";
	if( b.atom_style == block_data::MOLECULAR ) atom_style = "molecular";
	out << "Atoms # " << atom_style << "\n\n";


	for( int i = 0; i < b.N; ++i ){
		out << id[i];
		if( b.atom_style == block_data::MOLECULAR ){
			out << " " << mol[i];
		}
		out << " " << type[i] << " " << x[i]
		    << " " << y[i] << " " << z[i];

		if( has_image_flags ){
			out << " " << ix[i] << " " << iy[i] << " " << iz[i];
		}


		out << "\n";
	}
}


void block_to_lammps_dump( const std::string &fname,
                           const block_data &b, int fformat )
{
	switch(fformat){
		default:
			my_runtime_error(__FILE__, __LINE__,
			                 "Unknown file format!" );
		case lammps_tools::readers::PLAIN: {
			std::ofstream out( fname );
			block_to_lammps_dump_text( out, b );
			break;
		}
		case lammps_tools::readers::BIN:{
			std::ofstream out( fname, std::ios::binary );
			block_to_lammps_dump_bin( out, b );
			break;
		}
		case lammps_tools::readers::GZIP:{
#ifdef HAVE_BOOST_GZIP
			using namespace boost::iostreams;

			std::ofstream in(fname);
			filtering_stream<output> out;
			out.push(gzip_compressor());
			out.push(in);
			block_to_lammps_dump_text( out, b );
			break;
#else
			my_runtime_error( __FILE__, __LINE__,
			                  "Not compiled with boost support! "
			                  "Cannot write to GZIP!" );
#endif
		}
	}
}

void block_to_lammps_dump( std::ostream &out,
                          const block_data &b, int fformat )
{
	switch(fformat){
		default:
			my_runtime_error(__FILE__, __LINE__,
			                 "Unknown file format!" );
		case lammps_tools::readers::PLAIN:{
			block_to_lammps_dump_text( out, b );
			break;
		}
		case lammps_tools::readers::BIN:{
			block_to_lammps_dump_bin( out, b );
			break;
		}
		case lammps_tools::readers::GZIP:{
			block_to_lammps_dump_text( out, b );
			break;
		}
	}
}

std::vector<std::string> sort_names( const block_data &b )
{
	std::vector<std::string> headers( b.n_data_fields() );
	std::vector<std::string> data_names = b.get_data_names();
	int current = 0;
	headers[current++] = b.get_special_field_name( block_data::ID );
	if( b.atom_style == block_data::MOLECULAR ){
		headers[current++] = b.get_special_field_name( block_data::MOL );
	}
	headers[current++] = b.get_special_field_name( block_data::TYPE );
	headers[current++] = b.get_special_field_name( block_data::X );
	headers[current++] = b.get_special_field_name( block_data::Y );
	headers[current++] = b.get_special_field_name( block_data::Z );

	for( const std::string &name : data_names ){
		if( util::contains( headers, name ) ) continue;

		headers[current++] = name;
	}

	std::cerr << "Names in order are:";
	for( const std::string &name : headers ){
		std::cerr << " " << name;
	}
	std::cerr << "\n";

	return headers;
}

void output_atom_line( const block_data &b, std::ostream &out,
                       const std::vector<std::string> &data_order, int i )
{
	using dfd = data_field_double;
	using dfi = data_field_int;

	for( const std::string &name : data_order ){
		const data_field *df = b.get_data( name );
		int type = df->type();
		if( type == data_field::INT ){
			const dfi *di = static_cast<const dfi*>(df);
			out << (*di)[i] << " ";
		}else if( type == data_field::DOUBLE ){
			const dfd *dd = static_cast<const dfd*>(df);
			out << (*dd)[i] << " ";
		}
	}
	out << "\n";
}

void block_to_lammps_dump_text( std::ostream &out, const block_data &b )
{
	std::string boxline = "ITEM: BOX BOUNDS";
	if( b.dom.periodic & domain::BIT_X ) boxline += " pp";
	else                                 boxline += " ff";
	if( b.dom.periodic & domain::BIT_Y ) boxline += " pp";
	else                                 boxline += " ff";
	if( b.dom.periodic & domain::BIT_Z ) boxline += " pp";
	else                                 boxline += " ff";


	out<< "ITEM: TIMESTEP\n" << b.tstep << "\nITEM: NUMBER OF ATOMS\n";
	out << b.N << "\n" << boxline << "\n";
	out << b.dom.xlo[0] << " " << b.dom.xhi[0] << "\n";
	out << b.dom.xlo[1] << " " << b.dom.xhi[1] << "\n";
	out << b.dom.xlo[2] << " " << b.dom.xhi[2] << "\n";

	//std::vector<std::string> data_order = sort_names( b );
	std::vector<std::string> data_order = b.get_data_names();
	std::string header_line = "ITEM: ATOMS";
	for( const std::string &header : data_order ){
		header_line += " " + header;
	}

	out << header_line << "\n";

	for( int i = 0; i < b.N; ++i ){
		output_atom_line( b, out, data_order, i );
	}
}



template <typename T>
void write_bin( std::ostream &out, const T &val, int size )
{
	const void *address = &val;
	const char *buf = static_cast<const char *>( address );
	out.write( buf, size );
}

template <typename T>
void write_bin( std::ostream &out, const T &val )
{
	write_bin( out, val, sizeof(T) );
}



// void block_to_lammps_dump_bin( std::ostream &out, const block_data &b )
void block_to_lammps_dump_bin( std::ostream &out, const block_data &b )
{
	// Some defaults we don't use yet:
	int triclinic = 0;
	int xy = 0;
	int xz = 0;
	int yz = 0;
	int boundary[3][2];
	int size_one = b.n_data_fields();

	if( b.dom.periodic & domain::BIT_X ){
		boundary[0][0] = boundary[0][1] = 0;
	}else{
		boundary[0][0] = boundary[0][1] = 1;
	}
	if( b.dom.periodic & domain::BIT_Y ){
		boundary[1][0] = boundary[1][1] = 0;
	}else{
		boundary[1][0] = boundary[1][1] = 1;
	}
	if( b.dom.periodic & domain::BIT_Z ){
		boundary[2][0] = boundary[2][1] = 0;
	}else{
		boundary[2][0] = boundary[2][1] = 1;
	}
	write_bin( out, b.tstep );
	write_bin( out, b.N );
	write_bin( out, triclinic );
	write_bin( out, boundary[0][0], 6*sizeof(int) );
	write_bin( out, b.dom.xlo[0] );
	write_bin( out, b.dom.xhi[0] );
	write_bin( out, b.dom.xlo[1] );
	write_bin( out, b.dom.xhi[1] );
	write_bin( out, b.dom.xlo[2] );
	write_bin( out, b.dom.xhi[2] );

	if( triclinic ){
		write_bin( out, xy );
		write_bin( out, xz );
		write_bin( out, yz );
	}
	write_bin( out, size_one );

	// In multi-core writes or if the file is split in chunks this changes:
	int data_size = b.N * b.n_data_fields();
	int nchunk = 1;
	write_bin( out, nchunk );
	write_bin( out, data_size );

	std::vector<std::string> data_names = b.get_data_names();
	std::vector<int> data_types( size_one );
	for( int j = 0; j < size_one; ++j ){
		const data_field *df = b.get_data( data_names[j] );
		data_types[j] = df->type();
	}

	using dfd = lammps_tools::data_field_double;
	using dfi = lammps_tools::data_field_int;

	for( int j = 0; j < b.N; ++j ){
		for( int k = 0; k < size_one; ++k ){
			const data_field *df = b.get_data( data_names[k] );
			if( data_types[k] == data_field::INT ){
				const dfi* di = static_cast<const dfi*>( df );
			        double value = (*di)[j]; // Implicitly cast.
			        write_bin( out, value );
			}else if( data_types[k] == data_field::DOUBLE ){
				const dfd* dd = static_cast<const dfd*>( df );
				write_bin( out, (*dd)[j] );
			}
		}
	}
}


} // namespace writers

} // lammps_tools
