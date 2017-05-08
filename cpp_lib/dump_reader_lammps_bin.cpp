#include "dump_reader_lammps_bin.hpp"

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

using namespace lammps_tools;
using namespace readers;

namespace lammps_tools {

namespace readers {

dump_reader_lammps_bin::dump_reader_lammps_bin( const std::string &fname )
	: in( nullptr )
{
	in = fopen( fname.c_str(), "rb" );
	my_assert( __FILE__, __LINE__, in, "Failed to open dump file!" );
}

dump_reader_lammps_bin::dump_reader_lammps_bin( const std::string &fname,
                                                std::vector<std::string> h )
	: in( nullptr )
{
	in = fopen( fname.c_str(), "rb" );
	my_assert( __FILE__, __LINE__, in, "Failed to open dump file!" );
}

dump_reader_lammps_bin::~dump_reader_lammps_bin()
{
	if( in ) fclose( in );
}



int dump_reader_lammps_bin::get_next_block( block_data &block )
{
	if( !in ) return -1;

	if( feof(in) ){
		return 1;
	}
	int size_one, nchunk;
	block_data tmp;
	int status = next_block_meta( tmp, size_one, nchunk );
	if( status > 0 ){
		if( status == 1 && !quiet ) std::cerr << "EOF reached.\n";
		return status;
	}else if( status < 0 ){
		std::cerr << "Failed to get meta!\n";
		return -1;
	}
	// status was 0, so meta was gotten properly.
	status = next_block_body( tmp, size_one, nchunk );
	if( status ){
		std::cerr << "Failed to get body!\n";
		return status;
	}
	block = std::move( tmp );
	return 0;
}


int dump_reader_lammps_bin::next_block_meta( block_data &block,
                                             int &size_one, int &nchunk )
{
	bigint ntimestep, natoms;
	double xlo[3], xhi[3];
	int boundary[3][2];
	char boundstr[9];
	int triclinic;
	int xy,xz,yz;

	while( true ){
		fread( &ntimestep, sizeof(bigint), 1, in );
		if( !in ) return -1;

		if( feof(in) ){
			return 1;
		}

		fread(&natoms,sizeof(bigint),1,in);
		fread(&triclinic,sizeof(int),1,in);
		fread(&boundary[0][0],6*sizeof(int),1,in);
		fread(xlo  ,sizeof(double),1,in);
		fread(xhi  ,sizeof(double),1,in);
		fread(xlo+1,sizeof(double),1,in);
		fread(xhi+1,sizeof(double),1,in);
		fread(xlo+2,sizeof(double),1,in);
		fread(xhi+2,sizeof(double),1,in);
		if (triclinic) {
			fread(&xy,sizeof(double),1,in);
			fread(&xz,sizeof(double),1,in);
			fread(&yz,sizeof(double),1,in);
		}
		fread(&size_one,sizeof(int),1,in);
		fread(&nchunk,sizeof(int),1,in);

		block.N = natoms;
		block.tstep = ntimestep;
		block.dom.xlo[0] = xlo[0];
		block.dom.xlo[1] = xlo[1];
		block.dom.xlo[2] = xlo[2];

		block.dom.xhi[0] = xhi[0];
		block.dom.xhi[1] = xhi[1];
		block.dom.xhi[2] = xhi[2];

		int m = 0;
		for (int idim = 0; idim < 3; idim++) {
			for (int iside = 0; iside < 2; iside++) {
				if (boundary[idim][iside] == 0){
					boundstr[m++] = 'p';
				}else if (boundary[idim][iside] == 1){
					boundstr[m++] = 'f';
				}else if (boundary[idim][iside] == 2){
					boundstr[m++] = 's';
				}else if (boundary[idim][iside] == 3){
					boundstr[m++] = 'm';
				}
			}
			boundstr[m++] = ' ';
		}
		boundstr[8] = '\0';

		block.dom.periodic = 0;
		if( (boundary[0][0] == 0) && (boundary[0][1] == 0) ){
			block.dom.periodic += domain::BIT_X;
		}
		if( (boundary[1][0] == 0) && (boundary[1][1] == 0) ){
			block.dom.periodic += domain::BIT_Y;
		}
		if( (boundary[2][0] == 0) && (boundary[2][1] == 0) ){
			block.dom.periodic += domain::BIT_Z;
		}

		return 0;
	}
}

int dump_reader_lammps_bin::next_block_body( block_data &block,
                                             int size_one, int nchunk )
{
	using dfi = data_field_int;
	using dfd = data_field_double;

	int maxbuf = 0;
	double *buf = nullptr;
	int n = 0;
	block.set_natoms( block.N );
	const std::vector<std::string> &headers = get_column_headers();
	std::vector<data_field*> data_fields(size_one);


	my_assert( __FILE__, __LINE__, !headers.empty(),
	           "Column headers required for binary LAMMPS dump files!" );


	my_assert( __FILE__, __LINE__, headers.size() == size_one,
	           "Column number does not match number of headers!" );

	for( int i = 0; i < nchunk; i++ ){
		fread(&n,sizeof(int),1,in);

		// extend buffer to fit chunk size

		if( n > maxbuf ){
			if (buf) delete [] buf;
			buf = new double[n];
			maxbuf = n;
		}

		// read chunk and write as size_one values per line

		fread(buf,sizeof(double),n,in);
		n /= size_one;


		int m = 0, j = 0;
		for( std::string h : headers ){
			// Depending on the keyword, you want to take
			// either an int or a double.
			if( is_int_data_field( h ) ){
				dfi *new_field = new dfi( h, block.N );
				data_fields[i] = new_field;
			}else{
				// Assume it's a double.
				dfd *new_field = new dfd( h, block.N );
				data_fields[i] = new_field;
			}
			++i;
		}

		for( j = 0; j < n; j++ ){
			for( int k = 0; k < size_one; k++ ){
				int type = data_fields[k]->type();
				data_field *df = data_fields[k];
				if( type == data_field::INT ){
					dfi *field = static_cast<dfi*>( df );
					(*field)[j] = buf[m++];
				}else if( type == data_field::DOUBLE ){
					dfd *field = static_cast<dfd*>( df );
					(*field)[j] = buf[m++];
				}else{
					// No clue.
				}
			}
		}
	}

	// Copy all data fields into the block data.
	add_custom_data_fields( data_fields, block );

	if( buf ) delete [] buf;

	return 0;
}

bool dump_reader_lammps_bin::check_eof()  const
{
	return std::feof(in);
}

bool dump_reader_lammps_bin::check_good() const
{
	return !std::ferror(in);
}

} // namespace readers

} // namespace lammps_tools
