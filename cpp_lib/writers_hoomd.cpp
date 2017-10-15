#include "writers_hoomd.hpp"

#ifdef HAVE_GSD
#include <gsd.h>
#endif // HAVE_GSD

#include <iostream>

#include "block_data.hpp"
#include "block_data_access.hpp"
#include "enums.hpp"
#include "types.hpp"
#include "my_assert.hpp"

namespace lammps_tools {

namespace writers {


int block_to_hoomd_gsd( const std::string &fname, const block_data &b,
                         const std::string &write_mode )
{
#ifndef HAVE_GSD
	my_runtime_error( __FILE__, __LINE__,
	                  "Lammps-tools not compiled with GSD support! "
	                  "Cannot write to GSD format!" );
	return -1;

#else

	gsd_handle gh;
	uint32_t schema_version = gsd_make_version( 1, 1 );

	gsd_open_flag flag = GSD_OPEN_READWRITE;
	bool got_append = false;
	bool got_write = false;
	bool got_binary = false;
	for( char c : write_mode ){
		switch(c){
			case 'a':
				got_append = true;
				break;
			case 'w':
				got_write = true;
				break;
			case 'b':
				got_binary = true;
				break;
			case '+':
				std::cerr << "Ignoring '+' in write mode.\n";
				break;
			default:
				std::cerr << "Unrecognized char '" << c
				          << "' in write mode!\n";
				break;
		}
	}

	if( got_append && got_write ){
		std::cerr << "Ambigous write mode passed ('"
		          << write_mode << "')!\n";
		return -2;
	}
	if( !(got_append || got_write) ){
		std::cerr << "No write mode passed ('"
		          << write_mode << "')!\n";
		return -4;
	}
	if( !got_binary ){
		std::cerr << "GSD can only be written to binary but got '"
		          << write_mode << "'!\n";
		return -3;
	}

	if( got_append ){
		flag = GSD_OPEN_APPEND;
	}else{
		flag = GSD_OPEN_READWRITE;
	}

	// Check if file already exists:
	int status;
	if( !util::file_exists( fname ) || got_write ){
		status = gsd_create_and_open( &gh, fname.c_str(),
		                              "lammps-tools", "hoomd",
		                              schema_version, flag, 0 );
	}else{
		status = gsd_open( &gh, fname.c_str(), flag );
	}

	if( status ){
		std::cerr << "An error occured creating GSD file!\n";
		return status;
	}else{
		status = block_to_hoomd_gsd( &gh, b );
	}
	gsd_close( &gh );
	if( status ){
		std::cerr << "An error occured writing to GSD file!\n";
		return status;
	}

	return 0;

#endif // HAVE_GSD

}


/**
   \brief reconstructs GSD-style buffer from vectors of data_fields.

   Assumes dest is allocated and the right size.
*/
template <int data_field_type, typename T_to>
int reconstruct_fields_as( T_to *dest, int n_arr, const block_data &b,
                           std::vector<std::string> field_names )
{
	using dfi = data_field_int;
	using dfd = data_field_double;

	auto is_null = []( const void *ptr ){ return ptr == nullptr; };



	std::size_t M = field_names.size();
	std::size_t N = b.N;
	my_assert( __FILE__, __LINE__, n_arr == M,
	           "Array count mismatch!" );



	if( data_field_type == data_field::DOUBLE ){
		std::vector<const dfd *> arr_from(n_arr, nullptr);
		for( std::size_t i = 0; i < b.n_data_fields(); ++i ){
			for( int name_idx = 0; name_idx < n_arr; ++name_idx ){
				if( b[i].name == field_names[name_idx] ){
					const dfd *ptr = static_cast<const dfd *>( &b[i] );
					arr_from[name_idx] = ptr;
				}
			}
		}
		if( std::any_of( arr_from.begin(), arr_from.end(), is_null ) ){
			// Abort!
			return -1;
		}

		for( std::size_t i = 0; i < N; ++i ){
			for( int n = 0; n < n_arr; ++n ){
				int idx = n_arr * i + n;
				const dfd *ptr = arr_from[n];
				my_assert( __FILE__, __LINE__, ptr,
				           "wtf? ptr is null!" );
				const dfd &data = *ptr;
				dest[idx] = data[i];
			}
		}
	}else if( data_field_type == data_field::INT ){
		std::vector<const dfi *> arr_from(n_arr, nullptr);
		for( std::size_t i = 0; i < b.n_data_fields(); ++i ){
			for( int name_idx = 0; name_idx < n_arr; ++name_idx ){
				if( b[i].name == field_names[name_idx] ){
					const dfi *ptr = static_cast<const dfi *>( &b[i] );
					arr_from[name_idx] = ptr;
				}
			}
		}
		bool everything_ok = std::any_of( arr_from.begin(),
		                                  arr_from.end(), is_null );
		my_assert( __FILE__, __LINE__, everything_ok,
		           "Failed to properly map arrays!" );

		for( std::size_t i = 0; i < N; ++i ){
			for( int n = 0; n < n_arr; ++n ){
				int idx = n_arr * i + n;
				const dfi &data = (*arr_from[n]);
				dest[idx] = data[i];
			}
		}
	}

	return 0;
}

int block_to_hoomd_gsd( gsd_handle *gh, const block_data &b )
{
#ifndef HAVE_GSD
	my_runtime_error( __FILE__, __LINE__,
	                  "Lammps-tools not compiled with GSD support! "
	                  "Cannot write to GSD format!" );
	return -1;
#else

	int status;
	uint64_t step = b.tstep;
	uint8_t  dims = 3;
	uint32_t N    = b.N;
	double L[3];
	float box[6];
	L[0] = b.dom.xhi[0] - b.dom.xlo[0];
	L[1] = b.dom.xhi[1] - b.dom.xlo[1];
	L[2] = b.dom.xhi[2] - b.dom.xlo[2];
	box[3] = box[4] = box[5] = 0.0;

	box[0] = L[0];
	box[1] = L[1];
	box[2] = L[2];


	status = gsd_write_chunk( gh, "configuration/step", GSD_TYPE_UINT64,
	                          1, 1, 0, &step );
	my_assert( __FILE__, __LINE__, status == 0,
	           "Failed to write time step" );
	status = gsd_write_chunk( gh, "configuration/dimensions", GSD_TYPE_UINT8,
	                          1, 1, 0, &dims );
	my_assert( __FILE__, __LINE__, status == 0,
	           "Failed to write box dimensions" );
	status = gsd_write_chunk( gh, "configuration/box", GSD_TYPE_FLOAT,
	                          6, 1, 0, box );
	my_assert( __FILE__, __LINE__, status == 0,
	           "Failed to write domain info" );
	status = gsd_write_chunk( gh, "particles/N", GSD_TYPE_UINT32,
	                          1, 1, 0, &N );
	my_assert( __FILE__, __LINE__, status == 0,
	           "Failed to write particle number" );

	float    *xx    = new float[3*N];
	uint32_t *types	= new uint32_t[N];
	int n_types = 0;

	const std::vector<int> &id = get_id(b);
	const std::vector<int> &type = get_type(b);

	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);


	for( std::size_t i = 0; i < N; ++i ){
		// Sort them along id:
		int j = id[i] - 1;

		// Remap the positions to -0.5L and 0.5L.
		double xi[3];
		for( int d = 0; d < 3; ++d ){
			double xd = 0;
			switch(d){
				case 0:
					xd = x[i];
					break;
				case 1:
					xd = y[i];
					break;
				case 2:
					xd = z[i];
					break;
			}
			xi[d] = xd - b.dom.xlo[d];
			xi[d] -= 0.5*L[d];

			// Check box bounds:
			if( xi[d] > 0.5*L[d] || xi[d] < -0.5*L[d] ){
				std::cerr << "Particle " << id[i] << " is "
				          << "out of box bound in dim " << d
				          << " ( " << -0.5*L[d] << ", "
				          << 0.5*L[d] << " ) with x = "
				          << xi[d] << ".\n";
				my_runtime_error( __FILE__, __LINE__,
				                  "Invalid particle position" );
			}

			xx[3*j+d] = xi[d];
		}
		int current_type = type[i];
		types[j] = current_type - 1;
		if( current_type > n_types ) n_types = current_type;

	}
	status = gsd_write_chunk( gh, "particles/position", GSD_TYPE_FLOAT,
	                          N, 3, 0, xx );
	my_assert( __FILE__, __LINE__, status == 0,
	           "Failed to write particle positions" );
	status = gsd_write_chunk( gh, "particles/typeid", GSD_TYPE_UINT32,
	                          N, 1, 0, types );
	my_assert( __FILE__, __LINE__, status == 0,
	           "Failed to write particle ids" );


	// Write the actual type names.
	const uint buff_size = gsd::TYPE_BUFFER_SIZE;
	char *type_names = new char[buff_size * n_types];

	for( int t = 0; t < n_types; ++t ){
		// Typenames are now stored:
		char *current_name = type_names + t*buff_size;
		// Remember the +1, lammps-tools doesn't use 0-indexed one.
		std::string tname = b.ati.type_names[t+1];

		std::size_t idx = 0;
		for( idx = 0; idx < tname.length(); ++idx ){
			current_name[idx] = tname[idx];
		}
		current_name[idx] = '\0';
	}
	status = gsd_write_chunk( gh, "particles/types", GSD_TYPE_UINT8,
	                          n_types, buff_size, 0, type_names );
	my_assert( __FILE__, __LINE__, status == 0,
	           "Failed to write particle types" );

	float *orientation = new float[4*b.N];

	constexpr const int double_type = data_field::DOUBLE;
	// constexpr const int int_type = data_field::INT;

	reconstruct_fields_as<double_type, float>( orientation, 4, b,
	                                          {"orientation.x", "orientation.y",
	                                           "orientation.z", "orientation.w"} );

	status = gsd_write_chunk( gh, "particles/orientation",
	                          GSD_TYPE_FLOAT, b.N, 4, 0, orientation );
	delete [] orientation;



	status = gsd_end_frame( gh );

	delete [] xx;
	delete [] types;
	delete [] type_names;

	return status;

#endif // HAVE_GSD

}

} // namespace writers

} // namespace lammps_tools
