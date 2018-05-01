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
                        const std::string &write_mode, uint props )
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
		status = block_to_hoomd_gsd( &gh, b, props );
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
int reconstruct_fields_as_gsd( T_to *dest, int n_arr, const block_data &b,
                               const std::vector<std::string> &field_names )
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
		bool everything_ok = !std::any_of( arr_from.begin(),
		                                   arr_from.end(), is_null );
		std::cerr << "names (ptr) are:";
		for( int i = 0; i < n_arr; ++i ){
			const std::string &n = field_names[i];
			std::cerr << " " << n << "(" << arr_from[i] << ")";
		}
		std::cerr << "\n";

		my_assert( __FILE__, __LINE__, everything_ok,
		           "Failed to properly map arrays!" );

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
		bool everything_ok = !std::any_of( arr_from.begin(),
		                                   arr_from.end(), is_null );
		std::cerr << "names (ptr) are:";
		for( int i = 0; i < n_arr; ++i ){
			const std::string &n = field_names[i];
			std::cerr << " " << n << "(" << arr_from[i] << ")";
		}
		std::cerr << "\n";
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






int block_to_hoomd_gsd( gsd_handle *gh, const block_data &b, uint props )
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

	if( props & TIME_STEP ){
		status = gsd_write_chunk( gh, "configuration/step",
		                          GSD_TYPE_UINT64, 1, 1, 0, &step );
		my_assert( __FILE__, __LINE__, status == 0,
		           "Failed to write time step" );
	}

	if( props & DIMENSIONS ){
		status = gsd_write_chunk( gh, "configuration/dimensions",
		                          GSD_TYPE_UINT8, 1, 1, 0, &dims );
		my_assert( __FILE__, __LINE__, status == 0,
		           "Failed to write box dimensions" );
	}

	if( props & BOX ){
		status = gsd_write_chunk( gh, "configuration/box",
		                          GSD_TYPE_FLOAT, 6, 1, 0, box );
		my_assert( __FILE__, __LINE__, status == 0,
		           "Failed to write domain info" );
	}

	if( props & PARTICLE_NUMBER ){
		status = gsd_write_chunk( gh, "particles/N",
		                          GSD_TYPE_UINT32, 1, 1, 0, &N );
		my_assert( __FILE__, __LINE__, status == 0,
		           "Failed to write particle number" );
	}

	const std::vector<int> &id   = get_id(b);
	const std::vector<int> &type = get_type(b);

	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);

	float *xx = nullptr;
	uint32_t *types = nullptr;

	if( props & POSITIONS ){
		xx = new float[3*N];
	}

	types = new uint32_t[N];
	int n_types = 0;

	for( std::size_t i = 0; i < N; ++i ){
		// Sort them along id:
		int j = id[i] - 1;

		// NOTE: n_types might be needed for types instead of typeid,
		// so always generate it
		int current_type = type[i];
		types[j] = current_type - 1;
		if( current_type > n_types ) n_types = current_type;

		if( !(props & POSITIONS) ){
			continue;
		}

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
				std::cerr << "Particle " << id[i]
				          << " is out of box bound in "
				          << "dim " << d << " ( "
				          << -0.5*L[d] << ", "
				          << 0.5*L[d] << " ) with x = "
				          << xi[d] << ".\n";
				const char *msg = "Invalid particle position";
				my_runtime_error( __FILE__, __LINE__,
				                  msg );
			}
			xx[3*j+d] = xi[d];
		}
	}


	if( props & TYPES ){
		// Write the actual type names.
		// const uint buff_size = n_types*gsd::TYPE_BUFFER_SIZE;
		std::size_t longest_name = 0;
		std::cerr << "Writing types...\n";
		for( int t = 0; t < n_types; ++t ){
			std::string name = b.ati.type_names[t+1];
			longest_name = std::max( longest_name, name.length() );
		}
		std::cerr << "Longest name is " << longest_name << " long.\n";
		std::size_t stride = longest_name + 1;

		const uint buff_size = n_types * stride;
		char *type_names = new char[buff_size]();

		std::cerr << "The type names are:";

		for( int t = 0; t < n_types; ++t ){
			// Typenames are now stored:
			char *current_name = type_names + t*stride;
			// Remember the +1, lammps-tools doesn't use 0-indexed one.
			std::string tname = b.ati.type_names[t+1];
			std::cerr << " " << tname;
			std::size_t idx = 0;
			for( idx = 0; idx < tname.length(); ++idx ){
				current_name[idx] = tname[idx];
			}
			current_name[idx] = '\0';
		}

		std::cerr << "\nWriting type names to a buffer of size "
		          << n_types << " times " << stride << "\n";
		/*
		std::cerr << "That buffer is:\n";
		for( int i = 0; i < buff_size; ++i ){
			std::cerr << type_names[i] << "\n";
		}
		*/
		status = gsd_write_chunk( gh, "particles/types", GSD_TYPE_INT8,
		                          n_types, stride, 0, type_names );

		my_assert( __FILE__, __LINE__, status == 0,
		           "Failed to write particle types" );
		delete [] type_names;
	}


	if( props & TYPEID ){
		std::cerr << "Writing types to a buffer of size " << N << "\n";
		status = gsd_write_chunk( gh, "particles/typeid",
		                          GSD_TYPE_UINT32, N, 1, 0, types );
		my_assert( __FILE__, __LINE__, status == 0,
		           "Failed to write particle ids" );


		delete [] types;
	}



	if( props & BODY ){
		using namespace lammps_tools;
		const data_field_int *mol_field;
		const data_field *ptr;
		ptr = b.get_special_field( block_data::MOL );
		if( !ptr ){
			std::cerr << "Cannot write MOL because block does not "
			          << "have MOL field!\n";
		}else{
			mol_field = static_cast<const data_field_int *>(ptr);
			int32_t *body = new int32_t[b.N];
			constexpr const int int_type = data_field::INT;
			for( std::size_t i = 0; i < b.N; ++i ){
				int j = id[i] - 1;
				body[j] = (*mol_field)[i];
			}

			status = gsd_write_chunk( gh, "particles/body",
			                          GSD_TYPE_INT32, b.N, 1, 0,
			                          body );
			delete [] body;

		}
	}

	if( props & MOM_INERTIA ){

		std::cerr << "Writing moment of inertia!\n";

		float *mom_inertia = new float[3*b.N];
		constexpr const int double_type = data_field::DOUBLE;

		reconstruct_fields_as_gsd<double_type, float>( mom_inertia, 3, b,
		                                               {"mom_inertia.x",
		                                                "mom_inertia.y",
		                                                "mom_inertia.z"} );

		status = gsd_write_chunk( gh, "particles/mom_inertia",
		                          GSD_TYPE_FLOAT, b.N, 3, 0, mom_inertia );
		delete [] mom_inertia;


	}

	if( props & POSITIONS ){
		status = gsd_write_chunk( gh, "particles/position",
		                          GSD_TYPE_FLOAT, N, 3, 0, xx );
		my_assert( __FILE__, __LINE__, status == 0,
		           "Failed to write particle positions" );

		delete [] xx;
	}



	if( props & ORIENTATION ) {
		float *orient = new float[4*b.N];
		constexpr const int double_type = data_field::DOUBLE;

		reconstruct_fields_as_gsd<double_type, float>( orient, 4, b,
		                                           {"orientation.x",
				                            "orientation.y",
				                            "orientation.z",
				                            "orientation.w"} );

		status = gsd_write_chunk( gh, "particles/orientation",
		                          GSD_TYPE_FLOAT, b.N, 4, 0, orient );
		delete [] orient;
	}




	status = gsd_end_frame( gh );


	return status;

#endif // HAVE_GSD
}

} // namespace writers


} // namespace lammps_tools
