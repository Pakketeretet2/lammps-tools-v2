#ifndef DUMP_READER_HOOMD_GSD_HPP
#define DUMP_READER_HOOMD_GSD_HPP

/**
   \file dump_reader_hoomd_gsd.hpp

   Functionality to read in HOOMD gsd files.
*/

#include "enums.hpp"
#include "dump_reader.hpp"

#include <bitset>
#include <iosfwd>
#include <memory>

/// Declare a gsd_handle, but don't include yet.
struct gsd_handle;

namespace lammps_tools {

class block_data;

/// Contains functions and classes that are related to reading dump files.
namespace readers {


/// A generic class for reading in dump files.
class dump_reader_hoomd_gsd : public dump_reader
{
public:
	/// Constructor for opening from file name:
	dump_reader_hoomd_gsd( const std::string &fname );

	/// Destructor:
	virtual ~dump_reader_hoomd_gsd();

private:
	virtual int  get_next_block( block_data &block );
	virtual bool check_eof()  const;
	virtual bool check_good() const;

	/// Reads in type names, and does a bunch of corner case checking.
	int get_type_names( std::vector<std::string> &type_names );

	/**
	   \brief extracts named chunk and stores it in dest.

	   On error, dest and store shall be unchanged.
	   The return values are:
	      0: Succesfully read data and stored fresh data in store.
	      1: Encountered EOF while reading.
	      2: Could not find chunk, defaulted to store instead.
	     -1: I/O failure
	     -2: Invalid input
	     -3: Invalid data file
	     -4: A generic error.

	   \param name  The name of the chunk to read
	   \param dest  The destination to write to
	   \param size  The size of the destination
	   \param store If the chunk could not be found, default dest to this

	   \returns a status code.
	*/
	template <typename T>
	int get_chunk_data( const char *name, T *dest, uint size, T *store );

	/**
	   \brief overloads get_chunk_data for ref types.

	   This variant is more useful for containers.
	*/
	template <typename T>
	int get_chunk_data( const char *name, T &dest, T &store );

	/**
	   \brief adds additional data fields that are not mandatory to block

	   Returns:
	      1 if the data field does not exist anywhere in dump
	      0 on success
	     -1 on some generic failure.

	     \param b          The block to add data to
	     \param data       The data to add
	     \param store_data
	     \param names

	   \returns a status code
	*/
	template <int data_type, typename T>
	int add_optional_data( block_data &b, T &data, uint n_fields,
	                       const std::vector<std::string> &names );


	int status;              ///< Keeps track of return codes.
	gsd_handle *gh;          ///< Handle to open file.
	uint64_t current_frame;  ///< Current frame in file
	uint64_t max_frame;      ///< Total number of frames at time of opening
	bool eof_, good_;        ///< Flags for file status.

	// These variables are updated whenever a new value of them is
	// encountered, and simultaneously get_next_block relies on these to
	// get default values for data chunks that are absent from the dump.
	uint64_t store_tstep;   ///< Current time step
	uint8_t store_dims;     ///< Dimensions of system
	float store_box[6];     ///< Box dimensions
	uint32_t store_N;       ///< Number of particles

	std::vector<uint32_t>    store_type_ids;  ///< Types.
	std::vector<std::string> store_typenames; ///< Type names.
	std::vector<float>       store_x;         ///< Positions.

	// Enumerate optional fields as bitset.
	enum optional_data {
		BODY = 0,
		MASS,
		CHARGE,
		DIAMETER,
		MOMENT_INERTIA,
		ORIENTATION,
		VELOCITY,
		ANGMOM,
		IMAGE,

		// Dummy that is the number of total optional data entries:
		N_OPTIONAL_DATA_FIELDS
	};

	std::bitset<N_OPTIONAL_DATA_FIELDS> optional_data_found;

	std::vector<int32_t>     store_body;      ///< Body.

	// There's a whole bunch of optional stuff that could be in data file:
	std::vector<float>   store_v;              ///< Velocities.
	std::vector<float>   store_orient;         ///< Orientation (quaternion)
	std::vector<float>   store_moment_inertia; ///< Moment of inertia.


};


} // namespace readers

} // lammps_tools

#endif //DUMP_READER_HOOMD_GSD_HPP
