#include "lt_transformations.h"
#include "../cpp_lib/transformations.hpp"
#include "../cpp_lib/block_data_access.hpp"

#include <iostream>
#include <vector>

/**
   \brief Rotates all particles in block_data about an axis in given point.

   \warning This function _mutates_ the block data!

   \param[in/out] bdh    The block data to rotate
   \param[in]     axis   The axis of rotation
   \param[in]     origin The origin of the axis of rotation
   \param[in]     angle  The angle with which to rotate, in radians.
*/
void lt_transformations_rotate_all( lt_block_data_handle bdh, const double *axis,
                                    const double *origin, double angle )
{
	const std::vector<double> &x = lammps_tools::get_x( *bdh.bd );
	std::cerr << "x before is " << x[0] << "\n";
	lammps_tools::block_data tmp;
	tmp = lammps_tools::transformations::rotate_all( *bdh.bd, axis,
	                                                 origin, angle );
	*bdh.bd = tmp;
	std::cerr << "x after  is " << x[0] << "\n";
}

/**
   \brief Shifts all particles in block_data over given distance.

   \warning This function _mutates_ the block data!

   \param[in/out] bdh    The block data to rotate
   \param[in]     delta  The distance over which to shift.
*/
void lt_transformations_shift_all( lt_block_data_handle bdh, const double *delta )
{
	*bdh.bd = lammps_tools::transformations::shift_all( *bdh.bd, delta );
}

/**
   \brief Centers the box and all particles on given origin.

   \warning This function _mutates_ the block data!

   \param[in/out] bdh  The block data to rotate
   \param[in]     o    The new origin of the data
*/
void lt_transformations_center_box_on( lt_block_data_handle bdh, const double *o )
{
	*bdh.bd = lammps_tools::transformations::center_box_on( *bdh.bd, o );
}

/**
   \brief Shifts particles with given ids over given distance.

   \warning This function _mutates_ the block data!

   \param[in/out] bdh      The block data to rotate
   \param[in]     delta    The distances over which to shift
   \param[in]     id_size  The size of the id array
   \param[in]     ids      Array containing the ids of particles to shift
*/
void lt_transformations_shift( lt_block_data_handle bdh, const double *delta,
                               size_t id_size, const int *ids )
{
	std::vector<int> v_id( ids, ids + id_size );
	*bdh.bd = lammps_tools::transformations::shift( *bdh.bd, delta, v_id );
}

/**
   \brief Shifts particles with given ids over given distance.

   \warning This function _mutates_ the block data!

   \param[in/out] bdh      The block data to rotate
   \param[in]     axis     The axis of rotation
   \param[in]     origin   The origin of the axis of rotation
   \param[in]     angle	   The angle with which to rotate, in radians
   \param[in]     id_size  The size of the id array
   \param[in]     ids	   Array containing the ids of particles to shift
*/
void lt_transformations_rotate( lt_block_data_handle bdh, const double *axis,
                                const double *origin, double angle,
                                size_t id_size, const int *ids )
{
	std::vector<int> v_id( ids, ids + id_size );
	*bdh.bd = lammps_tools::transformations::rotate( *bdh.bd, axis, origin,
	                                                 angle, v_id );
}
