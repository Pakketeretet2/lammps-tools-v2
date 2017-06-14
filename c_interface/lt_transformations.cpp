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
void lt_transformations_rotate_all( lt_block_data_handle *bdh, const double *axis,
                                    const double *origin, double angle )
{
	lammps_tools::block_data *b = bdh->bd;
	lammps_tools::transformations::rotate_all( b, axis,
	                                           origin, angle );
}

/**
   \brief Shifts all particles in block_data over given distance.

   \warning This function _mutates_ the block data!

   \param[in/out] bdh    The block data to rotate
   \param[in]     delta  The distance over which to shift.
*/
void lt_transformations_shift_all( lt_block_data_handle *bdh, const double *delta )
{
	lammps_tools::block_data *b = bdh->bd;
	lammps_tools::transformations::shift_all( b, delta );
}

/**
   \brief Centers the box and all particles on given origin.

   \warning This function _mutates_ the block data!

   \param[in/out] bdh  The block data to rotate
   \param[in]     o    The new origin of the data
*/
void lt_transformations_center_box_on( lt_block_data_handle *bdh, const double *o )
{
	lammps_tools::block_data *b = bdh->bd;
	lammps_tools::transformations::center_box_on( b, o );
}

/**
   \brief Shifts particles with given ids over given distance.

   \warning This function _mutates_ the block data!

   \param[in/out] bdh      The block data to rotate
   \param[in]     delta    The distances over which to shift
   \param[in]     id_size  The size of the id array
   \param[in]     ids      Array containing the ids of particles to shift
*/
void lt_transformations_shift( lt_block_data_handle *bdh, const double *delta,
                               size_t id_size, const int *ids )
{
	lammps_tools::block_data *b = bdh->bd;
	std::vector<int> v_id( ids, ids + id_size );
	lammps_tools::transformations::shift( b, delta, v_id );
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
void lt_transformations_rotate( lt_block_data_handle *bdh, const double *axis,
                                const double *origin, double angle,
                                size_t id_size, const int *ids )
{
	lammps_tools::block_data *b = bdh->bd;
	std::vector<int> v_id( ids, ids + id_size );
	lammps_tools::transformations::rotate( b, axis, origin,
	                                       angle, v_id );
}


/**
   \brief Unfolds atoms belonging to same molecule over PBCs

   \warning The block_data needs to be of ATOM_TYPE MOLECULAR

   \param[in/out] bdh  The block data to unfold.
*/
void lt_transformations_unfold_mols( lt_block_data_handle *bdh )
{
	lammps_tools::transformations::unfold_mols( bdh->bd );
}
