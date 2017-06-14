#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP

#include "block_data.hpp"
#include "geometry.hpp"

namespace lammps_tools {

namespace transformations {

void rotate_all( block_data *b, point axis, point origin, double angle );
void shift_all( block_data *b, point delta );

void rotate( block_data *b, point axis, point origin,
             double angle, const std::vector<int> &ids );
void shift( block_data *b, point delta,
            const std::vector<int> &ids );

void center_box_on( block_data *b, point origin, point &shift );
void center_box_on( block_data *b, point origin );
void shift_box( block_data *b, point delta );
void unfold_mols( block_data *b );


// Non-destructive counterparts:
inline
block_data rotate_all( block_data b, point axis, point origin, double angle )
{
	rotate_all( &b, axis, origin, angle );
	return b;
}

inline
block_data shift_all( block_data b, point delta )
{
	shift_all( &b, delta );
	return b;
}

inline
block_data rotate( block_data b, point axis, point origin,
                   double angle, const std::vector<int> &ids )
{
	rotate( &b, axis, origin, angle, ids );
	return b;
}

inline
block_data shift( block_data b, point delta,
                  const std::vector<int> &ids )
{
	shift( &b, delta, ids );
	return b;
}

inline
block_data center_box_on( block_data b, point origin )
{
	center_box_on( &b, origin );
	return b;
}

inline
block_data shift_box( block_data b, point delta )
{
	shift_box( &b, delta );
	return b;
}





} // namespace transformations

} // namespace lammps_tools

#endif // TRANSFORMATIONS_HPP
