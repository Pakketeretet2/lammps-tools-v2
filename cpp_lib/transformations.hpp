#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP

#include "block_data.hpp"
#include "geometry.hpp"

namespace lammps_tools {

namespace transformations {

block_data rotate_all( block_data b, point axis, point origin, double angle );
block_data shift_all( block_data b, point delta );

block_data rotate( block_data b, point axis, point origin,
                   double angle, const std::vector<int> &ids );
block_data shift( block_data b, point delta,
                  const std::vector<int> &ids );

block_data center_box_on( block_data b, point origin );
block_data shift_box( block_data b, point delta );


// Interface:



} // namespace transformations

} // namespace lammps_tools

#endif // TRANSFORMATIONS_HPP
