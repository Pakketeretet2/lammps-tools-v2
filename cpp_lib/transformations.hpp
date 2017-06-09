#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP

#include "block_data.hpp"
#include "quaternion.hpp"

namespace lammps_tools {

namespace transformations {

block_data rotate( block_data b, point axis, point origin );
block_data shift( block_data b, point delta );



} // namespace transformations

} // namespace lammps_tools

#endif // TRANSFORMATIONS_HPP
