#include "block_data.hpp"
#include "block_data_access.hpp"
#include "id_map.hpp"
#include "transformations.hpp"

namespace lammps_tools {

namespace transformations {

block_data rotate_all( block_data b, point axis, point origin )
{
	return rotate( b, axis, origin, get_id(b) );
}

block_data rotate( block_data b, point axis, point origin,
                   const std::vector<int> &ids )
{
	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);

	return b;
}

block_data shift_all( block_data b, point delta )
{
	return shift( b, delta, get_id( b ) );
}


block_data shift( block_data b, point delta,
                  const std::vector<int> &ids )
{
	std::vector<double> &x = get_x_rw(b);
	std::vector<double> &y = get_y_rw(b);
	std::vector<double> &z = get_z_rw(b);

	data_field_int *ix = static_cast<data_field_int*>(
		b.get_special_field_rw( block_data::IX ) );
	data_field_int *iy = static_cast<data_field_int*>(
		b.get_special_field_rw( block_data::IY ) );
	data_field_int *iz = static_cast<data_field_int*>(
		b.get_special_field_rw( block_data::IZ ) );
	bool image_flags = ix && iy && iz;

	id_map im( get_id(b) );

	for( int idi : ids ){
		int i = im[idi];
		double xx[3] = { x[i] + delta[0],
		                 y[i] + delta[1],
		                 z[i] + delta[2] };

		// Reapply periodic boundaries:
		if( image_flags ){
			int ixx[3] = { (*ix)[i], (*iy)[i], (*iz)[i] };
			b.dom.rewrap(xx, ixx);
			(*ix)[i] = ixx[0];
			(*iy)[i] = ixx[1];
			(*iz)[i] = ixx[2];
		}else{
			b.dom.rewrap(xx);
		}

		x[i] = xx[0];
		y[i] = xx[1];
		z[i] = xx[2];
	}


	return b;
}

} // namespace transformations

} // namespace lammps_tools
