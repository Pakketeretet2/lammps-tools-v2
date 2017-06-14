#include "block_data.hpp"
#include "block_data_access.hpp"
#include "id_map.hpp"
#include "transformations.hpp"

namespace lammps_tools {

namespace transformations {

block_data rotate_all( block_data b, point axis, point origin, double angle )
{
	return rotate( b, axis, origin, angle, get_id(b) );
}

block_data rotate( block_data b, point axis, point origin,
                   double angle, const std::vector<int> &ids )
{
	// To rotate around an origin, center the entire box on origin.
	b = center_box_on( b, origin );
	std::vector<double> &x = get_x_rw(b);
	std::vector<double> &y = get_y_rw(b);
	std::vector<double> &z = get_z_rw(b);

	double n2 = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
	double inv_na = 1.0 / std::sqrt( n2 );
	axis[0] *= inv_na;
	axis[1] *= inv_na;
	axis[2] *= inv_na;

	double t_h = 0.5*angle;
	double c = std::cos(t_h);
	double s = std::sin(t_h);
	quat rot( c, axis[0]*s, axis[1]*s, axis[2]*s );

	for( int i = 0; i < b.N; ++i ){
		// Standard quaternion rotation:
		quat pos( 0.0, x[i], y[i], z[i] );
		// Conjugate is fine because rot is normalised.
		quat res = (rot*pos)*(rot.conj());
		x[i] = res[1];
		y[i] = res[2];
		z[i] = res[3];
	}

	// Rotate box too?

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
			b.dom.rewrap_position(xx, ixx);
			(*ix)[i] = ixx[0];
			(*iy)[i] = ixx[1];
			(*iz)[i] = ixx[2];
		}else{
			b.dom.rewrap_position(xx);
		}

		x[i] = xx[0];
		y[i] = xx[1];
		z[i] = xx[2];
	}


	return b;
}

block_data shift_box( block_data b, point delta )
{
	std::vector<double> &x = get_x_rw(b);
	std::vector<double> &y = get_y_rw(b);
	std::vector<double> &z = get_z_rw(b);

	b.dom.xlo[0] += delta[0];
	b.dom.xlo[1] += delta[1];
	b.dom.xlo[2] += delta[2];

	b.dom.xhi[0] += delta[0];
	b.dom.xhi[1] += delta[1];
	b.dom.xhi[2] += delta[2];

	for( int i = 0; i < b.N; ++i ){
		x[i] += delta[0];
		y[i] += delta[1];
		z[i] += delta[2];
	}

	return b;
}


block_data center_box_on( block_data b, point origin )
{
	double *xlo = b.dom.xlo;
	double *xhi = b.dom.xhi;
	double L[3];
	L[0] = xhi[0] - xlo[0];
	L[1] = xhi[1] - xlo[1];
	L[2] = xhi[2] - xlo[2];

	double lambda[3];
	lambda[0] = (origin[0] - xlo[0]) / L[0];
	lambda[1] = (origin[1] - xlo[1]) / L[1];
	lambda[2] = (origin[2] - xlo[2]) / L[2];

	std::cerr << "In box [ " << xlo[0] << ", " << xlo[1] << ", " << xlo[2]
	          << " ] x [ " << xhi[0] << ", " << xhi[1] << ", " << xhi[2]
	          << " ] the point ( " << origin[0] << ", " << origin[1] << ", "
	          << origin[2] << " ) has lambda coord ( " << lambda[0] << ", "
	          << lambda[1] << ", " << lambda[2] << " ).\n";
	double shift[3] = { -lambda[0]*L[0], -lambda[1]*L[1], -lambda[2]*L[2] };

	return shift_box( b, shift );
}



} // namespace transformations

} // namespace lammps_tools
