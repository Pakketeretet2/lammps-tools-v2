#include "block_data.hpp"
#include "block_data_access.hpp"
#include "id_map.hpp"
#include "transformations.hpp"

#include <algorithm>
#include <vector>


namespace lammps_tools {

namespace transformations {

void rotate_all( block_data *b, point axis, point origin, double angle )
{
	rotate( b, axis, origin, angle, get_id(*b) );
}

void rotate( block_data *b, point axis, point origin,
                   double angle, const std::vector<int> &ids )
{
	// To rotate around an origin, center the entire box on origin.
	point shift;
	center_box_on( b, origin, shift );

	std::vector<double> &x = get_x_rw(*b);
	std::vector<double> &y = get_y_rw(*b);
	std::vector<double> &z = get_z_rw(*b);

	double n2 = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
	double inv_na = 1.0 / std::sqrt( n2 );
	axis[0] *= inv_na;
	axis[1] *= inv_na;
	axis[2] *= inv_na;

	double t_h = 0.5*angle;
	double c = std::cos(t_h);
	double s = std::sin(t_h);
	quat rot( c, axis[0]*s, axis[1]*s, axis[2]*s );

	for( int i = 0; i < b->N; ++i ){
		// Standard quaternion rotation:
		quat pos( 0.0, x[i], y[i], z[i] );
		// Conjugate is fine because rot is normalised.
		quat res = (rot*pos)*(rot.conj());
		x[i] = res[1];
		y[i] = res[2];
		z[i] = res[3];
	}

	// Rotate box too?

	// After rotating, undo the centering!
	shift[0] *= -1;
	shift[1] *= -1;
	shift[2] *= -1;

	shift_box( b, shift );

}




void shift_all( block_data *b, point delta )
{
	shift( b, delta, get_id( *b ) );
}


void shift( block_data *b, point delta,
            const std::vector<int> &ids )
{
	std::vector<double> &x = get_x_rw(*b);
	std::vector<double> &y = get_y_rw(*b);
	std::vector<double> &z = get_z_rw(*b);

	data_field_int *ix = static_cast<data_field_int*>(
		b->get_special_field_rw( block_data::IX ) );
	data_field_int *iy = static_cast<data_field_int*>(
		b->get_special_field_rw( block_data::IY ) );
	data_field_int *iz = static_cast<data_field_int*>(
		b->get_special_field_rw( block_data::IZ ) );
	bool image_flags = ix && iy && iz;

	id_map im( get_id(*b) );

	for( int idi : ids ){
		int i = im[idi];
		double xx[3] = { x[i] + delta[0],
		                 y[i] + delta[1],
		                 z[i] + delta[2] };
		// Reapply periodic boundaries:
		if( image_flags ){
			int ixx[3] = { (*ix)[i], (*iy)[i], (*iz)[i] };
			b->dom.rewrap_position(xx, ixx);
			(*ix)[i] = ixx[0];
			(*iy)[i] = ixx[1];
			(*iz)[i] = ixx[2];
		}else{
			b->dom.rewrap_position(xx);
		}

		x[i] = xx[0];
		y[i] = xx[1];
		z[i] = xx[2];
	}
}

void shift_box( block_data *b, point delta )
{
	std::vector<double> &x = get_x_rw(*b);
	std::vector<double> &y = get_y_rw(*b);
	std::vector<double> &z = get_z_rw(*b);

	b->dom.xlo[0] += delta[0];
	b->dom.xlo[1] += delta[1];
	b->dom.xlo[2] += delta[2];

	b->dom.xhi[0] += delta[0];
	b->dom.xhi[1] += delta[1];
	b->dom.xhi[2] += delta[2];

	for( int i = 0; i < b->N; ++i ){
		x[i] += delta[0];
		y[i] += delta[1];
		z[i] += delta[2];
	}
}


void center_box_on( block_data *b, point origin )
{
	point dont_care;
	center_box_on( b, origin, dont_care );
}

void center_box_on( block_data *b, point origin, point &shift )
{
	double *xlo = b->dom.xlo;
	double *xhi = b->dom.xhi;
	double L[3];
	L[0] = xhi[0] - xlo[0];
	L[1] = xhi[1] - xlo[1];
	L[2] = xhi[2] - xlo[2];

	double lambda[3];
	lambda[0] = (origin[0] - xlo[0]) / L[0];
	lambda[1] = (origin[1] - xlo[1]) / L[1];
	lambda[2] = (origin[2] - xlo[2]) / L[2];

	/*
	std::cerr << "In box [ " << xlo[0] << ", " << xlo[1] << ", " << xlo[2]
	          << " ] x [ " << xhi[0] << ", " << xhi[1] << ", " << xhi[2]
	          << " ] the point ( " << origin[0] << ", " << origin[1] << ", "
	          << origin[2] << " ) has lambda coord ( " << lambda[0] << ", "
	          << lambda[1] << ", " << lambda[2] << " ).\n";
	*/
	shift[0] = -lambda[0]*L[0];
	shift[1] = -lambda[1]*L[1];
	shift[2] = -lambda[2]*L[2];

	shift_box( b, shift );
}


void unfold_mols( block_data *b )
{
	my_assert( __FILE__, __LINE__, b->atom_style == ATOM_STYLE_MOLECULAR,
	           "Cannot unfold molecules on data that is not molecular!" );

	std::vector<double> &x = get_x_rw(*b);
	std::vector<double> &y = get_y_rw(*b);
	std::vector<double> &z = get_z_rw(*b);
	const std::vector<int> &mol = get_mol(*b);

	double *xlo = b->dom.xlo;
	double *xhi = b->dom.xhi;
	double L[3];
	L[0] = xhi[0] - xlo[0];
	L[1] = xhi[1] - xlo[1];
	L[2] = xhi[2] - xlo[2];
	double Lt[3] = { L[0]/3.0, L[1]/3.0, L[2]/3.0 };

	// First find all molecules:
	int max_mol = *std::max_element( mol.begin(), mol.end() );
	// enumerate the atom indices per molecule:

	std::vector<std::vector<int> > id_per_mol( max_mol + 1 );
	for( int i = 0; i < b->N; ++i ){
		int mol_id = mol[i];
		id_per_mol[ mol_id ].push_back( i );
	}

	// 0 is left, 1 is center, 2 is right.
	std::vector<int> x_left_right_center( b->N, 0 );
	std::vector<int> y_left_right_center( b->N, 0 );
	std::vector<int> z_left_right_center( b->N, 0 );
	std::vector<bool> needs_remap( max_mol+1, false );
	int remap_count = 0;
	// Now you can iterate over the molecule atoms:
	for( int imol = 1; imol <= max_mol; ++imol ){
		int x_bins_got = 0;
		int y_bins_got = 0;
		int z_bins_got = 0;

		// Check if this needs unfolding or not:
		for( int idx : id_per_mol[ imol ] ){
			double xi = x[idx];
			double yi = y[idx];
			double zi = z[idx];

			int xbin = (xi - xlo[0]) / Lt[0];
			int ybin = (yi - xlo[1]) / Lt[1];
			int zbin = (zi - xlo[2]) / Lt[2];

			x_left_right_center[idx] = xbin;
			y_left_right_center[idx] = ybin;
			z_left_right_center[idx] = zbin;

			if( xbin == 0 ) x_bins_got |= 1;
			if( ybin == 0 ) y_bins_got |= 1;
			if( zbin == 0 ) z_bins_got |= 1;

			if( xbin == 1 ) x_bins_got |= 2;
			if( ybin == 1 ) y_bins_got |= 2;
			if( zbin == 1 ) z_bins_got |= 2;

			if( xbin == 2 ) x_bins_got |= 4;
			if( ybin == 2 ) y_bins_got |= 4;
			if( zbin == 2 ) z_bins_got |= 4;

		}
		if( x_bins_got == 5 || y_bins_got == 5 || z_bins_got == 5 ){
			needs_remap[imol] = true;
			++remap_count;
		}
	}
	std::cerr << remap_count << " of " << max_mol
	          << " mols need remapping.\n";
	int n_remapped = 0;

	// Do the remapping itself:
	for( int imol = 1; imol <= max_mol; ++imol ){
		// There is no unambiguous way to unfold. I just count the
		// numbers left and right, and unfold the least.
		if( !needs_remap[imol] ) continue;

		int n_left[3]  = {0,0,0};
		int n_right[3] = {0,0,0};

		for( int idx : id_per_mol[ imol ] ){
			int xbin = x_left_right_center[idx];
			int ybin = y_left_right_center[idx];
			int zbin = z_left_right_center[idx];

			if( xbin == 0 ) n_left[0]++;
			if( ybin == 0 ) n_left[1]++;
			if( zbin == 0 ) n_left[2]++;

			if( xbin == 2 ) n_right[0]++;
			if( ybin == 2 ) n_right[1]++;
			if( zbin == 2 ) n_right[2]++;
		}

		// Determine along which axes to unfold:
		if( n_left[0] > 0 && n_right[0] > 0 ){
			bool right_to_left = false;
			if( n_left[0] > n_right[0] ){
				// Fold left to right:
				right_to_left = true;
			}

			for( int idx : id_per_mol[ imol ] ){
				if( right_to_left &&
				    x_left_right_center[idx] == 2 ){
					x[idx] -= L[0];
					n_remapped++;
				}else if( !right_to_left &&
				          x_left_right_center[idx] == 0 ){
					x[idx] += L[0];
					n_remapped++;
				}
			}
		}
		if( n_left[1] > 0 && n_right[1] > 0 ){
			// unfold y
			bool right_to_left = false;
			if( n_left[1] > n_right[1] ){
				// Fold left to right:
				right_to_left = true;
			}

			for( int idx : id_per_mol[ imol ] ){
				if( right_to_left &&
				    y_left_right_center[idx] == 2 ){
					y[idx] -= L[1];
					n_remapped++;
				}else if( !right_to_left &&
				          y_left_right_center[idx] == 0 ){
					y[idx] += L[1];
					n_remapped++;
				}
			}
		}
		if( n_left[2] > 0 && n_right[2] > 0 ){
			// unfold z
			bool right_to_left = false;
			if( n_left[2] > n_right[2] ){
				// Fold left to right:
				right_to_left = true;
			}

			for( int idx : id_per_mol[ imol ] ){
				if( right_to_left &&
				    z_left_right_center[idx] == 2 ){
					z[idx] -= L[2];
					n_remapped++;
				}else if( !right_to_left &&
				          z_left_right_center[idx] == 0 ){
					z[idx] += L[2];
					n_remapped++;
				}
			}
		}
	}
	std::cerr << "Done remapping, remapped particles " << n_remapped
	          << " times in 3 dimensions.\n";

}




} // namespace transformations

} // namespace lammps_tools
