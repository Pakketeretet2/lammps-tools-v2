#include "block_data_access.hpp"
#include "triangulate.hpp"
#include "neighborize_bin.hpp"

#include "util.hpp"

#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>


namespace lammps_tools {

namespace triangulate {

using namespace lammps_tools;


double vec_dist( const double *x1, const double *x2 )
{
	double dx = x1[0] - x2[0];
	double dy = x1[1] - x2[1];
	double dz = x1[2] - x2[2];
	return std::sqrt( dx*dx + dy*dy + dz*dz );
}



int insert_triangle( const block_data &b, int i, int j, int k, int **out,
                     std::vector<triangle> &triangles, std::vector<int> &neighs )
{
	int added = 0;
	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);

	if( (i < k) && (j < k) ){
		if( util::contains( neighs, k) ){
			double xi[3], xj[3], xk[3];
			xi[0] = x[i];
			xi[1] = y[i];
			xi[2] = z[i];

			xj[0] = x[j];
			xj[1] = y[j];
			xj[2] = z[j];

			xk[0] = x[k];
			xk[1] = y[k];
			xk[2] = z[k];

			out[i][j]++;
			out[i][k]++;
			out[j][k]++;

			triangle tri( i, j, k, xi, xj, xk );
			triangles.push_back( tri );
			++added;
		}
	}
	return added ? 1 : 0;
}

void triangulate_block( const block_data &b, double rc, int periodic,
                        int dims, int method, std::vector<triangle> &triangles )
{
	// Make a neighbour list first:
	bigint N = b.N;

	const std::vector<int> &id = get_id( b );
	std::vector<int> all( id.begin(), id.end() );

	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);
	std::vector<std::vector<int> > neighs;

	neighborize::neighborizer_bin n( b, all, all, dims, rc );
	id_map im( id );

	int N2 = N*N;
	int *idx_out_  = new int[N2];
	int **idx_out  = new int*[N];


	for( bigint i = 0; i < N; ++i ){
		idx_out[i] = idx_out_ + N*i;
		for( bigint j = 0; j < N; ++j ){
			idx_out_[i + j*N] = 0;
		}
	}


	for( bigint ii = 0; ii < N; ++ii ){
		// loop over all neighs of i1:

		const std::vector<int> &l = neighs[ii];
		for( int j : l ){
			if( (ii < j) && (idx_out[ii][j] < 3) ){
				// Now you got two ids... ids[i] and j. All you need
				// to know now is if neighs[j] contains ids that are
				// also in neighs[i]! However, if i and j were already
				// added twice, you can safely skip them.
				for( int k : neighs[j] ){
					if( (ii < k) && (j < k) &&
					    util::contains( neighs[k], ii ) ){

						// Add triangle:
						double xi[3], xj[3], xk[3];
						xi[0] = x[ii];
						xi[1] = y[ii];
						xi[2] = z[ii];

						xj[0] = x[j];
						xj[1] = y[j];
						xj[2] = z[j];

						xk[0] = x[k];
						xk[1] = y[k];
						xk[2] = z[k];


						triangles.push_back(
							triangle(ii,j,k,xi,
							         xj, xk ) );
					}
				}
			}
		}
	}


	delete [] idx_out_;
	delete [] idx_out;

}


double triangulation_area( block_data &b, std::vector<triangle> &triangles )
{
	double A = 0.0;
	for( const triangle &t : triangles ){
		A += t.area();
	}
	return A;
}


} // namespace triangulate

} // namespace lammps_tools
