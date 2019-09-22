#include "domain.hpp"

#include "block_data_access.hpp"
#include "data_field.hpp"

#include <algorithm>

using namespace lammps_tools;

namespace lammps_tools {

namespace /* anonymous */ {
	
struct vec3{

	vec3() : x(0), y(0), z(0) {}
	vec3(const vec3& o) : x(o.x), y(o.y), z(o.z) {}
	
	double x, y, z;
	vec3 &operator+=(const vec3 &o) {
		x += o.x;
		y += o.y;
		z += o.z;
		
		return *this;
	}
	vec3 &operator*=(double scalar) {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	vec3 operator*(double scalar) {
		vec3 o(*this);
		o *= scalar;
		return o;
	}
};

} // anonymous namespace	

void swap( domain &f, domain &s )
{
	using std::swap;

	swap( f.xlo, s.xlo );
	swap( f.xhi, s.xhi );
	
	swap( f.periodic, s.periodic );
}


domain::domain( const domain &o )
	: periodic( o.periodic )
{
	std::copy( o.xlo, o.xlo + 3, xlo );
	std::copy( o.xhi, o.xhi + 3, xhi );
	
}


void domain::reconstruct_image_flags(const block_data &b,
                                     std::vector<int> &image_x,
                                     std::vector<int> &image_y,
                                     std::vector<int> &image_z) const
{
	if (image_x.empty() || image_x.size() < b.N)  image_x.resize(b.N);
	if (image_y.empty() || image_y.size() < b.N)  image_y.resize(b.N);
	if (image_z.empty() || image_z.size() < b.N)  image_z.resize(b.N);
	
	
	long max_mol = 0;
	const std::vector<int> &mol = get_mol(b);
	const std::vector<int> &type = get_type(b);

	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);

	
	for (int i = 0; i < mol.size(); ++i) {
		if (mol[i] > max_mol) {
			max_mol = mol[i];
		}
	}
	long n_mols = max_mol+1;
	std::vector<int> atom2mol(b.N);
	std::vector<std::vector<int> > mol2atom(n_mols);
	
	for (int i = 0; i < b.N; ++i) {
		mol2atom[mol[i]].push_back(i);
		atom2mol[i] = mol[i];
	}

	// Now we know which molecule each atom belongs to.
	// We construct the image flags as follows:
	// For the first atom position, if it is inside the box,
	// we use that position. Otherwise, we adjust the image flags so
	// that it is whithin the box.

	
	for (int imol = 0; imol < n_mols; ++imol) {
		if (mol2atom[imol].empty()) continue;

		double xp[3];
		double Lx = b.dom.xhi[0] - b.dom.xlo[0];
		double Lz = b.dom.xhi[1] - b.dom.xlo[1];
		double Ly = b.dom.xhi[2] - b.dom.xlo[2];

		{
			// Corner case for first atom in mol:
			int i = 0;
			// Wrap the first particle inside the box:

			int imx = 0, imy = 0, imz = 0;
			int atom_idx = mol2atom[imol][i];
			image_x[atom_idx] = imx;
		        image_y[atom_idx] = imy;
			image_z[atom_idx] = imz;

			xp[0] = x[atom_idx];
			xp[1] = y[atom_idx];
			xp[2] = z[atom_idx];
			b.dom.rewrap_vector(xp);
			
		}
		
		for (int j = 1; j < mol2atom[imol].size(); ++j) {
			int i = mol2atom[imol][j];
			double xx[3];
			
			xx[0] = x[i];
			xx[1] = y[i];
			xx[2] = z[i];

			// For each of these, we find image flags that minimize
			// the distance to the previous atom in the molecule.
			// This is fine if the molecules are much smaller than
			// the box.
			
			int imx = 0, imy = 0, imz = 0;
			int imx_m = 0, imy_m = 0, imz_m = 0;
			// each r2 is guaranteed to be smaller than this squared
			double r2_min = Lx*Lx + Ly*Ly + Lz*Lz;

			double r[3];
			int ximgs[27] = {0, 1, -1, 0, 1, -1,  0,  1, -1,
			                 0, 1, -1, 0, 1, -1,  0,  1, -1,
			                 0, 1, -1, 0, 1, -1,  0,  1, -1};
			int yimgs[27] = {0, 0,  0, 1, 1,  1, -1, -1, -1,
			                 0, 0,  0, 1, 1,  1, -1, -1, -1,
			                 0, 0,  0, 1, 1,  1, -1, -1, -1};
			int zimgs[27] = {0, 0,  0, 0, 0,  0,  0,  0,  0,
			                 1 ,1,  1, 1, 1,  1,  1,  1,  1,
			                 -1, -1, -1, -1, -1, -1, -1, -1, -1};
			             
			for (int flag_i = 0; flag_i < 27; ++flag_i) {
				double xdum3[3];
				xdum3[0] = xx[0];
				xdum3[1] = xx[1];
				xdum3[2] = xx[2];

				xdum3[0] += Lx*ximgs[flag_i];
				xdum3[1] += Ly*yimgs[flag_i];
				xdum3[2] += Lz*zimgs[flag_i];
				double r2 = b.dom.dist_2(xp, xdum3, r, false);
				if (r2 < r2_min) {
					r2_min = r2;
					imx_m = ximgs[flag_i];
					imy_m = yimgs[flag_i];
					imz_m = zimgs[flag_i];
				}
			}
			
			xp[0] = Lx*imx_m + xx[0];
			xp[1] = Ly*imy_m + xx[1];
			xp[2] = Lz*imz_m + xx[2];

			image_x[i] = imx_m;
		        image_y[i] = imy_m;
			image_z[i] = imz_m;
		}
	}
		
}


} // namespace lammps_tools

