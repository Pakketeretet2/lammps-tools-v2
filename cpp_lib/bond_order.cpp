#include "bond_order.hpp"
#include "block_data_access.hpp"
#include "constants.hpp"
#include "fast_math.hpp"
#include "neighborize.hpp"

namespace lammps_tools {

namespace order_parameters {

double angle_2pi( const point &v1, const point &v2 )
{
	double d = util::dot( v1, v2 );
	point c = util::cross( v1, v2 );
	double v1n = util::dot( v1, v1 );
	double v2n = util::dot( v2, v2 );

	bool v1_unit = fabs( v1n - 1.0 ) < 1e-12;
	bool v2_unit = fabs( v2n - 1.0 ) < 1e-12;

	if( !v1_unit ){
		std::cerr << "v1 is not a unit vector but ( "
		          << v1[0] << " " << v1[1] << " " << v1[2] << ")!\n";
	}
	if( !v2_unit ){
		std::cerr << "v2 is not a unit vector but ( "
		          << v2[0] << " " << v2[1] << " " << v2[2] << ")!\n";
	}

	my_assert( __FILE__, __LINE__, v1_unit && v2_unit,
	           "v1 and/or v2 is not a unit vector!" );

	double a = std::acos(d);

	if( c.z >= 0 ){
		// Positive angle:
		return a;
	}else{
		// Negative angle:
		return fast_math::pi2 - a;
	}
}

double compute_psi_n( const block_data &b,
                      const neighborize::neigh_list &neighs,
                      int n, const point &axis,
                      std::vector<double> &psi_n_real,
                      std::vector<double> &psi_n_imag )
{
	// 1. Calculate all bonds and bond angles.
	// 2. For each bond, determine psi_n and add to psi_n of atoms involved.
	// 3. Average for each atom.

	std::vector<double> angles;
	std::vector<bond> bonds_from = neighborize::neigh_list_to_bonds( b, neighs );
	const std::vector<bond> &bonds = bonds_from;
	const std::vector<int> &id = get_id(b);

	// 1
	relative_bond_angles( b, neighs, axis, bonds, angles );
	std::size_t N = b.N;
	if( psi_n_real.size() != N ) psi_n_real.resize( N );
	if( psi_n_imag.size() != N ) psi_n_imag.resize( N );

	std::vector<double> bond_counts( N, 0.0 );

	// 2
	for( std::size_t ii = 0; ii < bonds.size(); ++ii ){
		const bond &bb = bonds[ii];
		double angle = angles[ii];

		int i = bb.particle1;
		int j = bb.particle2;

		double psi_n_r, psi_n_i;
		// fast_math::sincos( n*angle, psi_n_i, psi_n_r );
		double f_angle = n*angle;
		psi_n_i = std::sin( f_angle );
		psi_n_r = std::cos( f_angle );

		psi_n_real[i] += psi_n_r;
		psi_n_imag[i] += psi_n_i;

		psi_n_real[j] += psi_n_r;
		psi_n_imag[j] += psi_n_i;

		bond_counts[i] += 1.0;
		bond_counts[j] += 1.0;
	}

	// 3:
	double psi_imag_avg = 0.0;
	double psi_real_avg = 0.0;
	for( std::size_t i = 0; i < N; ++i ){
		if( bond_counts[i] == 0 ) continue;

		psi_n_real[i] /= bond_counts[i];
		psi_n_imag[i] /= bond_counts[i];

		psi_imag_avg += psi_n_imag[i];
		psi_real_avg += psi_n_real[i];
	}

	double inv_N = 1.0 / N;
	psi_real_avg *= inv_N;
	psi_imag_avg *= inv_N;

	double abs_psi = psi_real_avg*psi_real_avg + psi_imag_avg*psi_imag_avg;

	return std::sqrt( abs_psi );
}

void relative_bond_angles( const block_data &b,
                           const neighborize::neigh_list &neighs,
                           point axis,
                           const std::vector<bond> &bonds,
                           std::vector<double> &angles )
{
	// Loop over all bonds:
	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);


	double ax_norm2 = util::dot( axis, axis );
	double inv_ax_norm  = 1.0 / std::sqrt( ax_norm2 );
	axis.x *= inv_ax_norm;
	axis.y *= inv_ax_norm;
	axis.z *= inv_ax_norm;

	for( const bond &bb : bonds ){
		int i = bb.particle1;
		int j = bb.particle2;
		// Force z == 0.
		my_assert( __FILE__, __LINE__, std::fabs( z[i] ) < 1e-12,
		           "Z-coefficient of particle i not 0!" );
		my_assert( __FILE__, __LINE__, std::fabs( z[j] ) < 1e-12,
		           "Z-coefficient of particle j not 0!" );

		double xi[3] = { x[i], y[i], 0.0 };
		double xj[3] = { x[j], y[j], 0.0 };
		double rij[3];
		double r2 = b.dom.dist_2( xj, xi, rij );

		double inv_rr = 1.0 / std::sqrt(r2);
		point rr_axis( rij );

		rr_axis.x *= inv_rr;
		rr_axis.y *= inv_rr;
		rr_axis.z *= inv_rr;

		double angle = angle_2pi( rr_axis, axis );
		angles.push_back( angle );
	}
	my_assert( __FILE__, __LINE__, bonds.size() == angles.size(),
	           "Size of angles and bonds mismatches!" );
}


void relative_bond_angles_( const block_data &b,
                            const neighborize::neigh_list &neighs,
                            const std::vector<double> &axis,
                            std::vector<double> &angles,
                            std::vector<bond> &bonds )
{
	point aaxis( axis[0], axis[1], axis[2] );
	relative_bond_angles( b, neighs, aaxis, bonds, angles );
}

} // order_parameters

} // lammps_tools
