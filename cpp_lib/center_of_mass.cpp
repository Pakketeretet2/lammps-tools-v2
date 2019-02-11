#include "block_data.hpp"
#include "block_data_access.hpp"
#include "center_of_mass.hpp"
#include "id_map.hpp"

#include <list>
#include <vector>

namespace lammps_tools {

point center_of_mass( const block_data &b )
{
	std::vector<int> s;
	all_ids(b, s);
	return center_of_mass( b, s.begin(), s.end() );
}



point geometric_center( const block_data &b )
{
	std::vector<int> s;
	all_ids(b, s);
	return geometric_center( b, s.begin(), s.end() );
}


template <typename iter>
point center_of_mass( const block_data &b, iter it, iter end )
{
	auto w = []( const block_data &b, int idx )
	{ int ttype = get_type(b)[idx];
	  double m = b.ati.mass[ttype];
	  return m; };

	return weighted_pos_avg( b, it, end, w );
}

template <typename iter>
point geometric_center( const block_data &b, iter it, iter end )
{
	auto w = []( const block_data &b, int idx )
	{ return 1.0; };
	return weighted_pos_avg( b, it, end, w );
}



template <typename iter, typename functor>
point weighted_pos_avg( const block_data &b, iter it, iter end, functor w )
{
	point p;

	id_map im( get_id(b) );
	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);

	// Check if the block has image flags:
	const std::vector<int> &ix = get_ix(b);
	const std::vector<int> &iy = get_iy(b);
	const std::vector<int> &iz = get_iz(b);


	double c = 0.0;
	for( iter i = it; i != end; ++i ){
		int idi = *i;
		int idx = im[idi];

		// Convert position with the image flags.

		double ww = w( b, idx );

		p.x += ww*x[ idx ];
		p.y += ww*y[ idx ];
		p.z += ww*z[ idx ];

		c += ww;
	}
	double inv_c = 1.0 / c;
	p.x *= inv_c;
	p.y *= inv_c;
	p.z *= inv_c;

	return p;
}

using v_c_iter = std::vector<int>::const_iterator;
using l_c_iter = std::vector<int>::const_iterator;

template <>
point center_of_mass<v_c_iter>(
	const block_data &b, v_c_iter it, v_c_iter end );
template <>
point center_of_mass<l_c_iter>(
	const block_data &b, l_c_iter it, l_c_iter end );

template <>
point geometric_center<v_c_iter>(
	const block_data &b, l_c_iter it, v_c_iter end );
template <>
point geometric_center<l_c_iter>(
	const block_data &b, l_c_iter it, l_c_iter end );




} // lammps_tools
