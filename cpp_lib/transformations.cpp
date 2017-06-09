#include "transformations.hpp"

block_data rotate( block_data b, point axis, point origin )
{
	std::vector<double> &x = get_x(b);
	std::vector<double> &y = get_y(b);
	std::vector<double> &z = get_z(b);

	return b;
}

block_data shift( block_data b, point delta )
{
	std::vector<double> &x = get_x(b);
	std::vector<double> &y = get_y(b);
	std::vector<double> &z = get_z(b);

	return b;
}
