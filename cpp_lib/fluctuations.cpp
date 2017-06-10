#include "block_data.hpp"
#include "block_data_access.hpp"
#include "fluctuations.hpp"
#include "id_map.hpp"

#include <cmath>
#include <memory>

namespace lammps_tools {

namespace fluctuations {


double rmsd( readers::dump_reader *reader, std::list<double> &rmsds,
             std::list<bigint> &time, const std::list<int> &ids, bool got_ids )
{
	const std::list<int> *my_ids = nullptr;
	block_data b;
	std::cerr << "Grabbing first block.\n";
	int status = reader->next_block(b);
	std::cerr << "Done! Got " << b.N << " atoms and " << b.n_data_fields()
	          << " data fields.\n";
	my_assert( __FILE__, __LINE__, status == 0,
	           "Failed to get block from dump reader!" );
	if( !got_ids ){
		std::list<int> *all = new std::list<int>;
		const std::vector<int> &ids = get_id( b );
		for( int i = 0; i < b.N; ++i ){
			all->push_back( ids[i] );
		}
		my_ids = all;
	}else{
		my_ids = &ids;
	}


	std::cerr << "Fields:";
	for( std::size_t i = 0; i < b.n_data_fields(); ++i ){
		std::cerr << " " << b[i].name;
	}
	std::cerr << "\n";

	std::vector<double> x_avg = get_x(b);
	std::vector<double> y_avg = get_y(b);
	std::vector<double> z_avg = get_z(b);
	std::vector<vec4> rmsd_curr(b.N, {0.0,0.0,0.0,0.0});
	id_map im0( get_id( b ) );
	double alpha = 0.95;
	double total_rmsd = 0.0;
	while( (status = reader->next_block(b)) == 0 ){
		id_map im( get_id(b) );
		const std::vector<double> &x = get_x(b);
		const std::vector<double> &y = get_y(b);
		const std::vector<double> &z = get_z(b);
		double trmsd = 0.0;
		for( int idi : *my_ids ){
			int i0 = im0[idi];
			int ii = im[idi];

			double dx = x_avg[i0] - x[ii];
			double dy = y_avg[i0] - y[ii];
			double dz = z_avg[i0] - z[ii];
			double x2 = dx*dx;
			double y2 = dy*dy;
			double z2 = dz*dz;

			rmsd_curr[i0][0] = x2;
			rmsd_curr[i0][1] = y2;
			rmsd_curr[i0][2] = z2;
			rmsd_curr[i0][3] = x2 + y2 + z2;

			x_avg[i0] = (1-alpha)*x_avg[i0] + alpha*x[ii];
			y_avg[i0] = (1-alpha)*y_avg[i0] + alpha*y[ii];
			z_avg[i0] = (1-alpha)*z_avg[i0] + alpha*z[ii];

			trmsd += rmsd_curr[i0][3];
		}
		rmsds.push_back( std::sqrt( trmsd ) );
		time.push_back( b.tstep );

	}
	if( got_ids ) delete my_ids;
	return total_rmsd;
}


double rmsd( readers::dump_reader *reader, std::list<double> &rmsds,
             std::list<bigint> &time )
{
	std::list<int> dummy;
	return rmsd( reader, rmsds, time, dummy, false );
}


} // namespace fluctuations

} // namespace lammps_tools
