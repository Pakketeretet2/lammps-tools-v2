#include "block_data_access.hpp"
#include "skeletonize.hpp"
#include "triangulate.hpp"
#include "neighborize.hpp"
#include <cmath>
#include <array>

namespace lammps_tools {

namespace skeletonize {

struct graph_vertex
{
	graph_vertex() : i(-1), pt{0,0,0}
	{}
	graph_vertex(int i, const std::array<double,3> &p) : i(i), pt(p)
	{}
	graph_vertex(int i, double x, double y, double z) : i(i), pt{x,y,z}
	{}

	std::size_t i;
	std::array<double,3> pt;
};



template <typename container>
void mean_and_var( const container &c, double &m, double &v )
{
	m = v = 0.0;
	double n = 0.0;
	for( double value : c ){
		m += value;
		n += 1.0;
	}
	m /= n;
	for( double value : c ){
		double dx = value - m;
		v += dx*dx;
	}
	v /= (n-1);
}


double sphere_dist( const std::array<double,3> &xi,
                    const std::array<double,3> &xj, double R )
{
	double R2 = R * R;
	double dot = xi[0]*xj[0] + xi[1]*xj[1] + xi[2]*xj[2];
	double theta = std::acos( dot / R2 );
	return theta*R;

}

double sphere_dist( const double *xi, const double *xj, double R )
{
	double R2 = R * R;
	double dot = xi[0]*xj[0] + xi[1]*xj[1] + xi[2]*xj[2];
	double theta = std::acos( dot / R2 );
	return theta*R;

}




std::vector<double> get_insideness( const class block_data &b,
                                    const neighborize::neigh_list &neighs )
{

	std::vector<double> insideness( b.N );
	for( int i = 0; i < b.N; ++i ){
		insideness[i] = -1.0;
	}


	for( int i = 0; i < b.N; ++i ){
		if( neighs[i].size() < 6 ){
			insideness[i] = 0.0;
		}
	}

	const std::vector<int> &ids = get_id( b );

	// Stategy: Recursively loop over all points in network. Those that
	// are not yet determined to be anywhere are assigned the lowest value
	// of their neighboring insideness plus one.
	id_map im( ids );
	bool assigned_one = false;
	double current_val = 0.0;
	// std::cerr << "Finding insideness. At loop 0...";
	do{
		assigned_one = false;
		for( int i = 0; i < b.N; ++i ){
			if( insideness[i] >= 0 ){
				continue;
			}

			const std::vector<int> &ni = neighs[i];
			bool has_current_val = false;
			for( int idx : ni ){
				if( insideness[idx] == current_val ){
					has_current_val = true;
					break;
				}
			}
			if( has_current_val ){
				insideness[i] = current_val + 1.0;
				assigned_one = true;
			}
		}
		current_val += 1.0;
		// std::cerr << " " << current_val << "...";
	}while( assigned_one );

	// std::cerr << "Done assigning...\n";


	return insideness;
}


void get_edge( const block_data &b, const std::vector<double> &insideness,
               std::list<int> &edge )
{
	for( int i = 0; i < b.N; ++i ){
		if( insideness[i] == 0 ){
			edge.push_back(i);
		}
	}
}


std::vector<double> euclidian_distance_transform( const class block_data &b,
                                                  const std::vector<double> &insideness,
                                                  double R )
{
	std::vector<double> edt( b.N );
	std::list<int> edge;
	get_edge( b, insideness, edge );

	const std::vector<double> &x = data_as<double>(
		b.get_special_field( block_data::X ) );
	const std::vector<double> &y = data_as<double>(
		b.get_special_field( block_data::Y ) );
	const std::vector<double> &z = data_as<double>(
		b.get_special_field( block_data::Z ) );

	for( int i = 0; i < b.N; ++i ){
		if( insideness[i] == 0 ){
			edt[i] = 0.0;
		}else{
			double xi[3];
			xi[0] = x[i];
			xi[1] = y[i];
			xi[2] = z[i];

			double R2 = R*R;
			double min_dist = 16*R;

			// Find the closest edge:
			for( int j : edge ){
				double xx[3];
				xx[0] = x[j];
				xx[1] = y[j];
				xx[2] = z[j];
				double dot = xi[0]*xx[0] + xi[1]*xx[1] + xi[2]*xx[2];
				double theta = std::acos( dot / R2 );
				min_dist = std::min( min_dist, theta*R );
			}
			edt[i] = min_dist;
		}
	}

	return edt;
}



void get_local_maxima( const std::vector<std::vector<int> > & neighs,
                       const std::vector<double> &field,
                       std::vector<int> &max_indices,
                       std::function< double(double) > f )
{
	// d_out << "Indices of maxima:";
	for( std::size_t idx = 0; idx < neighs.size(); ++idx ){
		double vi = f( field[idx] );
		bool largest = true;
		for( std::size_t idj : neighs[idx] ){
			if( f( field[idj] ) >= vi ){
				largest = false;
				break;
			}
		}
		if( neighs[idx].empty() ){
			largest = false;
		}
		if( largest ){
			// d_out << " " << idx;
			max_indices.push_back( idx );
		}
	}
	// d_out << "\n";
}

void get_local_maxima( const std::vector<std::vector<int> > &neighs,
                       const std::vector<double> &field,
                       std::vector<int> &max_indices )
{
	auto unit_operator = []( double x ){ return x; };
	get_local_maxima( neighs, field, max_indices, unit_operator );
}



void get_stats( const std::vector<double> &values,
                const std::vector<int> &idx,
                double &avg, double &var, double &min_val, double &max_val,
		std::vector<double> &largest )
{
	avg = var = min_val = max_val = 0.0;
	std::vector<double> sorted( idx.size() );
	std::size_t j = 0;
	for( std::size_t i : idx ){
		sorted[j] = values[i];
		j++;
	}

	mean_and_var( sorted, avg, var );

	std::sort( sorted.begin(), sorted.end() );
	min_val = sorted[0];
	max_val = sorted[sorted.size()-1];
	for( std::size_t i = 0; i < largest.size(); ++i ){
		largest[i] = sorted[sorted.size()-1-i];
	}
}


// There is no get_ribbon_data yet.


void skeletonize_edt( const class block_data &b, std::vector<int> &skeleton,
                      const std::vector<double> &insideness, double R )
{
	std::vector<double> edt = euclidian_distance_transform( b, insideness, R );
	//skeletonize( b, skeleton, edt, R );
}


void skeletonize_alg1( const block_data &b, std::vector<int> &skeleton,
                       const std::vector<double> &distance_map, double R )
{
	// Algorithm:
	// 0. Find point with largest distance map, i
	// 1. Find neighbor of that point j for which
	//    (DM(i) - DM(j)) / (xi - xj) is minimal.


	std::vector<bool> out( distance_map.size() );
	for( std::size_t i = 0; i < out.size(); ++i ){
		out[i] = false;
	}
	int i = 0;
	double d_max = 0.0;
	for( std::size_t j = 0; j < distance_map.size(); ++j ){
		if( distance_map[j] > d_max ){
			d_max = distance_map[j];
			i = j;
		}
	}
	std::cerr << "max dist = " << d_max << " for particle " << i << ".\n";
	out[i] = true;

	std::vector<std::vector<int> > neighs;
	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);

	int dims = 3;
	double rc = 1.3;

	neighborize::make_list_dist( neighs, b, 0, 0,
	                             neighborize::DIST_BIN, dims, rc );

	skeleton.push_back(i);

	bool got_one_in = false;
	do{
		// Check the particles in the neighborhood of i
		// that are not yet out, and find the one that has
		// the most positive gradient in EDT.
		for( int j : neighs[i] ){
			double max_grad = -1000;
			int max_j = j;
			got_one_in = false;
			if( !out[j] ){
				// Compute the gradient in EDT:
				double d_edt = distance_map[j] - distance_map[i];
				double xi[3], xj[3];
				xi[0] = x[i];
				xi[1] = y[i];
				xi[2] = z[i];

				xj[0] = x[j];
				xj[1] = y[j];
				xj[2] = z[j];

				double d_x   = sphere_dist( xi, xj, R );
				double grad = d_edt / d_x;
				std::cerr << "Gradient between " << i << " and "
				          << j << " = " << grad << ".\n";
				got_one_in = true;
				if( grad > max_grad ){
					max_grad = grad;
					max_j = j;
				}
			}
			// Select j as the new i.
			skeleton.push_back(max_j);
			out[j] = true;
			i = j;
		}
	}while( got_one_in );
}



void skeletonize( const block_data &b, std::vector<int> &skeleton,
                  const std::vector<double> &distance_map, double R )
{
	skeletonize_alg1( b, skeleton, distance_map, R );
}


std::vector<double> get_ribbon_widths( const block_data &b,
                                       const neighborize::neigh_list &neighs,
                                       const std::vector<double> &edt,
                                       const std::vector<double> &insideness )
{
	std::vector<int> maxima;
	std::vector<double> widths;
	get_local_maxima( neighs, edt, maxima );

	std::list<int> edge;
	for( int i = 0; i < b.N; ++i ){
		if( insideness[i] == 0.0 ) edge.push_back(i);
	}

	const std::vector<double> &x = get_x(b);
	const std::vector<double> &y = get_y(b);
	const std::vector<double> &z = get_z(b);


	// Get the distance to the edge.
	for( int i : maxima ){
		// Get the shortest distance to any edge particle:
		double xi[3] = { x[i], y[i], z[i] };
		double shortest_dist2 = 100000000.0;

		for( int j : edge ){
			double xj[3] = { x[j], y[j], z[j] };
			double rr[3];
			double r2 = b.dom.dist_2( xi, xj, rr );
			if( r2 < shortest_dist2 ){
				shortest_dist2 = r2;
			}
		}
		/*
		std::cerr << "Closest atom to " << i << "( id = " << id[i]
		          << " ) is " << shortest_idx << " ( id = "
		          << id[shortest_idx] << " )\n";
		*/
		widths.push_back( std::sqrt( shortest_dist2 ) );
	}

	return widths;
}


std::vector<double> neighbor_strain( block_data &b, double r0,
                                     int itype, int jtype, int method,
                                     int dims, double rc )
{
	const std::vector<double> &x = get_x( b );
	const std::vector<double> &y = get_y( b );
	const std::vector<double> &z = get_z( b );
	std::vector<double> strain( b.N );

	neighborize::neigh_list nl;
	neighborize::make_list_dist( nl, b, itype, jtype, method, dims, rc );

	for( bigint i = 0; i < b.N; ++i ){
		int nc = nl[i].size();
		double s_avg = 0.0;
		double xi[3] = { x[i], y[i], z[i] };
		double c = 1.0 / nc;
		for( int j : nl[i] ){
			double xj[3] = { x[j], y[j], z[j] };
			double rr[3];
			double r2 = b.dom.dist_2( xi, xj, rr );
			double r = std::sqrt(r2);
			s_avg += (r - r0) * c;
		}
		strain[i] = s_avg;
	}

	return strain;
}



} // namespace skeletonize

} // namespace lammps_tools
