#include "block_data.hpp"
#include "block_data_access.hpp"

#include "cluster_finder.hpp"
#include "markov_state_capsid.hpp"
#include "my_assert.hpp"
#include "neighborize.hpp"

#include <iostream>

namespace lammps_tools {

namespace msm {

std::ostream &operator<<( std::ostream &o,
                          const msm_id_capsid::edge &e )
{
	int smaller = std::min( e.i, e.j );
	int larger  = std::max( e.i, e.j );
	o << "(" << smaller << ", " << larger << ")";

	return o;
}



bool is_in( int idx, const std::vector<int> &cluster )
{
	return std::find( cluster.begin(), cluster.end(), idx )
		!= cluster.end();
}




void msm_id_capsid::cluster2edges( std::vector<edge> &edges,
                                   const std::vector<int> &cluster,
                                   const neigh_list &mol_nlist )
{
	for( int idx : cluster ){
		for( std::size_t jdx : mol_nlist[idx] ){
			bool add = true;

			// Make sure to add edge only if it is unique:
			edge ee( idx, jdx );
			for( const auto &e : edges ){
				if( e == ee ){
					add = false;
					break;
				}
			}
			if( add ) edges.push_back( ee );
		}
	}
}



double msm_id_capsid::shortest_cluster_dist(
	const block_data &b,
	const std::vector<int> &neighc,
	const std::vector<int> &largest_cluster,
	const std::vector<int> &not_in_cluster,
	const std::vector<int> &mol2com )
{
	double dist = b.dom.xhi[0] - b.dom.xlo[0];

	// You only need to consider the largest cluster, and of those, only
	// the particles inside the cluster with less than 3 bonds.

	std::vector<int> consider_these;
	consider_these.reserve( largest_cluster.size() );
	for( int idx : largest_cluster ){
		if( neighc[idx] < 3 ){
			consider_these.push_back(idx);
		}
	}

	const std::vector<double> &x = get_x( b );
	const std::vector<double> &y = get_y( b );
	const std::vector<double> &z = get_z( b );

	double min_r2 = dist*dist;
	for( int idx : consider_these ){
		int i = mol2com[idx];
		double xi[3] = { x[i], y[i], z[i] };
		for( int jdx : not_in_cluster ){
			int j = mol2com[jdx];
			double r[3];
			double xj[3] = { x[j], y[j], z[j] };
			double r2 = b.dom.dist_2( xi, xj, r );
			if( r2 < min_r2 ){
				min_r2 = r2;
			}
		}
	}

	dist = std::sqrt(min_r2);
	return dist;
}


int msm_id_capsid::prune_edges( std::vector<edge> &edges )
{
	int n_removed = 0;
	std::map< int, int > vertex_counts;

	for( const auto &ei : edges ){
		int idx = ei.i;
		int jdx = ei.j;
		vertex_counts[idx]++;
		vertex_counts[jdx]++;
	}

	for( const auto vc : vertex_counts ){
		if( vc.second == 1 ){
			// Remove all edges with this vertex in them:
			int target = vc.first;
			for( auto ei = edges.begin(); ei != edges.end();  ){
				if( ei->i == target || ei->j == target ){
					ei = edges.erase( ei );
					n_removed++;
				}else{
					++ei;
				}
			}
			return n_removed;
		}
	}

	// No premature return meanse all done:
	return 0;
}


// Calculates the number of different vertices in the given graph
std::size_t msm_id_capsid::vertex_count( const std::vector<edge> &ed )
{
	// std::cerr << "    Size of { ";
	std::size_t size = 0;
	std::map<int,bool> out;
	for( auto ei = ed.begin(); ei != ed.end(); ++ei ){
		//std::cerr << *ei << " ";
		if( !out.count( ei->i ) ){
			++size;
			out[ei->i] = true;
		}
		if( !out.count( ei->j ) ){
			++size;
			out[ei->j] = true;
		}
	}
	//std::cerr << "} = " << size << "\n";
	return size;

}




void msm_id_capsid::mol_cluster2graph( std::vector<edge> &edges,
                                       const std::vector<int> &cluster,
                                       const neigh_list &mol_nlist )
{
	/*
	std::cerr << "  Analyzing cluster with entries";
	for( int idx : cluster ){
		std::cerr << " " << idx;
	}
	std::cerr << "\n";
	*/

	cluster2edges( edges, cluster, mol_nlist );
	/*
	std::cerr << "  Cluster has " << edges.size()
	          << " edges before pruning:";
	for( const auto &e : edges ){
		std::cerr << " " << e;
	}
	std::cerr << "\n";
	*/

	if( prune ){
		int n_removed = 0;
		int n_pruned = 0;
		do{
			n_removed = prune_edges( edges );
			n_pruned += n_removed;
		}while( n_removed > 0 );
	}
}


int msm_id_capsid::to_markov_state( const lammps_tools::block_data &b )
{
	using namespace lammps_tools;

	// Types to consider for binding.
	std::vector<int> type_i = { 3, 5, 7 };
	std::vector<int> type_j = { 4, 6, 8 };
	int type_com = 1;


	// First we need to find the largest cluster.

	const std::vector<int> &mols = get_mol( b );
	const std::vector<int> &types = get_type( b );
	std::size_t Nmols = *std::max_element( mols.begin(), mols.end() );
	std::vector<int> mol2com( Nmols + 1 );

	int nmols = *std::max_element( mols.begin(), mols.end() );
	int atoms_per_mol = b.N / nmols;
	int target_atom = 25180;
	int target_mol = target_atom / atoms_per_mol;

	/*
	  Analysis scheme:

	  1.  Construct molecular network.
	  2.  Count clusters.
	  3.  Find largest cluster
	  4a. If largest cluster is 1, find shortest distrance sd
	  between triangle COMS. If sd < thresh, state = 1, else state = 0
	  4b. Else, state is cluster size
	*/
	std::vector<int> ilist, jlist;

	for( std::size_t idx = 0; idx < b.N; ++idx ){
		int ttype = types[idx];

		if( ttype == type_i[0] ||
		    ttype == type_i[1] ||
		    ttype == type_i[2] ){
			ilist.push_back(idx);
		}
		if( ttype == type_j[0] ||
		    ttype == type_j[1] ||
		    ttype == type_j[2] ){
			jlist.push_back(idx);
		}

		if( ttype == type_com ){
			int moli = mols[idx];
			mol2com[moli] = idx;
		}
	}

	my_assert( ilist.size() == jlist.size() &&
	           "Binding sites not assigned properly!" );

	double scale_factor = static_cast<double>(b.N) / nmols;
	neighborize::neigh_list nl;
	double rc = 1.25;
	double avg_neighs = neighborize::make_list_dist_indexed( nl, b,
	                                                         ilist, jlist,
	                                                         neighborize::DIST_BIN, 3,
	                                                         rc );

	// Step 1:  Construct a list of molecular clusters
	auto mol_nlist = neighborize::get_molecular_connections( b, nl, false );

	// Step 2:  Now we need to convert the molecular connections
	//          to actual clusters.
	auto mol_clusters = neighborize::neigh_list_to_clusters( mol_nlist );


	// At this point, we have the following variables at our disposal:
	//
	// - mol_nlist:    Contains a per-molecule neighbor list, i.e., which
	//                 molecule is less than rc removed from another, in
	//                 terms of the bonding particles' distance.
	//
	// - mol_clusters: Contains a list of molecule clusters. That is, each
	//                 entry is a vector containing molecule IDs that are in
	//                 the same cluster.
	//
	// From this we will construct a bunch of possibly useful criteria:
	//
	// - min_dist:   The shortest distance between a particle in the
	//               largest cluster that is not fully bonded and a
	//               free particle.
	//
	// - mol_neighc: Contains the number of neighbors of each molecule
	//
	// - gamma:      The gamma parameter of Perkett & Hagan, J Chem Phys
	//               140, 214101, 2014.
	//
	// - edges:      Contains the bond topology of the cluster.
	//               It is a vector of graph_edge structs that
	//               encode bonded molecules.
	//
	// - max_bonds:  Number of bonds in largest cluster
	//
	// - max_size:   Number of particles in largest cluster.
	//


	// Step 3: Find largest cluster. Also construct a list of all particles
	//         that are and are not in clusters.


	std::size_t max_i    = 0;
	std::size_t max_size = 0;
	double max_bond_avg = 0;
	std::size_t max_bonds = 0;
	int max_gamma = 0;

	std::vector<int> in_cluster, not_in_cluster;
	std::vector<int> mol_neighc( nmols + 1, 0 );

	// Determine the cluster with the best bonds/particle size.

	for( std::size_t i = 0; i < mol_clusters.size(); ++i ){

		// Construct in_cluster and not_in_cluster vectors:
		const auto &cluster = mol_clusters[i];

		if( cluster.size() > 1 ){
			// Identify all these particles as part of a cluster:
			for( int i : cluster ) in_cluster.push_back(i);
		}else{
			for( int i : cluster ) not_in_cluster.push_back(i);
		}


		// Construct cluster bond topology:
		std::vector<edge> edges;
		mol_cluster2graph( edges, cluster, mol_nlist );
		std::size_t size = cluster.size();

		if( size > 1 ){
			my_assert( size == vertex_count(edges) &&
			           "Cluster size != number of vertices!" );
		}

		std::size_t N_bonds = edges.size();
		double bond_avg = N_bonds / static_cast<double>( size );

		if( size > 20 ){
			// Possibly malformed, ignore it.
			std::cerr << "Ignoring possibly malformed capsid of "
			          << "size " << size << "!\n";
			continue;
		}

		if( bond_avg > max_bond_avg ){

			max_i = i;
			max_size = size;
			max_bonds = N_bonds;
			max_bond_avg = bond_avg;
			max_gamma = 10*( N_bonds - size ) + N_bonds;
		}

		for( const auto &ge : edges ){
			mol_neighc[ge.i]++;
			mol_neighc[ge.j]++;
		}
	}


	// Note: if your bond_average is 1.5 then technically
	// this step is not needed:
	double shortest_dist = 0.0;
	if( max_bonds == 30 && max_size == 20 ){
		shortest_dist = shortest_cluster_dist( b, mol_neighc,
		                                       mol_clusters[max_i],
		                                       not_in_cluster,
		                                       mol2com );
	}

	bool nearby = false;
	if( shortest_dist < dist_thresh ){
		// It is a "nearby" state.
		nearby = true;
	}else{
		// It is a "far" state.
	}

	// There must be some unique way of encoding this...

	// This is at least a unique number for each state:
	int most_bonds_possible = 30;
	int most_mols_possible  = 20;

	int markov_state = max_bonds + most_bonds_possible * max_size;
	if( nearby ){
		markov_state += most_bonds_possible*most_mols_possible;
	}

	return markov_state;
}

} // namespace msm

} // namespace lammps_tools
