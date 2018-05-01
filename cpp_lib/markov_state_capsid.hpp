#ifndef MARKOV_STATE_CAPSID_HPP
#define MARKOV_STATE_CAPSID_HPP

/**
   \file markov_state_capsid.hpp

   \brief Contains a class specialized to perform markov state modelling
          on capsid-like structures.
*/

#include "markov_state_model.hpp"

#include "neighborize.hpp"

#include <vector>


namespace lammps_tools {

namespace msm {

// We use a general class to identify the Markov states,
// so that it is more easily extended. All user classes should
// satisfy this interface:

class msm_id_capsid : public msm_identifier
{

public:
	/**
	   \brief A simple struct that represents an edge in a graph.
	*/
	struct edge
	{
		edge() : i(0), j(0) {}
		edge( int i, int j ) : i(i), j(j) {}
		int i, j;

		bool operator==( const edge &o ) const
		{
			return (i == o.i && j == o.j) ||
				(i == o.j && j == o.i);
		}

		int operator<( const edge &e2 )
		{
			return ( i == e2.i ) ? j < e2.j : i < e2.i;
		}
	};


	using neigh_list = lammps_tools::neighborize::neigh_list;


	/**
	   \brief Public constructor

	   \param dist_thresh  Distance criterion for determining whether the
	                       closest free molecule is "nearby" or "far".
	   \param prune        If true, dangling edges and vertices are removed
	                       from the graph before categorizing.
	*/
	msm_id_capsid( double dist_thresh, bool prune = false )
		: dist_thresh(dist_thresh), prune(prune)  {}

	/**
	   \brief Convert a given block to a state following
	          the Perkett & Hagan approach.

	   \param[in] b the block to categorize.

	   \returns an integer that classifies the current state.
	*/
	virtual int to_markov_state( const lammps_tools::block_data &b );






private:


	double dist_thresh;
	bool prune;



	/**
	   \brief Converts a given cluster, that is, a vector of molecule or
	   particle IDs, into a set of edges by using the given neighbor list

	   \param edges   Will contain all edges in the cluster's graph
	   \param cluster The cluster to determine the edges for
	   \param n_list  Neighbor list for the data set. That is, n_list[i]
	                  will contain all neigbors of the molecule or particle
	                  _indexed_ by i.
	*/
	void cluster2edges( std::vector<edge> &edges,
	                    const std::vector<int> &cluster,
	                    const neigh_list &mol_nlist );


	/**
	   \brief Calculate the shortest distance between a free particle and
	          a particle in the cluster that is not fully bonded yet.

           \param b                The block_data to analyze
           \param neigc            The number of neighbors of the molecules
           \param largest_cluster  The largest cluster in b
           \param not_in_cluster   All molecules not in a cluster
           \param mol2com          Converts a molecule index to the index of
                                   the center of mass of that molecule.

	   \return the shortest distance between a free molecule and
	           the largest cluster.
	*/
	double shortest_cluster_dist( const lammps_tools::block_data &b,
	                              const std::vector<int> &neighc,
	                              const std::vector<int> &largest_cluster,
	                              const std::vector<int> &not_in_cluster,
	                              const std::vector<int> &mol2com );


	/**
	   \brief "prunes" the given graph, that is, removes edges that
	          connect only a single vertex to the rest.

            \param[in/out] edges  The graph to prune

            \returns the number of edges pruned.
	*/
	int prune_edges( std::vector<edge> &edges );


	/**
	   \brief Count the number of vertices in the graph.

	   \param[in] ed  The vector of edges (or graph).

	   \returns the number of vertices in the graph.
	*/
	std::size_t vertex_count( const std::vector<edge> &ed );


	/**
	   \brief Convert a molecular cluster to a graph. The graph is a
	   collection of edges that encode which molecules in the given
	   cluster are directly connected to each other

	   \param[out] edges      The graph to construct for given cluster
	   \param      cluster    The cluster to construct the graph for
	   \param      mol_nlist  The molecular neighbor list.
	*/
	void mol_cluster2graph( std::vector<edge> &edges,
	                        const std::vector<int> &cluster,
	                        const neigh_list &mol_nlist );
};


} // namespace msm

} // namespace lammps_tools


#endif // MARKOV_STATE_CAPSID_HPP
