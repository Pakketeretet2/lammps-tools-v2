#ifndef MARKOV_STATE_MODEL_HPP
#define MARKOV_STATE_MODEL_HPP

/**
   \file markov_state_model.hpp

   \brief Contains files relating to Markov state modelling

*/

namespace lammps_tools {

// Forward-decl
class block_data;

namespace msm {


// We use a general class to identify the Markov states,
// so that it is more easily extended. All user classes should
// satisfy this interface:

class msm_identifier
{
public:

	msm_identifier() {}
	virtual ~msm_identifier() {}

	// Convert a given block to a state.
	virtual int to_markov_state( const block_data & ) const = 0;

};



} // namespace msm

} // namespace lammps_tools



#endif // MARKOV_STATE_MODEL_HPP
