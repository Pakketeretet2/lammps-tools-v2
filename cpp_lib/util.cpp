#include "util.hpp"


#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <unistd.h>

#include "gsd.h"

using namespace lammps_tools;

namespace lammps_tools {

namespace util {

std::vector<std::string> split( std::string line )
{
	std::stringstream ss(line);
	std::vector<std::string> words;
	words.reserve( 10 ); // Should be plenty.
	std::string word;
	while( ss >> word ) words.push_back( word );

	return words;
}




} // namespace util

} // namespace lammps_tools
