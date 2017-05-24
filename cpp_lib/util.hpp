#ifndef UTIL_HPP
#define UTIL_HPP

/**
   \file util.hpp

   Declarations/definitions of various helper functions.
*/

#include <algorithm>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <sstream>
#include <vector>

#include "zip.hpp"

namespace lammps_tools {


/**
   \brief Wraps useful functions in a namespace to prevent name clashes.
*/
namespace util {

/**
   Checks if string starts with substring.

   \param s      String to check.
   \param begin  Beginning to check for.

   \returns true if s starts with begin, false otherwise.
*/
inline bool starts_with( const std::string &s, const std::string &begin )
{
	return s.compare(0, begin.length(), begin ) == 0;
}

/**
   Checks if string ends with substring.

   \param s      String to check.
   \param end    Ending to check for.

   \returns true if s ends with end, false otherwise.
*/
inline bool ends_with( const std::string &s, const std::string &end )
{
	if( end.length() > s.length() ){
		return false;
	}else{
		return s.compare(s.length() - end.length(), s.length(), end ) == 0;
	}
}


/**
   Strips string starting at given delimiter.

   \param s      String to strip.
   \param delim  Delimiter to strip at.

   \returns s without the given delimiter and everything following it.
            s is not modified if delim could not be found.
*/
inline std::string rstrip( std::string s, char delim )
{
	std::size_t idx = s.find( delim );
	if( idx == s.length() ) return s;

	s = s.substr( 0, idx );
	std::size_t top = s.find_last_not_of( " \t" ) + 1;
	return s.substr( 0, top );
}


/**
   Splits line into a list of words, splitting on any whitespace.

   \param    line The line to split.
   \returns       A vector of words.
 */
std::vector<std::string> split( std::string line );


/**
   Counts the number of words in given string.

   \param s   String whose words to count.

   \returns   The number of words in given string.
*/
inline int word_count( const std::string &s )
{
	std::stringstream ss(s);
	std::string w;
	int count = 0;
	while( ss >> w ) ++count;
	return count;
}


/**
   Checks if container contains value.

   \param l    Container to check.
   \param val  Value to look for.

   \returns true if container contains val, false otherwise.
*/
template <typename container_type, typename val_type>
inline bool contains( const container_type &l, const val_type &val )
{
	return std::find( l.begin(), l.end(), val ) != l.end();
}

/**
   Specialises contains for strings.

   \overloads contains

   \param s   String to check.
   \param sub Substring to find.

   \returns true if container contains val, false otherwise.
*/
inline bool str_contains( const std::string &s, const std::string &sub )
{
	return s.find(sub) != std::string::npos;
}


/**
   Returns maximum of \p a and \p b.

   \param a
   \param b
*/
template <typename T1, typename T2> inline
T1 max( T1 a, T2 b )
{
	return a > b ? a : b;
}

/**
   Returns minimum of \p a and \p b.

   \param a
   \param b
*/
template <typename T1, typename T2> inline
T1 min( T1 a, T2 b )
{
	return a > b ? b : a;
}


/**
   Calculates average of all entries in container.

   \param c  Container whose values to average.

   \returns  The average of all values in container.
*/
template <typename container_type> inline
double mean( const container_type &c )
{
	double s = std::accumulate( c.begin(), c.end(), 0.0 );
	return s / static_cast<double>( c.size() );
}

/**
   Calculates variance of all entries in container.

   \param c  Container whose values to calculate the variance for.

   \returns  The variance of all values in container.
*/
template <typename container_type> inline
double var( const container_type &c )
{
	double m = mean(c);
	double s = std::inner_product( c.begin(), c.end(), c.begin(), 0.0 );
	return s / static_cast<double>( c.size() - 1 ) - m*m;
}


/**
   Sets bit \p bit to on.

   \param x  Variable whose bit to set.
*/
template <unsigned int bit, typename T> inline
void set_bit( T& x )
{
	x |= (1u << (bit-1) );
}

/**
   Sets bit \p bit to off.

   \param x  Variable whose bit to clear.
*/
template <unsigned int bit, typename T> inline
void clear_bit( T& x )
{
	x &= ~(1u << (bit-1) );
}

/**
   Bitwise if.

   \param x  Variable whose bit to check
*/

template <unsigned int bit, typename T> inline
bool is_bit( T& x )
{
	return x & (1u << (bit-1) );
}


/**
   Checks if string is a non-negative number.

   \param s  String to check.

   \returns  True if s is non-negative, false otherwise.
*/
inline bool is_non_negative( const std::string &s )
{
	int i = 0;
	while( s[i] ){
		if( !std::isdigit( s[i] ) ){
			return false;
		}
	}
	return true;
}

/**
   Checks if file name exists.

   \param  fname File name to check.

   \returns True if file with name exists, false otherwise.
*/
inline bool file_exists( const std::string &fname )
{
	return std::ifstream(fname).good();
}

/**
   Checks if entry is unique in container.

   \param c container to check.
   \param v value to check.

   \returns True if v appears once in container, false otherwise.
*/
template <typename container, typename val>
inline bool is_unique( const container &c, const val &v )
{
	return std::count( c.begin(), c.end(), v ) == 1;
}

/**
   Code to generate a permutation that will sort given vector according to comp.

   Code taking from Timothy Shields at "http://stackoverflow.com/questions/
   	17074324/how-can-i-sort-two-vectors-in-the-same-way-with-
	criteria-that-uses-only-one-of"

   \param vec  The vector to generate the sorting permutation for.
   \param comp The comparator to use.

   \returns a permutation vector p with the sorting permutation, such that
            for( std::size_t i : p ) vec[ i ] is sorted.
*/
template <typename T, typename comparator>
std::vector<std::size_t> sort_permutation( const std::vector<T> &vec,
                                           const comparator &comp )
{
	std::vector<std::size_t> p(vec.size(), 0);
	std::iota(p.begin(), p.end(), 0);
	auto apply_comp = [&](std::size_t i, std::size_t j)
		{ return comp( vec[i], vec[j] ); };
	std::sort( p.begin(), p.end(), apply_comp );
	return p;
}


/**
   Code to apply a permutation in place to a given vector.

   Code taking from Timothy Shields at "http://stackoverflow.com/questions/
   	17074324/how-can-i-sort-two-vectors-in-the-same-way-with-
	criteria-that-uses-only-one-of"

   \param vec  The vector to sort according to permutation p
   \param p    The permutation to use.
*/
template <typename T>
void apply_permutation( std::vector<T>& vec, const std::vector<std::size_t>& p )
{
	std::vector<bool> done(vec.size());
	for (std::size_t i = 0; i < vec.size(); ++i) {
		if (done[i]) continue;

		done[i] = true;
		std::size_t prev_j = i;
		std::size_t j = p[i];

		while( i != j ){
			std::swap(vec[prev_j], vec[j]);
			done[j] = true;
			prev_j = j;
			j = p[j];
		}
	}
}

/**
   Removes double entries from given container in O(n log(n)) complexity.

   \param c  reference to container whose elements to remove.
*/
template <typename container>
void remove_doubles( container &c )
{
	std::sort( c.begin(), c.end() );
	c.erase( std::unique( c.begin(), c.end() ), c.end() );
}


/**
   \brief custom implementation of make_unique for
   backwards compatibility with C++11.
*/
template <typename T>
std::unique_ptr<T> make_unique( T *t )
{
	return std::unique_ptr<T>(t);
}


// Iterates

} // namespace util

} // namespace lammps_tools

#endif // UTIL_HPP
