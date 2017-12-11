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

#include "types.hpp"
#include "zip.hpp"


// Forward-declare for some utils.
struct gsd_handle;


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
#ifndef _WIN32
	int i = 0;
	while( s[i] ){
		if( !std::isdigit( s[i] ) ){
			return false;
		}
	}
	return true;
#else
	// Windows sucks:
	int i = 0;
	while( s[i] ){
		if( !iswdigit( s[i] ) ){
			return false;
		}
	}
	return true;
#endif // _WIN32
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

/**
   \brief norm of fixed-length array or vector
*/
template <typename C, std::size_t N> inline
double norm2( const C &c )
{
	double nn2 = 0.0;
	for( std::size_t i = 0; i < N; ++i ){
		nn2 += c[i]*c[i];
	}
	return nn2;
}

template <std::size_t N> inline
double norm2( const double *c )
{
	double nn2 = 0.0;
	for( std::size_t i = 0; i < N; ++i ){
		nn2 += c[i]*c[i];
	}
	return nn2;
}



/**
   \brief dot of fixed-length array or vector
*/
template <typename C, std::size_t N> inline
double dot( const C &c1, const C &c2 )
{
	double ddot = 0.0;
	for( std::size_t i = 0; i < N; ++i ){
		ddot += c1[i]*c2[i];
	}
	return ddot;
}

template <std::size_t N> inline
double dot( const double *c1, const double *c2 )
{
	double ddot = 0.0;
	for( std::size_t i = 0; i < N; ++i ){
		ddot += c1[i]*c2[i];
	}
	return ddot;
}


/**
   \brief extracts next UTF-16 char from char stream

   \param stream The stream to read
   \param c      The extracted UTF character

   \returns the number of bytes the character actually is, or 0 on error.
*/
inline
int next_utf16_char( const char *stream, utf16_char &c, int size )
{
	utf16_char tmp;
	for( int i = 0; i < 4 && i < size; ++i ){
		tmp.c[i] = stream[i];
	}
	c.c[0] = c.c[1] = c.c[2] = c.c[3] = 0;

	if( tmp.c[0] < 0x80 ){
		// One bit suffices, this is a normal ASCII char.
		c.c[0] = tmp.c[0];
		return 1;
	}

	// If still here, if might be UTF-8.
	utf8_char utf8;
	utf8.c[0] = tmp.c[0];
	utf8.c[1] = tmp.c[1];
	if( utf8.d < 0x800 ){
		// It was indeed UTF-8.
		c.c[0] = utf8.c[0];
		c.c[1] = utf8.c[1];
		return 2;
	}

	// If still here, it might be UTF-16.
	utf16_char utf16;
	utf16.c[0] = tmp.c[0];
	utf16.c[1] = tmp.c[1];
	utf16.c[2] = tmp.c[1];

	if( utf16.d < 0x10000 ){
		// It was indeed UTF-8.
		c.c[0] = utf16.c[0];
		c.c[1] = utf16.c[1];
		c.c[2] = utf16.c[2];
		return 3;
	}

	// From here on out I don't support yet.
	return 0;
}


/**
   \brief "downcasts" a utf16 to a utf8 character.

   \param src  The char to downcast
   \param dest The destination to cast to.

   If the downcast could not be made without information loss,
   dest shall be unchanged.

   \returns 0 if the utf16 could be cast without info loss, -1 otherwise.
*/
inline
int down_cast_utf16_char( const utf16_char &src, utf8_char &dest )
{
	if( src.d < 0x800 ){
		// Was utf.
		dest.c[0] = src.c[0];
		dest.c[1] = src.c[1];
		return 0;
	}else{
		return -1;
	}
}

/**
   \brief merges two std containers.
*/
template <typename container1, typename container2> inline
container1 merge( const container1 &c1, const container2 &c2 )
{
	container1 c3;
	c3.reserve( c1.size() + c2.size() );
	c3.insert( c3.end(), c1.begin(), c1.end() );
	c3.insert( c3.end(), c2.begin(), c2.end() );
	return c3;
}


/**
   \brief inserts element in such a way that the container remains sorted.
*/
template <typename container, typename val_type, typename comp_func> inline
typename container::iterator insert_sorted( container &c,
                                            typename container::iterator it,
                                            typename container::iterator end,
                                            val_type v, const comp_func &comp )
{
	my_assert( __FILE__, __LINE__, it <= end, "Invalid begin iterator!" );
	while( it != end && comp( *it, v ) ){
		++it;
	}
	return c.insert( it, v );
}

/**
   \brief inserts element in such a way that the container remains sorted and
          only if element is not in container yet

   This defaults to < as comparator.
   \overload insert_sorted.
*/
template <typename container, typename val_type, typename comp_func> inline
typename container::iterator insert_sorted_unique( container &c,
                                                   typename container::iterator it,
                                                   typename container::iterator end,
                                                   val_type v, const comp_func &comp )
{
	my_assert( __FILE__, __LINE__, it <= end, "Invalid begin iterator!" );
	while( it != end && comp( *it, v ) ){
		++it;
	}
	if( it == end || *it != v ){
		c.insert( it, v );
	}
	return it;
}




/**
   \brief inserts element in such a way that the container remains sorted.
*/
template <typename container, typename val_type, typename comp_func> inline
typename container::iterator insert_sorted( container &c, val_type val,
                                            const comp_func &comp )
{
	return insert_sorted( c, c.begin(), c.end(), val, comp );
}


/**
   \brief inserts element in such a way that the container remains sorted and
          only if element is not in container yet

   This defaults to < as comparator.
   \overload insert_sorted.
*/
template <typename container, typename val_type, typename comp_func> inline
typename container::iterator insert_sorted_unique( container &c, val_type val,
                                                   const comp_func &comp )
{
	return insert_sorted_unique( c, c.begin(), c.end(), val, comp );
}



/**
   \brief inserts element in such a way that the container remains sorted.

   This defaults to < as comparator.
   \overload insert_sorted.
*/
template <typename container, typename val_type> inline
typename container::iterator insert_sorted( container &c, val_type val )
{
	auto comp = []( val_type v1, val_type v2 ){ return v1 < v2; };
	return insert_sorted( c, val, comp );
}


/**
   \brief inserts element in such a way that the container remains sorted and
          only if element is not in container yet

   This defaults to < as comparator.
   \overload insert_sorted_unique.
*/
template <typename container, typename val_type> inline
typename container::iterator insert_sorted_unique( container &c, val_type val )
{
	auto comp = []( val_type v1, val_type v2 ){ return v1 < v2; };
	return insert_sorted_unique( c, val, comp );
}


/**
   \brief finds the arg max.
*/
template <typename container> inline
typename container::const_iterator arg_max( const container &c )
{
	auto it = c.begin();
	auto max_it = it;
	while( it != c.end() ){
		if( *it > *max_it ) max_it = it;
		++it;
	}
	return max_it;
}

template <typename container> inline
typename container::const_iterator arg_min( const container &c )
{
	auto it = c.begin();
	auto min_it = it;
	while( it != c.end() ){
		if( *it < *min_it ) min_it = it;
		++it;
	}
	return min_it;
}

} // namespace util

} // namespace lammps_tools

#endif // UTIL_HPP
