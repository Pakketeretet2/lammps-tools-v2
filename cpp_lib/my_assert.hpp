#ifndef MY_ASSERT_HPP
#define MY_ASSERT_HPP

/**
   \file my_assert.hpp

   Contains some custom macros for error handling and assertions.
*/

#include <iostream>
#include <string>
#include <sstream>
#include <exception>

namespace lammps_tools {

#ifdef NDEBUG
constexpr const bool no_debug = true;
#else
constexpr const bool no_debug = false;
#endif // NDEBUG


#ifdef USE_EXCEPTIONS
constexpr const bool use_exceptions = true;
#else
constexpr const bool use_exceptions = false;
#endif // USE_EXCEPTIONS


// Implementation of assertions and the like in terms of std::terminate:
inline void my_assert_terminate( const std::string &file, int line, bool test,
                                 const std::string &msg )
{
	if( !(test) ){
		std::cerr << file << " ( " << line << " ): Assertion failed: "
		          << msg << "\n";
		std::terminate();
	}
}

inline void my_logic_error_terminate( const std::string &file, int line,
                                      const std::string &msg )
{
	std::cerr << file << " ( " << line << " ): Logic error: "
	          << msg << "\n";
	std::terminate();
}
inline void my_runtime_error_terminate( const std::string &file, int line,
                                        const std::string &msg )
{
	std::cerr << file << " ( " << line << " ): Runtime error: "
	          << msg << "\n";
	std::terminate();
}

#ifndef LEGACY_COMPILER
// Implementation of assertions and the like that throw an exception:
inline void my_assert_except( const std::string &file, int line, bool test,
                              const std::string &msg )
{
	if( !test ){
		std::stringstream ss;
		ss << file << " ( " << line << " ): " << msg;
		throw std::runtime_error( ss.str() );
	}
}
inline void my_logic_error_except( const std::string &file, int line,
                                   const std::string &msg )
{
	std::stringstream ss;
	ss << file << " ( " << line << " ): " << msg;
	throw std::logic_error( ss.str() );
}

inline void my_runtime_error_except( const std::string &file, int line,
                                     const std::string &msg )
{
	std::stringstream ss;
	ss << file << " ( " << line << " ): " << msg;
	throw std::runtime_error( ss.str() );
}
#endif // #ifndef LEGACY_COMPILER

/**
   \brief asserts that test is true, and if not, terminates.
          Does nothing if NDEBUG is defined.

   \param file  File name (use __FILE__ macro)
   \param line  Line number (use __LINE__ macro(
   \param test  The test to check.
   \param msg   Error message to print upon failing of test.
*/
inline void my_assert( const std::string &file, int line, bool test,
                       const std::string &msg )
{
	if( !no_debug ){
#ifdef LEGACY_COMPILER
		my_assert_terminate( file, line, test, msg );
#else
		if( use_exceptions ){
			my_assert_except( file, line, test, msg );
		}else{
			my_assert_terminate( file, line, test, msg );
		}
#endif
	}
}


/**
   \brief Prints useful error message and then terminates.

   \param file  File name (use __FILE__ macro)
   \param line  Line number (use __LINE__ macro(
   \param msg   Error message to print before termination.
*/
inline void my_logic_error( const std::string &file, int line,
                            const std::string &msg )
{
	if( !no_debug ){
#ifdef LEGACY_COMPILER
		my_logic_error_terminate( file, line, msg );
#else
		if( use_exceptions ){
			my_logic_error_except( file, line, msg );
		}else{
			my_logic_error_terminate( file, line, msg );
		}
#endif
	}
}

/**
   \brief Prints useful error message and then terminates.

   \param file  File name (use __FILE__ macro)
   \param line  Line number (use __LINE__ macro(
   \param msg   Error message to print before termination.
*/
inline void my_runtime_error( const std::string &file, int line,
                              const std::string &msg )
{
	if( !no_debug ){
#ifdef LEGACY_COMPILER
		my_runtime_error_terminate( file, line, msg );
#else
		if( use_exceptions ){
			my_runtime_error_except( file, line, msg );
		}else{
			my_runtime_error_terminate( file, line, msg );
		}
#endif
	}
}

/**
   \brief Prints a warning.

   \param file   File name (use __FILE__ macro)
   \param line   Line number (use __LINE__ macro)
   \param msg    The warning message to print.
   \param check  Warn if this is true.
*/
inline void my_warning( const std::string &file, int line, const std::string &msg )
{
	std::cerr << file << " ( " << line << " ): Warning: "
	          << msg << "\n";
}


/**
   \brief Prints a warning if assertion is not met.

   \param file   File name (use __FILE__ macro)
   \param line   Line number (use __LINE__ macro)
   \param msg    The warning message to print.
   \param check  Warn if this is true.
*/
inline void my_warning_if( const std::string &file, int line,  bool check,
                           const std::string &msg )
{
	if( check ){
		my_warning( file, line, msg );
	}
}





} // namespace lammps_tools



#endif // MY_ASSERT_HPP
