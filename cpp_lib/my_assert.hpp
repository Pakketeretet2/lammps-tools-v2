#ifndef MY_ASSERT_HPP
#define MY_ASSERT_HPP

/**
   \file my_assert.hpp

   Contains some custom macros for error handling and assertions.
*/

#include <iostream>
#include <string>

namespace lammps_tools {

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
#ifndef NDEBUG
	if( !(test) ){
		std::cerr << file << " ( " << line << " ): Assertion failed: "
		          << msg << "\n";
		std::terminate();
	}
#endif
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
	std::cerr << file << " ( " << line << " ): Logic error: "
	          << msg << "\n";
	std::terminate();
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
	std::cerr << file << " ( " << line << " ): Runtime error: "
	          << msg << "\n";
	std::terminate();
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
		std::cerr << file << " ( " << line << " ): Warning: "
		          << msg << "\n";
	}
}



} // namespace lammps_tools



#endif // MY_ASSERT_HPP
