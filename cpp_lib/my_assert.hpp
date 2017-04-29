#ifndef MY_ASSERT_HPP
#define MY_ASSERT_HPP

/**
   \file my_assert.hpp
   
   Contains some custom macros for error handling and assertions.
*/

#include <iostream>

namespace lammps_tools {

#define FILELINE __FILE__,__LINE__


#ifdef NDEBUG
#define my_assert( file, line, test, msg )  \
	do {} while(0)
#else
/**
   \brief asserts that test is true, and if not, terminates.
          Does nothing if NDEBUG is defined.
   
   \param file  File name (use __FILE__ macro)
   \param line  Line number (use __LINE__ macro(
   \param test  The test to check.
   \param msg   Error message to print upon failing of test.
*/
#define my_assert( file, line, test, msg )  \
	do { \
	if( !(test) ){ \
		std::cerr << file << " ( " << line << " ): Assertion failed: " \
		          << msg << "\n"; \
		std::terminate(); \
	} \
	} while(0)
#endif

/**
   \brief Prints useful error message and then terminates.
   
   \param file  File name (use __FILE__ macro)
   \param line  Line number (use __LINE__ macro(
   \param msg   Error message to print before termination.
*/
#define my_logic_error( file, line, msg )	  \
	do { \
		std::cerr << file << " ( " << line << " ): Logic error: " \
		          << msg << "\n"; \
		std::terminate(); \
	}while(0)

/**
   \brief Prints useful error message and then terminates.
   
   \param file  File name (use __FILE__ macro)
   \param line  Line number (use __LINE__ macro(
   \param msg   Error message to print before termination.
*/
#define my_runtime_error( file, line, msg ) \
	do { \
		std::cerr << file << " ( " << line << " ): Runtime error: " \
		          << msg << "\n"; \
		std::terminate(); \
	}while(0)


} // namespace lammps_tools



#endif // MY_ASSERT_HPP
