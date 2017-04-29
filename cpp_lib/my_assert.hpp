#ifndef MY_ASSERT_HPP
#define MY_ASSERT_HPP

#include <iostream>

#define FILELINE __FILE__,__LINE__


#ifdef NDEBUG
#define my_assert( file, line, test, msg )  \
	do {} while(0)
#else
#define my_assert( file, line, test, msg )  \
	do { \
	if( !(test) ){ \
		std::cerr << file << " ( " << line << " ): Assertion failed: " \
		          << msg << "\n"; \
		std::terminate(); \
	} \
	} while(0)
#endif


#define my_logic_error( file, line, msg )	  \
	do { \
		std::cerr << file << " ( " << line << " ): Logic error: " \
		          << msg << "\n"; \
	}while(0)

#define my_runtime_error( file, line, msg ) \
	do { \
		std::cerr << file << " ( " << line << " ): Runtime error: " \
		          << msg << "\n"; \
	}while(0)






#endif // MY_ASSERT_HPP
