#ifndef MY_ASSERT_HPP
#define MY_ASSERT_HPP

#include <iostream>

#define FILELINE __FILE__,__LINE__


#ifdef NDEBUG
#define my_assert( test, msg, file, line ) \
	do {} while(0)
#else
#define my_assert( test, msg, file, line )  \
	do { \
	if( !(test) ){ \
		std::cerr << file << " ( " << line << " ): Assertion failed: " \
		          << msg << "\n"; \
		std::terminate(); \
	} \
	} while(0)
#endif


#define logic_error( msg, file, line ) \
	do { \
		std::cerr << file << " ( " << line << " ): Logic error: " \
		          << msg << "\n"; \
	}while(0)






#endif // MY_ASSERT_HPP
