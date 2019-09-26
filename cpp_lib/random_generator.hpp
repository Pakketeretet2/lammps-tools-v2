/**
   \file random_generator.h

   Defines some classes for pseudo-random number generation.
*/

#ifndef RANDOM_GENERATOR_HPP
#define RANDOM_GENERATOR_HPP

#include <cmath>
#include <random>


/**
   \brief Class that wraps around the C++11 Mersenne twister.
*/
class RanMT {
public:
	RanMT( int );
	~RanMT();

	double uniform();
	unsigned int random_int();
private:
	std::mt19937_64 mt_rand;
	std::uniform_real_distribution<double> u_dist;
};


/**
   \brief This class provides functionality for drawing random numbers.
*/
class RanGen {
public:
	/// Constructor, takes random seed and ziggurat size.
	RanGen(int, int = 256);
	~RanGen();

	/// Random unsigned int
	uint random_int();

	/// A random sign (i.e., either +1 or -1 with 50% chance)
	int random_sign();

	/// Unform [0,1) interval
	double uniform();
	/// Gaussian with mean 0 and variance 1
	double gaussian();
	/// Exponential with mean = variance = 1
	double exponential();

	/// Gaussian through the Marsaglia polar method
	double gaussian_polar();
	/// Gaussian from the internal ziggurat
	double gaussian_ziggurat();
	/// Exponential through inversion (i.e., -log(uniform()))
	double exponential_inversion();
	/// Exponential from the internal ziggurat.
	double exponential_ziggurat();

	/// Supported ziggurat distributions:
	enum distributions { GAUSS, EXPONENTIAL };

private:

	/// Helper function to get a random number from a table.
	double draw_from_table( int distribution );

	/// Internal ziggurat structure.
	struct ziggurat {
		ziggurat( int N );
		~ziggurat();

		const double pi    = 3.14159265358979323846;
		const double sq2   = 1.41421356237309504880;
		const double sqhpi = sqrt( 0.5 * pi );

		int N;
		double *xi_gauss;
		double *yi_gauss;

		double *xi_exponential;
		double *yi_exponential;

		long int n_fallbacks_exp;
		long int n_fallbacks_gauss;

		long int n_exp;
		long int n_gauss;

		/// Generates exponential table
		void generate_table_exponential( double *table_x, double *table_y,
		                                 double x0 );
		/// Generates gaussian table
		void generate_table_gauss( double *table_x, double *table_y, double x0 );

	};


	ziggurat internal_z;
	RanMT rng;

	bool save;
	double second;
};


#endif // RANDOM_GENERATOR_HPP
