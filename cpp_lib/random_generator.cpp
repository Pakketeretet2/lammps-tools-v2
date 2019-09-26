#include <cmath>
#include <cassert>
#include <iostream>

#include "random_generator.hpp"


// Some compilers do not define uint:
typedef unsigned int uint;


/**
   \brief Initialize the random generator.
*/
RanGen::RanGen(int seed, int ziggurat_size)
	: internal_z(ziggurat_size), rng(seed), save(false), second(0)
{
}

/**
   \brief Outout some info about the random numbers used.
*/
RanGen::~RanGen()
{ }


double RanGen::uniform()
{
	return rng.uniform();
}


uint RanGen::random_int()
{
	return rng.random_int();
}


double RanGen::gaussian()
{
	return gaussian_polar();
}


double RanGen::gaussian_polar()
{
	double first;
	if (!save) {
		double v1,v2,rsq,fac;
		do {
			v1 = 2.0*uniform()-1.0;
			v2 = 2.0*uniform()-1.0;
			rsq = v1*v1 + v2*v2;
		} while ((rsq >= 1.0) || (rsq == 0.0));
		fac = sqrt(-2.0*log(rsq)/rsq);
		second = v1*fac;
		first = v2*fac;
		save = true;
	} else {
		first = second;
		save = false;
	}
	return first;
}



double RanGen::exponential()
{
	return exponential_inversion();
}



double RanGen::gaussian_ziggurat()
{
	internal_z.n_gauss++;
	// The ziggurat draws only from the positive size so
	// multiply by a random sign.
	double g = draw_from_table( GAUSS );

	// Branch:
	return g * random_sign();
}


int RanGen::random_sign()
{
	// It seems like such a waste to discard a whole random int
	// just to check one bit...
	return -1 + 2*( random_int() & 0x00000000000000000000000000000001 );
}


double RanGen::exponential_ziggurat()
{
	internal_z.n_exp++;
	return draw_from_table( EXPONENTIAL );
}


// pdf must match table!
double RanGen::draw_from_table( int distribution )
{
	int layer = internal_z.N*uniform();
	double *table_x = nullptr, *table_y = nullptr;

	switch(distribution){
		case GAUSS:
			table_x = internal_z.xi_gauss;
			table_y = internal_z.yi_gauss;
			break;
		case EXPONENTIAL:
			table_x = internal_z.xi_exponential;
			table_y = internal_z.yi_exponential;
			break;
		default:
			// This should error.
			break;
	}

	// Branch:
	if( layer == 0 ){
		// Fallback for the tail (depends on distribution!)
		switch(distribution){
			case GAUSS:
				internal_z.n_fallbacks_gauss++;
				return table_x[0] + gaussian_polar();
			case EXPONENTIAL:
				internal_z.n_fallbacks_exp++;
				return table_x[0] + exponential_inversion();
			default:
				break;
		}
	}

	double x = uniform()*table_x[layer];
	// Branch:
	if( (layer < internal_z.N - 1) && (x < table_x[layer+1]) ){
		return x;
	}

	double y = table_y[layer] + uniform()*( table_y[layer+1] - table_y[layer] );
	double pdf_x = 0;
	switch(distribution){
		case GAUSS:
			pdf_x = exp(-0.5*x*x);
			break;
		case EXPONENTIAL:
			pdf_x = exp(-x);
			break;
		default:
			break;
	}
	// Branch:
	if( y < pdf_x ){
		return x;
	}

	return draw_from_table( distribution );
}



double RanGen::exponential_inversion()
{
	return -log( uniform() );
}


/**
   \brief Initialize ziggurat struct of N layers.

   N must be 32, 64, 128, 256 or 512.
*/
RanGen::ziggurat::ziggurat( int N )
	: N(N), xi_gauss(NULL), yi_gauss(NULL),
	  xi_exponential(NULL), yi_exponential(NULL),
	  n_fallbacks_exp(0), n_fallbacks_gauss(0)
{
	xi_gauss = new double[N];
	xi_exponential = new double[N];

	yi_gauss = new double[N];
	yi_exponential = new double[N];



	if( ! (xi_gauss && xi_exponential &&
	       yi_gauss && yi_exponential) ){
		// Error.
	}

	// These numbers for the first layer are pre-computed.
	double x0_gauss = 3.654152885360915;
	double x0_exponential = 7.697117470131029;

	switch( N ){
		case 32:
			x0_exponential = 5.220000000000102;
			x0_gauss = 2.961300121264017;
			break;
		case 64:
			x0_exponential = 6.070000000000084;
			x0_gauss = 3.213657627158888;
			break;
		case 128:
			x0_exponential = 6.898315116615643;
			x0_gauss = 3.442619855896625;
			break;
		case 256:
			x0_exponential = 7.697117470131029;
			x0_gauss = 3.654152885360915;
			break;
		case 512:
			x0_exponential = 8.481739963222662;
			x0_gauss = 3.852046150367979;
			break;
		default:
			// Not supported!
			std::cerr << "N = " << N
			          << " not supported by ziggurat!\n";
	}

	generate_table_exponential( xi_exponential, yi_exponential,
	                            x0_exponential );
	generate_table_gauss( xi_gauss, yi_gauss, x0_gauss );
}


void RanGen::ziggurat::generate_table_exponential( double *xi, double *yi,
                                                    double x0 )
{
	xi[0] = x0;
	auto pdf = []( double x ) { return exp(-x); };
	auto cdf = []( double x ) { return 1 - exp(-x); };
	auto ipdf = []( double x ) { return -log(x); };

	double f_im = pdf(x0);
	yi[0] = f_im;
	double V0 = x0*f_im + (1.0 - cdf(x0));

	for( int i = 1; i < N; ++i ){
		xi[i] = ipdf( f_im + V0 / xi[i-1] );
		double f_i = pdf(xi[i]);
		yi[i] = f_i;
		V0 = xi[i-1]*(f_i - f_im);
		f_im = f_i;
	}

	// Just set the final x to 0?
	xi[N-1] = 0.0;
}


void RanGen::ziggurat::generate_table_gauss( double *xi, double *yi,
                                              double x0 )
{
	xi[0] = x0;

	double pi_const = pi;

	// Note that these lambdas here _need_ to be normalized for
	// the equations to work, even though the table makes no use of them!

	auto pdf  = [pi_const]( double x )
		{ return (x>=0)*sqrt(2.0/pi_const)*exp(-0.5*x*x); };
	auto cdf  = []( double x )
		{ return erf(x/sqrt(2.0)); };
	auto ipdf = [pi_const]( double x )
		{ return sqrt( -2 * log(x*sqrt(pi_const/2.0)) ); };

	double f_im = pdf(x0);
	yi[0] = f_im;
	double V0 = x0*f_im + (1.0 - cdf(x0));
	for( int i = 1; i < N; ++i ){
		xi[i] = ipdf( f_im + V0 / xi[i-1] );
		double f_i = pdf(xi[i]);
		yi[i] = f_i;
		V0 = xi[i-1]*(f_i - f_im);
		f_im = f_i;
	}

	// Just set the final x to 0?
	xi[N-1] = 0.0;
}



RanGen::ziggurat::~ziggurat()
{
	if( xi_gauss ) delete [] xi_gauss;
	if( yi_gauss ) delete [] yi_gauss;

	if( xi_exponential ) delete [] xi_exponential;
	if( yi_exponential ) delete [] yi_exponential;
}


/*** Other RNG stuff: ***/
RanMT::RanMT( int seed ) : mt_rand(seed), u_dist(0,1)
{}


RanMT::~RanMT()
{}


uint RanMT::random_int()
{
	auto max_int = mt_rand.max();
	std::uniform_int_distribution<uint> ui(0, max_int);
	return ui(mt_rand);
}

double RanMT::uniform()
{
	return u_dist(mt_rand);
}
