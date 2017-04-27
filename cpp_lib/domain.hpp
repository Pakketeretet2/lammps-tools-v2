#ifndef DOMAIN_HPP
#define DOMAIN_HPP

struct domain {

	enum periodic_bits { BIT_X = 1,
	                     BIT_Y = 2,
	                     BIT_Z = 4 };
	
	double xlo[3], xhi[3];
	int periodic;
};

#endif // DOMAIN_HPP
