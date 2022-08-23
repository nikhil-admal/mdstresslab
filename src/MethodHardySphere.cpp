#include "MethodHardySphere.h"
#include <math.h>

// weighting function:

//w(r)= c if r<= R/2,
//	= 2c(1-r/R) if R/2<=r<=R
//where c= 8/(5\pi R^3)
//
MethodHardySphere::MethodHardySphere(double averagingDomainSize) : Sphere(averagingDomainSize)
{
    normalizer= 8/(5*M_PI*pow(averagingDomainSize,3));
	// constant polynomial
	std::vector<double> constant{normalizer};
	Polynomial constantPolynomial{constant};
	//linear polynomial
	std::vector<double> linear{2*normalizer,-2*normalizer/averagingDomainSize};
	Polynomial linearPolynomial{linear};

	std::pair<double,double> interval;

	// construct piecewise polynomial
	interval.first= 0;
	interval.second= averagingDomainSize/2;
	piecewisePolynomial[interval]= constantPolynomial;

	interval.first= averagingDomainSize/2;
	interval.second= averagingDomainSize;
	piecewisePolynomial[interval]= linearPolynomial;
}
MethodHardySphere::MethodHardySphere(const MethodHardySphere& _hardy) : Sphere(_hardy)
{
	*this= _hardy;
}

MethodHardySphere::~MethodHardySphere() {
	// TODO Auto-generated destructor stub
}



