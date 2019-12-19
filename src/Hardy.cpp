/*
 * Hardy.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#include "Hardy.h"
#include <math.h>
#include <iostream>
#include <limits>
#include "typedef.h"


Hardy::Hardy(double averagingDomainSize):averagingDomainSize(averagingDomainSize),normalizer(8/(5*M_PI*pow(averagingDomainSize,3)))
{ }
Hardy::Hardy(const Hardy& _hardy)
{
	*this= _hardy;
}

Hardy::~Hardy() {
	// TODO Auto-generated destructor stub
}

// weighting function:

//w(r)= c if r<= R/2,
//	= 2c(1-r/R) if R/2<=r<=R
//where c= 8/(5\pi R^3)

double Hardy::integrateConstant(const double& a,
							    const double& b,
							    const double& r1,
								const double& r2)
{
	assert(r1<=r2);
	assert (a*pow(r2,2)+b+epsilon>0 && a*pow(r1,2)+b+epsilon>0);
	return (std::sqrt(a*pow(r2,2)+b+epsilon) - std::sqrt(a*pow(r1,2)+b+epsilon))/a;
}

double Hardy::integrateLinear(const double& a,
							  const double& b,
							  const double& r1,
							  const double& r2)
{
	assert(r1<=r2);
	assert(a*r1*r1+b+epsilon>0 && a*r2*r2+b+epsilon>0 && a>0);
	double temp1= std::sqrt(a)*std::sqrt(a*r1*r1 + b + epsilon);
	double temp2= std::sqrt(a)*std::sqrt(a*r2*r2 + b + epsilon);

	return ((r2*temp2-b*log(temp2+a*r2))-
	   	    (r1*temp1-b*log(temp1+a*r1)))/(2*pow(a,1.5));
}
double Hardy::operator()(const Vector3d& vec1, const Vector3d& vec2)
{
	double r1= vec1.norm();
	double r2= vec2.norm();

	double a= 4*(vec1-vec2).squaredNorm();
	double b= 4*pow(vec2.dot(vec1)-vec1.dot(vec1),2)-a*r1*r1;

	double sPerp= (vec1.dot(vec1)-vec2.dot(vec1))/(vec1-vec2).squaredNorm();

	Vector3d vecPerp;
	vecPerp= (1-sPerp)*vec1+sPerp*vec2;
	double rPerp= vecPerp.norm();
	assert(rPerp<std::min(r1,r2)+epsilon);

	if (rPerp>=averagingDomainSize) return 0;

	r1= std::min(r1,averagingDomainSize);
	r2= std::min(r2,averagingDomainSize);

	double rmin;
	double rmax;
	double result= 0;
	if (sPerp<0 || sPerp>1)
	{
		rmin= std::min(r1,r2);
		rmax= std::max(r1,r2);

		if (abs(rmax-rmin)>epsilon) return integrate(a,b,rmin,rmax);
		else return 0;
	}
	else
	{
		rmin= rPerp;
		if (abs(r1-rmin)>epsilon) result= integrate(a,b,rmin,r1);
		if (abs(r2-rmin)>epsilon) result= result+integrate(a,b,rmin,r2);

		return result;
	}
	assert(0);
}

double Hardy::integrate(const double& a,
					    const double& b,
					    const double& rmin,
					    const double& rmax)
{
	assert(rmin<=rmax);
	if (rmin >= averagingDomainSize/2)
		return 2*normalizer*2*(integrateConstant(a,b,rmin,rmax) -
							 integrateLinear(a,b,rmin,rmax)/averagingDomainSize);

	else if (rmax > averagingDomainSize/2)
		return 2*normalizer*(integrateConstant(a,b,rmin,averagingDomainSize/2) +
						   2*integrateConstant(a,b,averagingDomainSize/2,rmax) - 2*integrateLinear(a,b,averagingDomainSize/2,rmax)/averagingDomainSize);
	else
		return 2*normalizer*integrateConstant(a,b,rmin,rmax);

	assert(0);
}
