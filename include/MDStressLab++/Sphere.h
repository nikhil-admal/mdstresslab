/*
 * Sphere.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include "typedef.h"
#include "Polynomial.h"
#include <map>

class Sphere
{
public:
	Sphere(double);
	Sphere(const Sphere&);
	virtual ~Sphere();

	std::map<std::pair<double,double>,Polynomial> piecewisePolynomial;


	double averagingDomainSize;
	double operator()(const Vector3d& vec);
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2);



	double normalizer;
	double integratePolynomial(const int&, const double&,const double&,const double&,const double&);
	double integrate(const double&,const double&,const double&,const double&);
};

#endif /* SPHERE_H_ */
