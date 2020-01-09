/*
 * Ldad.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#include <math.h>
#include <Ldad.h>
#include <iostream>
#include <limits>
#include "typedef.h"


Ldad::Ldad(const Matrix3d& ldadVectors):ldadVectors(ldadVectors)
{
	// initialize the normalizer here
	normalizer= 0;

	// initialize the averagingDomainSize
	averagingDomainSize= 0;
}

Ldad::Ldad(const Ldad& _Ldad)
{
	*this= _Ldad;
}

Ldad::~Ldad() {
	// TODO Auto-generated destructor stub
}

double Ldad::operator()(const Vector3d& vec)
{
	// computes the weighting function
	return 0;
}

double Ldad::bondFunction(const Vector3d& vec1, const Vector3d& vec2)
{
	// use integrate(t1,t2) -> helper function;
	return 0;
}

double Ldad::integrate(const double& t1, const double& t2)
{
	return 0;
}
