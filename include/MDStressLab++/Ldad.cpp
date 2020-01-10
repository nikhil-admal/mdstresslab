/*
 * Ldad.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#include <Ldad.h>
#include <math.h>
#include <iostream>
#include <limits>
#include "typedef.h"


template<typename T>
Ldad<T>::Ldad(const Matrix3d& ldadVectors):ldadVectors(ldadVectors)
{
	// initialize the normalizer here using oneDFunction.integrate(-1,1)
	normalizer= 0;

	// initialize the averagingDomainSize
	averagingDomainSize= 0;
}

template<typename T>
Ldad<T>::Ldad(const Ldad& _Ldad)
{
	*this= _Ldad;
}

template<typename T>
Ldad<T>::~Ldad() {
	// TODO Auto-generated destructor stub
}

template<typename T>
double Ldad<T>::operator()(const Vector3d& vec)
{
	double t1,t2,t3;
	// compute t1, t2, t3
	return oneDFunction(t1)*oneDFunction(t2)*oneDFunction(t3);
}

template<typename T>
double Ldad<T>::bondFunction(const Vector3d& vec1, const Vector3d& vec2)
{
	// use oneDFunction.integrate(t1,t2) -> helper function;
	return 0;
}

