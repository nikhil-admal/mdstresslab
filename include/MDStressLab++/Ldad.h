/*
 * Ldad.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef LDADBASE_H_
#define LDADBASE_H_

#include "typedef.h"
#include <map>

template<typename T>
class Ldad
{
public:
	Ldad(const Matrix3d& ldadVectors);
	Ldad(const Ldad&);
	virtual ~Ldad();

	double averagingDomainSize;
	double operator()(const Vector3d& vec);
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2);


protected:
	double normalizer;
	T oneDFunction;
	Matrix3d ldadVectors, InverseldadVectors;
};

#include "Ldad.cpp"
#endif /* LDADBASE_H_ */
