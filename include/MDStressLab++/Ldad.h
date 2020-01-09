/*
 * Ldad.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef LDAD_H_
#define LDAD_H_

#include "typedef.h"
#include <map>

class Ldad
{
public:
	Ldad(const Matrix3d& ldadVectors);
	Ldad(const Ldad&);
	virtual ~Ldad();

	double averagingDomainSize;
	double operator()(const Vector3d& vec);
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2);


private:
	double normalizer;
	Matrix3d ldadVectors;
};

#endif /* LDAD_H_ */
