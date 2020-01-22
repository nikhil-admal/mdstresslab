/*
 * Trigonometric.h
 *
 *  Created on: Jan 9, 2020
 *      Author: Nikhil
 */

#ifndef TRIGONOMETRIC_H_
#define TRIGONOMETRIC_H_

#include "typedef.h"

class Trigonometric{
public:
	Trigonometric();
	virtual ~Trigonometric();

	double operator()(const double& t);
	double integrate(const Vector3d& vec1_pull_seg, const Vector3d& vec2_pull_seg);
};



#endif /* TRIGONOMETRIC_H_ */
