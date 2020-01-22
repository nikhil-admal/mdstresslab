/*
 * Constant.h
 *
 *  Created on: Jan 9, 2020
 *      Author: Nikhil
 */

#ifndef CONSTANT_H_
#define CONSTANT_H_

#include "typedef.h"

class Constant {
public:
	Constant();
	virtual ~Constant();
	double operator()(const double& t);
	double integrate(const Vector3d& vec1_pull_seg, const Vector3d& vec2_pull_seg);
};

#endif /* CONSTANT_H_ */
