/*
 * Constant.h
 *
 *  Created on: Jan 9, 2020
 *      Author: Nikhil
 */

#ifndef CONSTANT_H_
#define CONSTANT_H_

class Constant {
public:
	Constant();
	virtual ~Constant();
	double operator()(const double& t);
	double integrate(const double& t1, const double& t2);
};

#endif /* CONSTANT_H_ */
