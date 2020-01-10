/*
 * Trigonometric.h
 *
 *  Created on: Jan 9, 2020
 *      Author: Nikhil
 */

#ifndef TRIGONOMETRIC_H_
#define TRIGONOMETRIC_H_


class Trigonometric{
public:
	Trigonometric();
	virtual ~Trigonometric();

	double operator()(const double& t);
	double integrate(const double& t1, const double& t2);
};



#endif /* TRIGONOMETRIC_H_ */
