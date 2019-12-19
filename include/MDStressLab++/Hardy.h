/*
 * Hardy.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef HARDY_H_
#define HARDY_H_

#include "typedef.h"

class Hardy
{
public:
	Hardy(double);
	Hardy(const Hardy&);
	virtual ~Hardy();

	double averagingDomainSize;
	double operator()(const Vector3d& vec1, const Vector3d& array2);



	double normalizer;
	double integrateConstant(const double&,const double&,const double&,const double&);
	double integrateLinear(const double&,const double&,const double&,const double&);
	double integrate(const double&,const double&,const double&,const double&);
};

#endif /* HARDY_H_ */
