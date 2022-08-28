/*
 * MethodUser.h
 *
 *  Created on: Aug 26, 2022
 *      Author: Nikhil
 */

#ifndef METHODUSER_H_
#define METHODUSER_H_

#include "Method.h"
#include "typedef.h"

class MethodUser : public Method<MethodUser>
{
public:
	MethodUser(double averagingDomainSize):Method<MethodUser>(averagingDomainSize)
    { }
	virtual ~MethodUser()=default;

	virtual double operator()(const Vector3d& vec) const= 0;
	virtual double bondFunction(const Vector3d& vec1, const Vector3d& vec2) const= 0;
};

#endif /* METHODUSER_H_ */
