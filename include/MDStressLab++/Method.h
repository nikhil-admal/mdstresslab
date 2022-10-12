/*
 * Method.h
 *
 *  Created on: Aug 26, 2022
 *      Author: Nikhil
 */

#ifndef METHOD_H_
#define METHOD_H_

#include "typedef.h"

template<typename TMethod>
class Method
{
public:
	virtual ~Method();

	double operator()(const Vector3d& vec) const;
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2) const;
    double getAveragingDomainSize() const;

protected:
	double averagingDomainSize;
private:
    Method();
    Method(double);
    Method(const Method<TMethod>& method);
    friend TMethod;
};

#include "Method.cpp"
#endif /* METHOD_H_ */
