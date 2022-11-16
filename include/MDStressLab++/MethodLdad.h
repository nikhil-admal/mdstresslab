/*
 * MethodLdad.h
 *
 *  Created on: Nov 5, 2019
 *  Modified on: Aug 27, 2022
 *      Author: Nikhil
 */

#ifndef METHODLDADBASE_H_
#define METHODLDADBASE_H_

#include "Method.h"
#include "typedef.h"
#include "Trigonometric.h"
#include "Constant.h"
#include <map>

template<typename T>
class MethodLdad : public Method<MethodLdad<T>>
{
    friend class Method<MethodLdad<T>>;
public:
	MethodLdad(const Matrix3d& ldadVectors);
	MethodLdad(const MethodLdad&);
	virtual ~MethodLdad();

protected:
	double operator()(const Vector3d& vec) const;
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2) const;

private:
	double normalizer;
	T oneDFunction;
	Matrix3d ldadVectors, inverseLdadVectors;
};

/*!
 * @example TestLdadLJ.cpp
 * @example TestLdadSW.cpp
 */
#include "MethodLdad.cpp"

// instantiate MethodLdad templates
//
template class MethodLdad<Constant>;
template class MethodLdad<Trigonometric>;
typedef MethodLdad<Constant> MethodLdadConstant;
typedef MethodLdad<Trigonometric> MethodLdadTrigonometric;

#endif /* METHODLDADBASE_H_ */
