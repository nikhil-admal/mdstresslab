/*
 * MethodSphere.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef METHODSPHERE_H_
#define METHODSPHERE_H_

#include "Method.h"
#include "typedef.h"
#include "Polynomial.h"
#include <map>

class MethodSphere : public Method<MethodSphere>
{          
    friend class Method<MethodSphere>;

public:
	MethodSphere(double, std::string);
	MethodSphere(double, std::map<double,double>);
	MethodSphere(const MethodSphere&);
	virtual ~MethodSphere();

protected:
	double operator()(const Vector3d& vec) const;
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2) const;


private:
	std::map<std::pair<double,double>,Polynomial> piecewisePolynomial;
	double normalizer;
	double integratePolynomial(const int&, const double&,const double&,const double&,const double&) const;
	double integrate(const double&,const double&,const double&,const double&) const;
};

#endif /* METHODSPHERE_H_ */
