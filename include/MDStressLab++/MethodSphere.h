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

/*!
 * @example testLJ.cpp
 * An example demonstrating the computation of Piola and Cauchy stress tensors
 * using various spherical averaging domains and grids
 *
 * -# Open the configuration file (*.data in MDStressLab format) and read the number of particles
 * @snippet{lineno} testLJ.cpp Read
 *
 * -# Initialize the BoxConfiguration and read the reference and current atomic coordinates
 *  from the configuration file. The reference configuration is an fcc Ar crystal in the relaxed state.
 *  The deformed/current configuration is the reference crystal strained in the \f$y\f$-direction.
 * @snippet{lineno} testLJ.cpp Configuration
 *
 * -# Initialize the Kim model
 * @snippet{lineno} testLJ.cpp Model
 *
 * -# Create grids
 * @snippet{lineno} testLJ.cpp Grid
 *
 * -# Construct MethodSphere objects of various averaging domain diameters
 * @snippet{lineno} testLJ.cpp Method
 *
 * -# Construct Stress objects using the MethodSphere and Grid objects
 * @snippet{lineno} testLJ.cpp Stress
 *
 * -# Calculate stresses. The following snippet shows that stresses can be calculated all at once or
 * with various combinations.
 * @snippet{lineno} testLJ.cpp Calculate
 *
 * -# Write stresses
 * @snippet{lineno} testLJ.cpp Write
 *
 * -# We compare our results with the exact results for unit testing purposes.
 * @snippet{lineno} testLJ.cpp Compare
 *
 * Full code:
 */


#endif /* METHODSPHERE_H_ */
