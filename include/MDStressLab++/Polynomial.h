/*
 * Polynomial.h
 *
 *  Created on: Dec 21, 2019
 *      Author: Nikhil
 */

#ifndef SRC_POLYNOMIAL_H_
#define SRC_POLYNOMIAL_H_

#include <vector>
#include <cmath>

class Polynomial {
public:
	Polynomial(const std::vector<double>& coefficients):coefficients(coefficients)
	{}
	Polynomial()=default;
	virtual ~Polynomial(){}
	std::vector<double> coefficients;

	double operator()(const double& argument) const
	{
		int degree= 0;
		double value= 0;
		for(const auto& coefficient : coefficients)
		{
			value= value + coefficient*pow(argument,degree);
			degree++;
		}

		return value;
	}
};

#endif /* SRC_POLYNOMIAL_H_ */
