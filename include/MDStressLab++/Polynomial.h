/*
 * Polynomial.h
 *
 *  Created on: Dec 21, 2019
 *      Author: Nikhil
 */

#ifndef SRC_POLYNOMIAL_H_
#define SRC_POLYNOMIAL_H_

#include <deque>
#include <cmath>

class Polynomial {
public:
	Polynomial(const std::deque<double>& coefficients):coefficients(coefficients)
	{}
	Polynomial()=default;
	virtual ~Polynomial(){}
	std::deque<double> coefficients;

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

    Polynomial integrate() const
    {
        Polynomial p;
		int degree= 0;
        p.coefficients.push_back(0);
		for(const auto& coefficient : coefficients)
		{
            p.coefficients.push_back(coefficient/(degree+1));
			degree++;
		}
        return p;
    }
};

#endif /* SRC_POLYNOMIAL_H_ */
