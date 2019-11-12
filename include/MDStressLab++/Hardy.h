/*
 * Hardy.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef HARDY_H_
#define HARDY_H_

#include <array>
#include "Stress.h"
#include "typedef.h"

class Hardy : public Stress{
public:
	Hardy(const bool&,const bool&);
	virtual ~Hardy();
	double avgDomainRadius;
	virtual double bondFunction(const std::array<double,DIM>& , const std::array<double,DIM>& );
};

#endif /* HARDY_H_ */
