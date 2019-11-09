/*
 * Stress.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef STRESS_H_
#define STRESS_H_

#include "Grid.h"
#include "typedef.h"
#include <functional>

class Stress {
public:
	Stress(const bool&, const bool&);
	Stress();
	virtual ~Stress();
	bool requestPiola,requestCauchy;
	MatrixXd* gridCoordinates;
	std::shared_ptr<Matrix3d[]> field;
	virtual double bondFunction(const std::array<double,DIM>&, const std::array<double,DIM>&)=0;
};

#endif /* STRESS_H_ */
