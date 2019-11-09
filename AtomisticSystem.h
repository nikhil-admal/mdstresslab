/*
 * AtomisticSystem.h
 *
 *  Created on: Nov 4, 2019
 *      Author: Nikhil
 */

#ifndef ATOMISTICSYSTEM_H_
#define ATOMISTICSYSTEM_H_

#include <iostream>
#include "typedef.h"

class AtomisticSystem {
public:
	AtomisticSystem();
	virtual ~AtomisticSystem();
	int numberOfParticles;
	VectorXi speciesCode;
	MatrixXd coordinates;
};

#endif /* ATOMISTICSYSTEM_H_ */
