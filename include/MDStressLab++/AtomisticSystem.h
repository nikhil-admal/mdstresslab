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
	AtomisticSystem(int);
	virtual ~AtomisticSystem();
	int numberOfParticles;
	VectorXi speciesCode;
	Vector3d sampleSize;
	MatrixXd coordinates, reference_coordinates;
};

#endif /* ATOMISTICSYSTEM_H_ */
