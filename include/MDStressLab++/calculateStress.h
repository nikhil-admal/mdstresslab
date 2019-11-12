/*
 * calculateStress.h
 *
 *  Created on: Nov 7, 2019
 *      Author: Nikhil
 */

#ifndef CALCULATESTRESS_H_
#define CALCULATESTRESS_H_

#include <vector>
#include "AtomisticSystem.h"
#include "Hardy.h"
#include "typedef.h"

void calculateStress(const AtomisticSystem&,const MatrixXd& gridCoordinates,std::vector<Stress*>);

#endif /* CALCULATESTRESS_H_ */
