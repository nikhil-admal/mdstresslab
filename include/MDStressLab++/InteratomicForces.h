/*
 * InteratomicForces.h
 *
 *  Created on: Nov 21, 2019
 *      Author: Nikhil
 */

#ifndef SRC_INTERATOMICFORCES_H_
#define SRC_INTERATOMICFORCES_H_
#include "neighbor_list.h"

class InteratomicForces {
public:
	InteratomicForces(NeighList*);
	virtual ~InteratomicForces();
	const NeighListOne* nlOne_ptr;
	std::vector<double> fij;
};

#endif /* SRC_INTERATOMICFORCES_H_ */
