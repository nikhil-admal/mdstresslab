/*
 * Grid.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef GRID_H_
#define GRID_H_
#include <iostream>
#include "typedef.h"

class Grid {
public:
	Grid(const int&, const int&, const int&);
	virtual ~Grid();
	int ngrid;
	MatrixXd coordinates;
};

#endif /* GRID_H_ */
