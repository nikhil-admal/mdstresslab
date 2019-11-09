/*
 * Grid.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#include "Grid.h"

Grid::Grid(const int& nx, const int& ny, const int& nz) {
	this->ngrid = nx*ny*nz;
	this->coordinates.resize(ngrid,3);
}

Grid::~Grid() {
	// TODO Auto-generated destructor stub
}

