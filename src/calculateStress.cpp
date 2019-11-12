/*
 * calculateStress.cpp
 *
 *  Created on: Nov 7, 2019
 *      Author: Nikhil
 */

#include "calculateStress.h"
void calculateStress(const AtomisticSystem& body,
					 const MatrixXd& gridCoordinates,
					 std::vector<Stress*> pStress)
{
	std::array<double,DIM> vec1,vec2;
	std::cout << (*(pStress[0])).bondFunction(vec1,vec2);
}


