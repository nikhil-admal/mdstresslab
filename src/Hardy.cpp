/*
 * Hardy.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#include "Hardy.h"

Hardy::Hardy(const bool& reqPiola, const bool& reqCauchy) : Stress(reqPiola,reqCauchy){
	// TODO Auto-generated constructor stub

}

Hardy::~Hardy() {
	// TODO Auto-generated destructor stub
}

double Hardy::bondFunction(const std::array<double,DIM>& vec1, const std::array<double,DIM>& vec2)
{
	return 0.0;
}
