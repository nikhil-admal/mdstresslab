/*
 * Stress.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#include "Stress.h"

Stress::Stress(const bool& reqPiola, const bool& reqCauchy):gridCoordinates(nullptr){
	this->requestPiola= reqPiola;
	this->requestCauchy= reqCauchy;
}

Stress::Stress():gridCoordinates(nullptr),requestPiola(false),requestCauchy(false)
{
}

Stress::~Stress() {
	// TODO Auto-generated destructor stub
}

