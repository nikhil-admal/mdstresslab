/*
 * AtomisticSystem.cpp
 *
 *  Created on: Nov 4, 2019
 *      Author: Nikhil
 */

#include "AtomisticSystem.h"
#include "Exception.h"

AtomisticSystem::AtomisticSystem(int _numberOfParticles):numberOfParticles(_numberOfParticles)
{
	if (numberOfParticles < 0) throw Exception("Error: Negative number of particles.\n");
	speciesCode.resize(numberOfParticles);
	coordinates.resize(numberOfParticles,3);
	reference_coordinates.resize(numberOfParticles,3);
}

AtomisticSystem::~AtomisticSystem() {
	// TODO Auto-generated destructor stub
}

