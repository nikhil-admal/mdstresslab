/*
 * Configuration.h
 *
 *  Created on: Nov 25, 2019
 *      Author: Nikhil
 */

#ifndef INCLUDE_MDSTRESSLAB___CONFIGURATION_H_
#define INCLUDE_MDSTRESSLAB___CONFIGURATION_H_

#include <vector>
#include <string>
#include "typedef.h"
#include <iostream>
#include <set>
#include <map>

/*!
 * Configuration class describes the properties (position, velocity, species) of atoms.
 * It is the base class of the BoxConfiguration class.
 */
class Configuration{
public:
	Configuration(int,int);
	virtual ~Configuration();
	int numberOfParticles;
	std::vector<std::string> species;
	std::map<ConfigType,MatrixXd> coordinates;
	MatrixXd velocities;
	Configuration* getLocalConfiguration(const std::set<int>& localParticleList) const;

protected:
	void allocate(int _numberOfParticles, int _referenceAndFinal);
};


#endif /* INCLUDE_MDSTRESSLAB___CONFIGURATION_H_ */
