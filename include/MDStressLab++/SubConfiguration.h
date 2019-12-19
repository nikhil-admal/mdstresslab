/*
 * SubConfiguration.h
 *
 *  Created on: Dec 16, 2019
 *      Author: Nikhil
 */

#ifndef SRC_SUBCONFIGURATION_H_
#define SRC_SUBCONFIGURATION_H_

#include "Configuration.h"
#include "Stencil.h"

class SubConfiguration : public Configuration
{
public:
	const Configuration& parent;
	std::map<int,int> globalLocalMap;
	VectorXi particleContributing;
	SubConfiguration(const Stencil& stencil);
	~SubConfiguration();
};

#endif /* SRC_SUBCONFIGURATION_H_ */
