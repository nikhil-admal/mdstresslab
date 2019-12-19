/*
 * Configuration.cpp
 *
 *  Created on: Nov 25, 2019
 *      Author: Nikhil
 */

#include "Configuration.h"
#include <iostream>

Configuration::Configuration(int numberOfParticles, int referenceAndFinal): numberOfParticles(numberOfParticles)
{
	assert(numberOfParticles>=0);
	// Set all particles as non-contributing by default
	coordinates[Current]= MatrixXd::Zero(numberOfParticles,DIM);
	coordinates[Reference]= MatrixXd::Zero(0,0);
	if(referenceAndFinal) coordinates[Reference].resize(numberOfParticles,DIM);
	velocities.resize(numberOfParticles,DIM);
	species.reserve(numberOfParticles);
}

Configuration* Configuration::getLocalConfiguration(const std::set<int>& localParticleList) const
{
	// Ensure the local particle list is within the range
	assert( std::all_of(localParticleList.cbegin(),
						localParticleList.cend(),
						[this](int particle){return (particle >= 0 && particle <= numberOfParticles); }) &&
						"Local particle list is out of range");

	int referenceAndFinal= (coordinates.at(Reference).rows()>0);
	int localNumberOfParticles= localParticleList.size();
	std::cout << "Creating a local configuration of size = " << localNumberOfParticles << std::endl;
	Configuration* plocalConfiguration= new Configuration{localNumberOfParticles,referenceAndFinal};

	int i_localParticle=0;
	for(const auto& localParticle : localParticleList)
	{
		plocalConfiguration->coordinates.at(Current).row(i_localParticle)= coordinates.at(Current).row(localParticle);
		if (referenceAndFinal) plocalConfiguration->coordinates.at(Reference).row(i_localParticle)= coordinates.at(Reference).row(localParticle);
		plocalConfiguration->velocities.row(i_localParticle)= velocities.row(localParticle);
		plocalConfiguration->species.push_back(species[localParticle]);
		i_localParticle++;
	}
	std::cout << "Local configuration creation successful" << std::endl;
	return plocalConfiguration;
}

Configuration::~Configuration() {
	// TODO Auto-generated destructor stub
}




