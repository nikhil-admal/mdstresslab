/*
 * SubConfiguration.cpp
 *
 *  Created on: Dec 16, 2019
 *      Author: Nikhil
 */

#include "SubConfiguration.h"

SubConfiguration::SubConfiguration(const Stencil& stencil) :
	Configuration(stencil.particleContributingMap.size(), stencil.parent.coordinates.at(Reference).rows()>0),
	parent(stencil.parent)
{
	int numberOfParticlesInParent= parent.numberOfParticles;

	particleContributing.resize(numberOfParticles);
	particleContributing.setZero();


	int referenceAndFinal= (coordinates.at(Reference).rows()>0);

	int i_localParticle=0;
	for(const auto& particleContributingPair : stencil.particleContributingMap)
	{
		int particle= particleContributingPair.first;
		int contributing= particleContributingPair.second;

		auto result= globalLocalMap.insert({particle,i_localParticle});
		assert(result.second && "Insertion failed in the creation of a glocalLocalMap");
        result= localGlobalMap.insert({i_localParticle,particle});
        assert(result.second && "Insertion failed in the creation of a localGlobalMap");
		particleContributing[i_localParticle]= contributing;

		this->coordinates.at(Current).row(i_localParticle)= parent.coordinates.at(Current).row(particle);
		if (referenceAndFinal) this->coordinates.at(Reference).row(i_localParticle)= parent.coordinates.at(Reference).row(particle);
		this->velocities.row(i_localParticle)= parent.velocities.row(particle);
		this->species.push_back(parent.species[particle]);
		i_localParticle++;
	}

}

SubConfiguration::~SubConfiguration() {
	// TODO Auto-generated destructor stub
}

