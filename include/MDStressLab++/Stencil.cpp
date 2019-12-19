/*
 * Stencil.cpp
 *
 *  Created on: Dec 13, 2019
 *      Author: Nikhil
 */

#include "Stencil.h"
#include "SpatialHash.h"
#include "Configuration.h"

Stencil::Stencil(const Configuration& parent) : parent(parent) { }

template<ConfigType configType>
void Stencil::expandStencil(const Grid<configType>* pgrid, const double& contributingNeighborhoodSize, const double& noncontributingNeighborhoodSize)
{
		double padding= contributingNeighborhoodSize + noncontributingNeighborhoodSize;

		Vector3d origin,step;
		origin.setConstant(0.0);
		step.setConstant(padding);

		// Hash the coordinates
		ConstSpatialHash hashParticles(origin,step,parent.coordinates.at(configType));
		// Hash the grid points
		ConstSpatialHash hashGrid(origin,step,pgrid->coordinates);

		int i_gridPoint= 0;
		// for each grid point
		for(const auto& gridPoint : pgrid->coordinates)
		{
			Triplet bin= hashGrid.hashFunction(i_gridPoint);
			// for each neighboring bin of a grid point's bin
			for (const auto& neighborBin : bin.neighborList())
			{
				std::vector<int>& particleList= hashParticles.hashTable[neighborBin];
				// for each particle in a neighboring bin
				for(const auto& particle : particleList)
				{
					double distanceSquared= (parent.coordinates.at(configType).row(particle)-gridPoint).squaredNorm();

					if (distanceSquared <= pow(contributingNeighborhoodSize,2))
//						particleContributingMap.insert_or_assign({particle,1});
						particleContributingMap[particle]= 1;
					else if (distanceSquared<=pow(padding,2))
						particleContributingMap.insert({particle,0});
				}
			}
			i_gridPoint++;
		}
}

void Stencil::emptyStencil()
{
	particleContributingMap.clear();
}
Stencil::~Stencil() {
	// TODO Auto-generated destructor stub
}

