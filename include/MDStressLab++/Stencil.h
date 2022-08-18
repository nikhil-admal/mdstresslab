/*
 * Stencil.h
 *
 *  Created on: Dec 13, 2019
 *      Author: Nikhil
 */

#ifndef SRC_STENCIL_H_
#define SRC_STENCIL_H_

#include "SpatialHash.h"
#include "Configuration.h"
#include "typedef.h"
#include <map>

class Configuration;

template<ConfigType T>
class Grid;
class GridBase;

class Stencil {
public:
	const Configuration& parent;
	std::map<int,int> particleContributingMap;

	//template<ConfigType configType>
	//void expandStencil(const Grid<configType>* pgrid, const double&, const double&);

    template<ConfigType configType>
    void expandStencil(const Grid<configType>* pgrid, const double& contributingNeighborhoodSize, const double& noncontributingNeighborhoodSize)
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
						//particleContributingMap.insert_or_assign({particle,1});
						particleContributingMap[particle]= 1;
					else if (distanceSquared<=pow(padding,2))
						particleContributingMap.insert({particle,0});
				}
			}
			i_gridPoint++;
		}
    }


	void emptyStencil();
	Stencil(const Configuration&);
	virtual ~Stencil();

};

//#include "Stencil.cpp"

#endif /* SRC_STENCIL_H_ */
