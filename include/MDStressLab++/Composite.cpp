/*
 * Composite.cpp
 *
 *  Created on: Dec 2, 2019
 *      Author: Nikhil
 */

#include <iostream>
#include "Composite.h"
#include "SpatialHash.h"
#include <map>

// Typename T refers to 'Reference' / 'Current'

template<ConfigType T>
Composite<T>::Composite(const Configuration& config,
					    const std::map<const Grid<T>*, double> gridNeighborhoodSizeMap,
					    const double& noncontributingNeighborhoodSize)
{
//	localConfigList: A list of global particle numbers representing a local configuration.
//	localContributingConfigList: A subset of localConfigList representing contributing partilces.
	std::set<int> localConfigList, localContributingConfigList;

//	NeighborLists: An array of sets, where each set describes the contributing neighbor
//					   list of global particles
	std::map<const Grid<T>*,std::vector<std::set<int>>> gridNeighborListsMap;

//	Loop over each grid
	for(const auto& gridNeighborhoodSize : gridNeighborhoodSizeMap)
	{
		double contributingNeighborhoodSize= gridNeighborhoodSize.second;
		const Grid<T>* pgrid= gridNeighborhoodSize.first;
		double padding= contributingNeighborhoodSize + noncontributingNeighborhoodSize;

		Vector3d origin,step;
		origin.setConstant(0.0);
		step.setConstant(padding);

		// Hash the coordinates
		ConstSpatialHash hashParticles(origin,step,config.coordinates.at(T));
		// Hash the grid points
		ConstSpatialHash hashGrid(origin,step,pgrid->coordinates);



		int i_gridPoint= 0;
		// for each grid point
		for(const auto& gridPoint : pgrid->coordinates)
		{

//		Variables starting with 'grid' correspond to one grid point
//		gridPaddingList: A set of particles in the neighborhood of radius=
//						 averagingDomainSize+2*influence_distance of a grid point
//		gridContributingList: A set of particles in the neighborhood of radius=
//						 	  averagingDomainSize+influence_distance of a grid point

//		By construction gridContributingList \subset gridPaddingList

			std::set<int> gridPaddingList,gridContributingList;
			Triplet bin= hashGrid.hashFunction(i_gridPoint);
			Triplet neighborBin;
			// for each neighboring bin of a grid point's bin
			for (const auto& neighborBin : bin.neighborList(bin))
			{
				std::vector<int>& particleList= hashParticles.hashTable[neighborBin];
				// for each particle in a neighboring bin
				for(const auto& particle : particleList)
				{
					double distance= (config.coordinates.at(T).row(particle)-gridPoint).squaredNorm();
					if (distance<=pow(contributingNeighborhoodSize,2))
					{
						gridContributingList.insert(particle);
					}
					if (distance<=pow(padding,2)) gridPaddingList.insert(particle);
				}
			}
			localConfigList.insert(gridPaddingList.begin(),gridPaddingList.end());
			localContributingConfigList.insert(gridContributingList.begin(),gridContributingList.end());

			gridNeighborListsMap[pgrid].push_back(gridContributingList);

			i_gridPoint++;
		}
	}

	plocalConfiguration.reset(config.getLocalConfiguration(localConfigList));
	particleContributing.resize(localConfigList.size());
	particleContributing.setZero();

//	Construct a global to local map, and assign particleContributing using localContributingConfigList
	std::map<int,int> globalLocalMap;
	int i_localParticle= 0;
	for(const auto& particle : localConfigList)
	{
		globalLocalMap[particle]= i_localParticle;
		if (localContributingConfigList.find(particle) != localContributingConfigList.end())
			particleContributing[i_localParticle] = 1;
		i_localParticle++;
	}

//	for each (grid,NeighborLists) pair
	for(const auto& gridNeighborLists : gridNeighborListsMap)
	{
		const Grid<T>* pgrid= gridNeighborLists.first;
//		for each list in the NeighborLists
		for (const auto& gridContributingList : gridNeighborLists.second)
		{
			std::set<int> gridLocalContributingList;
//			for each particle in the list
			for (const auto& particle : gridContributingList)
			{
				auto it= globalLocalMap.find(particle);
				assert(it != globalLocalMap.end() && "The contributing list is not a subset of the local configuration");
				gridLocalContributingList.insert(it->second);
			}
			gridLocalNeighborListsMap[pgrid].push_back(gridLocalContributingList);
		}
	}
}


template<ConfigType T>
Composite<T>::~Composite() {
	// TODO Auto-generated destructor stub
}

