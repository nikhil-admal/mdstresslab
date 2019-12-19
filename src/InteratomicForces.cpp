/*
 * InteratomicForces.cpp
 *
 *  Created on: Nov 21, 2019
 *      Author: Nikhil
 */

#include "InteratomicForces.h"
#include<iostream>

// Instantiate interatomic forces from neighbor list
InteratomicForces::InteratomicForces(NeighList* _nl_ptr)
{
	int numberOfNeighborLists= _nl_ptr->numberOfNeighborLists;
	double maxCutoff= 0.0;
	int chosenListIndex= 0;
	for(int i=0; i<numberOfNeighborLists; i++)
	{
		double cutoff= _nl_ptr->lists[i].cutoff;
		if ( cutoff > maxCutoff)
		{
			maxCutoff= cutoff;
			chosenListIndex= i;
		}
	}

	int size=0;
	int numberOfParticles= _nl_ptr->lists[chosenListIndex].numberOfParticles;
	for(int i=0; i<numberOfParticles;i++) size+= _nl_ptr->lists[chosenListIndex].Nneighbors[i];
	fij.resize(size,0.0);
	nlOne_ptr= _nl_ptr->lists+chosenListIndex;
}

InteratomicForces::~InteratomicForces() {
	// TODO Auto-generated destructor stub
}

