/*
 * main.cpp
 *
 *  Created on: Nov 3, 2019
 *      Author: Nikhil Admal
 */
#include <iostream>
#include <vector>
#include "AtomisticSystem.h"
#include "Hardy.h"
#include "Grid.h"
#include "typedef.h"
#include "calculateStress.h"

enum StressDefinitions
{
	HARDY=  1<<0,
	VIRIAL= 1<<2,
	TSAI=	1<<3
};


int main()
{
	int stressFlag= {HARDY|VIRIAL};

//	-------------------------------------------------------------------
// Input configuration and potential
	AtomisticSystem body;
	int numberOfParticles= 10;
	body.numberOfParticles= numberOfParticles;
	body.speciesCode.resize(numberOfParticles);
	body.coordinates.resize(numberOfParticles,3);


// Create grid
	int ngridx,ngridy,ngridz;
	ngridx= ngridy= ngridz= 10;
	Grid grid(ngridx,ngridy,ngridz);

	if (stressFlag & HARDY)
	{
		// Create hardyStress object
		bool requestCauchy= true;
		bool requestPiola= true;
		Hardy hardyStress(requestPiola,requestCauchy);
		// Calculate stress
		std::vector<Stress*> pStressFields;
		pStressFields.push_back(&hardyStress);
		calculateStress(body,grid.coordinates,pStressFields);
	}
	return 0;
}


