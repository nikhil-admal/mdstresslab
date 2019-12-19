/*
 * main.cpp
 *
 *  Created on: Nov 3, 2019
 *      Author: Nikhil Admal
 */
#include <string>
#include <iostream>
#include <tuple>
#include <fstream>
#include "BoxConfiguration.h"
#include "Hardy.h"
#include "calculateStress.h"
#include "Grid.h"
#include "typedef.h"


int main()
{
	int numberOfParticles;
	int referenceAndFinal= true;
	std::string configFileName= "config.data";
	std::string modelname= "SW_StillingerWeber_1985_Si__MO_405512056662_005";

//	-------------------------------------------------------------------
// Input configuration and potential
//	-------------------------------------------------------------------

	std::ifstream file(configFileName);
	if(!file) MY_ERROR("ERROR: config.dat could not be opened for reading!");

	file >> numberOfParticles;
	if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");

	BoxConfiguration body{numberOfParticles,referenceAndFinal};
	body.read(configFileName,referenceAndFinal);

	Kim kim(modelname);

//	-------------------------------------------------------------------
// Create grid
//	-------------------------------------------------------------------
	int ngrid=100;
	Grid<Current> randomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid);
	Grid<Reference> referenceRandomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid);

	ngrid= 20;
	Grid<Reference> reference_grid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid);

	ngrid= 125;
	Grid<Current> gridFromFile(ngrid);
	gridFromFile.read("grid_cauchy.data");


//	-------------------------------------------------------------------
// Calculate stress on the grid
//	-------------------------------------------------------------------
	// Create hardyStress object

	// Hardy stress 1
	Hardy hardy1(5.29216036151419);

	//TODO The bond function should be accepted as a reference
	Stress<Hardy,Cauchy> hardyStress1("hardy1",hardy1,&gridFromFile);

	// Hardy stress 2
	Hardy hardy2(20);
	Stress<Hardy,Cauchy> hardyStress2("hardy2",hardy2,&gridFromFile);

	// Hardy stress 3
	Hardy hardy3(5);
	Stress<Hardy,Piola> hardyStress3("hardy3",hardy3,&reference_grid);

	// Hardy stress 4
	Hardy hardy4(7);
	Stress<Hardy,Piola> hardyStress4("hardy4",hardy4,&reference_grid);

	Hardy hardyRandom(9);
	Stress<Hardy,Cauchy> hardyStressRandomCauchy("hardyRandomCauchy",hardyRandom,&randomGrid);

	Stress<Hardy,Piola> hardyStressRandomPiola("hardyRandomPiola",hardyRandom,&referenceRandomGrid);

//	Calculate only  Piola
	calculateStress(body,kim,
					std::tie(hardyStressRandomCauchy));

//	Calculate only Cauchy
	calculateStress(body,kim,
					std::tie(hardyStressRandomPiola));

//  Calculate none
	calculateStress(body,kim,
					std::tie(),
					std::tie());

//  Calculate all
	calculateStress(body,kim,
					std::tie(hardyStress3,hardyStress4),
					std::tie(hardyStress1));

	calculateStress(body,kim,
					std::tie(),
					std::tie(hardyStress2));


	return 0;
}


