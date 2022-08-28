/*
 * main.cpp
 *
 *  Created on: Nov 3, 2019
 *      Author: Nikhil Admal
 */
#include "MethodSphere.h"
#include <string>
#include <iostream>
#include <tuple>
#include <fstream>
#include "BoxConfiguration.h"
#include "calculateStress.h"
#include "Grid.h"
#include "typedef.h"
#include "../compareStress.cpp"


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
	int ngrid=5;
	Grid<Current> randomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid,ngrid);
	Grid<Reference> referenceRandomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid,ngrid);

	ngrid= 3;
	Grid<Reference> reference_grid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid,ngrid);

	Grid<Current> gridFromFile("grid_cauchy.data");


//	-------------------------------------------------------------------
// Calculate stress on the grid
//	-------------------------------------------------------------------
	// Create hardyStress object

	// MethodSphere stress 1
	MethodSphere hardy1(5.29216036151419,"hardy");

	Stress<MethodSphere,Cauchy> hardyStress1("hardy1",hardy1,&gridFromFile);

	// MethodSphere stress 2
	MethodSphere hardy2(20,"hardy");
	Stress<MethodSphere,Cauchy> hardyStress2("hardy2",hardy2,&gridFromFile);

	// MethodSphere stress 3
	MethodSphere hardy3(5,"hardy");
	Stress<MethodSphere,Piola> hardyStress3("hardy3",hardy3,&reference_grid);

	// MethodSphere stress 4
	MethodSphere hardy4(7,"hardy");
	Stress<MethodSphere,Piola> hardyStress4("hardy4",hardy4,&reference_grid);

	MethodSphere hardyRandom(9,"hardy");
	Stress<MethodSphere,Cauchy> hardyStressRandomCauchy("hardyRandomCauchy",hardyRandom,&randomGrid);

	Stress<MethodSphere,Piola> hardyStressRandomPiola("hardyRandomPiola",hardyRandom,&referenceRandomGrid);

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

    hardyStress1.write();
    hardyStress2.write();
    hardyStress3.write();
    hardyStress4.write();
    hardyStressRandomPiola.write();
    hardyStressRandomCauchy.write();

	compareStress("hardy1");
	compareStress("hardy2");
	compareStress("hardy3");
	compareStress("hardy4");
	compareStress("hardyRandomCauchy");
	compareStress("hardyRandomPiola");

	return 0;
}


