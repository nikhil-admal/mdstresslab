/*
 * main.cpp
 *
 *  Created on: Nov 3, 2019
 *      Author: Nikhil Admal
 */
#include "Sphere.h"
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
	std::string modelname= "LJ_Smoothed_Bernardes_1958_Ar__MO_764178710049_001";

//	-------------------------------------------------------------------
// Input configuration and potential
//	-------------------------------------------------------------------

	std::ifstream file(configFileName);
	if(!file) MY_ERROR("ERROR: config.data could not be opened for reading!\n");

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

	// Sphere stress 1
	Sphere hardy1(5.29216036151419);

	//TODO The bond function should be accepted as a reference
	Stress<Sphere,Cauchy> hardyStress1("hardy1",hardy1,&gridFromFile);

	// Sphere stress 2
	Sphere hardy2(20);
	Stress<Sphere,Cauchy> hardyStress2("hardy2",hardy2,&gridFromFile);

	// Sphere stress 3
	Sphere hardy3(5);
	Stress<Sphere,Piola> hardyStress3("hardy3",hardy3,&reference_grid);

	// Sphere stress 4
	Sphere hardy4(7);
	Stress<Sphere,Piola> hardyStress4("hardy4",hardy4,&reference_grid);

	Sphere hardyRandom(9);
	Stress<Sphere,Cauchy> hardyStressRandomCauchy("hardyRandomCauchy",hardyRandom,&randomGrid);

	Stress<Sphere,Piola> hardyStressRandomPiola("hardyRandomPiola",hardyRandom,&referenceRandomGrid);


//  Calculate none
	calculateStress(body,kim,
					std::tie(),
					std::tie());

//  Calculate all
	calculateStress(body,kim,
					std::tie(hardyStress3,hardyStress4,hardyStressRandomPiola),
					std::tie(hardyStress1,hardyStress2,hardyStressRandomCauchy));

    hardyStress1.write();
    hardyStress2.write();
    hardyStress3.write();
    hardyStress4.write();
    hardyStressRandomPiola.write();
    hardyStressRandomCauchy.write();

	compareStress("hardy1");
	compareStress("hardy3");
	compareStress("hardy4");
	compareStress("hardy2");
	compareStress("hardyRandomCauchy");
	compareStress("hardyRandomPiola");


	return 0;
}


