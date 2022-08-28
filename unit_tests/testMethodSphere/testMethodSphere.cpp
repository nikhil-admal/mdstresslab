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
#include <map>
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
	int ngrid=10;
	Grid<Current> randomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid);
	Grid<Reference> referenceRandomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid);

	ngrid= 20;
	Grid<Reference> reference_grid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid);

	ngrid= 125;
	Grid<Current> gridFromFile(ngrid);
	gridFromFile.read("grid_cauchy.data");


//	-------------------------------------------------------------------
// Calculate stress on the grid
//	-------------------------------------------------------------------
	// Create hardyStress object

	// MethodSphere stress 1
    double averagingDomainSize= 5.29216036151419;
	MethodSphere hardy1(averagingDomainSize,{{0,1},{1./3,1},{1./3,1},{1./2,1},{1,0}});
	Stress<MethodSphere,Cauchy> hardyStress1("hardy1",hardy1,&gridFromFile);

	// MethodSphere stress 2
    averagingDomainSize= 20;
	MethodSphere hardy2(averagingDomainSize,{{0,1},{1./2,1},{2./3,2./3},{1,0}});
	Stress<MethodSphere,Cauchy> hardyStress2("hardy2",hardy2,&gridFromFile);

	// MethodSphere stress 3
    averagingDomainSize= 5;
	MethodSphere hardy3(averagingDomainSize,{{0,1},{1./2,1},{1,0}});
	Stress<MethodSphere,Piola> hardyStress3("hardy3",hardy3,&reference_grid);

	// MethodSphere stress 4
    averagingDomainSize= 7;
	MethodSphere hardy4(averagingDomainSize,{{0,1},{1./2,1},{1,0}});
	Stress<MethodSphere,Piola> hardyStress4("hardy4",hardy4,&reference_grid);

    averagingDomainSize= 9;
	MethodSphere hardyRandom(averagingDomainSize,{{0,1},{1./2,1},{1,0}});
	Stress<MethodSphere,Cauchy> hardyStressRandomCauchy("hardyRandomCauchy",hardyRandom,&randomGrid);

	Stress<MethodSphere,Piola> hardyStressRandomPiola("hardyRandomPiola",hardyRandom,&referenceRandomGrid);


//  Calculate none
	calculateStress(body,kim,
					std::tie(),
					std::tie());

//  Calculate all
	calculateStress(body,kim,
					std::tie(hardyStress3,hardyStress4,hardyStressRandomPiola),
					std::tie(hardyStress1,hardyStress2,hardyStressRandomCauchy));

	calculateStress(body,kim,
					std::tie(hardyStress1));

	calculateStress(body,kim,
                    std::tie(hardyStress3),
					std::tie(hardyStress1));

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


