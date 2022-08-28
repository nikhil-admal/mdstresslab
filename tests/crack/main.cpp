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


int main()
{
	int numberOfParticles;
	int referenceAndFinal= true;
	std::string configFileName= "config.data";
	std::string modelname= "SW_BalamaneHauchShi_2017Brittle_Si__MO_381114941873_002";

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
	int ngrid;
	ngrid= 90601;
	Grid<Current> gridFromFile(ngrid);
	gridFromFile.read("grid_cauchy.data");


//	-------------------------------------------------------------------
// Calculate stress on the grid
//	-------------------------------------------------------------------
	// Create hardyStress object

	// Hardy stress
	MethodSphere hardy5(5,"hardy");
	MethodSphere hardy10(10,"hardy");
	MethodSphere hardy15(15,"hardy");
	MethodSphere hardy20(20,"hardy");

	Stress<MethodSphere,Cauchy> hardyStress5("hardy5",hardy5,&gridFromFile);
	Stress<MethodSphere,Cauchy> hardyStress10("hardy10",hardy10,&gridFromFile);
	Stress<MethodSphere,Cauchy> hardyStress15("hardy15",hardy15,&gridFromFile);
	Stress<MethodSphere,Cauchy> hardyStress20("hardy20",hardy20,&gridFromFile);

	calculateStress(body,kim,
					std::tie(),
					std::tie(hardyStress5,hardyStress10,hardyStress15,hardyStress20));
    hardyStress5.write();
    hardyStress10.write();
    hardyStress15.write();
    hardyStress20.write();


	return 0;
}


