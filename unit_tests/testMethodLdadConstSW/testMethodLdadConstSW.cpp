/*
 * main.cpp
 *
 *  Created on: Nov 3, 2019
 *      Author: Nikhil Admal
 */
#include "MethodLdadConst.h"
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

	int ngrid = 125;
	Grid<Reference> gridFromFile_ref(ngrid);
	gridFromFile_ref.read("grid_pk1.data");

    Grid<Current> gridFromFile_def(ngrid);
	gridFromFile_def.read("grid_cauchy.data");

//	-------------------------------------------------------------------
// Calculate stress on the grid
//	-------------------------------------------------------------------
	// Create hardyStress object

    Matrix3d ldadVectors_ref;

    ldadVectors_ref << 5.43094977840521, 0.0, 0.0, 
                       0.0, 5.43094977840521, 0.0,
                       0.0,  0.0, 5.43094977840521;   

	// Ldad stress ref
	MethodLdadConst ldad_Constant_ref(ldadVectors_ref);

	//TODO The bond function should be accepted as a reference
	Stress<MethodLdadConst,Piola> ldad_Constant_Stress_ref("ldad_Constant_ref",ldad_Constant_ref,&gridFromFile_ref);

	calculateStress(body,kim,std::tie(ldad_Constant_Stress_ref));

    Matrix3d ldadVectors_def;

    ldadVectors_def << 5.43094977840521, 0.0, 0.0,
	                   0.0, 5.4852592761892621, 0.0,
					   0.0, 0.0, 5.43094977840521;

	// Ldad stress def
	MethodLdadConst ldad_Constant_def(ldadVectors_def);

	//TODO The bond function should be accepted as a reference
	Stress<MethodLdadConst,Cauchy> ldad_Constant_Stress_def("ldad_Constant_def",ldad_Constant_def,&gridFromFile_def);

	calculateStress(body,kim,std::tie(ldad_Constant_Stress_def));

	return 0;
}


