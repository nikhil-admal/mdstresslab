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
#include "MethodHardySphere.h"
#include "calculateStress.h"
#include "Mls.h"
#include "Grid.h"
#include "typedef.h"

int main()
{
	int numberOfParticles;
	int referenceAndFinal= true;
	std::string configFileName= "config_T.data";
	std::string modelname= "SW_BalamaneHauchShi_2017Brittle_Si__MO_381114941873_002";

    //	-------------------------------------------------------------------
    //  Input configuration and potential
    //	-------------------------------------------------------------------

	std::ifstream file(configFileName);
	if(!file) MY_ERROR("ERROR: config_T.data could not be opened for reading!");

	file >> numberOfParticles;
	if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");

	BoxConfiguration body{numberOfParticles,referenceAndFinal};
	body.read(configFileName,referenceAndFinal);

	Kim kim(modelname);

	int ngrid = 10;
    Grid<Reference> gridFromFile(ngrid);
	gridFromFile.read("grid_ten.data");


	// MethodHardySphere stress 3
	MethodHardySphere hardy3(5);
	Stress<MethodHardySphere,Piola> hardyStress3("hardy3",hardy3,&gridFromFile);

    calculateStress(body,kim,
					std::tie(hardyStress3));

    // MLS
    double MlsRadius = 16.29285;
	Mls testMls(body,&gridFromFile,MlsRadius,"hardyStress3");
	std::vector<Matrix3d> cauchyPushedField3;
    testMls.pushToCauchy(hardyStress3.field,cauchyPushedField3);
    testMls.writeDeformationGradient();
	testMls.writeGridPushed();
	testMls.writePushedCauchyStress(cauchyPushedField3);

	return 0;
}


