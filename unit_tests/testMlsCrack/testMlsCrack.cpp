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
#include "MethodSphere.h"
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

	std::ifstream file(configFileName);
	if(!file) MY_ERROR("ERROR: config_T.data could not be opened for reading!");

	file >> numberOfParticles;
	if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");

	BoxConfiguration body{numberOfParticles,referenceAndFinal};
	body.read(configFileName,referenceAndFinal);

	Kim kim(modelname);

	int ngrid = 301;
    Grid<Reference> grid({9950,1870,13.577375}, {10050,1970.0001,13.577376},ngrid,ngrid);
    grid.write("referenceGrid");


	MethodSphere hardy(5,"hardy");
	Stress<MethodSphere,Piola> hardyPiolaStress("hardy",hardy,&grid);
    calculateStress(body,kim,
					std::tie(hardyPiolaStress));
    hardyPiolaStress.write("hardyPiola");


    double MlsRadius = 16.29285;
	Mls testMls(body,&grid,MlsRadius,"mls_crack");
	std::vector<Matrix3d> cauchyPushedField;
    testMls.pushToCauchy(hardyPiolaStress.field,cauchyPushedField);
    testMls.writeDeformationGradient();
	testMls.writeGridPushed();
	testMls.writePushedCauchyStress(cauchyPushedField);

	return 0;
}


