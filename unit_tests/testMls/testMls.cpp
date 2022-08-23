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
#include "Ldad.h"
#include "Constant.h"
#include "Trigonometric.h"
#include "BoxConfiguration.h"
#include "Sphere.h"
#include "calculateStress.h"
#include "Mls.h"
#include "Grid.h"
#include "typedef.h"

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
	if(!file) MY_ERROR("ERROR: config.data could not be opened for reading!");

	file >> numberOfParticles;
	if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");

	BoxConfiguration body{numberOfParticles,referenceAndFinal};
	body.read(configFileName,referenceAndFinal);

	Kim kim(modelname);

	int ngrid = 125;
	Grid<Reference> gridFromFile_ref(ngrid);
	gridFromFile_ref.read("grid_pk1.data");

    Grid<Current> gridFromFile_def(ngrid);
	gridFromFile_def.read("grid_cauchy.data");


    Matrix3d ldadVectors_ref;

    ldadVectors_ref << 5.29216036151419, 0.0, 0.0, 
                       0.0, 5.29216036151419, 0.0,
                       0.0, 0.0, 5.29216036151419;   

	// Ldad stress 1
	Ldad<Constant> ldad_Constant_ref(ldadVectors_ref);

	//TODO The bond function should be accepted as a reference
	Stress<Ldad<Constant>,Piola> ldad_Constant_Stress_ref("ldad_Constant_ref",ldad_Constant_ref,&gridFromFile_ref);

	calculateStress(body,kim,std::tie(ldad_Constant_Stress_ref));


    double MlsRadius = 15.87648108454257;
    //Mls testMls(body.coordinates.at(Reference),body.coordinates.at(Current),gridFromFile.coordinates,MlsRadius,"hardyStress3");
    //double MlsRadius = 16.29285;
	Mls testMls(body,&gridFromFile_ref,MlsRadius,"ldad_Constant_Stress_ref");
	std::vector<Matrix3d> cauchyPushedField;
    testMls.pushToCauchy(ldad_Constant_Stress_ref.field,cauchyPushedField);
    testMls.writeDeformationGradient();
	testMls.writeGridPushed();
	testMls.writePushedCauchyStress(cauchyPushedField);

	return 0;
}


