/*
 * main.cpp
 *
 *  Created on: Nov 9, 2020
 *      Author: Min Shi
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
#include "../compareStress.cpp"

int main()
{

//	-------------------------------------------------------------------
//	Input configuration and potential
//	-------------------------------------------------------------------

	int numberOfParticles;
	int referenceAndFinal=true;
	std::string configFileName="config.data";
	std::string modelname=\
	"SW_StillingerWeber_1985_Si__MO_405512056662_005";

	std::ifstream file(configFileName);
	if(!file) \
	MY_ERROR("ERROR: config.data could not be opened for reading!");
	
	file >> numberOfParticles;
	if (numberOfParticles < 0) \
	MY_ERROR("Error: Negative number of particles.\n");

	BoxConfiguration body{numberOfParticles,referenceAndFinal};
	body.read(configFileName,referenceAndFinal);

	Kim kim(modelname);

//	-------------------------------------------------------------------
//	Create grid
//	-------------------------------------------------------------------

	int ngrid_pk1 = 125;
	Grid<Reference> gridFromFile_ref(ngrid_pk1);
	gridFromFile_ref.read("grid_pk1.data");

	int ngrid_cauchy = 125;
	Grid<Current> gridFromFile_def(ngrid_cauchy);
	gridFromFile_def.read("grid_cauchy.data");

//	-------------------------------------------------------------------
//	Calculate stress on the grid
//	-------------------------------------------------------------------

	Sphere hardy1(20);
	Stress<Sphere,Cauchy> \
	hardyStress1("hardy1",hardy1,&gridFromFile_def);

	calculateStress(body,kim,std::tie(),std::tie(hardyStress1));

	Matrix3d ldadVectors_ref;

	ldadVectors_ref << 5.430949777364731, 0.0, 0.0, 
			0.0, 5.430949777364731, 0.0,
			0.0,  0.0, 5.430949777364731;   

	// Ldad stress ref
	Ldad<Constant> ldad_Constant_ref(ldadVectors_ref);
	Ldad<Trigonometric> ldad_Trigonometric_ref(ldadVectors_ref);

	//TODO The bond function should be accepted as a reference
	Stress<Ldad<Constant>,Piola> \
	ldad_Constant_Stress_ref("ldad_Constant_ref",\
	ldad_Constant_ref,&gridFromFile_ref);
	Stress<Ldad<Trigonometric>,Piola> \
	ldad_Trigonometric_Stress_ref("ldad_Trigonometric_ref",\
	ldad_Trigonometric_ref,&gridFromFile_ref);

	calculateStress(body,kim,std::tie(ldad_Constant_Stress_ref));
	calculateStress(body,kim,std::tie(ldad_Trigonometric_Stress_ref));

//	-------------------------------------------------------------------
//	Perform MLS
//	-------------------------------------------------------------------

	double MlsRadius = 5.430949777364731 * 3.0;
	Mls testMls(body,&gridFromFile_ref,MlsRadius,\
	"ldad_Constant_Stress_ref");
	std::vector<Matrix3d> cauchyPushedField;
	testMls.pushToCauchy(ldad_Constant_Stress_ref.field,\
	cauchyPushedField);
	testMls.writeDeformationGradient();
	testMls.writeGridPushed();
	testMls.writePushedCauchyStress(cauchyPushedField);

	return 0;
}

