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
	if(!file) MY_ERROR("ERROR: config.data could not be opened for reading!");

	file >> numberOfParticles;
	if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");

	BoxConfiguration body{numberOfParticles,referenceAndFinal};
	body.read(configFileName,referenceAndFinal);

	Kim kim(modelname);

//	-------------------------------------------------------------------
// Create grid
//	-------------------------------------------------------------------
	//int ngrid=100;
	//Grid<Current> randomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid);
	//Grid<Reference> referenceRandomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid);
/*
	int ngrid = 1;
    Grid<Reference> gridFromFile_ref(ngrid);
	gridFromFile_ref.read("grid_pk1_one.data");
*/
	int ngrid = 125;
	Grid<Reference> gridFromFile_ref(ngrid);
	gridFromFile_ref.read("grid_pk1.data");

    Grid<Current> gridFromFile_def(ngrid);
	gridFromFile_def.read("grid_cauchy.data");

	//ngrid= 125;
	//Grid<Current> gridFromFile(ngrid);
	//gridFromFile.read("grid_cauchy.data");


//	-------------------------------------------------------------------
// Calculate stress on the grid
//	-------------------------------------------------------------------
	// Create hardyStress object

	// Sphere stress 1
	//Sphere hardy1(5.29216036151419);

	//TODO The bond function should be accepted as a reference
	//Stress<Sphere,Cauchy> hardyStress1("hardy1",hardy1,&gridFromFile);

	// Sphere stress 2
	//Sphere hardy2(20);
	//Stress<Sphere,Cauchy> hardyStress2("hardy2",hardy2,&gridFromFile);

	// Sphere stress 3
	//Sphere hardy3(5);
	//Stress<Sphere,Piola> hardyStress3("hardy3",hardy3,&gridFromFile);

	// Sphere stress 4
	//Sphere hardy4(7);
	//Stress<Sphere,Piola> hardyStress4("hardy4",hardy4,&reference_grid);

	//Sphere hardyRandom(9);
	//Stress<Sphere,Cauchy> hardyStressRandomCauchy("hardyRandomCauchy",hardyRandom,&randomGrid);

	//Stress<Sphere,Piola> hardyStressRandomPiola("hardyRandomPiola",hardyRandom,&referenceRandomGrid);

//	Calculate only  Piola
	//calculateStress(body,kim,
	//				std::tie(hardyStressRandomCauchy));

//	Calculate only Cauchy
	//calculateStress(body,kim,
	//				std::tie(hardyStressRandomPiola));

//  Calculate none
	//calculateStress(body,kim,
	//				std::tie(),
	//				std::tie());

//  Calculate all
	//calculateStress(body,kim,
	//				std::tie(hardyStress3,hardyStress4),
	//				std::tie(hardyStress1));
    //calculateStress(body,kim,
	//				std::tie(hardyStress3));
	//calculateStress(body,kim,
	//				std::tie(),
	//				std::tie(hardyStress2));

/*
	compareStress("hardy1");
	compareStress("hardy2");
	compareStress("hardy3");
	compareStress("hardy4");
	compareStress("hardyRandomCauchy");
	compareStress("hardyRandomPiola");
*/

    // MLS

//	-------------------------------------------------------------------
// Calculate stress on the grid
//	-------------------------------------------------------------------
	// Create hardyStress object

    Matrix3d ldadVectors_ref;

    ldadVectors_ref << 5.29216036151419, 0.0, 0.0, 
                       0.0, 5.29216036151419, 0.0,
                       0.0, 0.0, 5.29216036151419;   

	// Ldad stress 1
	Ldad<Constant> ldad_Constant_ref(ldadVectors_ref);
    //Ldad<Trigonometric> ldad_Trigonometric_ref(ldadVectors_ref);

	//TODO The bond function should be accepted as a reference
	Stress<Ldad<Constant>,Piola> ldad_Constant_Stress_ref("ldad_Constant_ref",ldad_Constant_ref,&gridFromFile_ref);
    //Stress<Ldad<Trigonometric>,Piola> ldad_Trigonometric_Stress_ref("ldad_Trigonometric_ref",ldad_Trigonometric_ref,&gridFromFile_ref);

	calculateStress(body,kim,std::tie(ldad_Constant_Stress_ref));
	//calculateStress(body,kim,std::tie(ldad_Trigonometric_Stress_ref));

/*

    Matrix3d ldadVectors_def;

    ldadVectors_def << 5.29216036151419, 0.0, 0.0, 
                       0.0, 5.3450819651293319, 0.0,
                       0.0, 0.0, 5.29216036151419;   

	// Ldad stress 1
	Ldad<Constant> ldad_Constant_def(ldadVectors_def);
    Ldad<Trigonometric> ldad_Trigonometric_def(ldadVectors_def);

	//TODO The bond function should be accepted as a reference
	Stress<Ldad<Constant>,Cauchy> ldad_Constant_Stress_def("ldad_Constant_def",ldad_Constant_def,&gridFromFile_def);
    Stress<Ldad<Trigonometric>,Cauchy> ldad_Trigonometric_Stress_def("ldad_Trigonometric_def",ldad_Trigonometric_def,&gridFromFile_def);

	calculateStress(body,kim,std::tie(ldad_Constant_Stress_def));
	calculateStress(body,kim,std::tie(ldad_Trigonometric_Stress_def));
*/

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


