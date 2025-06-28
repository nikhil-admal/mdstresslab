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
#include "MethodLdad.h"
#include "BoxConfiguration.h"
#include "calculateStress.h"
#include "Mls.h"
#include "Grid.h"
#include "typedef.h"
#include "../compareStress.cpp"



int main()
{
    /*! [Read] */
	int referenceAndFinal= true;
	std::string configFileName= "config.data";
	std::ifstream file(configFileName);
	if(!file) MY_ERROR("ERROR: config.data could not be opened for reading!");
	int numberOfParticles;
	file >> numberOfParticles;
	if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");
    /*! [Read] */

    /*! [Configuration] */
	BoxConfiguration body{numberOfParticles,referenceAndFinal};
	body.read(configFileName,referenceAndFinal);
    /*! [Configuration] */

    /*! [Model] */
	std::string modelname= "LJ_Smoothed_Bernardes_1958_Ar__MO_764178710049_001";
	Kim kim(modelname);
    /*! [Model] */

    /*! [Grid] */
	int ngrid = 125;
	Grid<Reference> gridFromFile_ref(ngrid);
	gridFromFile_ref.read("reference.grid");
    /*! [Grid] */


    /*! [LDAD] */
    Matrix3d ldadVectors;
    ldadVectors << 5.29216036151419, 0.0, 0.0, 
                   0.0, 5.29216036151419, 0.0,
                   0.0, 0.0, 5.29216036151419;
	MethodLdadConstant ldad(ldadVectors);
    /*! [LDAD] */

    /*! [Stress] */
	//TODO The bond function should be accepted as a reference
	Stress<MethodLdadConstant,Piola> ldad_constant_stress_ref("ldad",ldad,&gridFromFile_ref);
    /*! [Stress] */


    /*! [Calculate] */
	calculateStress(body,kim,std::tie(ldad_constant_stress_ref));
    /*! [Calculate] */


    /*! [MLS] */
    double MlsRadius = 15.87648108454257;
	Mls testMls(body,&gridFromFile_ref,MlsRadius,"mls");
    /*! [MLS] */


    /*! [Output] */
	std::vector<Matrix3d> cauchyPushedField;
    testMls.pushToCauchy(ldad_constant_stress_ref.field,cauchyPushedField);
    testMls.writeDeformationGradient();
	testMls.writeGridPushed();
	testMls.writePushedCauchyStress(cauchyPushedField);
    /*! [Output] */

    /*! [Compare] */
    compareStress("mlsPushed");
    /*! [Compare] */


	return 0;
}


