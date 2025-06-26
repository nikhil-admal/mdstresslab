/*
 * main.cpp
 *
 *  Created on: Nov 3, 2019
 *      Author: Nikhil Admal
 */
#include "MethodLdad.h"
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
    /*! [Read] */
	int numberOfParticles;
	std::string configFileName= "config.data";
	std::ifstream file(configFileName);
	if(!file) MY_ERROR("ERROR: config.dat could not be opened for reading!");
	file >> numberOfParticles;
	if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");
    /*! [Read] */

    /*! [Configuration] */
    int referenceAndFinal= true;
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
    Grid<Current> gridFromFile_def(ngrid);
	gridFromFile_ref.read("grid_pk1.data");
	gridFromFile_def.read("grid_cauchy.data");
    /*! [Grid] */


    /*! [LDAD] */
    Matrix3d ldadVectors_ref;
    ldadVectors_ref << 5.29216036151419, 0.0, 0.0,
                       0.0, 5.29216036151419, 0.0,
                       0.0, 0.0, 5.29216036151419;
	MethodLdadConstant ldad_constant_ref(ldadVectors_ref);
    MethodLdadTrigonometric ldad_trigonometric_ref(ldadVectors_ref);
    /*! [LDAD] */

    /*! [Stress] */
	//TODO The bond function should be accepted as a reference
	Stress<MethodLdadConstant,Piola> ldad_constant_stress_ref("ldad_constant_ref",ldad_constant_ref,&gridFromFile_ref);
    Stress<MethodLdadTrigonometric,Piola> ldad_trigonometric_stress_ref("ldad_trigonometric_ref",ldad_trigonometric_ref,&gridFromFile_ref);
    /*! [Stress] */

    /*! [Calculate] */
	calculateStress(body,kim,std::tie(ldad_constant_stress_ref));
    ldad_constant_stress_ref.write();
	calculateStress(body,kim,std::tie(ldad_trigonometric_stress_ref));
    ldad_trigonometric_stress_ref.write();
    /*! [Calculate] */

    /*! [Cauchy] */
    Matrix3d ldadVectors_def;
    ldadVectors_def << 5.29216036151419, 0.0, 0.0,
                       0.0, 5.3450819651293319, 0.0,
                       0.0, 0.0, 5.29216036151419;   
	MethodLdadConstant ldad_constant_def(ldadVectors_def);
    MethodLdadTrigonometric ldad_trigonometric_def(ldadVectors_def);
	Stress<MethodLdadConstant,Cauchy> ldad_constant_stress_def("ldad_constant_def",ldad_constant_def,&gridFromFile_def);
    Stress<MethodLdadTrigonometric,Cauchy> ldad_trigonometric_stress_def("ldad_trigonometric_def",ldad_trigonometric_def,&gridFromFile_def);
	calculateStress(body,kim,std::tie(ldad_constant_stress_def));
    ldad_constant_stress_def.write();
	calculateStress(body,kim,std::tie(ldad_trigonometric_stress_def));
    ldad_trigonometric_stress_def.write();
    /*! [Cauchy] */

    /*! [Compare] */
	compareStress("ldad_constant_ref");
	compareStress("ldad_constant_def");
	compareStress("ldad_trigonometric_ref");
	compareStress("ldad_trigonometric_def");
    /*! [Compare] */

	return 0;
}


