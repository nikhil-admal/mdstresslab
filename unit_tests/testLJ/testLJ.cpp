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
#include "../compareStress.cpp"

/*!
 * @example testLJ.cpp
 * An example demonstrating the computation of Piola and Cauchy stress tensors
 * using various spherical averaging domains and grids
 *
 * Read the configuration of particles
 * @snippet{lineno} testLJ.cpp Read
 *
 * Link to the Kim model
 * @snippet{lineno} testLJ.cpp Model
 *
 * Create grids
 * @snippet{lineno} testLJ.cpp Grid
 *
 * Create Stress objects
 * @snippet{lineno} testLJ.cpp Stress
 *
 * Stresses can be calculated all at once or with any combinations
 * @snippet{lineno} testLJ.cpp Calculate
 *
 * Write stresses
 * @snippet{lineno} testLJ.cpp Write
 *
 * Compare stresses
 * @snippet{lineno} testLJ.cpp Compare
 *
 * Full code:
 */

int main()
{
    /*! [Read] */
	int numberOfParticles;
	int referenceAndFinal= true; // Does the input file include reference and final configurations
	std::string configFileName= "config.data";
    std::ifstream file(configFileName);
    if(!file) MY_ERROR("ERROR: config.data could not be opened for reading!\n");

    file >> numberOfParticles;
    if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");

	BoxConfiguration body{numberOfParticles,referenceAndFinal};
	body.read(configFileName,referenceAndFinal);
    /*! [Read] */

    /*! [Model] */
    std::string modelname= "LJ_Smoothed_Bernardes_1958_Ar__MO_764178710049_001";
	Kim kim(modelname);
    /*! [Model] */

    /*! [Grid] */
	int ngrid=10;
	Grid<Current> randomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid);
	Grid<Reference> referenceRandomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid);

	ngrid= 20;
	Grid<Reference> reference_grid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid);

	ngrid= 125;
	Grid<Current> gridFromFile(ngrid);
	gridFromFile.read("grid_cauchy.data");
    /*! [Grid] */

    /*! [Stress] */
	// MethodSphere stress 1
	MethodSphere hardy1(5.29216036151419,"hardy");
	Stress<MethodSphere,Cauchy> hardyStress1("hardy1",hardy1,&gridFromFile);

	// MethodSphere stress 2
	MethodSphere hardy2(20,"hardy");
	Stress<MethodSphere,Cauchy> hardyStress2("hardy2",hardy2,&gridFromFile);

	// MethodSphere stress 3
	MethodSphere hardy3(5,"hardy");
	Stress<MethodSphere,Piola> hardyStress3("hardy3",hardy3,&reference_grid);

	// MethodSphere stress 4
	MethodSphere hardy4(7,"hardy");
	Stress<MethodSphere,Piola> hardyStress4("hardy4",hardy4,&reference_grid);

	MethodSphere hardyRandom(9,"hardy");
	Stress<MethodSphere,Cauchy> hardyStressRandomCauchy("hardyRandomCauchy",hardyRandom,&randomGrid);

	Stress<MethodSphere,Piola> hardyStressRandomPiola("hardyRandomPiola",hardyRandom,&referenceRandomGrid);
    /*! [Stress] */


    /*! [Calculate] */
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
    /*! [Calculate] */

    /*! [Write] */
    hardyStress1.write();
    hardyStress2.write();
    hardyStress3.write();
    hardyStress4.write();
    hardyStressRandomPiola.write();
    hardyStressRandomCauchy.write();
    /*! [Write] */

    /*! [Compare] */
	compareStress("hardy1");
	compareStress("hardy3");
	compareStress("hardy4");
	compareStress("hardy2");
	compareStress("hardyRandomCauchy");
	compareStress("hardyRandomPiola");
    /*! [Compare] */


	return 0;
}


