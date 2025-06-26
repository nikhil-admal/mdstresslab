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
	int ngrid=10;
	Grid<Current> randomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid);
	Grid<Reference> referenceRandomGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid);

	ngrid= 20;
	Grid<Reference> reference_grid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid);

	ngrid= 125;
	Grid<Current> gridFromFile(ngrid);
	gridFromFile.read("grid_cauchy.data");
    /*! [Grid] */

    /*! [Method] */
	MethodSphere hardy1(5.29216036151419,"hardy");
	MethodSphere hardy2(20,"hardy");
    MethodSphere hardy3(5,"hardy");
    MethodSphere hardy4(7,"hardy");
    MethodSphere hardyRandom(9,"hardy");
    /*! [Method] */

    /*! [Stress] */
    Stress<MethodSphere,Cauchy> hardyStress1("hardy1",hardy1,&gridFromFile);
	Stress<MethodSphere,Cauchy> hardyStress2("hardy2",hardy2,&gridFromFile);
	Stress<MethodSphere,Piola> hardyStress3("hardy3",hardy3,&reference_grid);
	Stress<MethodSphere,Piola> hardyStress4("hardy4",hardy4,&reference_grid);
	Stress<MethodSphere,Cauchy> hardyStressRandomCauchy("hardyRandomCauchy",hardyRandom,&randomGrid);
	Stress<MethodSphere,Piola> hardyStressRandomPiola("hardyRandomPiola",hardyRandom,&referenceRandomGrid);
    /*! [Stress] */


    /*! [Calculate] */
	calculateStress(body,kim,
					std::tie(),
					std::tie());

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


