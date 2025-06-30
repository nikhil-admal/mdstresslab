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
	int ngrid=200;
	Grid<Current> current2DGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid);
	Grid<Reference> reference2DGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid,ngrid);

	ngrid= 500;
	Grid<Reference> reference1DGrid(Vector3d(0,0,0),Vector3d(60,60,60),ngrid);
    /*! [Grid] */

    /*! [Method] */
	MethodSphere hardy5(5.0,"hardy");
    MethodSphere virial5(5.0,"virial");
    /*! [Method] */

    /*! [Stress] */
    Stress<MethodSphere,Cauchy> hardyCauchy("hardyCauchy",hardy5,&current2DGrid);
    Stress<MethodSphere,Cauchy> virialCauchy("virialCauchy",virial5,&current2DGrid);

	Stress<MethodSphere,Piola> hardyPiola("hardyPiola",hardy5,&reference1DGrid);
    Stress<MethodSphere,Piola> virialPiola("virialPiola",virial5,&reference1DGrid);
    /*! [Stress] */


    /*! [Calculate] */
	calculateStress(body,kim,
					std::tie(),
					std::tie());

	calculateStress(body,kim,
					std::tie(hardyCauchy));

	calculateStress(body,kim,
                    std::tie(virialPiola));

    calculateStress(body,kim,
                    std::tie(hardyPiola,virialPiola),
                    std::tie(hardyCauchy,virialCauchy));

    /*! [Calculate] */

    /*! [Write] */
    hardyCauchy.write();
    hardyPiola.write();
    virialCauchy.write();
    virialPiola.write();
    /*! [Write] */

    /*! [Compare] */
	compareStress("hardyCauchy");
	compareStress("hardyPiola");
	compareStress("virialCauchy");
	compareStress("virialPiola");
    /*! [Compare] */


	return 0;
}


