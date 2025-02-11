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
#include <regex>


int main()
{
	int numberOfParticles;
	int referenceAndFinal= true;
	std::string prefix= "defANDundefSystem_";
    std::vector<std::string> modelnames={
            "MEAM_LAMMPS_DuLenoskyHennig_2011_Si__MO_883726743759_002",
            "EDIP_JustoBazantKaxiras_1998_Si__MO_958932894036_002",
            "Tersoff_LAMMPS_Tersoff_1988T3_Si__MO_186459956893_004",   //  - has processdEdr
            "ThreeBodyCluster_Gong_Gong_1993_Si__MO_407755720412_000",
            "SW_BalamaneHauchShi_2017Brittle_Si__MO_381114941873_002"};


           /*
            "MEAM_LAMMPS_DuLenoskyHennig_2011_Si__MO_883726743759_002",
            "SNAP_ZuoChenLi_2019_Si__MO_869330304805_000",
            "SNAP_ZuoChenLi_2019quadratic_Si__MO_721469752060_000",
            "ThreeBodyBondOrder_PPM_PurjaPunMishin_2017_Si__MO_566683736730_000",
            "EDIP_JustoBazantKaxiras_1998_Si__MO_958932894036_002",
            "EDIP_LAMMPS_JustoBazantKaxiras_1998_Si__MO_315965276297_000",
            "Tersoff_LAMMPS_Tersoff_1988T2_Si__MO_245095684871_004",   //  - has processdEdr
            "Tersoff_LAMMPS_Tersoff_1988T3_Si__MO_186459956893_004",   //  - has processdEdr
            "ThreeBodyBondOrder_KDS_KhorDasSarma_1988_Si__MO_722489435928_000",
            "ThreeBodyBondOrder_PPM_PurjaPunMishin_2017_Si__MO_566683736730_000",
            "ThreeBodyBondOrder_WR_WangRockett_1991_Si__MO_081872846741_000",
            "ThreeBodyCluster_BH_BiswasHamann_1987_Si__MO_019616213550_000",
            "ThreeBodyCluster_Gong_Gong_1993_Si__MO_407755720412_000",
            "ThreeBodyCluster_KP_KaxirasPandey_1988_Si__MO_072486242437_000",
            "SW_BalamaneHauchShi_2017Brittle_Si__MO_381114941873_002",// - has processdEdr
            "ThreeBodyCluster_SRS_StephensonRadnySmith_1996_Si__MO_604248666067_000"
            */

//	-------------------------------------------------------------------
// Create grid
//	-------------------------------------------------------------------
    int ngrid;
    ngrid= 250;
    Vector3d lowerLimit(20*5.5,20*5.5,2*5.5);
    Vector3d upperLimit(80*5.5,80*5.5,10*5.5);
    Grid<Current> gridFromFile(lowerLimit,upperLimit,ngrid,ngrid,1);

//	-------------------------------------------------------------------
// Calculate stress on the grid
//	-------------------------------------------------------------------
    // Create hardyStress object

    // Hardy stress
    MethodSphere hardy5(2,"hardy");
    MethodSphere hardy10(4,"hardy");
    MethodSphere hardy15(6,"hardy");
    MethodSphere hardy20(8,"hardy");

    for (const auto modelname : modelnames)
    {
        //	-------------------------------------------------------------------
        // Input configuration and potential
        //	-------------------------------------------------------------------
        Kim kim(modelname);

        std::string configFileName= prefix+modelname+".data";
        std::ifstream file(configFileName);
        if(!file) MY_ERROR("ERROR: config.dat could not be opened for reading!");
        file >> numberOfParticles;
        if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");
        BoxConfiguration body{numberOfParticles,referenceAndFinal};
        body.read(configFileName,referenceAndFinal);



        // stress calculation using projected forces
        try {
            Stress<MethodSphere,Cauchy> hardyStress5 (hardy5,&gridFromFile);
            Stress<MethodSphere,Cauchy> hardyStress10(hardy10,&gridFromFile);
            Stress<MethodSphere,Cauchy> hardyStress15(hardy15,&gridFromFile);
            Stress<MethodSphere,Cauchy> hardyStress20(hardy20,&gridFromFile);

            // Calculate stress using the projected forces
            calculateStress(body, kim,
                            std::tie(),
                            std::tie(hardyStress5, hardyStress10, hardyStress15, hardyStress20), true);
            hardyStress5.write("project_hardy5_" + modelname);
            hardyStress10.write("project_hardy10_" + modelname);
            hardyStress15.write("project_hardy15_" + modelname);
            hardyStress20.write("project_hardy20_" + modelname);
        }
        catch(const std::runtime_error& e){
            std::cout << e.what() << std::endl;
            std::cout << "Compute stress with projected forces failed. Moving on" << std::endl;
        }

        // stress calculation using process_dedr
        try{
            Stress<MethodSphere,Cauchy> hardyStress5 (hardy5,&gridFromFile);
            Stress<MethodSphere,Cauchy> hardyStress10(hardy10,&gridFromFile);
            Stress<MethodSphere,Cauchy> hardyStress15(hardy15,&gridFromFile);
            Stress<MethodSphere,Cauchy> hardyStress20(hardy20,&gridFromFile);

            // Calculate stress using the process_dedr, if possible
            calculateStress(body, kim,
                            std::tie(),
                            std::tie(hardyStress5, hardyStress10, hardyStress15, hardyStress20));
            hardyStress5.write("hardy5_" + modelname);
            hardyStress10.write("hardy10_" + modelname);
            hardyStress15.write("hardy15_" + modelname);
            hardyStress20.write("hardy20_" + modelname);
        }
        catch(const std::runtime_error& e){
            std::cout << e.what() << std::endl;
            std::cout << "Compute stress with process_dedr failed. Moving on" << std::endl;
        }
    }

	return 0;
}


