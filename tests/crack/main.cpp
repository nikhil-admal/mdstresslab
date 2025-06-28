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
#include "kim_query.cpp"
#include <regex>


int main()
{
    /*![QueryModels]*/
    std::vector<std::string> modelnames={
            "Tersoff_LAMMPS_Tersoff_1988T3_Si__MO_186459956893_004"
            };

    std::cout << "Querying models." << std::endl;
    for (int i = 0; i < modelnames.size(); ++i) {
        std::cout << "Model: " << modelnames[i] << std::endl;
        std::regex rgx("MO_\\d+_\\d{3}");

        std::smatch match;
        if (!std::regex_search(modelnames[i], match, rgx))
            std::cout << "Regex error!" << std::endl;


        std::string latticeConstant= kim_query(
                "get_lattice_constant_cubic",
                {
                        "model=[\"" + match.str()+ "\"]",
                        "crystal=[\"diamond\"]",
                        "species=[\"Si\"]",
                        "units=[\"angstrom\"]"
                }
        );
        std::cout << "Lattice constant = " << latticeConstant << std::endl;
        std::string elasticConstants= kim_query(
                "get_elastic_constants_isothermal_cubic",
                {
                        "model=[\"" + match.str()+ "\"]",
                        "crystal=[\"diamond\"]",
                        "species=[\"Si\"]",
                        "units=[\"eV/angstrom^3\"]"
                }
        );
        std::cout << "Elastic constants = " << elasticConstants << std::endl;
        std::cout << "------------------------------------------ " << std::endl;
    }
    /*![QueryModels]*/

    /*![ReadConfiguration]*/
    int numberOfParticles;
    int referenceAndFinal= true;
    std::string configFileName= "config.data";
    std::ifstream file(configFileName);
    if(!file) MY_ERROR("ERROR: config.dat could not be opened for reading!");

    file >> numberOfParticles;
    if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");

    BoxConfiguration body{numberOfParticles,referenceAndFinal};
    body.read(configFileName,referenceAndFinal);
    /*![ReadConfiguration]*/

    /*![LoadGrid]*/
    int ngrid;
    ngrid= 90601;
    Grid<Current> gridFromFile(ngrid);
    gridFromFile.read("grid_cauchy.data");
    /*![LoadGrid]*/

    /*![ComputeStress]*/
    MethodSphere hardy(6,"hardy");
    for (const auto modelname : modelnames)
    {
        Kim kim(modelname);
        try {
            Stress<MethodSphere,Cauchy> hardyStress(hardy,&gridFromFile);

            calculateStress(body, kim,
                            std::tie(),
                            std::tie(hardyStress), true);
            hardyStress.write("project_hardy_" + modelname);
        }
        catch(const std::runtime_error& e){
            std::cout << e.what() << std::endl;
            std::cout << "Compute stress with projected forces failed. Moving on" << std::endl;
        }

        try{
            Stress<MethodSphere,Cauchy> hardyStress(hardy,&gridFromFile);

            // Calculate stress using the process_dedr, if possible
            calculateStress(body, kim,
                            std::tie(),
                            std::tie(hardyStress));
            hardyStress.write("hardy_" + modelname);
        }
        catch(const std::runtime_error& e){
            std::cout << e.what() << std::endl;
            std::cout << "Compute stress with process_dedr failed. Moving on" << std::endl;
        }
    }
    /*![ComputeStress]*/

	return 0;
}


