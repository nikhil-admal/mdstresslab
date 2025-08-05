/*
 * main.cpp
 *
 *  Created on: Aug 5, 2025
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
  std::vector<std::string> modelnames={
    "EAM_Dynamo_RaveloGermannGuerrero_2013Ta1_Ta__MO_816821594689_001"
  };

  /*![ReadConfiguration]*/
  int numberOfParticles;
  int referenceAndFinal= true;
  std::string configFileName= "output.xyz";
  std::ifstream file(configFileName);
  if(!file) MY_ERROR("ERROR: output.xyz could not be opened for reading!");

  file >> numberOfParticles;
  if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");


  BoxConfiguration body{numberOfParticles,referenceAndFinal};
  body.read(configFileName,referenceAndFinal);
  /*![ReadConfiguration]*/

  //	-------------------------------------------------------------------
  // Create grid
  //	-------------------------------------------------------------------
  int nx,ny,nz;
  nx = 80;
  ny = 1;
  nz = 80;  
  Vector3d lowerLimit(0.,0.,0.);
  Vector3d upperLimit(217.8,16.5,198.0);
  Grid<Current> gridFromFile(lowerLimit,upperLimit,nx,ny,nz);
  
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
        hardyStress.write_voxel_grid("hardy_" + modelname,nx,ny,nz,lowerLimit,upperLimit);
      }
      catch(const std::runtime_error& e){
        std::cout << e.what() << std::endl;
        std::cout << "Compute stress with process_dedr failed. Moving on" << std::endl;
      }
    }
  /*![ComputeStress]*/

	return 0;
}


