/*
 * main.cpp
 *
 *  Created on: Nov 3, 2019
 *      Author: Nikhil Admal
 */
#include "Mls.h"
#include "MethodLdad.h"
#include "MethodSphere.h"
#include <string>
#include <iostream>
#include <tuple>
#include <cfloat>
#include <fstream>
#include "BoxConfiguration.h"
#include "calculateStress.h"
#include "Grid.h"
#include "typedef.h"
#include <regex>


int main()
{
    std::string line;
    auto next_line = [&]() -> std::string {
        while (std::getline(std::cin, line)) {
            if (line.empty() || line[0] == '#' || line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
            return line;
        }
        return "";  // EOF
    };

    // Read LAMMPS data file
    std::istringstream configFileStream(next_line());

    std::string currentConfigFileName, referenceConfigFileName;
    std::ifstream referenceConfigFile, currentConfigFile;
    if (configFileStream >> currentConfigFileName) {
        std::cout << "LAMMPS data file for current configuration: " << currentConfigFileName << std::endl;
        currentConfigFile.open(currentConfigFileName);
        if(!currentConfigFile) MY_ERROR("ERROR: current configuration file could not be opened for reading!");
    }
    if (configFileStream >> referenceConfigFileName) {
        std::cout << "LAMMPS data file for reference configuration: " << referenceConfigFileName << std::endl;
        referenceConfigFile.open(referenceConfigFileName);
        if(!referenceConfigFile) MY_ERROR("ERROR: reference configuration file could not be opened for reading!");
    }
    if(referenceConfigFileName.empty())
        referenceConfigFileName= currentConfigFileName;


    // get number of particles
    int numberOfParticles;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    while (std::getline(currentConfigFile, line)) {
        std::string loweredLine = line;
        std::transform(loweredLine.begin(), loweredLine.end(), loweredLine.begin(), ::tolower);

        if (loweredLine.find("atoms") != std::string::npos && (std::stringstream(line) >> numberOfParticles) ) {
        }
        else if (loweredLine.find("xlo xhi") != std::string::npos) {
            std::istringstream ss(line);
            ss >> xlo >> xhi;
        } else if (loweredLine.find("ylo yhi") != std::string::npos) {
            std::istringstream ss(line);
            ss >> ylo >> yhi;
        } else if (loweredLine.find("zlo zhi") != std::string::npos) {
            std::istringstream ss(line);
            ss >> zlo >> zhi;
            break;
        }
    }

    if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");

    BoxConfiguration body{numberOfParticles, 1};
    body.readLMP(currentConfigFileName,referenceConfigFileName);

    // Read periodic boundaries
    Vector3i pbc;
    std::istringstream pbStream(next_line());
    pbStream >> pbc(0) >> pbc(1) >> pbc(2);
    std::cout << "PBCs: " << pbc << std::endl;
    body.pbc= pbc;

    std::ofstream writeConfigReference("reference.txt");
    std::ofstream writeConfigCurrent("current.txt");
    for (int i=0; i<body.numberOfParticles; ++i)
    {
        writeConfigReference << body.species[i] << " " << body.coordinates[Reference].row(i) << std::endl;
        writeConfigCurrent << body.species[i] << " " << body.coordinates[Current].row(i) << std::endl;
    }

    // Read KIM ID
    std::string kimID = next_line();
    std::cout << "KIM ID: " << kimID << std::endl;
    Kim kim(kimID);

    // Read atom types
    std::istringstream atomStream(next_line());
    std::vector<std::string> atomTypes;
    std::string atom;
    while (atomStream >> atom) atomTypes.push_back(atom);
    std::cout << "Atom types: ";
    for (const auto& a : atomTypes) std::cout << a << " ";
    std::cout << std::endl;

    // Number of grids
    int numGrids = std::stoi(next_line());
    std::cout << "Number of grids: " << numGrids << "\n";
    Vector3d lowerLimit, upperLimit;
    int ngridx, ngridy, ngridz;
    double deltax, deltay, deltaz;
    std::string outPrefix;

    for (int i = 0; i < numGrids; ++i) {
        //-------- Read stress method
        MY_HEADING("Reading stress and grid parameters from input file");
        bool project= false;
        std::istringstream stressStream(next_line());
        std::string stressMethod, averagingDomain, averagingDomainParameter;
        double averagingDomainSize;

        stressStream >> stressMethod >> averagingDomain ;
        std::cout << "stress method: " << stressMethod << "\n";
        std::cout << "averaging domain: " << averagingDomain << "\n";
        if (stressMethod == "project")
            project= true;

        Vector3d xdir, ydir, zdir;
        Matrix3d ldadVectors;

        if (averagingDomain=="ldad")
        {
            stressStream >> averagingDomainParameter;
            std::cout << "averaging domain parameter: " << averagingDomainParameter << "\n";

            // read crystallographic directions and domain size to calculate ldadVectors
            if (averagingDomainParameter=="bcc" || averagingDomainParameter=="fcc")
            {
                stressStream >> averagingDomainSize >> xdir(0) >> xdir(1) >> xdir(2) >>
                                                       ydir(0) >> ydir(1) >> ydir(2) >>
                                                       zdir(0) >> zdir(1) >> zdir(2);
                std::cout << "averaging domain size : " << averagingDomainSize << "\n";

                Vector3d basis1, basis2, basis3;
                if (averagingDomainParameter=="bcc"){
                    basis1 << -0.5,0.5,0.5; basis2 << 0.5,-0.5,0.5; basis3 << 0.5,0.5,-0.5;
                }
                else if (averagingDomainParameter== "fcc"){
                    basis1 << 0,0.5,0.5; basis2 << 0.5,0,0.5; basis3 << 0.5,0.5,0;
                }
                Matrix3d rotation;

                // calculate lattice vectors
                rotation.row(0)= xdir.normalized();
                rotation.row(1)= ydir.normalized();
                rotation.row(2)= zdir.normalized();
                assert((rotation.transpose()*rotation).isIdentity());
                ldadVectors.col(0)= rotation*basis1.transpose();
                ldadVectors.col(1)= rotation*basis2.transpose();
                ldadVectors.col(2)= rotation*basis3.transpose();
                ldadVectors= ldadVectors* averagingDomainSize;
            }
            // instead directly read ldad lattice vectors
            else if (averagingDomainParameter=="lat")
                stressStream >> ldadVectors(0,0) >> ldadVectors(1,0) >> ldadVectors(2,0) >>
                                ldadVectors(0,1) >> ldadVectors(1,1) >> ldadVectors(2,1) >>
                                ldadVectors(0,2) >> ldadVectors(1,2) >> ldadVectors(2,2);
        }
        else if (averagingDomain=="sphere"){
                stressStream >> averagingDomainSize;
                std::cout << "averaging domain size : " << averagingDomainSize << "\n";
        }

        //------------ Read grid
        std::istringstream gridStream(next_line());
        std::string token[9];
        for (auto & j : token) {
            gridStream >> j;
        }
        // Lambda to convert token or use fallback
        auto parseOrFallback = [](const std::string& s, double fallback) -> double {
            return (s == "*") ? fallback : std::stod(s);
        };
        // Use fallback values if needed
        //lowerLimit[0]= parseOrFallback(token[0], xlo)-FLT_EPSILON;
        //lowerLimit[1]= parseOrFallback(token[1], ylo)-FLT_EPSILON;
        //lowerLimit[2]= parseOrFallback(token[2], zlo)-FLT_EPSILON;
        lowerLimit[0]= parseOrFallback(token[0], xlo);
        lowerLimit[1]= parseOrFallback(token[1], ylo);
        lowerLimit[2]= parseOrFallback(token[2], zlo);

        upperLimit[0]= parseOrFallback(token[3], xhi);
        upperLimit[1]= parseOrFallback(token[4], yhi);
        upperLimit[2]= parseOrFallback(token[5], zhi);

        // Parse delta values and diameter directly
        deltax = std::stod(token[6]);
        deltay = std::stod(token[7]);
        deltaz = std::stod(token[8]);

        ngridx= (abs(deltax) > FLT_EPSILON) ? floor((upperLimit(0)-lowerLimit(0)+FLT_EPSILON)/deltax) : 1;
        ngridy= (abs(deltay) > FLT_EPSILON) ? floor((upperLimit(1)-lowerLimit(1)+FLT_EPSILON)/deltay) : 1;
        ngridz= (abs(deltaz) > FLT_EPSILON) ? floor((upperLimit(2)-lowerLimit(2)+FLT_EPSILON)/deltaz) : 1;

        std::cout << "Grid Limits: ";
        std::cout << lowerLimit << " " << upperLimit << std::endl;
        std::cout << "Number of grid points: " << ngridx << " " << ngridy << " " << ngridz << std::endl;

        // ----------- Output prefix
        outPrefix= next_line();
        std::cout << "Output file: " << outPrefix << "\n";


        if (averagingDomain=="ldad") {
            try {
                Grid<Reference> grid(lowerLimit, upperLimit, ngridx, ngridy, ngridz);
                //Matrix3d ldadVectors= averagingDomainSize*Matrix3d::Identity();

                MethodLdadTrigonometric ldadDomain(ldadVectors);
                Stress<MethodLdadTrigonometric, Piola> ldadTrigonometricStress(ldadDomain, &grid);
                calculateStress(body, kim,
                                std::tie(ldadTrigonometricStress),
                                std::tie(), project);

                double mlsRadius= 10.0;
                Mls mls(body,&grid,mlsRadius,outPrefix);
                std::vector<Matrix3d> cauchyPushedField;
                mls.pushToCauchy(ldadTrigonometricStress.field,cauchyPushedField);
                mls.writePushedCauchyStress(cauchyPushedField);
                mls.writeDeformationGradient();
            }
            catch (const std::runtime_error &e) {
                std::cout << e.what() << std::endl;
                std::cout << "Compute stress with projected forces failed. Moving on" << std::endl;
            }
        }
        else if(averagingDomain=="sphere"){
            Grid<Current> grid(lowerLimit, upperLimit, ngridx, ngridy, ngridz);
            MethodSphere hardy(averagingDomainSize, "hardy");

            // stress calculation using projected forces
            try {
                Stress<MethodSphere, Cauchy> hardyStress(hardy, &grid);

                calculateStress(body, kim,
                                std::tie(),
                                std::tie(hardyStress), project);
                hardyStress.write(outPrefix);
            }
            catch (const std::runtime_error &e) {
                std::cout << e.what() << std::endl;
                std::cout << "Compute stress with projected forces failed. Moving on" << std::endl;
            }
        }
        else
            MY_ERROR("Unknown stress method: " + stressMethod);
    }

	return 0;
}


