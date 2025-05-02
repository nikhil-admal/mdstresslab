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
    std::string configFileName= next_line();
    std::cout << "LAMMPS data file: " << configFileName << std::endl;
    std::ifstream configFile(configFileName);
    if(!configFile) MY_ERROR("ERROR: configuration file could not be opened for reading!");

    // get number of particles
    int numberOfParticles;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    while (std::getline(configFile, line)) {
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
    BoxConfiguration body{numberOfParticles,false};
    body.readLMP(configFileName);

    // Read periodic boundaries
    Vector3i pbc;
    std::istringstream pbStream(next_line());
    pbStream >> pbc(0) >> pbc(1) >> pbc(2);
    std::cout << "PBCs: " << pbc << std::endl;

    // Read KIM ID
    std::string kimID = next_line();
    std::cout << "KIM ID: " << kimID << "\n";
    Kim kim(kimID);

    // Read atom types
    std::istringstream atomStream(next_line());
    std::vector<std::string> atomTypes;
    std::string atom;
    while (atomStream >> atom) atomTypes.push_back(atom);
    std::cout << "Atom types: ";
    for (const auto& a : atomTypes) std::cout << a << " " << std::endl;

    // Read stress method
    bool project= false;
    std::string stressMethod = next_line();
    std::cout << "Stress Method: " << stressMethod << "\n";
    if (stressMethod == "project")
        project= true;

    // Number of grids
    int numGrids = std::stoi(next_line());
    std::cout << "Number of grids: " << numGrids << "\n";

    std::vector<Vector3d> lowerLimit(numGrids), upperLimit(numGrids);
    std::vector<int> ngridx(numGrids), ngridy(numGrids), ngridz(numGrids);
    std::vector<double> deltax(numGrids), deltay(numGrids), deltaz(numGrids);
    std::vector<double> d(numGrids);
    std::vector<std::string> outPrefix(numGrids);
    // Grid info
    for (int i = 0; i < numGrids; ++i) {
        std::istringstream gridStream(next_line());
        std::string token[10];
        for (auto & j : token) {
            gridStream >> j;
        }
        // Lambda to convert token or use fallback
        auto parseOrFallback = [](const std::string& s, double fallback) -> double {
            return (s == "*") ? fallback : std::stod(s);
        };

// Use fallback values if needed
        lowerLimit[i][0]= parseOrFallback(token[0], xlo);
        lowerLimit[i][1]= parseOrFallback(token[1], ylo);
        lowerLimit[i][2]= parseOrFallback(token[2], zlo);

        upperLimit[i][0]= parseOrFallback(token[3], xhi);
        upperLimit[i][1]= parseOrFallback(token[4], yhi);
        upperLimit[i][2]= parseOrFallback(token[5], zhi);

// Parse delta values and diameter directly
        deltax[i] = std::stod(token[6]);
        deltay[i] = std::stod(token[7]);
        deltaz[i] = std::stod(token[8]);
        d[i]      = std::stod(token[9]);

        ngridx[i]= (abs(deltax[i]) > FLT_EPSILON) ? floor((upperLimit[i](0)-lowerLimit[i](0))/deltax[i]) : 0;
        ngridy[i]= (abs(deltay[i]) > FLT_EPSILON) ? floor((upperLimit[i](1)-lowerLimit[i](1))/deltay[i]) : 1;
        ngridz[i]= (abs(deltaz[i]) > FLT_EPSILON) ? floor((upperLimit[i](2)-lowerLimit[i](2))/deltaz[i]) : 0;

        std::cout << "Grid " << i << ":" << std::endl;
        std::cout << "  Limits: ";
        std::cout << lowerLimit[i] << " " << upperLimit[i] << std::endl;
        std::cout << "  Grid points: " << ngridx[i] << " " << ngridy[i] << " " << ngridz[i] << std::endl;
        std::cout << "  Diameter: " << d[i] << std::endl;

        // Output prefix
        outPrefix[i]= next_line();
        std::cout << "Output prefix: " << outPrefix[i] << "\n";

    }

    for(int i=0; i<numGrids; ++i) {
        Grid<Current> grid(lowerLimit[i], upperLimit[i], ngridx[i], ngridy[i], ngridz[i]);
        MethodSphere hardy(d[i], "hardy");


        // stress calculation using projected forces
        try {
            Stress<MethodSphere, Cauchy> hardyStress(hardy, &grid);

            // Calculate stress using the projected forces
            calculateStress(body, kim,
                            std::tie(),
                            std::tie(hardyStress), project);
            hardyStress.write(outPrefix[i] + "_grid" + std::to_string(i));
        }
        catch (const std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            std::cout << "Compute stress with projected forces failed. Moving on" << std::endl;
        }
    }

	return 0;
}


