/*
 * BoxConfiguration.cpp
 *
 *  Created on: Nov 4, 2019
 *      Author: Nikhil
 */

#include <fstream>
#include "BoxConfiguration.h"
#include "Configuration.h"
#include "typedef.h"
#include "neighbor_list.h"
#include "helper.hpp"


BoxConfiguration::BoxConfiguration(int numberOfParticles, int referenceAndFinal):
		Configuration(numberOfParticles,referenceAndFinal)
{
	MY_HEADING("Initializing a box configuration");
	if (referenceAndFinal)
		std::cout << "Creating a configuration of " << numberOfParticles << " particles" <<
					 " along with reference coordinates in a box of size zero" << std::endl;
	else
		std::cout << "Creating a configuration of " << numberOfParticles << " particles" <<
					 " in a box of size zero" << std::endl;
	// Default to zero box sizes and no pbc
	box.setZero();
	reference_box.setZero();
	pbc.setZero();
}


void BoxConfiguration::read(std::string configFileName, int referenceAndFinal)
{
	// Read in the atomistic system
	std::cout << "Reading the box configuration from file " << configFileName << std::endl;
	std::ifstream file(configFileName);
	if(!file)
	{
	// Print an error and exit
		std::cerr << "ERROR: config.dat could not be opened for reading!" << std::endl;
		exit(1);
	}

	int numberOfParticlesInFile;
	file >> numberOfParticlesInFile;
	if (numberOfParticles != numberOfParticlesInFile)
		MY_ERROR("Error: Number of particles in file does not equal to that of BoxConfiguration");

	for(int i=0;i<DIM*DIM;++i)
		if(!(file >> reference_box(i)))  MY_ERROR("ERROR: Reference box size.");
	for(int i=0;i<DIM*DIM;++i)
		if(!(file >> box(i))) 			 MY_ERROR("ERROR: Box size.");
	for(int i=0;i<DIM;++i)
		if(!(file >> pbc(i))) 			 MY_ERROR("ERROR: PBC.");


	std::string speciesTemp;
	for(int i=0;i<numberOfParticles;++i)
	{
		if(!(file >> speciesTemp)) 			 MY_ERROR("ERROR: Species code of particle " + std::to_string(i));
		species.push_back(speciesTemp);
		for(int j=0;j<DIM;++j)
			if(!(file >> coordinates[Current](i,j)))  MY_ERROR("ERROR: Coordinate of particle " + std::to_string(i));
		for(int j=0;j<DIM;++j)
			if(!(file >> velocities(i,j)))	 MY_ERROR("ERROR: Velocity of particle " + std::to_string(i));
		if (referenceAndFinal == true)
		{
			for(int j=0;j<DIM;++j)
				if(!(file >> coordinates[Reference](i,j)))
					MY_ERROR("ERROR: Reference coordinate of particle " + std::to_string(i) + "\n");
		}
		else
		{
			file.ignore(32767, '\n');
		}
	}
	std::cout << std::endl;
	std::cout << "Box size = " << std::endl;
	std::cout << box << std::endl;
	std::cout << std::endl;
	std::cout << "Reference box size = " << std::endl;
	std::cout << reference_box << std::endl;
	std::cout << std::endl;
	std::cout << "Periodic boundary conditions = " << pbc << std::endl;


}

Configuration* BoxConfiguration::getConfiguration(double padding) const
{
	// Build padding atoms
	int numberOfPaddings{0};
	std::vector<double> reference_coordinatesOfPaddings,coordinatesOfPaddings;
	std::vector<std::string> speciesOfPaddings;
	std::vector<int> masterOfPaddings;
	int referenceAndFinal= (coordinates.at(Reference).rows()>0);

	nbl_create_paddings(numberOfParticles,
						padding,
						reference_box.data(),
						box.data(),
						pbc.data(),
						coordinates.at(Reference).data(),
						coordinates.at(Current).data(),
						species,
						numberOfPaddings,
						reference_coordinatesOfPaddings,
						coordinatesOfPaddings,
						speciesOfPaddings,
						masterOfPaddings,
						referenceAndFinal);

	int total= numberOfParticles + numberOfPaddings;

	Configuration* config_ptr(new Configuration{total,referenceAndFinal});

	// copy the coordinates, particleContributing and species
	// of contributing atoms from BoxConfiguration to Configuration
	if (referenceAndFinal) (config_ptr->coordinates.at(Reference)).topRows(numberOfParticles)= coordinates.at(Reference);
	(config_ptr->coordinates.at(Current)).topRows(numberOfParticles)= coordinates.at(Current);
	for (auto it= species.begin();it!= species.end();it++)
		config_ptr->species.push_back(*it);

	if (numberOfPaddings)
	{
		if (referenceAndFinal) config_ptr->coordinates.at(Reference).bottomRows(numberOfPaddings)=
		*new Eigen::Map<MatrixXd> (reference_coordinatesOfPaddings.data(),numberOfPaddings,DIM);
		config_ptr->coordinates.at(Current).bottomRows(numberOfPaddings)=
		*new Eigen::Map<MatrixXd> (coordinatesOfPaddings.data(),numberOfPaddings,DIM);
		for (auto it= speciesOfPaddings.begin();it!= speciesOfPaddings.end();it++)
			config_ptr->species.push_back(*it);
	}

	return config_ptr;
}
BoxConfiguration::~BoxConfiguration() {
	// TODO Auto-generated destructor stub
}

