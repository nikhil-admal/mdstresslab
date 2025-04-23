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
		std::cerr << "ERROR: " << configFileName << " could not be opened for reading!" << std::endl;
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

void BoxConfiguration::readLMP(const std::string& configFileName,
                               const ConfigType& configType){
    // Read in the atomistic system
    if (configType==Current)
        std::cout << "Reading the current box configuration from lammps data file " << configFileName << std::endl;
    else if (configType==Current)
        std::cout << "Reading the reference box configuration from lammps data file " << configFileName << std::endl;
    std::ifstream file(configFileName);
    if(!file)
    {
        std::cerr << "ERROR: " << configFileName << " could not be opened for reading!" << std::endl;
        exit(1);
    }

    // Lambda version of hasEnding
    auto hasEnding = [](const std::string& fullString, const std::string& ending) -> bool {
        return fullString.size() >= ending.size() &&
               fullString.compare(fullString.size() - ending.size(), ending.size(), ending) == 0;
    };

    if (hasEnding(configFileName, ".lmp")) {
        lmpParser(file,configType);
    } else {
        std::cerr << "ERROR: Expecting file with extension .lmp!" << std::endl;
        exit(1);
    }
}

void BoxConfiguration::readLMP(const std::string& referenceConfigFileName,
                               const std::string& currentConfigFileName){
    readLMP(currentConfigFileName,Current);
    if(coordinates[Reference].rows()>0)
        readLMP(referenceConfigFileName,Reference);
    else
        MY_ERROR("Error: Memory not assigned to store reference configuration.");
}

void BoxConfiguration::lmpParser(std::ifstream& file, const ConfigType& configType)
{
    std::string line;
    int numAtoms = 0;
    int numberOfAtomTypes = 0;
    std::unordered_map<int, std::string> typeToSpecies;

    while (std::getline(file, line)) {
        // Normalize to lowercase for keyword checks (optional but helpful)
        std::string loweredLine = line;
        std::transform(loweredLine.begin(), loweredLine.end(), loweredLine.begin(), ::tolower);

        // Parse total number of atoms
        if (loweredLine.find("atoms") != std::string::npos && (std::stringstream(line) >> numAtoms) ) {
            //std::istringstream ss(line);
            //ss >> numAtoms;
            if (numberOfParticles != numAtoms)
                MY_ERROR("Error: Number of particles in file does not equal to that of BoxConfiguration");
        }
            // Parse number of atom types
        else if (loweredLine.find("atom types") != std::string::npos) {
            std::istringstream ss(line);
            ss >> numberOfAtomTypes;
        }
            // Parse box dimensions
        else if (loweredLine.find("xlo xhi") != std::string::npos) {
            std::istringstream ss(line);
            double xlo, xhi;
            ss >> xlo >> xhi;
            box(0, 0) = xhi - xlo;
        } else if (loweredLine.find("ylo yhi") != std::string::npos) {
            std::istringstream ss(line);
            double ylo, yhi;
            ss >> ylo >> yhi;
            box(1, 1) = yhi - ylo;
        } else if (loweredLine.find("zlo zhi") != std::string::npos) {
            std::istringstream ss(line);
            double zlo, zhi;
            ss >> zlo >> zhi;
            box(2, 2) = zhi - zlo;
        }

            // Process Masses section
        else if (loweredLine.find("masses") != std::string::npos) {
            int massLinesRead = 0;

            // Read lines until we've read all the atom types
            while (std::getline(file, line)) {
                // Skip blank or whitespace-only lines
                if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos)
                    continue;

                // Skip comment lines
                if (line[0] == '#')
                    continue;

                // Read: <type> <mass> # optional comment with species name
                std::istringstream ss(line);
                int type;
                double mass;
                std::string comment;

                ss >> type >> mass;
                std::getline(ss, comment); // grab remainder of line (comment)

                std::string speciesName = "Unknown";

                // Extract species name from comment if present
                size_t hashPos = comment.find('#');
                if (hashPos != std::string::npos) {
                    speciesName = comment.substr(hashPos + 1);
                    // Trim whitespace
                    speciesName.erase(0, speciesName.find_first_not_of(" \t"));
                    speciesName.erase(speciesName.find_last_not_of(" \t\r\n") + 1);
                }

                typeToSpecies[type] = speciesName;

                if (++massLinesRead >= numberOfAtomTypes)
                    break; // We've read all the expected mass lines
            }
        }

            // Process Atoms section
        else if (loweredLine.find("atoms") != std::string::npos) {
            // Skip lines until we reach actual atom data
            while (std::getline(file, line)) {
                if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos)
                    continue;
                if (line[0] == '#')
                    continue;
                break; // first data line found
            }

            // Read atom lines
            for (int i = 0; i < numberOfParticles; ++i) {
                std::istringstream ss(line); // first valid line
                int id, type;
                double x, y, z;
                if (!(ss >> id >> type >> x >> y >> z))
                    MY_ERROR("ERROR: Coordinate of particle " + std::to_string(i));

                species.push_back(typeToSpecies[type]);
                coordinates[configType](i, 0) = x;
                coordinates[configType](i, 1) = y;
                coordinates[configType](i, 2) = z;

                if (i < numberOfParticles - 1) {
                    std::getline(file, line); // read next line
                }
            }

            break; // done reading atom data
        }
    }

    // Set reference box same as current box
    reference_box = box;

    // Assume no periodicity for now
    pbc = Eigen::Vector3i::Zero();
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

