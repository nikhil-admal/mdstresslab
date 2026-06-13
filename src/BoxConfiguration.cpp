/*
 * BoxConfiguration.cpp
 *
 *  Created on: Nov 4, 2019
 *      Author: Nikhil
 */

#include <fstream>
#include <algorithm>
#include <sstream>
#include <unordered_map>
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
	box_origin.setZero();
	reference_box_origin.setZero();
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

	file.ignore(32767, '\n');
	auto nextDataLine = [&]() {
		std::string dataLine;
		do
		{
			if (!std::getline(file,dataLine))
				MY_ERROR("ERROR: Unexpected end of file while reading box configuration.");
		}
		while (dataLine.empty() ||
			   dataLine.find_first_not_of(" \t\r\n") == std::string::npos ||
			   dataLine[dataLine.find_first_not_of(" \t\r\n")] == '#');
		return dataLine;
	};

	std::string firstBoxLine= nextDataLine();
	std::istringstream keywordStream(firstBoxLine);
	std::string keyword;
	keywordStream >> keyword;
	std::transform(keyword.begin(),keyword.end(),keyword.begin(),::tolower);
	if (keyword=="origin")
	{
		Vector3d origin;
		if (!(keywordStream >> origin(0) >> origin(1) >> origin(2)))
			MY_ERROR("ERROR: Expected three coordinates after origin.");
		reference_box_origin= origin;
		box_origin= origin;
		firstBoxLine= nextDataLine();
	}

	std::istringstream firstBoxLineStream(firstBoxLine);
	auto readBoxScalar = [&](double& value, const std::string& errorMessage) {
		if (firstBoxLineStream >> value)
			return;
		if (!(file >> value))
			MY_ERROR(errorMessage);
	};

	for(int i=0;i<DIM*DIM;++i)
		readBoxScalar(reference_box(i),"ERROR: Reference box size.");
	for(int i=0;i<DIM*DIM;++i)
		readBoxScalar(box(i),"ERROR: Box size.");
	for(int i=0;i<DIM;++i)
	{
		if (firstBoxLineStream >> pbc(i))
			continue;
		if (!(file >> pbc(i))) 			 MY_ERROR("ERROR: PBC.");
	}

	file.ignore(32767, '\n');
	std::string speciesMassLine;
	do
	{
		if (!std::getline(file, speciesMassLine))
			MY_ERROR("ERROR: Expected species-mass line after PBC.");
	}
	while (speciesMassLine.empty() || speciesMassLine.find_first_not_of(" \t\r\n") == std::string::npos);

	std::map<std::string,double> speciesMasses;
	std::istringstream speciesMassStream(speciesMassLine);
	std::string speciesName;
	double speciesMass;
	while (speciesMassStream >> speciesName)
	{
		if (!(speciesMassStream >> speciesMass))
			MY_ERROR("ERROR: Expected mass after species " + speciesName + " in species-mass line.");
		speciesMasses[speciesName]= speciesMass;
	}
	if (speciesMasses.empty())
		MY_ERROR("ERROR: No species masses found after PBC.");

	std::string speciesTemp;
	for(int i=0;i<numberOfParticles;++i)
	{
		if(!(file >> speciesTemp)) 			 MY_ERROR("ERROR: Species code of particle " + std::to_string(i));
		if (speciesMasses.find(speciesTemp) == speciesMasses.end())
			MY_ERROR("ERROR: Species " + speciesTemp + " has no mass in the species-mass line.");
		species.push_back(speciesTemp);
		masses(i)= speciesMasses.at(speciesTemp);
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
    else if (configType==Reference)
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

void BoxConfiguration::readLMP(const std::string& currentConfigFileName,
                               const std::string& referenceConfigFileName){
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
    std::unordered_map<int, double> typeToMass;
    double xlo=0.0, xhi=0.0, ylo=0.0, yhi=0.0, zlo=0.0, zhi=0.0;
    double xy=0.0, xz=0.0, yz=0.0;
    auto updateBox = [&]() {
        double xloTrue= xlo - std::min({0.0,xy,xz,xy+xz});
        double xhiTrue= xhi - std::max({0.0,xy,xz,xy+xz});
        double yloTrue= ylo - std::min(0.0,yz);
        double yhiTrue= yhi - std::max(0.0,yz);
        Matrix3d lammpsBox= Matrix3d::Zero();
        lammpsBox.col(0)= Vector3d(xhiTrue-xloTrue,0.0,0.0);
        lammpsBox.col(1)= Vector3d(xy,yhiTrue-yloTrue,0.0);
        lammpsBox.col(2)= Vector3d(xz,yz,zhi-zlo);
        if (configType==Current)
        {
            box= lammpsBox;
            box_origin= Vector3d(xloTrue,yloTrue,zlo);
        }
        else
        {
            reference_box= lammpsBox;
            reference_box_origin= Vector3d(xloTrue,yloTrue,zlo);
        }
    };

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
            ss >> xlo >> xhi;
            updateBox();
        } else if (loweredLine.find("ylo yhi") != std::string::npos) {
            std::istringstream ss(line);
            ss >> ylo >> yhi;
            updateBox();
        } else if (loweredLine.find("zlo zhi") != std::string::npos) {
            std::istringstream ss(line);
            ss >> zlo >> zhi;
            updateBox();
        } else if (loweredLine.find("xy xz yz") != std::string::npos) {
            std::istringstream ss(line);
            ss >> xy >> xz >> yz;
            updateBox();
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
                typeToMass[type] = mass;

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

            if(configType==Current) species.resize(numberOfParticles);
            // Read atom lines
            for (int i = 0; i < numberOfParticles; ++i) {
                std::istringstream ss(line); // first valid line
                int id, type;
                double x, y, z;
                if (!(ss >> id >> type >> x >> y >> z))
                    MY_ERROR("ERROR: Coordinate of particle " + std::to_string(i));
                if (typeToSpecies.find(type) == typeToSpecies.end() ||
                    typeToMass.find(type) == typeToMass.end())
                    MY_ERROR("ERROR: Atom type " + std::to_string(type) + " is missing from the Masses section.");

                int idxFlagx= 0; int idxFlagy= 0; int idxFlagz= 0;
                int tmpx, tmpy, tmpz;
                if (ss >> tmpx >> tmpy >> tmpz) {
                    idxFlagx = tmpx;
                    idxFlagy = tmpy;
                    idxFlagz = tmpz;
                }

                //if(configType==Current) species.push_back(typeToSpecies[type]);
                if(configType==Current) {
                    species[id-1]= typeToSpecies.at(type);
                    masses(id-1)= typeToMass.at(type);
                }
                if(configType==Reference) {
                    assert(species[id-1] == typeToSpecies.at(type) &&
                           "Species in the reference configuration do not match with those in the "
                           "current configuration");
                    assert(masses(id-1) == typeToMass.at(type) &&
                           "Masses in the reference configuration do not match with those in the "
                           "current configuration");
                }
                const Matrix3d& boxMatrix= (configType==Reference) ? reference_box : box;
                Vector3d imageShift= idxFlagx*boxMatrix.col(0).transpose() +
                                     idxFlagy*boxMatrix.col(1).transpose() +
                                     idxFlagz*boxMatrix.col(2).transpose();
                coordinates[configType].row(id - 1)= Vector3d(x,y,z) + imageShift;

                if (i < numberOfParticles - 1) {
                    std::getline(file, line); // read next line
                }
            }

        }

            // Process Velocities section
        else if (loweredLine.find("velocities") != std::string::npos) {
            // Skip lines until we reach actual velocity data
            while (std::getline(file, line)) {
                if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos)
                    continue;
                if (line[0] == '#')
                    continue;
                break; // first data line found
            }

            for (int i = 0; i < numberOfParticles; ++i) {
                std::istringstream ss(line);
                int id;
                double vx, vy, vz;
                if (!(ss >> id >> vx >> vy >> vz))
                    MY_ERROR("ERROR: Velocity of particle " + std::to_string(i));

                if(configType==Current) {
                    velocities(id - 1, 0) = vx;
                    velocities(id - 1, 1) = vy;
                    velocities(id - 1, 2) = vz;
                }

                if (i < numberOfParticles - 1) {
                    std::getline(file, line);
                }
            }
        }
    }

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
						reference_box_origin.data(),
						box_origin.data(),
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
	config_ptr->velocities.topRows(numberOfParticles)= velocities;
	config_ptr->masses.head(numberOfParticles)= masses;
	for (auto it= species.begin();it!= species.end();it++)
		config_ptr->species.push_back(*it);

	if (numberOfPaddings)
	{
		using RowMajorMatrixXd = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
		if (referenceAndFinal)
		{
			config_ptr->coordinates.at(Reference).bottomRows(numberOfPaddings)=
					Eigen::Map<const RowMajorMatrixXd> (reference_coordinatesOfPaddings.data(),numberOfPaddings,DIM);
		}
		config_ptr->coordinates.at(Current).bottomRows(numberOfPaddings)=
				Eigen::Map<const RowMajorMatrixXd> (coordinatesOfPaddings.data(),numberOfPaddings,DIM);
		for (int i_padding=0; i_padding<numberOfPaddings; ++i_padding)
		{
			int master= masterOfPaddings[i_padding];
			config_ptr->velocities.row(numberOfParticles+i_padding)= velocities.row(master);
			config_ptr->masses(numberOfParticles+i_padding)= masses(master);
		}
		for (auto it= speciesOfPaddings.begin();it!= speciesOfPaddings.end();it++)
			config_ptr->species.push_back(*it);
	}

	return config_ptr;
}
BoxConfiguration::~BoxConfiguration() {
	// TODO Auto-generated destructor stub
}
