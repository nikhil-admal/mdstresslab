/*
 * test_SpatialHash.cpp
 *
 *  Created on: Nov 30, 2019
 *      Author: Nikhil
 */




#include "typedef.h"
#include "SpatialHash.h"
#include <iostream>
#include <fstream>
#include <iomanip>

int main()
{
	const int numberOfPoints=  5000;
	MatrixXd configuration(numberOfPoints,DIM);
	configuration= MatrixXd::Random(numberOfPoints,DIM);

	// Wrap the configuration in a std::vector<Vector3d> to check an alternate constructor
	std::vector<Vector3d> configurationVector;
	for (int i_point=0; i_point<numberOfPoints; i_point++)
	{
		Eigen::Map<Vector3d> vector(&configuration(i_point,0),3);
		configurationVector.push_back(vector);
	}


	Vector3d origin,step;
	origin.setConstant(0);
	step.setConstant(1.5);

	ConstSpatialHash hash(origin,step,configuration);

	std::ofstream file("random_hashed.data");
	file << numberOfPoints << std::endl;
	file << std::endl;
	int binNumber= 0;

	for(auto pair : hash.hashTable)
	{
		for (auto i_particle : pair.second)
		{
			file << std::setw(5) << binNumber << std::setw(5) << i_particle << std::setw(5) << pair.first << std::setw(15) << configuration.row(i_particle) << std::endl;
		}
		binNumber++;

	}

	ConstSpatialHash hash_vector(origin,step,configurationVector);
	std::ofstream file_map("random_hashed_vector.data");
	file_map << numberOfPoints << std::endl;
	file_map << std::endl;
	int binNumberVector= 0;

	for(auto pair : hash.hashTable)
	{
		for (auto i_particle : pair.second)
		{
			file_map << std::setw(5) << binNumberVector << std::setw(5) << i_particle << std::setw(5) << pair.first << std::setw(15) << configurationVector[i_particle] << std::endl;
		}
		binNumberVector++;

	}
}
