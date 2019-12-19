/*
 * grid_test.cpp
 *
 *  Created on: Nov 29, 2019
 *      Author: Nikhil
 */
#include "Grid.h"
#include <iostream>
#include "typedef.h"
#include <fstream>

int main()
{
	Vector3d lowerLimit(1,2,3);
	Vector3d upperLimit(2,4,6);
	int ngrid= 1000;
	Grid<Reference> grid1(lowerLimit,upperLimit,ngrid);

	// write to file
	std::ofstream output_file("random.dat");
	output_file << ngrid << std::endl;
	for (const auto& position : grid1.coordinates)
		output_file << position.transpose() << std::endl;


	// read from file
	Grid<Current> grid2(ngrid);
	grid2.read("random.dat");

	std::ofstream read_file("random_read.dat");
	read_file << ngrid << std::endl;
	for (const auto& position : grid2.coordinates)
		read_file << position.transpose() << std::endl;
}





