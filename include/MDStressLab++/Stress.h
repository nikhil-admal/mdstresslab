/*
 * Stress.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef STRESS_H_
#define STRESS_H_

#include <vector>
#include "Grid.h"
#include "typedef.h"
#include "SpatialHash.h"
#include <string>
#include <fstream>
#include <iostream>

template<typename TStress,
		 StressType stressType,
		 typename TGrid = typename std::conditional<stressType==Piola,Grid<Reference>,Grid<Current>>::type>
class Stress {
public:
	std::vector<Matrix3d> field;
	TGrid* pgrid;
	TStress bondFunction;
	std::string name;

	Stress(std::string name,
		   const TStress& bondFunction,
		   TGrid* pgrid): name(name),pgrid(pgrid),bondFunction(bondFunction)
	{
		field.resize(pgrid->ngrid);
		for(auto& matrix : field)
			matrix= Matrix3d::Zero();
	}

	void write()
	{
		std::ofstream file(name+".stress");

		file << field.size() << "\n";
		file << "\n";
		for (auto& stress : field)
		{
			Eigen::Map<Eigen::Matrix<double,1,DIM*DIM>> stressRow(stress.data(), stress.size());
			Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "      ", "\n", "", "", "", "");
			file << stressRow.format(fmt) << std::endl;
		}
	}

	~Stress()
	{
		// TODO Auto-generated destructor stub
	}


};

#endif /* STRESS_H_ */
