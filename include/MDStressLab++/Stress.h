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
	TStress weightFunction;
	std::string name;

	Stress(std::string name,
		   const TStress& weightFunction,
		   TGrid* pgrid): name(name),pgrid(pgrid),weightFunction(weightFunction)
	{
		field.resize(pgrid->ngrid);
		for(auto& matrix : field)
			matrix= Matrix3d::Zero();
	}
	Stress(const TStress& weightFunction,
		   TGrid* pgrid): pgrid(pgrid),weightFunction(weightFunction)
	{
		field.resize(pgrid->ngrid);
		for(auto& matrix : field)
			matrix= Matrix3d::Zero();
	}

	void write()
	{
        if (name.empty())
            MY_ERROR("Stress object created without specifying a name. Use write(filename) instead of write()");
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

    void write(const std::string& filename)
    {
        if (name.empty())
            name= filename;
        else
            std::cout << "Stress object created with name " << name << ". Ignoring the filename: " << filename << "." << std::endl; 
        write();
    }

	~Stress()
	{
		// TODO Auto-generated destructor stub
	}


};

#endif /* STRESS_H_ */
