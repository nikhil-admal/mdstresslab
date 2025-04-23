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
#include "Method.h"
#include <string>
#include <fstream>
#include <iostream>

template<typename TMethod,
		 StressType stressType,
		 typename TGrid = typename std::conditional<stressType==Piola,Grid<Reference>,Grid<Current>>::type>
class Stress {
public:
	std::vector<Matrix3d> field;
	TGrid* pgrid;
	const Method<TMethod>& method;
	std::string name;

	Stress(std::string name,
		   const Method<TMethod>& method,
		   TGrid* pgrid): name(name),pgrid(pgrid),method(method)
	{
		field.resize(pgrid->ngrid);
		for(auto& matrix : field)
			matrix= Matrix3d::Zero();
	}
	Stress(const Method<TMethod>& method,
		   TGrid* pgrid): pgrid(pgrid),method(method)
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
		//file << "\n";
        int index= 0;
        //Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "      ", "\n", "", "", "");
        Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "      ", "\n", "", "", "");
        file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        file << "Properties=pos:R:3:stress:R:6" << std::endl;
        for (auto& stress : field)
		{
			//Eigen::Map<Eigen::Matrix<double,1,DIM*DIM>> stressRow(stress.data(), stress.size());
			//file << pgrid->coordinates[index].format(fmt) << std::setw(5) << stressRow.format(fmt) << std::endl;
            file << pgrid->coordinates[index].format(fmt)
                << std::setw(25) << stress(0,0)
                << std::setw(25) << stress(1,1)
                << std::setw(25) << stress(2,2)
                << std::setw(25) << stress(0,1)
                << std::setw(25) << stress(0,2)
                << std::setw(25) << stress(1,2)
                << std::endl;
            index++;
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
