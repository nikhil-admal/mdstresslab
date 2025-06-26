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

/*!Stress class template describes a three dimensional stress field
 * on a grid, computed using a prescribed averaging domain.
 * @tparam TMethod - Method template parameter. For example,
 *                   TMethod=MethodSphere for a spherical averaging domain
 *                   and TMethod=MethodLDAD for LDAD. For a user-defined averaging
 *                   domain, TMethod=MethodUser
 * @tparam stressType - Piola or Cauchy
 * @tparam TGrid - pointer to Grid<Reference> or Grid<Current> depending on whether
 *                 stressType is Piola or Cauchy, respectively
 *
 */
template<typename TMethod,
     StressType stressType,
     typename TGrid = typename std::conditional<stressType==Piola,Grid<Reference>,Grid<Current>>::type>
class Stress {
public:
/*!
 * \brief A three-dimensional stress field
 */
	std::vector<Matrix3d> field;

    /*!
     * \brief Pointer to the Grid on which the stress field is defined
     */
	TGrid* pgrid;

    /*!
     * \brief The method used to compute the stress field. The Method object
     * provides details about the weighting function and its support (averaging domain) and
     * includes the bond function to compute the stress field
     */
	const Method<TMethod>& method;

    /*!
     * \brief The prefix of the filename that will be outputted when the stress field
     * is written.
     */
	std::string name;

    /*!
     * \brief Constructs a Stress object
     */
	Stress(std::string name,
		   const Method<TMethod>& method,
		   TGrid* pgrid): name(name),pgrid(pgrid),method(method)
	{
		field.resize(pgrid->ngrid);
		for(auto& matrix : field)
			matrix= Matrix3d::Zero();
	}

    /*!
     * \brief Constructs a Stress object
     */
	Stress(const Method<TMethod>& method,
		   TGrid* pgrid): pgrid(pgrid),method(method)
	{
		field.resize(pgrid->ngrid);
		for(auto& matrix : field)
			matrix= Matrix3d::Zero();
	}

    /*!
     * This function writes the stress field to a filename with extension
     * .stress and prefix [name]. The output file is in a OVITO-readable format
     * with nine columns. The first three columns represent the three coordinates of
     * the grid points, and the last six columns are the six components of the stress
     * field - \f$\sigma_{xx}\f$, \f$\sigma_{yy}\f$, \f$\sigma_{zz}\f$,
     * \f$\sigma_{xy}\f$, \f$\sigma_{xz}\f$, and \f$\sigma_{yz}\f$
     */
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
