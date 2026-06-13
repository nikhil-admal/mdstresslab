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

/*! \brief Three-dimensional stress field on a grid.
 *
 * The stress field is computed using a prescribed averaging method.  For
 * `Cauchy` stresses, the object also stores continuum momentum density,
 * mass density, and the internal continuum velocity used to form the kinetic
 * stress from relative velocities.  These continuum fields are not defined for
 * `Piola` stresses; the kinetic contribution to Piola stress is taken to be
 * zero.
 *
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

	/*! \brief Cauchy-grid momentum density \f$\mathbf p(\mathbf x)\f$.
	 *
	 * Computed as \f$\sum_i m_i \mathbf v_i w(\mathbf x-\mathbf x_i)\f$ by
	 * `calculateKineticStress()`. Empty for Piola stress objects.
	 */
	std::vector<Vector3d> momentumDensityField;

	/*! \brief Cauchy-grid mass density \f$\rho(\mathbf x)\f$.
	 *
	 * Computed as \f$\sum_i m_i w(\mathbf x-\mathbf x_i)\f$ by
	 * `calculateKineticStress()`. Empty for Piola stress objects.
	 */
	std::vector<double> massDensityField;

	/*! \brief Internal Cauchy-grid continuum velocity.
	 *
	 * This is \f$\mathbf v(\mathbf x)=\mathbf p(\mathbf x)/\rho(\mathbf x)\f$
	 * where \f$\rho>0\f$, and zero otherwise.  It is used internally to form
	 * the kinetic Cauchy stress from relative velocities and is not written by
	 * `write()`.
	 */
	std::vector<Vector3d> velocityField;

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
		if constexpr (stressType==Cauchy)
		{
			momentumDensityField.resize(pgrid->ngrid,Vector3d::Zero());
			massDensityField.resize(pgrid->ngrid,0.0);
			velocityField.resize(pgrid->ngrid,Vector3d::Zero());
		}
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
		if constexpr (stressType==Cauchy)
		{
			momentumDensityField.resize(pgrid->ngrid,Vector3d::Zero());
			massDensityField.resize(pgrid->ngrid,0.0);
			velocityField.resize(pgrid->ngrid,Vector3d::Zero());
		}
	}

    /*!
     * \brief Write stress and, for Cauchy stress, density fields.
     *
     * This function writes the stress field to `[name].stress`. The stress file
     * is OVITO-readable and contains nine columns: grid coordinates followed by
     * \f$\sigma_{xx}\f$, \f$\sigma_{yy}\f$, \f$\sigma_{zz}\f$,
     * \f$\sigma_{xy}\f$, \f$\sigma_{xz}\f$, and \f$\sigma_{yz}\f$.
     *
     * For `Cauchy` stress objects, two additional OVITO-readable files are
     * written:
     * - `[name].momentum_density`, with grid coordinates and
     *   \f$\mathbf p(\mathbf x)\f$.
     * - `[name].mass_density`, with grid coordinates and
     *   \f$\rho(\mathbf x)\f$.
     *
     * The continuum velocity field is an internal intermediate and is not
     * written.
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

		if constexpr (stressType==Cauchy)
		{
			std::ofstream momentumDensityFile(name+".momentum_density");
			std::ofstream massDensityFile(name+".mass_density");
			momentumDensityFile << momentumDensityField.size() << "\n";
			massDensityFile << massDensityField.size() << "\n";
			momentumDensityFile << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
			massDensityFile << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
			momentumDensityFile << "Properties=pos:R:3:momentum_density:R:3" << std::endl;
			massDensityFile << "Properties=pos:R:3:mass_density:R:1" << std::endl;
			for (int i_grid=0; i_grid<pgrid->coordinates.size(); ++i_grid)
			{
				momentumDensityFile << pgrid->coordinates[i_grid].format(fmt)
									<< std::setw(25) << momentumDensityField[i_grid](0)
									<< std::setw(25) << momentumDensityField[i_grid](1)
									<< std::setw(25) << momentumDensityField[i_grid](2)
									<< std::endl;
				massDensityFile << pgrid->coordinates[i_grid].format(fmt)
								<< std::setw(25) << massDensityField[i_grid]
								<< std::endl;
			}
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
