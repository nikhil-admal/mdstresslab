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
#include <utility>

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
     *
     * Use `write_voxel_grid()` for optional structured-grid output in
     * LAMMPS dump-grid format.
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

	/*!
	 * \brief Write structured grid fields in LAMMPS dump-grid format.
	 *
	 * This output is intended for direct visualization of structured stress
	 * grids in OVITO. The caller must provide the grid dimensions and
	 * orthogonal bounding box used to create the grid. The grid dimensions are
	 * checked against `field.size()` before writing.
	 *
	 * The dump-grid `DIMENSION` is inferred from the number of grid counts
	 * greater than one. Bounds are shifted by half a grid spacing along active
	 * directions so that voxel cell centers coincide with the MDStressLab grid
	 * coordinates; inactive directions are written as a collapsed plane.
	 *
	 * The stress tensor is written to `[name].voxel_grid_stress` with component
	 * order `SXX SYY SZZ SYZ SXZ SXY`, matching LAMMPS/OVITO dump-grid
	 * conventions. For `Cauchy` stress objects, this function also writes
	 * `[name].voxel_grid_momentum_density` and `[name].voxel_grid_mass_density`.
	 *
	 * The values are written in the same order as `pgrid->coordinates`.
	 */
	void write_voxel_grid(const int nx,
	                      const int ny,
	                      const int nz,
	                      const Vector3d& lowerLimit,
	                      const Vector3d& upperLimit)
	{
        if (name.empty())
            MY_ERROR("Stress object created without specifying a name. Use write_voxel_grid(filename,...) instead of write_voxel_grid(...)");
		validateVoxelGridDimensions(nx,ny,nz);

		std::ofstream stressFile(name+".voxel_grid_stress");
		writeVoxelGridHeader(stressFile,nx,ny,nz,lowerLimit,upperLimit,"SXX SYY SZZ SYZ SXZ SXY");
		for (const auto& stress : field)
		{
			stressFile << std::setw(25) << stress(0,0)
			           << std::setw(25) << stress(1,1)
			           << std::setw(25) << stress(2,2)
			           << std::setw(25) << stress(1,2)
			           << std::setw(25) << stress(0,2)
			           << std::setw(25) << stress(0,1)
			           << std::endl;
		}

		if constexpr (stressType==Cauchy)
		{
			std::ofstream momentumDensityFile(name+".voxel_grid_momentum_density");
			writeVoxelGridHeader(momentumDensityFile,nx,ny,nz,lowerLimit,upperLimit,"PX PY PZ");
			for (const auto& momentumDensity : momentumDensityField)
			{
				momentumDensityFile << std::setw(25) << momentumDensity(0)
				                    << std::setw(25) << momentumDensity(1)
				                    << std::setw(25) << momentumDensity(2)
				                    << std::endl;
			}

			std::ofstream massDensityFile(name+".voxel_grid_mass_density");
			writeVoxelGridHeader(massDensityFile,nx,ny,nz,lowerLimit,upperLimit,"RHO");
			for (const auto& massDensity : massDensityField)
				massDensityFile << std::setw(25) << massDensity << std::endl;
		}
	}

	void write_voxel_grid(const std::string& filename,
	                      const int nx,
	                      const int ny,
	                      const int nz,
	                      const Vector3d& lowerLimit,
	                      const Vector3d& upperLimit)
	{
        if (name.empty())
            name= filename;
        else
            std::cout << "Stress object created with name " << name << ". Ignoring the filename: " << filename << "." << std::endl;
        write_voxel_grid(nx,ny,nz,lowerLimit,upperLimit);
	}

	~Stress()
	{
		// TODO Auto-generated destructor stub
	}

private:
	void validateVoxelGridDimensions(const int nx,
	                                 const int ny,
	                                 const int nz) const
	{
		if (nx<=0 || ny<=0 || nz<=0)
			MY_ERROR("Voxel grid dimensions must be positive.");
		const auto numberOfGridPoints= static_cast<std::size_t>(nx)*
		                               static_cast<std::size_t>(ny)*
		                               static_cast<std::size_t>(nz);
		if (numberOfGridPoints != field.size())
			MY_ERROR("Voxel grid dimensions do not match the number of stress grid points.");
	}

	void writeVoxelGridHeader(std::ofstream& file,
	                          const int nx,
	                          const int ny,
	                          const int nz,
	                          const Vector3d& lowerLimit,
	                          const Vector3d& upperLimit,
	                          const std::string& columns) const
	{
		file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
		file << "ITEM: TIMESTEP\n";
		file << "0\n";
		const auto xBounds= voxelGridBounds(nx,lowerLimit(0),upperLimit(0));
		const auto yBounds= voxelGridBounds(ny,lowerLimit(1),upperLimit(1));
		const auto zBounds= voxelGridBounds(nz,lowerLimit(2),upperLimit(2));
		file << "ITEM: BOX BOUNDS pp pp pp\n";
		file << xBounds.first << " " << xBounds.second << "\n";
		file << yBounds.first << " " << yBounds.second << "\n";
		file << zBounds.first << " " << zBounds.second << "\n";
		file << "ITEM: DIMENSION\n";
		file << voxelGridDimension(nx,ny,nz) << "\n";
		file << "ITEM: GRID SIZE nx ny nz\n";
		file << nx << " " << ny << " " << nz << "\n";
		file << "ITEM: GRID CELLS " << columns << "\n";
	}

	int voxelGridDimension(const int nx,
	                       const int ny,
	                       const int nz) const
	{
		int dimension= 0;
		if (nx>1) ++dimension;
		if (ny>1) ++dimension;
		if (nz>1) ++dimension;
		return std::max(1,dimension);
	}

	std::pair<double,double> voxelGridBounds(const int n,
	                                        const double lowerLimit,
	                                        const double upperLimit) const
	{
		if (n>1)
		{
			const double spacing= (upperLimit-lowerLimit)/static_cast<double>(n);
			return {lowerLimit-0.5*spacing,upperLimit-0.5*spacing};
		}
		return {lowerLimit,lowerLimit};
	}


};

#endif /* STRESS_H_ */
