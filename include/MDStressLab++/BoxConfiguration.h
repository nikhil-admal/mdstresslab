/*
 * Material.h
 *
 *  Created on: Nov 4, 2019
 *      Author: Nikhil
 */

#ifndef ATOMISTICSYSTEM_H_
#define ATOMISTICSYSTEM_H_

#include "typedef.h"
#include "Configuration.h"

/*!
 * BoxConfiguration class describes the shape of boxes that
 * bound the current and reference configurations in addition to periodic
 * boundary conditions. Moreover, it includes member functions for reading
 * atomic configurations in either MDStressLab or Lammps dump file formats.
 */
class BoxConfiguration : public Configuration{
public:
    /*! \brief The current and reference box vectors stored as columns of respective
     * matrices.
     */
    Matrix3d box, reference_box;

    /*! \brief Periodic boundary conditions. \ref pbc=(1,0,1) implies periodicity along the \f$x\f$
     * and \f$z\f$-directions.
     */
	Vector3i pbc;

	Configuration* getConfiguration(double) const;
	void read(std::string,int);
    void readLMP(const std::string&,const ConfigType& configType);
    void readLMP(const std::string& currentConfigFileName,
                 const std::string& referenceConfigFileName);
    void lmpParser(std::ifstream&, const ConfigType&);

	BoxConfiguration(int,int);
	BoxConfiguration(std::string,int);
	virtual ~BoxConfiguration();
};

#endif /* ATOMISTICSYSTEM_H_ */
