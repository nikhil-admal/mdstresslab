/*
 * MethodLdad.h
 *
 *  Created on: Nov 5, 2019
 *  Modified on: Aug 27, 2022
 *      Author: Nikhil
 */

#ifndef METHODLDADBASE_H_
#define METHODLDADBASE_H_

#include "Method.h"
#include "typedef.h"
#include "Trigonometric.h"
#include "Constant.h"
#include <map>

template<typename T>
class MethodLdad : public Method<MethodLdad<T>>
{
    friend class Method<MethodLdad<T>>;
public:
	MethodLdad(const Matrix3d& ldadVectors);
	MethodLdad(const MethodLdad&);
	virtual ~MethodLdad();

protected:
	double operator()(const Vector3d& vec) const;
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2) const;

private:
	double normalizer;
	T oneDFunction;
	Matrix3d ldadVectors, inverseLdadVectors;
};

/*!
 * @example{lineno} testLDADLJ.cpp
 * -# Open the configuration file (in MDStressLab format) and read the number of particles
 * @snippet{lineno} testLDADLJ.cpp Read
 *
 * -# Initialize the BoxConfiguration and read the reference and current atomic coordinates
 *  from the configuration file. The reference configuration is an fcc Ar crystal in the relaxed state.
 *  The deformed/current configuration is the reference crystal strained in the \f$y\f$-direction.
 * @snippet{lineno} testLDADLJ.cpp Configuration
 *
 * -# Initialize the Kim model
 * @snippet{lineno} testLDADLJ.cpp Model
 *
 * -# Initialize reference and current 1D grids and read their
 *    coordinates from the respective grid files
 * @snippet{lineno} testLDADLJ.cpp Grid
 *
 * -# Initialize LDAD vectors and construct LDAD methods
 *    MethodLdadConstant and MethodLdadTrigonometric. The LDAD vectors (columns) are the three
 *    lattice vectors of the conventional unit cell.
 * @snippet{lineno} testLDADLJ.cpp LDAD
 *
 * -# Construct two Piola stress objects corresponding to the two methods
 * @snippet{lineno} testLDADLJ.cpp Stress
 *
 * -# Calculate the two Piola stress fields corresponding to the two methods
 * @snippet{lineno} testLDADLJ.cpp Calculate
 *
 * -# Repeat the previous three steps to calculate and output Cauchy stress fields. Note
 * that the LDAD stress was originally formulated as a Piola stress since the LDAD vectors should
 * be commensurate with a lattice, which is possible if the reference configuration is a single crystal.
 * Since the deformed configuration is also a single crystal in this example, we compute the LDAD Cauchy
 * stress using the appropriately deformed lattice vectors.
 * @snippet{lineno} testLDADLJ.cpp Calculate
 *
 * -# We compare our results with the exact results for unit testing purposes.
 * @snippet{lineno} testLDADLJ.cpp Compare
 *
 * Full code:
 */

/*!
 * @example TestLDADSW.cpp
 */


#include "MethodLdad.cpp"
// instantiate MethodLdad templates
//
template class MethodLdad<Constant>;
template class MethodLdad<Trigonometric>;
typedef MethodLdad<Constant> MethodLdadConstant;
typedef MethodLdad<Trigonometric> MethodLdadTrigonometric;

#endif /* METHODLDADBASE_H_ */
