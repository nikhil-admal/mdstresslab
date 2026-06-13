/*
 * Constant.h
 *
 *  Created on: Jan 9, 2020
 *      Author: Nikhil
 */

#ifndef CONSTANT_H_
#define CONSTANT_H_

#include "typedef.h"

class Constant {
public:
	Constant();
	virtual ~Constant();
	/*! \brief Evaluate the one-dimensional constant LDAD shape function. */
	double operator()(const double& t) const;
	/*! \brief Return \f$\int_{-1}^{1}\phi(t)\,dt=2\f$ for this shape. */
	double integral() const;
	/*! \brief Evaluate the segment integral used by the LDAD bond function. */
	double integrate(const Vector3d& vec1_pull_seg, const Vector3d& vec2_pull_seg) const;
};

#endif /* CONSTANT_H_ */
