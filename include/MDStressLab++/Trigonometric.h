/*
 * Trigonometric.h
 *
 *  Created on: Jan 9, 2020
 *      Author: Nikhil
 */

#ifndef TRIGONOMETRIC_H_
#define TRIGONOMETRIC_H_

#include "typedef.h"

class Trigonometric{
public:
	Trigonometric();
	virtual ~Trigonometric();

	/*! \brief Evaluate \f$\phi(t)=\frac{1}{2}(1+\cos(\pi t))\f$ on \f$(-1,1)\f$. */
	double operator()(const double& t) const;
	/*! \brief Return \f$\int_{-1}^{1}\phi(t)\,dt=1\f$ for this shape. */
	double integral() const;
	/*! \brief Evaluate the expanded segment integral used by the LDAD bond function. */
	double integrate(const Vector3d& vec1_pull_seg, const Vector3d& vec2_pull_seg) const;
};



#endif /* TRIGONOMETRIC_H_ */
