/*
 * MethodHardySphere.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef METHODHARDYSPHERE_H_
#define METHODHARDYSPHERE_H_

#include "Sphere.h"

class MethodHardySphere : public Sphere
{
public:
	MethodHardySphere(double);
	MethodHardySphere(const MethodHardySphere&);
	virtual ~MethodHardySphere();
};

#endif /* METHODHARDYSPHERE_H_ */
