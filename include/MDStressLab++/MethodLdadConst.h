/*
 * MethodLdadConst.h
 *
 *  Created on: Aug 22, 2022
 *      Author: Nikhil
 */

#ifndef METHODLDADCONST_H_
#define METHODLDADCONST_H_
#include "Ldad.h"
#include "Constant.h"

class MethodLdadConst : public Ldad<Constant>
{
public:
    MethodLdadConst(const Matrix3d& ldadVectors) : Ldad(ldadVectors)
    {
    }
	MethodLdadConst(const MethodLdadConst& method) : Ldad(method)
    {
        *this= method;
    }
	virtual ~MethodLdadConst(){}
};
#endif /* METHODLDADCONST_H_ */
