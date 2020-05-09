/*
 * Stencil.cpp
 *
 *  Created on: Dec 13, 2019
 *      Author: Nikhil
 */

#include "Stencil.h"

Stencil::Stencil(const Configuration& parent) : parent(parent) { }

void Stencil::emptyStencil()
{
	particleContributingMap.clear();
}
Stencil::~Stencil() {
	// TODO Auto-generated destructor stub
}

