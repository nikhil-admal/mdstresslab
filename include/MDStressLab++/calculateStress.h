/*
 * calculateStress.h
 *
 *  Created on: Nov 7, 2019
 *      Author: Nikhil
 */

#ifndef CALCULATESTRESS_H_
#define CALCULATESTRESS_H_

int process_DEDr(const void* dataObject, const double de, const double r, const double* const dx, const int i, const int j);
#include "calculateStress.cpp"



#endif /* CALCULATESTRESS_H_ */
