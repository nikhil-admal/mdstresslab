/*
 * kim.h
 *
 *  Created on: Nov 18, 2019
 *      Author: Nikhil
 */

#ifndef KIM_H_
#define KIM_H_

#include "KIM_SimulatorHeaders.hpp"
#include "KIM_SupportedExtensions.hpp"
#include "BoxConfiguration.h"
#include "neighbor_list.h"
#include "InteratomicForces.h"

class Kim
{
public:
	KIM::Model* kim_ptr;
	KIM::ComputeArguments* computeArguments;
	double influenceDistance;
	std::vector<int> speciesCode;

	Kim(std::string);
	Kim(){kim_ptr=nullptr; computeArguments=nullptr; influenceDistance=0;}
	virtual ~Kim();
	void queryModel();
	void broadcastToModel(const Configuration* config_ptr,
						  const VectorXi& particleContributing,
						  NeighList* nl_ptr,
						  KIM::Function* get_neigh_ptr,
						  InteratomicForces* bonds,
						  KIM::Function* processDEDr_ptr);
	void compute();
	const double* getCutoffs() const;
	int getNumberOfNeighborLists() const;
	static void processDEDr(){};

};
#endif /* KIM_H_ */
