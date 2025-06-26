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

/*!
 * The Kim class links MDStressLab to an openkim interatomic potential model.
 */
class Kim
{
public:
    std::vector<int> speciesCode;
    std::string modelname;
    double influenceDistance;
	KIM::Model* kim_ptr;
	KIM::ComputeArguments* computeArguments;

    /*!
     * Constructs the Kim object
     * @param modelname - KIM model name
     */
	Kim(const std::string& modelname);
	Kim(){kim_ptr=nullptr; computeArguments=nullptr; influenceDistance=0;}
	virtual ~Kim();
	void queryModel();


    /*!
     * This function broadcasts the atomistic system to the KIM model
     * @param config_ptr - pointer to the configuration of atoms
     * @param particleContributing -
     * @param forces_ptr -
     * @param nl_ptr -
     * @param get_neigh_ptr -
     * @param bonds -
     * @param processDEDr_ptr -
     */
	void broadcastToModel(const Configuration* config_ptr,
						  const VectorXi& particleContributing,
                          const MatrixXd* forces_ptr,
						  NeighList* nl_ptr,
						  KIM::Function* get_neigh_ptr,
						  InteratomicForces* bonds,
						  KIM::Function* processDEDr_ptr);
	void compute();
	const double* getCutoffs() const;
	int getNumberOfNeighborLists() const;
};
#endif /* KIM_H_ */
