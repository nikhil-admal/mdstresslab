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
     * @param particleContributing - an integer array of size equal to the number of atoms. If a
     *                               particle is contributing it is marked as 1, and 0 otherwise.
     * @param forces_ptr - pointer to a matrix of size [3 x numberOfParticle] describing atomic forces
     * @param nl_ptr - pointer to the neighbor list
     * @param get_neigh_ptr - a pointer to a function that return the neighbor list of an atom
     * @param bonds - pointer to a InteratomicForces object
     * @param processDEDr_ptr - a pointer to the processdEdr function
     */
	void broadcastToModel(const Configuration* config_ptr,
						  const VectorXi& particleContributing,
                          const MatrixXd* forces_ptr,
						  NeighList* nl_ptr,
						  KIM::Function* get_neigh_ptr,
						  InteratomicForces* bonds,
						  KIM::Function* processDEDr_ptr);


    /*!
     * This function call the model's compute routine to calculate the atomic forces.
     * In addition, if \ref processDEDr_ptr!=nullptr,
     * interatomic forces are calculated using the model's processdEdr functionality.
     */
    void compute();

    const double* getCutoffs() const;
    int getNumberOfNeighborLists() const;
};
#endif /* KIM_H_ */
