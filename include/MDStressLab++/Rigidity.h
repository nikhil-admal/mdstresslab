//
// Created by Nikhil Chandra Admal on 1/29/25.
//

#ifndef MDSTRESSLAB_RIGIDITY_H
#define MDSTRESSLAB_RIGIDITY_H

#include "Eigen/Eigen/Sparse"
#include "neighbor_list.h"

// input: coordinates, input_Nl.
// private: halfNl, rigidity matrices
// Output bonds;
class Rigidity {
public:
    Eigen::SparseMatrix<double> matrixR, matrixRTR;
    const MatrixXd& coordinates;
    const NeighListOne* p_inputNl;
    Rigidity(const MatrixXd & coordinates, const NeighListOne* p_inputNl);
    std::vector<double> project(const Eigen::VectorXd& gi) const;
    ~Rigidity();
private:
    NeighListOne halfNl;
    //NeighListOne static getHalfNeighborList(const NeighListOne* p_inputNl);
    NeighListOne static getHalfNeighborList(const MatrixXd& coordinates, const double& cutoff);
};


#endif //MDSTRESSLAB_RIGIDITY_H
