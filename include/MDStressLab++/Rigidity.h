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
    Eigen::MatrixXd matrixR, matrixRTR;
    const MatrixXd& coordinates;
    Rigidity(const MatrixXd & coordinates);
    std::vector<double> project(const Eigen::VectorXd& gi) const;
};


#endif //MDSTRESSLAB_RIGIDITY_H
