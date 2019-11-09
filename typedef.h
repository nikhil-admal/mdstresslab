/*
 * typedef.h
 *
 *  Created on: Nov 6, 2019
 *      Author: Nikhil
 */

#ifndef TYPEDEF_H_
#define TYPEDEF_H_

#include "Eigen/Eigen/Dense"
const int DIM= 3;
typedef std::unique_ptr<double[]> array_dptr;
typedef std::unique_ptr<int> int_ptr;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatrixXd;
typedef Eigen::Matrix<double,DIM,DIM,Eigen::RowMajor> Matrix3d;
typedef Eigen::VectorXi VectorXi;




#endif /* TYPEDEF_H_ */
