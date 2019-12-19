/*
 * typedef.h
 *
 *  Created on: Nov 6, 2019
 *      Author: Nikhil
 */

#ifndef TYPEDEF_H_
#define TYPEDEF_H_

#include <iomanip>
#include <string>
#include <iostream>
#include "Eigen/Eigen/Dense"
#include <memory>

#define MY_ERROR(message)                                                \
  {                                                                      \
    std::cout << "* Error : \"" << message << "\" : " << __LINE__ << ":" \
              << __FILE__ << std::endl;                                  \
    exit(1);                                                             \
  }

#define MY_WARNING(message)                                                \
  {                                                                        \
    std::cout << "* Warning : \"" << message << "\" : " << __LINE__ << ":" \
              << __FILE__ << std::endl;                                    \
  }

#define MY_BANNER(announcement)                                            \
  {                                                                        \
	std::cout << "--------------------------------------------------------------------------------" << "\n";\
	std::cout << std::setw(40+strlen(announcement)/2)<< announcement << "\n";\
	std::cout << "--------------------------------------------------------------------------------" << "\n";\
  }
#define MY_HEADING(heading)                                             \
{																		\
	std::cout << "\n";													\
	std::cout << heading << "\n";										\
	std::cout << std::string(strlen(heading),'-')<< "\n";				\
}
#define MY_LINE(message)                                             	\
{																		\
	std::cout << message << "\n";										\
}
#define MY_SUBLINE(message)                                             \
{																		\
	std::cout << "     " << message << "\n";							\
}

const int DIM= 3;
typedef std::unique_ptr<double[]> array_dptr;
typedef std::unique_ptr<int> int_ptr;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatrixXd;
typedef Eigen::Matrix<double,DIM,DIM,Eigen::RowMajor> Matrix3d;
typedef Eigen::Matrix<int,1,Eigen::Dynamic,Eigen::RowMajor> VectorXi;
typedef Eigen::Matrix<double,1,Eigen::Dynamic,Eigen::RowMajor> VectorXd;
typedef Eigen::Matrix<double,1,DIM,Eigen::RowMajor> Vector3d;
typedef Eigen::Matrix<int,1,DIM,Eigen::RowMajor> Vector3i;
enum StressType {
	Cauchy,
	Piola
};

enum ConfigType {
	Reference,
	Current
};

const double epsilon= 1e-8;




#endif /* TYPEDEF_H_ */
