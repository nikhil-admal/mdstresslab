/*
 * Mls.h
 *
 *  Created on: Feb 11, 2020
 *      Author: Min
 */

#ifndef MLS_H_
#define MLS_H_

#include <vector>
#include "typedef.h"
#include "BoxConfiguration.h"

class Mls {
    public:
    std::string name;
    double radiusMls;
    std::vector<Matrix3d> deformationGradient;
    std::vector<Vector3d> gridPushed;

    //Mls(const MatrixXd& referenceCoordinates, const MatrixXd& currentCoordinates, const std::vector<Vector3d>& gridCoordinates, double radiusMls, const std::string name);
    //Mls(const BoxConfiguration& body, const std::vector<Vector3d>& gridCoordinates, double radiusMls, const std::string name);
    Mls(const BoxConfiguration& body, const std::vector<Vector3d>& gridCoordinates, double radiusMls, const std::string name);
    ~Mls();

    void pushToCauchy(const std::vector<Matrix3d>& piolaStress,std::vector<Matrix3d>& cauchyStress);
    void writeDeformationGradient();
    void writeGridPushed();
    void writePushedCauchyStress(std::vector<Matrix3d>& cauchyStress);

};

#endif /* MLS_H_ */


