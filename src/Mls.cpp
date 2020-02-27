/*
 * Mls.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: Min
 */

#include <fstream>
#include <vector>
#include "Mls.h"
#include "typedef.h"

Mls::Mls(const MatrixXd& referenceCoordinates, const MatrixXd& currentCoordinates, const std::vector<Vector3d>& gridCoordinates, \
double radiusMls, const std::string name):radiusMls(radiusMls),name(name)
{
    MY_BANNER("Constructing Deformation Gradient using the Moving Least Square method!");
    deformationGradient.reserve(gridCoordinates.size());
    gridPushed.reserve(gridCoordinates.size());
    MatrixXd displacements(currentCoordinates.rows(),currentCoordinates.cols());
    displacements = currentCoordinates - referenceCoordinates;

    std::cout << currentCoordinates.rows() << std::endl;
    std::cout << currentCoordinates.cols() << std::endl;
    std::cout << "Reference Coordinates." << std::endl;
    std::cout << referenceCoordinates(0,0) << " " << referenceCoordinates(0,1) << " " << referenceCoordinates(0,2) << std::endl; 
    std::cout << "Current Coordinates." << std::endl;
    std::cout << currentCoordinates(0,0) << " " << currentCoordinates(0,1) << " " << currentCoordinates(0,2) << std::endl;

    double gridRadiusMls;
    double r2;
    Vector3d r;
    Vector3d rSanityI, rSanityJ, rSanityK, rSanityQ;
    int sanityI, sanityJ, sanityK, sanityQ;
    bool sanityCheck;
    int cycle_count;    

    Matrix3d tensorF;
    Vector3d gptIPushedF;
    Matrix4d sanity, A, inverseA, dAdx, dAdy, dAdz;
    Matrix4d inverseAXdAdx, inverseAXdAdy, inverseAXdAdz;
    Matrix4d inverseAXdAdxXInverseA, inverseAXdAdyXInverseA, inverseAXdAdzXInverseA; 
    MatrixXd MlsDisp(4,3), dMlsDispdx(4,3), dMlsDispdy(4,3), dMlsDispdz(4,3);
    MatrixXd fintmdx, fintmdy,fintmdz, sintmdx, sintmdy, sintmdz;
    MatrixXd B, inverseAXB, dBdx, dBdy, dBdz;
    MatrixXd exactDisplacements;
    MatrixXd P, dWdx;
    VectorXd W;
    
    std::vector<int> gridMlsAtomList;
    double ngridMLSAtomList;
    
    std::cout << referenceCoordinates.rows() << std::endl;
    std::cout << referenceCoordinates.cols() << std::endl;
    std::cout << "Displacements." << std::endl;
    std::cout << displacements(0,0) << " " << displacements(0,1) << " " << displacements(0,2) << std::endl;
    // debug only calculate one grid point;
    for (std::vector<Vector3d>::size_type iGrid = 0; iGrid != gridCoordinates.size(); iGrid++)
    //for (std::vector<Vector3d>::size_type iGrid = 0; iGrid != 1; iGrid++) // debug only calculate one point
    {
        gridRadiusMls = radiusMls;
        A = Matrix4d::Zero();
        cycle_count = 1;

        do
        {
            if (cycle_count != 1)
            {
                MY_WARNING("Something is wrong and need to double the MLS_Radius for grid point " + std::to_string(iGrid + 1) \
                         + ". Starting Cycle " + std::to_string(cycle_count) + ".")
            }

            if (ngridMLSAtomList == referenceCoordinates.rows())
            {
                MY_ERROR("Model degenerate to 2D. Cannot use MLS3D to generate deformation gradient!")
            }

            gridMlsAtomList.clear();
            ngridMLSAtomList = 0;

            for (int jRefAtoms = 0; jRefAtoms != referenceCoordinates.rows(); jRefAtoms++)
            {
                r(0) = referenceCoordinates(jRefAtoms,0) - gridCoordinates[iGrid](0);
                r(1) = referenceCoordinates(jRefAtoms,1) - gridCoordinates[iGrid](1);
                r(2) = referenceCoordinates(jRefAtoms,2) - gridCoordinates[iGrid](2);
                /*for (int k = 0; k != DIM; k++)
                [
                    r(k) = referenceCoordinates(jRefAtoms,k) - gridCoordinates[iGrid](k);
                ]
                */
                if (r.norm() <= gridRadiusMls)
                {
                    //std::cout << r.norm() << std::endl;
                    gridMlsAtomList.push_back(jRefAtoms);
                    ngridMLSAtomList++;
                }
            }
            // sanity check start!
            sanityCheck = false;

            if (ngridMLSAtomList <= 3)
            {
                sanityCheck = false;
                A = Matrix4d::Zero();
                MY_WARNING("Not enough atoms in the MLS radius " + std::to_string(gridRadiusMls) + ". Doubling it.");
                gridRadiusMls = 2.0 * gridRadiusMls;
                goto EndLoop;
            }
            else
            {
                for (int i = 0; i <= ngridMLSAtomList - 4; i++)
                {
                    sanityI = gridMlsAtomList[i];
                    rSanityI(0) = referenceCoordinates(sanityI,0) - gridCoordinates[iGrid](0);
                    rSanityI(1) = referenceCoordinates(sanityI,1) - gridCoordinates[iGrid](1);
                    rSanityI(2) = referenceCoordinates(sanityI,2) - gridCoordinates[iGrid](2);
                    for(int j = i + 1; j <= ngridMLSAtomList - 3; j++)
                    {
                        sanityJ = gridMlsAtomList[j];
                        rSanityJ(0) = referenceCoordinates(sanityJ,0) - gridCoordinates[iGrid](0);
                        rSanityJ(1) = referenceCoordinates(sanityJ,1) - gridCoordinates[iGrid](1);
                        rSanityJ(2) = referenceCoordinates(sanityJ,2) - gridCoordinates[iGrid](2);
                        for (int k = j + 1; k <= ngridMLSAtomList - 2; k++)
                        {
                            sanityK = gridMlsAtomList[k];
                            rSanityK(0) = referenceCoordinates(sanityK,0) - gridCoordinates[iGrid](0);
                            rSanityK(1) = referenceCoordinates(sanityK,1) - gridCoordinates[iGrid](1);
                            rSanityK(2) = referenceCoordinates(sanityK,2) - gridCoordinates[iGrid](2);
                            for (int q = k + 1; q <= ngridMLSAtomList - 1; q++)
                            {
                                sanityQ = gridMlsAtomList[q];
                                rSanityQ(0) = referenceCoordinates(sanityQ,0) - gridCoordinates[iGrid](0);
                                rSanityQ(1) = referenceCoordinates(sanityQ,1) - gridCoordinates[iGrid](1);
                                rSanityQ(2) = referenceCoordinates(sanityQ,2) - gridCoordinates[iGrid](2);

                                sanity = Matrix4d::Constant(1.0);
                                sanity(0,0) = rSanityI(0);
                                sanity(0,1) = rSanityI(1);
                                sanity(0,2) = rSanityI(2);
                                sanity(1,0) = rSanityJ(0);
                                sanity(1,1) = rSanityJ(1);
                                sanity(1,2) = rSanityJ(2);
                                sanity(2,0) = rSanityK(0);
                                sanity(2,1) = rSanityK(1);
                                sanity(2,2) = rSanityK(2);
                                sanity(3,0) = rSanityQ(0);
                                sanity(3,1) = rSanityQ(1);
                                sanity(3,2) = rSanityQ(2);
                                
                                if (sanity.determinant() > epsilon)
                                {
                                    sanityCheck = true;
                                    //inverseA = Matrix4d::Identity();
                                    //std::cout << "Passed Sanity Check!" << std::endl;
                                    goto SanityDone;
                                }
                            }
                        }
                    }
                }
                if (!sanityCheck)
                {
                    A = Matrix4d::Zero();
                    MY_WARNING("All atoms lying on a surface for the MLS radius " + std::to_string(gridRadiusMls) + ". Doubling it.");
                    gridRadiusMls = 2.0 * gridRadiusMls;
                    goto EndLoop;
                }
            }

            SanityDone:;
            P.resize(ngridMLSAtomList,4);
            W.resize(ngridMLSAtomList);
            dWdx.resize(ngridMLSAtomList,3);
            B.resize(4,ngridMLSAtomList);
            dBdx.resize(4,ngridMLSAtomList);
            dBdy.resize(4,ngridMLSAtomList);
            dBdz.resize(4,ngridMLSAtomList);
            inverseAXB.resize(4,ngridMLSAtomList);
            exactDisplacements.resize(ngridMLSAtomList,3);
            fintmdx.resize(4,ngridMLSAtomList);
            fintmdy.resize(4,ngridMLSAtomList);
            fintmdz.resize(4,ngridMLSAtomList);
            sintmdx.resize(4,ngridMLSAtomList);
            sintmdy.resize(4,ngridMLSAtomList);
            sintmdz.resize(4,ngridMLSAtomList);

            tensorF = Matrix3d::Zero();
            gptIPushedF = Vector3d::Zero();
            P = MatrixXd::Zero(ngridMLSAtomList,4);
            W = VectorXd::Zero(ngridMLSAtomList);
            dWdx = MatrixXd::Zero(ngridMLSAtomList,3);
            A = Matrix4d::Zero();
            inverseA = Matrix4d::Zero();
            dAdx = Matrix4d::Zero();
            dAdy = Matrix4d::Zero();
            dAdz = Matrix4d::Zero();  
            B = MatrixXd::Zero(4,ngridMLSAtomList);  
            dBdx = MatrixXd::Zero(4,ngridMLSAtomList);
            dBdy = MatrixXd::Zero(4,ngridMLSAtomList);
            dBdz = MatrixXd::Zero(4,ngridMLSAtomList);
            inverseAXB = MatrixXd::Zero(4,ngridMLSAtomList);
            MlsDisp = MatrixXd::Zero(4,3);
            dMlsDispdx = MatrixXd::Zero(4,3);
            dMlsDispdy = MatrixXd::Zero(4,3);
            dMlsDispdz = MatrixXd::Zero(4,3);
            exactDisplacements = MatrixXd::Zero(ngridMLSAtomList,3);
            inverseAXdAdx = Matrix4d::Zero();
            inverseAXdAdy = Matrix4d::Zero();
            inverseAXdAdz = Matrix4d::Zero();
            inverseAXdAdxXInverseA = Matrix4d::Zero();
            inverseAXdAdyXInverseA = Matrix4d::Zero();
            inverseAXdAdzXInverseA = Matrix4d::Zero();
            fintmdx = MatrixXd::Zero(4,ngridMLSAtomList);
            fintmdy = MatrixXd::Zero(4,ngridMLSAtomList);
            fintmdz = MatrixXd::Zero(4,ngridMLSAtomList);
            sintmdx = MatrixXd::Zero(4,ngridMLSAtomList);
            sintmdy = MatrixXd::Zero(4,ngridMLSAtomList);
            sintmdz = MatrixXd::Zero(4,ngridMLSAtomList);            

            for (int i = 0; i < ngridMLSAtomList; i++)
            {
                r(0) = referenceCoordinates(gridMlsAtomList[i],0) - gridCoordinates[iGrid](0);
                r(1) = referenceCoordinates(gridMlsAtomList[i],1) - gridCoordinates[iGrid](1);
                r(2) = referenceCoordinates(gridMlsAtomList[i],2) - gridCoordinates[iGrid](2);
                // Form P
                P(i,0) = 1.0;
                P(i,1) = referenceCoordinates(gridMlsAtomList[i],0);
                P(i,2) = referenceCoordinates(gridMlsAtomList[i],1);
                P(i,3) = referenceCoordinates(gridMlsAtomList[i],2);
                // Form W and dWdx
                r2 = r.norm() / gridRadiusMls;
                if (r2 <= 0.5)
                {
                    W(i) = 2.0 / 3.0 - 4.0 * r2 * r2 + 4.0 * r2 * r2 * r2;
                    dWdx(i,0) = -8.0 * r(0) / gridRadiusMls / gridRadiusMls \
                              + 12.0 * r2 * r(0) / gridRadiusMls / gridRadiusMls;
                    dWdx(i,1) = -8.0 * r(1) / gridRadiusMls / gridRadiusMls \
                              + 12.0 * r2 * r(1) / gridRadiusMls / gridRadiusMls;
                    dWdx(i,2) = -8.0 * r(2) / gridRadiusMls / gridRadiusMls \
                              + 12.0 * r2 * r(2) / gridRadiusMls / gridRadiusMls;
                }
                else if (r2 >= 0.5 && r2 <= 1.0)
                {
                    W(i) = 4.0 / 3.0 - 4.0 * r2 + 4.0 * r2 * r2 - 4.0 / 3.0 * r2 * r2 * r2;
                    dWdx(i,0) = -4.0 * r(0) / gridRadiusMls / gridRadiusMls / r2 \
                              + 8.0 * r(0) / gridRadiusMls / gridRadiusMls \
                              - 4.0 * r2 * r(0) / gridRadiusMls / gridRadiusMls;
                    dWdx(i,1) = -4.0 * r(1) / gridRadiusMls / gridRadiusMls / r2 \
                              + 8.0 * r(1) / gridRadiusMls / gridRadiusMls \
                              - 4.0 * r2 * r(1) / gridRadiusMls / gridRadiusMls;
                    dWdx(i,2) = -4.0 * r(2) / gridRadiusMls / gridRadiusMls / r2 \
                              + 8.0 * r(2) / gridRadiusMls / gridRadiusMls \
                              - 4.0 * r2 * r(2) / gridRadiusMls / gridRadiusMls;
                }
                else
                {
                    W(i) = 0.0;
                    dWdx(i,0) = 0.0;
                    dWdx(i,1) = 0.0;
                    dWdx(i,2) = 0.0;
                }
                exactDisplacements(i,0) = displacements(gridMlsAtomList[i],0);
                exactDisplacements(i,1) = displacements(gridMlsAtomList[i],1);
                exactDisplacements(i,2) = displacements(gridMlsAtomList[i],2);              
            }

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < ngridMLSAtomList; j++)
                {
                    B(i,j) = P(j,i) * W(j);
                }
            }

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < ngridMLSAtomList; k++)
                    {
                        A(i,j) = A(i,j) + B(i,k) * P(k,j);
                    }
                }
            }

            inverseA = A.inverse();
            gridRadiusMls = 2.0 * gridRadiusMls;
/*
            std::cout << "P" << std::endl;
            std::cout << P << std::endl;
            std::cout << "W" << std::endl;
            std::cout << W << std::endl;

            std::cout << "B" << std::endl;
            std::cout << B << std::endl;
*/
            std::cout << "A" << std::endl;
            std::cout << A << std::endl;

            std::cout << "inverseA" << std::endl;
            std::cout << inverseA << std::endl;

            EndLoop:;
            cycle_count++;
        } while (A.determinant() < epsilon);

        gridRadiusMls = gridRadiusMls / 2.0;

        // inv_A * B * exact_disp
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < ngridMLSAtomList; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    inverseAXB(i,j) = inverseAXB(i,j) + inverseA(i,k) * B(k,j);
                }
            }
        }

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < ngridMLSAtomList; k++)
                {
                    MlsDisp(i,j) = MlsDisp(i,j) + inverseAXB(i,k) * exactDisplacements(k,j);
                }
            }
        }

        //std::cout << "MlsDisp: " << std::endl;
        //std::cout << MlsDisp << std::endl;

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < ngridMLSAtomList; j++)
            {
                dBdx(i,j) = P(j,i) * dWdx(j,0);
                dBdy(i,j) = P(j,i) * dWdx(j,1);
                dBdz(i,j) = P(j,i) * dWdx(j,2);
            }
        }
        
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < ngridMLSAtomList; k++)
                {
                    dAdx(i,j) = dAdx(i,j) + dBdx(i,k) * P(k,j);
                    dAdy(i,j) = dAdy(i,j) + dBdy(i,k) * P(k,j);
                    dAdz(i,j) = dAdz(i,j) + dBdz(i,k) * P(k,j);
                }
            }
        }

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    inverseAXdAdx(i,j) = inverseAXdAdx(i,j) + inverseA(i,k) * dAdx(k,j);
                    inverseAXdAdy(i,j) = inverseAXdAdy(i,j) + inverseA(i,k) * dAdy(k,j);
                    inverseAXdAdz(i,j) = inverseAXdAdz(i,j) + inverseA(i,k) * dAdz(k,j);
                }
            }
        }
        
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {      
                    inverseAXdAdxXInverseA(i,j) = inverseAXdAdxXInverseA(i,j) + inverseAXdAdx(i,k) * inverseA(k,j);
                    inverseAXdAdyXInverseA(i,j) = inverseAXdAdyXInverseA(i,j) + inverseAXdAdy(i,k) * inverseA(k,j);
                    inverseAXdAdzXInverseA(i,j) = inverseAXdAdzXInverseA(i,j) + inverseAXdAdz(i,k) * inverseA(k,j);
                }
            }
        }

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < ngridMLSAtomList; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                     fintmdx(i,j) = fintmdx(i,j) + inverseAXdAdxXInverseA(i,k) * B(k,j);
                     fintmdy(i,j) = fintmdy(i,j) + inverseAXdAdyXInverseA(i,k) * B(k,j);
                     fintmdz(i,j) = fintmdz(i,j) + inverseAXdAdzXInverseA(i,k) * B(k,j);
                }
            }
        }

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < ngridMLSAtomList; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    sintmdx(i,j) = sintmdx(i,j) + inverseA(i,k) * dBdx(k,j);
                    sintmdy(i,j) = sintmdy(i,j) + inverseA(i,k) * dBdy(k,j);
                    sintmdz(i,j) = sintmdz(i,j) + inverseA(i,k) * dBdz(k,j);
                }
            }
        }

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < ngridMLSAtomList; j++)
            {
                 fintmdx(i,j) = -fintmdx(i,j) + sintmdx(i,j);
                 fintmdy(i,j) = -fintmdy(i,j) + sintmdy(i,j);
                 fintmdz(i,j) = -fintmdz(i,j) + sintmdz(i,j);
            }
        }

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < ngridMLSAtomList; k++)
                {
                    dMlsDispdx(i,j) = dMlsDispdx(i,j) + fintmdx(i,k) * exactDisplacements(k,j);
                    dMlsDispdy(i,j) = dMlsDispdy(i,j) + fintmdy(i,k) * exactDisplacements(k,j);
                    dMlsDispdz(i,j) = dMlsDispdz(i,j) + fintmdz(i,k) * exactDisplacements(k,j);
                }
            }
        }

        tensorF(0,0) = 1.0 + MlsDisp(1,0) \
                              + 1.0 * dMlsDispdx(0,0) + gridCoordinates[iGrid](0) * dMlsDispdx(1,0) \
        + gridCoordinates[iGrid](1) * dMlsDispdx(2,0) + gridCoordinates[iGrid](2) * dMlsDispdx(3,0);
        tensorF(0,1) =       MlsDisp(2,0) \
                              + 1.0 * dMlsDispdy(0,0) + gridCoordinates[iGrid](0) * dMlsDispdy(1,0) \
        + gridCoordinates[iGrid](1) * dMlsDispdy(2,0) + gridCoordinates[iGrid](2) * dMlsDispdy(3,0);
        tensorF(0,2) =       MlsDisp(3,0) \
                              + 1.0 * dMlsDispdz(0,0) + gridCoordinates[iGrid](0) * dMlsDispdz(1,0) \
        + gridCoordinates[iGrid](1) * dMlsDispdz(2,0) + gridCoordinates[iGrid](2) * dMlsDispdz(3,0);
        tensorF(1,0) =       MlsDisp(1,1) \
                              + 1.0 * dMlsDispdx(0,1) + gridCoordinates[iGrid](0) * dMlsDispdx(1,1) \
        + gridCoordinates[iGrid](1) * dMlsDispdx(2,1) + gridCoordinates[iGrid](2) * dMlsDispdx(3,1);
        tensorF(1,1) = 1.0 + MlsDisp(2,1) \
                              + 1.0 * dMlsDispdy(0,1) + gridCoordinates[iGrid](0) * dMlsDispdy(1,1) \
        + gridCoordinates[iGrid](1) * dMlsDispdy(2,1) + gridCoordinates[iGrid](2) * dMlsDispdy(3,1);
        tensorF(1,2) =       MlsDisp(3,1) \
                              + 1.0 * dMlsDispdz(0,1) + gridCoordinates[iGrid](0) * dMlsDispdz(1,1) \
        + gridCoordinates[iGrid](1) * dMlsDispdz(2,1) + gridCoordinates[iGrid](2) * dMlsDispdz(3,1);
        tensorF(2,0) =       MlsDisp(1,2) \
                              + 1.0 * dMlsDispdx(0,2) + gridCoordinates[iGrid](0) * dMlsDispdx(1,2) \
        + gridCoordinates[iGrid](1) * dMlsDispdx(2,2) + gridCoordinates[iGrid](2) * dMlsDispdx(3,2);
        tensorF(2,1) =       MlsDisp(2,2) \
                              + 1.0 * dMlsDispdy(0,2) + gridCoordinates[iGrid](0) * dMlsDispdy(1,2) \
        + gridCoordinates[iGrid](1) * dMlsDispdy(2,2) + gridCoordinates[iGrid](2) * dMlsDispdy(3,2);
        tensorF(2,2) = 1.0 + MlsDisp(3,2) \
                              + 1.0 * dMlsDispdz(0,2) + gridCoordinates[iGrid](0) * dMlsDispdz(1,2) \
        + gridCoordinates[iGrid](1) * dMlsDispdz(2,2) + gridCoordinates[iGrid](2) * dMlsDispdz(3,2);

        gptIPushedF(0) = gridCoordinates[iGrid](0) \
        + MlsDisp(0,0) + MlsDisp(1,0) * gridCoordinates[iGrid](0) \
                       + MlsDisp(2,0) * gridCoordinates[iGrid](1) \
                       + MlsDisp(3,0) * gridCoordinates[iGrid](2);
        gptIPushedF(1) = gridCoordinates[iGrid](1) \
        + MlsDisp(0,1) + MlsDisp(1,1) * gridCoordinates[iGrid](0) \
                       + MlsDisp(2,1) * gridCoordinates[iGrid](1) \
                       + MlsDisp(3,1) * gridCoordinates[iGrid](2);
        gptIPushedF(2) = gridCoordinates[iGrid](2) \
        + MlsDisp(0,2) + MlsDisp(1,2) * gridCoordinates[iGrid](0) \
                       + MlsDisp(2,2) * gridCoordinates[iGrid](1) \
                       + MlsDisp(3,2) * gridCoordinates[iGrid](2);

/*
        std::cout << "tensorF: " << std::endl;
        std::cout << tensorF << std::endl;
        std::cout << "gptIPushedF: " << std::endl;
        std::cout << gptIPushedF << std::endl;
*/        
        deformationGradient.push_back(tensorF.transpose());
        gridPushed.push_back(gptIPushedF);
    }   
}

Mls::~Mls()
{
    // TODO Auto-generated destructor stub
}

void Mls::pushToCauchy(const std::vector<Matrix3d>& piolaStress,std::vector<Matrix3d>& cauchyStress)
{
    MY_HEADING("Pushing the Piola-Kirchhoff Stress to the Cauchy stress.");
    for (std::vector<Matrix3d>::size_type i = 0; i != piolaStress.size(); i++)
    {
        cauchyStress.push_back(piolaStress[i] * deformationGradient[i].transpose() / deformationGradient[i].determinant()); 
    }
}

void Mls::writeDeformationGradient()
{
    MY_HEADING("Writing Deformation Gradient.");
    std::ofstream file(name+".DeformationGradient");
	file << deformationGradient.size() << "\n";
	file << "\n";
	for (auto& F : deformationGradient)
	{
		Eigen::Map<Eigen::Matrix<double,1,DIM*DIM>> FRow(F.data(), F.size());
		Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "      ", "\n", "", "", "", "");
		file << FRow.format(fmt) << "      " << std::setprecision(16) << F.determinant() << std::endl;
	}
}

void Mls::writeGridPushed()
{
    MY_HEADING("Writing Pushed grid points.");
    std::ofstream file(name+".pushedGrids");
	file << gridPushed.size() << "\n";
	file << "\n";
	for (auto& grid : gridPushed)
	{
		Eigen::Map<Eigen::Matrix<double,1,DIM>> gridRow(grid.data(), grid.size());
		Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "      ", "\n", "", "", "", "");
		file << gridRow.format(fmt) << std::endl;
	}
}

void Mls::writePushedCauchyStress(std::vector<Matrix3d>& cauchyStress)
{
    MY_HEADING("Writing Pushed Cauchy Stress.");
    std::ofstream file(name+".CauchyPushed");
	file << cauchyStress.size() << "\n";
	file << "\n";
    for (auto& pushedStress : cauchyStress)
	{
        Eigen::Map<Eigen::Matrix<double,1,DIM*DIM>> pushedStressRow(pushedStress.data(), pushedStress.size());
		Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "      ", "\n", "", "", "", "");
		file << pushedStressRow.format(fmt) << std::endl;
	}
}




