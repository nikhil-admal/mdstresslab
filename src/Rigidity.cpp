//
// Created by Nikhil Chandra Admal on 1/29/25.
//

#include "Rigidity.h"
#include "typedef.h"
#include <set>
#include "neighbor_list.h"

Rigidity::Rigidity(const MatrixXd& coordinates)
    : coordinates(coordinates)
{
    std::cout << "Constructing rigidity matrix" << std::endl;
    size_t numberOfParticles= coordinates.rows();
    std::cout << "Number of particles = " << numberOfParticles << std::endl;
    int halfNlSize= (numberOfParticles*(numberOfParticles-1))/2;

    matrixR.resize(halfNlSize,DIM*numberOfParticles);
    matrixRTR.resize(DIM*numberOfParticles,DIM*numberOfParticles);
    matrixR.setZero();
    int index= -1;
    for(int i_particle1=0; i_particle1<numberOfParticles; ++i_particle1)
    {
        //	Loop through the neighbors of particle1 in the halfNl
        for(int i_particle2= i_particle1+1; i_particle2<numberOfParticles; i_particle2++)
        {
            index++;
            Vector3d rab= coordinates.row(i_particle1)-coordinates.row(i_particle2);
            for(size_t k=0; k<2*DIM; k++) {
                if (k < DIM)
                    matrixR(index,DIM * i_particle1 + k)= rab(k) / rab.norm();
                else
                    matrixR(index, DIM * i_particle2 + k - DIM)= -rab(k - DIM) / rab.norm();
            }
        }
    }
    matrixRTR= matrixR.transpose() * matrixR;
}

std::vector<double> Rigidity::project(const Eigen::VectorXd& b) const{
    VectorXd u,h;
    u.setZero();
    //Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
    //Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXd> cg;
    //u= cg.compute(matrixRTR).solve(b);
    //std::cout << "#iterations:     " << cg.iterations() << std::endl;
    //std::cout << "estimated error: " << cg.error()      << std::endl;


    Eigen::BDCSVD<Eigen::MatrixXd> cg;
    cg.compute(matrixRTR,Eigen::ComputeThinU | Eigen::ComputeThinV);
    u= cg.solve(b);
    double residualError= (matrixRTR*u.transpose()-b).norm()/b.norm();
    if (residualError>1e-5) {
        std::cout << "residual = " << residualError << ". Residual error too big." << std::endl;
    }

    h= matrixR*u.transpose();
    std::vector<double> hij(h.size(),0.0);

    for(int i=0; i<h.size(); ++i)
        hij[i]=h(i);

    return hij;
}


