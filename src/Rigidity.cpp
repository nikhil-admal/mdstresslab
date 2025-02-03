//
// Created by Nikhil Chandra Admal on 1/29/25.
//

#include "Rigidity.h"
#include "typedef.h"
#include <set>
#include "neighbor_list.h"

Rigidity::Rigidity(const MatrixXd& coordinates, const NeighListOne* p_inputNl)
    : coordinates(coordinates),
      p_inputNl(p_inputNl),
      halfNl(getHalfNeighborList(coordinates,p_inputNl->cutoff))
{
    MY_HEADING("Constructing rigidity");
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    int halfNlSize= 0;
    for(int i=0; i<halfNl.numberOfParticles;i++) halfNlSize+= halfNl.Nneighbors[i];
    std::cout << "half neighborlist size = " << halfNlSize << std::endl;
    int nnzmax= 2*DIM*halfNlSize;
    size_t numberOfParticles= halfNl.numberOfParticles;

    tripletList.reserve(nnzmax);

    for(size_t i_particle1=0; i_particle1<numberOfParticles; ++i_particle1)
    {
        int numberOfNeighborsOf1= halfNl.Nneighbors[i_particle1];

        //	Loop through the neighbors of particle1 in the halfNl
        for(size_t i_neighborOf1= 0; i_neighborOf1<numberOfNeighborsOf1; i_neighborOf1++)
        {
            size_t index= halfNl.beginIndex[i_particle1]+i_neighborOf1;
            size_t i_particle2= halfNl.neighborList[index];
            Vector3d rab= coordinates.row(i_particle1)-coordinates.row(i_particle2);
            for(size_t k=0; k<2*DIM; k++) {
                if (k < DIM)
                    tripletList.push_back(T(index, DIM * i_particle1 + k, rab(k) / rab.norm()));
                else
                    tripletList.push_back(T(index, DIM * i_particle2 + k - DIM, -rab(k - DIM) / rab.norm()));
            }
        }
    }
    matrixR.resize(halfNlSize,DIM*numberOfParticles);
    matrixR.setFromTriplets(tripletList.begin(), tripletList.end());

    matrixRTR.resize(DIM*numberOfParticles,DIM*numberOfParticles);
    matrixRTR= matrixR.transpose() * matrixR;
    std::cout << matrixRTR.rows() << "   " << matrixRTR.cols() << std::endl;
    std::cout << "Non-zeros in RTR = " << matrixRTR.nonZeros() << std::endl;
    //std::cout << "Diagonal of RTR = " << matrixRTR.diagonal() << std::endl;
    std::cout << "last row of RTR = " << matrixR.row(matrixR.rows()-1).norm() << std::endl;
}

std::vector<double> Rigidity::project(const Eigen::VectorXd& gi) const{
    VectorXd u,h;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
    //Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> cg;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cg;
    //cg.setMaxIterations(3*matrixRTR.cols());


    std::cout << "Norm of b :     " << gi.norm() << std::endl;
    u= cg.compute(matrixRTR).solve(gi);
    std::cout << "Norm of x :     " << u.norm() << std::endl;
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
    std::cout << "Residual error: " << (matrixRTR*u-gi).norm()/gi.norm()  << std::endl;
    h= matrixR*u.transpose();

    std::vector<double> hij;
    int inputNlSize= 0;
    for(int i=0; i<p_inputNl->numberOfParticles;i++) inputNlSize+= p_inputNl->Nneighbors[i];
    hij.resize(inputNlSize,0.0);

    size_t numberOfParticles= halfNl.numberOfParticles;
    for(size_t i_particle1=0; i_particle1<numberOfParticles; ++i_particle1) {
        int numberOfNeighborsOf1= p_inputNl->Nneighbors[i_particle1];

        //	Loop through the neighbors of particle 1 in the input NL
        for (size_t i_neighborOf1= 0; i_neighborOf1 < numberOfNeighborsOf1; i_neighborOf1++) {
            size_t i_Partcle2InParticle1Nl= p_inputNl->beginIndex[i_particle1] + i_neighborOf1;
            size_t i_particle2= p_inputNl->neighborList[i_Partcle2InParticle1Nl];

            int i_particlei= (i_particle1 < i_particle2) ? i_particle1 : i_particle2;
            int i_particlej= (i_particle1 < i_particle2) ? i_particle2 : i_particle1;

            // look for particlej in the half neighborlist of particlei
            for(int i_neighborOfi=0; i_neighborOfi<halfNl.Nneighbors[i_particlei]; ++i_neighborOfi) {
                int index= halfNl.beginIndex[i_particlei] + i_neighborOfi;
                int i_particlek = halfNl.neighborList[index];
                if (i_particlek == i_particlej) {
                    hij[i_Partcle2InParticle1Nl]= h(index);
                    break;
                }
            }
        }
    }
    return hij;
}

/*
NeighListOne Rigidity::getHalfNeighborList(const NeighListOne* p_inputNl)
{
    NeighListOne halfNl;
    int numberOfParticles= p_inputNl->numberOfParticles;
    halfNl.numberOfParticles= numberOfParticles;
    std::vector<std::set<int>> temp(numberOfParticles);

    for(size_t i_particle1=0; i_particle1<numberOfParticles; ++i_particle1) {
        int numberOfNeighborsOf1 = p_inputNl->Nneighbors[i_particle1];

        //	Loop over the neighbors of particle 1 in the input NL
        for (size_t i_neighborOf1 = 0; i_neighborOf1 < numberOfNeighborsOf1; i_neighborOf1++) {
            size_t index= p_inputNl->beginIndex[i_particle1] + i_neighborOf1;
            size_t i_particle2 = p_inputNl->neighborList[index];

            int i_particlei = (i_particle1 < i_particle2) ? i_particle1 : i_particle2;
            int i_particlej = (i_particle1 < i_particle2) ? i_particle2 : i_particle1;

            temp[i_particlei].insert(i_particlej);
        }
    }

    int halfNlSize= 0;
    for(int i=0; i<numberOfParticles;i++) halfNlSize+= temp[i].size();

    halfNl.Nneighbors= new int[numberOfParticles];
    halfNl.beginIndex= new int[numberOfParticles];
    halfNl.neighborList= new int[halfNlSize];

    int i_particle= 0;
    std::vector<int> tempNl;
    for(const auto neighbors : temp)
    {
        halfNl.Nneighbors[i_particle]= neighbors.size();
        halfNl.beginIndex[i_particle]= tempNl.size();
        // concaternate neighbors to the end of tempNl
        tempNl.insert(tempNl.end(),std::make_move_iterator(neighbors.begin()), std::make_move_iterator(neighbors.end()));
        i_particle++;
    }
    std::memcpy(halfNl.neighborList, tempNl.data(), sizeof(int) * tempNl.size());

    return halfNl;
}
 */
NeighListOne Rigidity::getHalfNeighborList(const MatrixXd& coordinates, const double& cutoff)
{
    NeighListOne halfNl;
    NeighListOne* fullNl;
    int numberOfParticles= coordinates.rows();
    VectorXi particleContributing(numberOfParticles);
    particleContributing.setOnes();


    MY_HEADING("Building neighbor list for rigidity calculation");
    NeighList* tempNl;
    nbl_initialize(&tempNl);
    nbl_build(tempNl,numberOfParticles,
              coordinates.data(),
              cutoff,
              1,
              &cutoff,
              particleContributing.data());
    fullNl= &(tempNl->lists[0]);
    int fullNlSize= 0;
    for(int i=0; i<numberOfParticles;i++) fullNlSize+= fullNl->Nneighbors[i];
    std::cout << "Size of neighbor list = " << fullNlSize << std::endl;


    halfNl.numberOfParticles= numberOfParticles;
    halfNl.Nneighbors= new int[numberOfParticles];
    halfNl.beginIndex= new int[numberOfParticles];
    halfNl.neighborList= new int[fullNlSize/2];

   int indexHalf= 0;
    for(size_t i_particle1=0; i_particle1<numberOfParticles; ++i_particle1)
    {
        halfNl.Nneighbors[i_particle1]= 0;
        int numberOfNeighborsOf1 = fullNl->Nneighbors[i_particle1];
        //	Loop over the neighbors of particle 1 in the input NL
        for (size_t i_neighborOf1 = 0; i_neighborOf1 < numberOfNeighborsOf1; i_neighborOf1++) {
            size_t indexFull= fullNl->beginIndex[i_particle1] + i_neighborOf1;
            size_t i_particle2 = fullNl->neighborList[indexFull];

            if (halfNl.Nneighbors[i_particle1]==0) halfNl.beginIndex[i_particle1]= indexHalf;
            if (i_particle1<i_particle2) {
                halfNl.neighborList[indexHalf] = i_particle2;
                halfNl.Nneighbors[i_particle1]++;
                indexHalf++;
            }
        }
    }
    nbl_clean(&tempNl);
    return halfNl;
}

Rigidity::~Rigidity()
{
    delete[] halfNl.Nneighbors;
    delete[] halfNl.beginIndex;
    delete[] halfNl.neighborList;
}
