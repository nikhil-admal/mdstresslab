#include "BoxConfiguration.h"
#include "Grid.h"
#include "MethodLdad.h"
#include "MethodSphere.h"
#include "Stress.h"
#include "calculateStress.h"
#include "typedef.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

namespace {
template<typename TMethod>
void validateIdealGasStress(const Stress<TMethod,Cauchy>& kineticStress,
                            const Matrix3d& expectedStress,
                            const Matrix3d& rawVelocityStress,
                            const double& expectedMassDensity,
                            const double& meanStressTolerance,
                            const double& pointwiseStressTolerance,
                            const double& densityTolerance)
{
    Matrix3d meanStress= Matrix3d::Zero();
    Matrix3d maxAbsDeviation= Matrix3d::Zero();
    double meanMassDensity= 0.0;

    for (int i_grid=0; i_grid<kineticStress.field.size(); ++i_grid)
    {
        const auto& stress= kineticStress.field[i_grid];
        meanStress+= stress;
        maxAbsDeviation= maxAbsDeviation.cwiseMax((stress-expectedStress).cwiseAbs());
        meanMassDensity+= kineticStress.massDensityField[i_grid];

        if (kineticStress.massDensityField[i_grid] < -epsilon)
            MY_ERROR("Mass density should be nonnegative.");

        if (kineticStress.massDensityField[i_grid] > epsilon)
        {
            Vector3d expectedVelocity=
                    kineticStress.momentumDensityField[i_grid]/kineticStress.massDensityField[i_grid];
            if ((kineticStress.velocityField[i_grid]-expectedVelocity).norm() > 1e-12)
                MY_ERROR("Continuum velocity does not equal momentum density divided by mass density.");
        }
        else
        {
            if (kineticStress.velocityField[i_grid].norm() > 1e-12)
                MY_ERROR("Continuum velocity should be zero where mass density is zero.");
        }
    }

    const double numberOfGridPoints= static_cast<double>(kineticStress.field.size());
    meanStress/= numberOfGridPoints;
    meanMassDensity/= numberOfGridPoints;

    if ((rawVelocityStress-expectedStress).cwiseAbs().maxCoeff() < 10.0*meanStressTolerance)
        MY_ERROR("Bulk velocity is too small to distinguish raw-velocity stress from relative-velocity stress.");

    if (std::abs(meanMassDensity-expectedMassDensity) > densityTolerance)
    {
        std::cout << "Expected mass density = " << expectedMassDensity << " amu/A^3" << std::endl;
        std::cout << "Mean computed mass density = " << meanMassDensity << " amu/A^3" << std::endl;
        MY_ERROR("Ideal gas mass density mean does not match the normalized-kernel value.");
    }

    if ((meanStress-expectedStress).cwiseAbs().maxCoeff() > meanStressTolerance)
    {
        std::cout << "Expected stress:\n" << expectedStress << std::endl;
        std::cout << "Mean computed stress:\n" << meanStress << std::endl;
        MY_ERROR("Ideal gas kinetic stress mean does not match the instantaneous ideal gas value.");
    }

    if (maxAbsDeviation.diagonal().maxCoeff() > pointwiseStressTolerance)
    {
        std::cout << "Expected stress:\n" << expectedStress << std::endl;
        std::cout << "Maximum pointwise absolute deviation:\n" << maxAbsDeviation << std::endl;
        MY_ERROR("Ideal gas kinetic stress field is too far from the expected uniform field.");
    }

    std::cout << kineticStress.name << " mean mass density = "
              << meanMassDensity << " amu/A^3" << std::endl;
    std::cout << kineticStress.name << " mean computed stress:\n" << meanStress << std::endl;
}
}

int main()
{
    const std::string configFileName= "idealGas.lmp";
    std::ifstream file(configFileName);
    if(!file) MY_ERROR("ERROR: idealGas.lmp could not be opened for reading.");

    int numberOfParticles= 0;
    std::string line;
    while (std::getline(file,line))
    {
        std::string loweredLine= line;
        std::transform(loweredLine.begin(),loweredLine.end(),loweredLine.begin(),::tolower);
        if (loweredLine.find("atoms") != std::string::npos && (std::stringstream(line) >> numberOfParticles))
            break;
    }
    if (numberOfParticles <= 0) MY_ERROR("ERROR: Could not read number of particles.");

    BoxConfiguration body{numberOfParticles,false};
    body.readLMP(configFileName,Current);
    body.pbc= Vector3i(1,1,1);

    const double volume= body.box.determinant();
    const double boltzmannConstantEvPerK= 8.617333262145e-5;
    double totalMass= 0.0;
    Vector3d totalMomentum= Vector3d::Zero();
    for (int i=0; i<body.numberOfParticles; ++i)
    {
        Vector3d velocity= body.velocities.row(i);
        totalMass+= body.masses(i);
        totalMomentum+= body.masses(i)*velocity;
    }
    Vector3d averageVelocity= totalMomentum/totalMass;

    Matrix3d expectedStress= Matrix3d::Zero();
    Matrix3d rawVelocityStress= Matrix3d::Zero();
    double thermalKineticEnergyFactor= 0.0;
    for (int i=0; i<body.numberOfParticles; ++i)
    {
        Vector3d velocity= body.velocities.row(i);
        Vector3d relativeVelocity= velocity - averageVelocity;
        expectedStress-= amuAngstromSquaredPerPicosecondSquaredToEv*
                        body.masses(i)*relativeVelocity.transpose()*relativeVelocity/volume;
        rawVelocityStress-= amuAngstromSquaredPerPicosecondSquaredToEv*
                           body.masses(i)*velocity.transpose()*velocity/volume;
        thermalKineticEnergyFactor+= body.masses(i)*relativeVelocity.squaredNorm();
    }
    const double instantaneousTemperature=
            amuAngstromSquaredPerPicosecondSquaredToEv*thermalKineticEnergyFactor/
            (3.0*body.numberOfParticles*boltzmannConstantEvPerK);

    Grid<Current> grid(Vector3d(0.0,0.0,30.0),Vector3d(60.0,60.0,31.0),12,12);
    MethodSphere virial(20.0,"virial");
    Matrix3d ldadVectors;
    ldadVectors << 20.0, 0.0, 0.0,
                   0.0, 20.0, 0.0,
                   0.0, 0.0, 20.0;
    MethodLdadConstant ldadConstant(ldadVectors);
    MethodLdadTrigonometric ldadTrigonometric(ldadVectors);

    Stress<MethodSphere,Cauchy> kineticStressSphere("idealGasKineticSphere",virial,&grid);
    Stress<MethodLdadConstant,Cauchy> kineticStressLdadConstant("idealGasKineticLdadConstant",ldadConstant,&grid);
    Stress<MethodLdadTrigonometric,Cauchy> kineticStressLdadTrigonometric("idealGasKineticLdadTrigonometric",ldadTrigonometric,&grid);

    calculateKineticStress(body,std::tie(kineticStressSphere,
                                         kineticStressLdadConstant,
                                         kineticStressLdadTrigonometric));
    kineticStressSphere.write();
    kineticStressLdadConstant.write();
    kineticStressLdadTrigonometric.write();

    const double pressure= -expectedStress.trace()/3.0;
    const double analyticalPressure=
            body.numberOfParticles*boltzmannConstantEvPerK*instantaneousTemperature/volume;
    const double expectedMassDensity= totalMass/volume;
    const double meanTolerance= 0.20*pressure;
    const double pointwiseTolerance= 0.75*pressure;
    const double densityTolerance= 0.20*expectedMassDensity;

    validateIdealGasStress(kineticStressSphere,expectedStress,rawVelocityStress,
                           expectedMassDensity,meanTolerance,pointwiseTolerance,densityTolerance);
    validateIdealGasStress(kineticStressLdadConstant,expectedStress,rawVelocityStress,
                           expectedMassDensity,meanTolerance,pointwiseTolerance,densityTolerance);
    validateIdealGasStress(kineticStressLdadTrigonometric,expectedStress,rawVelocityStress,
                           expectedMassDensity,meanTolerance,pointwiseTolerance,densityTolerance);

    std::cout << "Mass-weighted average velocity = " << averageVelocity << " A/ps" << std::endl;
    std::cout << "Instantaneous ideal-gas temperature = " << instantaneousTemperature << " K" << std::endl;
    std::cout << "Instantaneous ideal-gas pressure from stress = " << pressure << " eV/A^3" << std::endl;
    std::cout << "Analytical ideal-gas pressure NkBT/V = " << analyticalPressure << " eV/A^3" << std::endl;
    std::cout << "Raw-velocity stress before subtracting bulk motion:\n" << rawVelocityStress << std::endl;
    std::cout << "Expected stress:\n" << expectedStress << std::endl;

    return 0;
}
