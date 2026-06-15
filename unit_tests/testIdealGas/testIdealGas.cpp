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
void validateVoxelGridFile(const std::string& filename,
                           const std::string& expectedGridCellsHeader,
                           const int expectedNumberOfGridPoints)
{
    std::ifstream file(filename);
    if (!file) MY_ERROR("ERROR: " + filename + " could not be opened for reading.");

    std::string line;
    bool foundGridSize= false;
    bool foundGridCells= false;
    int dataLineCount= 0;
    while (std::getline(file,line))
    {
        if (line == "ITEM: GRID SIZE nx ny nz")
        {
            int nx,ny,nz;
            if (!(file >> nx >> ny >> nz))
                MY_ERROR("ERROR: Could not read voxel grid dimensions from " + filename);
            file.ignore(32767,'\n');
            if (nx*ny*nz != expectedNumberOfGridPoints)
                MY_ERROR("ERROR: Voxel grid dimensions do not match expected grid size in " + filename);
            foundGridSize= true;
        }
        else if (line == expectedGridCellsHeader)
        {
            foundGridCells= true;
            break;
        }
    }

    if (!foundGridSize)
        MY_ERROR("ERROR: Missing voxel grid size header in " + filename);
    if (!foundGridCells)
        MY_ERROR("ERROR: Missing expected voxel grid cells header in " + filename);

    while (std::getline(file,line))
        if (!line.empty())
            ++dataLineCount;

    if (dataLineCount != expectedNumberOfGridPoints)
        MY_ERROR("ERROR: Voxel grid data line count does not match grid size in " + filename);
}

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

    const int nx= 12;
    const int ny= 12;
    const int nz= 1;
    const Vector3d lowerLimit(0.0,0.0,30.0);
    const Vector3d upperLimit(60.0,60.0,31.0);
    Grid<Current> grid(lowerLimit,upperLimit,nx,ny,nz);
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
    kineticStressSphere.write_voxel_grid(nx,ny,nz,lowerLimit,upperLimit);
    kineticStressLdadConstant.write_voxel_grid(nx,ny,nz,lowerLimit,upperLimit);
    kineticStressLdadTrigonometric.write_voxel_grid(nx,ny,nz,lowerLimit,upperLimit);

    const int numberOfGridPoints= nx*ny*nz;
    validateVoxelGridFile("idealGasKineticSphere.voxel_grid_stress",
                          "ITEM: GRID CELLS SXX SYY SZZ SYZ SXZ SXY",
                          numberOfGridPoints);
    validateVoxelGridFile("idealGasKineticSphere.voxel_grid_momentum_density",
                          "ITEM: GRID CELLS PX PY PZ",
                          numberOfGridPoints);
    validateVoxelGridFile("idealGasKineticSphere.voxel_grid_mass_density",
                          "ITEM: GRID CELLS RHO",
                          numberOfGridPoints);
    validateVoxelGridFile("idealGasKineticLdadConstant.voxel_grid_stress",
                          "ITEM: GRID CELLS SXX SYY SZZ SYZ SXZ SXY",
                          numberOfGridPoints);
    validateVoxelGridFile("idealGasKineticLdadConstant.voxel_grid_momentum_density",
                          "ITEM: GRID CELLS PX PY PZ",
                          numberOfGridPoints);
    validateVoxelGridFile("idealGasKineticLdadConstant.voxel_grid_mass_density",
                          "ITEM: GRID CELLS RHO",
                          numberOfGridPoints);
    validateVoxelGridFile("idealGasKineticLdadTrigonometric.voxel_grid_stress",
                          "ITEM: GRID CELLS SXX SYY SZZ SYZ SXZ SXY",
                          numberOfGridPoints);
    validateVoxelGridFile("idealGasKineticLdadTrigonometric.voxel_grid_momentum_density",
                          "ITEM: GRID CELLS PX PY PZ",
                          numberOfGridPoints);
    validateVoxelGridFile("idealGasKineticLdadTrigonometric.voxel_grid_mass_density",
                          "ITEM: GRID CELLS RHO",
                          numberOfGridPoints);

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
