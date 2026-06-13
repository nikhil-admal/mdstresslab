/*
 * Regression test for LDAD stress under a non-orthogonal periodic basis.
 */
#include "MethodLdad.h"
#include <string>
#include <iostream>
#include <tuple>
#include <fstream>
#include "BoxConfiguration.h"
#include "calculateStress.h"
#include "Grid.h"
#include "typedef.h"

namespace {
void checkSkewedBox(const Matrix3d& box, const std::string& name)
{
    if (std::abs(box.col(0).dot(box.col(1))) < epsilon ||
        std::abs(box.col(0).dot(box.col(2))) < epsilon ||
        std::abs(box.col(1).dot(box.col(2))) < epsilon)
        MY_ERROR(name + " is not non-orthogonal.");
}

template<typename TMethod, StressType stressType>
void compareStressComponentsToReference(const Stress<TMethod,stressType>& stress,
                                        const std::string& referenceFilename)
{
    std::ifstream fileReference(referenceFilename);
    if(!fileReference) MY_ERROR("ERROR: " + referenceFilename + " could not be opened for reading!");

    int ngridReference;
    fileReference >> ngridReference;
    if(stress.field.size() != ngridReference)
        MY_ERROR("Test failed in " + referenceFilename + ". Number of grid points do not match.");

    std::string line;
    std::getline(fileReference,line);
    std::getline(fileReference,line);

    double maxDifference= 0.0;
    for (int i_point=0; i_point<ngridReference; ++i_point)
    {
        double x,y,z;
        double sxx,syy,szz,sxy,sxz,syz;
        fileReference >> x >> y >> z >> sxx >> syy >> szz >> sxy >> sxz >> syz;

        const Matrix3d& value= stress.field[i_point];
        maxDifference= std::max(maxDifference,std::abs(value(0,0)-sxx));
        maxDifference= std::max(maxDifference,std::abs(value(1,1)-syy));
        maxDifference= std::max(maxDifference,std::abs(value(2,2)-szz));
        maxDifference= std::max(maxDifference,std::abs(value(0,1)-sxy));
        maxDifference= std::max(maxDifference,std::abs(value(0,2)-sxz));
        maxDifference= std::max(maxDifference,std::abs(value(1,2)-syz));
    }

    const double tolerance= 1e-8;
    if (maxDifference > tolerance)
    {
        std::cout << "Maximum stress-component difference = " << maxDifference << std::endl;
        std::cout << "Tolerance = " << tolerance << std::endl;
        MY_ERROR("Non-orthogonal PBC stress regression failed against " + referenceFilename);
    }
}
}

int main()
{
    int numberOfParticles;
    int referenceAndFinal= true;
    std::string configFileName= "config.data";
    std::string modelname= "SW_StillingerWeber_1985_Si__MO_405512056662_005";

    std::ifstream file(configFileName);
    if(!file) MY_ERROR("ERROR: config.data could not be opened for reading!");

    file >> numberOfParticles;
    if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");

    BoxConfiguration body{numberOfParticles,referenceAndFinal};
    body.read(configFileName,referenceAndFinal);
    checkSkewedBox(body.reference_box,"Reference box");
    checkSkewedBox(body.box,"Current box");

    Kim kim(modelname);

    int ngrid = 125;
    Grid<Reference> gridFromFile_ref(ngrid);
    gridFromFile_ref.read("grid_pk1.data");

    Grid<Current> gridFromFile_def(ngrid);
    gridFromFile_def.read("grid_cauchy.data");

    Matrix3d ldadVectors_ref;
    ldadVectors_ref << 5.43094977840521, 0.0, 0.0,
                       0.0, 5.43094977840521, 0.0,
                       0.0, 0.0, 5.43094977840521;

    MethodLdadConstant ldad_constant_ref(ldadVectors_ref);
    MethodLdadTrigonometric ldad_trigonometric_ref(ldadVectors_ref);

    Stress<MethodLdadConstant,Piola> ldad_constant_stress_ref("ldad_constant_ref",ldad_constant_ref,&gridFromFile_ref);
    Stress<MethodLdadTrigonometric,Piola> ldad_trigonometric_stress_ref("ldad_trigonometric_ref",ldad_trigonometric_ref,&gridFromFile_ref);

    calculateStress(body,kim,std::tie(ldad_constant_stress_ref));
    ldad_constant_stress_ref.write();
    calculateStress(body,kim,std::tie(ldad_trigonometric_stress_ref));
    ldad_trigonometric_stress_ref.write();

    Matrix3d ldadVectors_def;
    ldadVectors_def << 5.43094977840521, 0.0, 0.0,
                       0.0, 5.4852592761892621, 0.0,
                       0.0, 0.0, 5.43094977840521;

    MethodLdadConstant ldad_constant_def(ldadVectors_def);
    MethodLdadTrigonometric ldad_trigonometric_def(ldadVectors_def);

    Stress<MethodLdadConstant,Cauchy> ldad_constant_stress_def("ldad_constant_def",ldad_constant_def,&gridFromFile_def);
    Stress<MethodLdadTrigonometric,Cauchy> ldad_trigonometric_stress_def("ldad_trigonometric_def",ldad_trigonometric_def,&gridFromFile_def);

    calculateStress(body,kim,std::tie(ldad_constant_stress_def));
    ldad_constant_stress_def.write();
    calculateStress(body,kim,std::tie(ldad_trigonometric_stress_def));
    ldad_trigonometric_stress_def.write();

    compareStressComponentsToReference(ldad_constant_stress_ref,"ldad_constant_refReference.stress");
    compareStressComponentsToReference(ldad_constant_stress_def,"ldad_constant_defReference.stress");
    compareStressComponentsToReference(ldad_trigonometric_stress_ref,"ldad_trigonometric_refReference.stress");
    compareStressComponentsToReference(ldad_trigonometric_stress_def,"ldad_trigonometric_defReference.stress");
    return 0;
}
