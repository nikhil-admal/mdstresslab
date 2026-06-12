#include "BoxConfiguration.h"
#include "Grid.h"
#include "typedef.h"
#include <cmath>

namespace {
void assertNear(double value, double expected, double tolerance, const std::string& message)
{
    if (std::abs(value-expected)>tolerance)
        MY_ERROR(message + ": expected " + std::to_string(expected) + ", got " + std::to_string(value));
}

void assertVectorNear(const Vector3d& value, const Vector3d& expected, double tolerance, const std::string& message)
{
    if ((value-expected).cwiseAbs().maxCoeff()>tolerance)
        MY_ERROR(message);
}

void assertMatrixNear(const Matrix3d& value, const Matrix3d& expected, double tolerance, const std::string& message)
{
    if ((value-expected).cwiseAbs().maxCoeff()>tolerance)
        MY_ERROR(message);
}
}

int main()
{
    BoxConfiguration body{2,false};
    body.readLMP("triclinic.lmp",Current);

    Matrix3d expectedBox= Matrix3d::Zero();
    expectedBox.col(0)= Vector3d(10.0,0.0,0.0);
    expectedBox.col(1)= Vector3d(2.0,20.0,0.0);
    expectedBox.col(2)= Vector3d(-1.0,3.0,30.0);
    assertMatrixNear(body.box,expectedBox,1e-12,"Triclinic cell matrix was parsed incorrectly.");
    assertVectorNear(body.box_origin,Vector3d(5.0,-4.0,2.0),1e-12,
                     "Triclinic cell origin was parsed incorrectly.");

    assertVectorNear(body.coordinates.at(Current).row(0),Vector3d(6.0,-2.0,5.0),1e-12,
                     "Atom 1 coordinate was parsed incorrectly.");
    assertVectorNear(body.coordinates.at(Current).row(1),Vector3d(15.0,-13.0,68.0),1e-12,
                     "Atom 2 image flags were not unwrapped using the full triclinic cell.");

    assertNear(body.masses(0),39.948,1e-12,"Atom 1 mass was parsed incorrectly.");
    assertNear(body.masses(1),39.948,1e-12,"Atom 2 mass was parsed incorrectly.");
    assertVectorNear(body.velocities.row(0),Vector3d(0.1,0.2,0.3),1e-12,
                     "Atom 1 velocity was parsed incorrectly.");
    assertVectorNear(body.velocities.row(1),Vector3d(-0.4,0.5,-0.6),1e-12,
                     "Atom 2 velocity was parsed incorrectly.");

    Vector3d upperCorner= body.box_origin + body.box.col(0).transpose() +
                          body.box.col(1).transpose() + body.box.col(2).transpose();
    Grid<Current> grid(body.box_origin,body.box,body.box_origin,upperCorner,2,2,2);
    assertVectorNear(grid.coordinates[0],Vector3d(5.0,-4.0,2.0),1e-12,
                     "Cell-aligned grid origin is incorrect.");
    assertVectorNear(grid.coordinates[4],Vector3d(10.0,-4.0,2.0),1e-12,
                     "Cell-aligned grid spacing along the first cell vector is incorrect.");
    assertVectorNear(grid.coordinates[2],Vector3d(6.0,6.0,2.0),1e-12,
                     "Cell-aligned grid spacing along the second cell vector is incorrect.");
    assertVectorNear(grid.coordinates[1],Vector3d(4.5,-2.5,17.0),1e-12,
                     "Cell-aligned grid spacing along the third cell vector is incorrect.");

    return 0;
}
