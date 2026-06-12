#include "BoxConfiguration.h"
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
    Matrix3d expectedReferenceBox= Matrix3d::Zero();
    expectedReferenceBox.col(0)= Vector3d(10.0,0.0,0.0);
    expectedReferenceBox.col(1)= Vector3d(0.0,20.0,0.0);
    expectedReferenceBox.col(2)= Vector3d(0.0,0.0,30.0);

    Matrix3d expectedCurrentBox= Matrix3d::Zero();
    expectedCurrentBox.col(0)= Vector3d(11.0,0.0,0.0);
    expectedCurrentBox.col(1)= Vector3d(0.0,21.0,0.0);
    expectedCurrentBox.col(2)= Vector3d(0.0,0.0,31.0);

    BoxConfiguration withOrigin{2,true};
    withOrigin.read("config_with_origin.data",true);
    assertVectorNear(withOrigin.reference_box_origin,Vector3d(5.0,-4.0,2.0),1e-12,
                     "Reference origin was not read from native config.");
    assertVectorNear(withOrigin.box_origin,Vector3d(5.0,-4.0,2.0),1e-12,
                     "Current origin was not read from native config.");
    assertMatrixNear(withOrigin.reference_box,expectedReferenceBox,1e-12,
                     "Reference box was parsed incorrectly after origin line.");
    assertMatrixNear(withOrigin.box,expectedCurrentBox,1e-12,
                     "Current box was parsed incorrectly after origin line.");
    assertVectorNear(withOrigin.coordinates.at(Current).row(0),Vector3d(6.0,-2.0,5.0),1e-12,
                     "Current coordinate was parsed incorrectly after origin line.");
    assertVectorNear(withOrigin.coordinates.at(Reference).row(1),Vector3d(7.0,-1.0,6.0),1e-12,
                     "Reference coordinate was parsed incorrectly after origin line.");
    assertNear(withOrigin.masses(0),39.948,1e-12,"Mass was parsed incorrectly after origin line.");

    BoxConfiguration legacy{1,true};
    legacy.read("config_legacy.data",true);
    assertVectorNear(legacy.reference_box_origin,Vector3d::Zero(),1e-12,
                     "Legacy native config should default reference origin to zero.");
    assertVectorNear(legacy.box_origin,Vector3d::Zero(),1e-12,
                     "Legacy native config should default current origin to zero.");
    assertMatrixNear(legacy.reference_box,expectedReferenceBox,1e-12,
                     "Legacy reference box was parsed incorrectly.");
    assertMatrixNear(legacy.box,expectedCurrentBox,1e-12,
                     "Legacy current box was parsed incorrectly.");

    return 0;
}
