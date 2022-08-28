#include "MethodUser.h"
#include "Polynomial.h"
#include <map>
#include "typedef.h"

class Custom: public MethodUser
{
public:
	Custom(double);
	double operator()(const Vector3d& vec) const override;
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2) const override;


private:
	std::map<std::pair<double,double>,Polynomial> piecewisePolynomial;
	double normalizer;
	double integratePolynomial(const int&, const double&,const double&,const double&,const double&) const;
	double integrate(const double&,const double&,const double&,const double&) const;
};
