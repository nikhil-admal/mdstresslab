#include "Custom.h"


Custom::Custom(double averagingDomainSize) : MethodUser(averagingDomainSize)
{ 
    
    // weighting function:

    //w(r)= c if r<= R/2,
    //	= 2c(1-r/R) if R/2<=r<=R
    //where c= 8/(5\pi R^3)

    normalizer= 8/(5*M_PI*pow(averagingDomainSize,3));
    // constant polynomial
    std::deque<double> constant{normalizer};
    Polynomial constantPolynomial{constant};
    //linear polynomial
    std::deque<double> linear{2*normalizer,-2*normalizer/averagingDomainSize};
    Polynomial linearPolynomial{linear};


    // construct piecewise polynomial
    std::pair<double,double> interval;
    interval.first= 0;
    interval.second= averagingDomainSize/2;
    piecewisePolynomial[interval]= constantPolynomial;

    interval.first= averagingDomainSize/2;
    interval.second= averagingDomainSize;
    piecewisePolynomial[interval]= linearPolynomial;
}


double Custom::integratePolynomial(const int& degree,
								   const double& a,
								   const double& b,
								   const double& r1,
								   const double& r2) const
{
	assert(r1<=r2);
	switch(degree)
	{
		case 0:
		{
			assert (a*pow(r2,2)+b+epsilon>0 && a*pow(r1,2)+b+epsilon>0);
			return (std::sqrt(a*pow(r2,2)+b+epsilon) - std::sqrt(a*pow(r1,2)+b+epsilon))/a;
		}

		case 1:
		{
			assert(a*r1*r1+b+epsilon>0 && a*r2*r2+b+epsilon>0 && a>0);
			double temp1= std::sqrt(a)*std::sqrt(a*r1*r1 + b + epsilon);
			double temp2= std::sqrt(a)*std::sqrt(a*r2*r2 + b + epsilon);

			return ((r2*temp2-b*log(temp2+a*r2))-
					   (r1*temp1-b*log(temp1+a*r1)))/(2*pow(a,1.5));
		}
		default:
			MY_ERROR("Attempt to integrate polynomial of degree not equal to 0 or 1.");
	}

}

double Custom::operator()(const Vector3d& vec) const
{
	double r= vec.norm();

    auto iter = std::find_if(piecewisePolynomial.cbegin(), piecewisePolynomial.cend(),
                          [=](const std::pair< std::pair<double,double>, Polynomial>& polynomial)
                          {
                             return r>=polynomial.first.first &&  r<polynomial.first.second;
                          });

    if (iter == piecewisePolynomial.end()) MY_ERROR("Argument out of range for the piecewise-defined weight function");
    return (iter->second)(r);
}

double Custom::bondFunction(const Vector3d& vec1, const Vector3d& vec2) const
{
	double r1= vec1.norm();
	double r2= vec2.norm();

	double a= 4*(vec1-vec2).squaredNorm();
	double b= 4*pow(vec2.dot(vec1)-vec1.dot(vec1),2)-a*r1*r1;

	double sPerp= (vec1.dot(vec1)-vec2.dot(vec1))/(vec1-vec2).squaredNorm();

	Vector3d vecPerp;
	vecPerp= (1-sPerp)*vec1+sPerp*vec2;
	double rPerp= vecPerp.norm();
	assert(rPerp<std::min(r1,r2)+epsilon);

    double averagingDomainSize= getAveragingDomainSize();
	if (rPerp>=averagingDomainSize) return 0;

	r1= std::min(r1,averagingDomainSize);
	r2= std::min(r2,averagingDomainSize);

	double rmin;
	double rmax;
	double result= 0;
	if (sPerp<0 || sPerp>1)
	{
		rmin= std::min(r1,r2);
		rmax= std::max(r1,r2);

		if (abs(rmax-rmin)>epsilon) return integrate(a,b,rmin,rmax);
		else return 0;
	}
	else
	{
		rmin= rPerp;
		if (abs(r1-rmin)>epsilon) result= integrate(a,b,rmin,r1);
		if (abs(r2-rmin)>epsilon) result= result+integrate(a,b,rmin,r2);

		return result;
	}
	assert(0);
}

double Custom::integrate(const double& a,
					    const double& b,
					    const double& rmin,
					    const double& rmax) const
{
    double averagingDomainSize= getAveragingDomainSize();
	assert(rmin<=rmax && rmax<=averagingDomainSize);
	double result= 0;
    auto iterMin = std::find_if(piecewisePolynomial.cbegin(), piecewisePolynomial.cend(),
                          [=](const std::pair< std::pair<double,double>, Polynomial>& polynomial)
                          {
                             return rmin>=polynomial.first.first &&  rmin<polynomial.first.second;
                          });
    if (iterMin == piecewisePolynomial.cend()) MY_ERROR("Lower integral limit out of range for the piecewise-defined weight function");
    auto iterMax = std::find_if(iterMin, piecewisePolynomial.cend(),
                          [=](const std::pair< std::pair<double,double>, Polynomial>& polynomial)
                          {
                             return rmax>=polynomial.first.first &&  rmax<polynomial.first.second;
                          });

	for(auto it=iterMin; it!=piecewisePolynomial.cend(); it++)
	{
		double lb,ub;
		if (it == iterMin) lb= rmin;
		else lb= it->first.first;
		if (it == iterMax) ub= rmax;
		else ub= it->first.second;

		int degree= 0;
		for(const auto& coefficient : it->second.coefficients)
		{
			result= result + 2*coefficient*integratePolynomial(degree,a,b,lb,ub);
			degree++;
		}

		if (it == iterMax) break;
	}
	return result;
}


