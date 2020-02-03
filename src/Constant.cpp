/*
 * Constant.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: Nikhil
 */

#include "Constant.h"
#include "typedef.h"

Constant::Constant(){}

Constant::~Constant(){}

double Constant::operator()(const double& t)
{
    // Constant weighting function
    if (t < -1.0 - epsilon || t > 1.0 + epsilon)
	{
        return 0.0;
	} 
	else if ((t >= -1.0 - epsilon && t <= -1.0 + epsilon) || (t <= 1.0 + epsilon && t >= 1.0 - epsilon))
    {
		return 0.5;
	} 
	else
    {
		return 1.0;
	}

}
double Constant::integrate(const Vector3d& vec1_pull_seg, const Vector3d& vec2_pull_seg)
{
	double result = 1.0;
	Vector3d vec_mid;
	vec_mid = (vec1_pull_seg + vec2_pull_seg) / 2.0;

	//debug
	//std::cout << "vec_mid: " << vec_mid << std::endl;

    for (int i = 0; i <= 2; i++)
	{
    	if (vec_mid(i) < -1.0 - epsilon || vec_mid(i) > 1.0 + epsilon)
		{
		    return 0.0;
	    } 
	    else if ((vec_mid(i) >= -1.0 - epsilon && vec_mid(i) <= -1.0 + epsilon) || (vec_mid(i) <= 1.0 + epsilon && vec_mid(i) >= 1.0 - epsilon))
        {
		    result = result * 0.5;
	    } 
	    else
        {
		    result = result * 1.0;
	    }
	}
	//std::cout << "result in oneDFunction: " << result << std::endl;
	return result;
}




