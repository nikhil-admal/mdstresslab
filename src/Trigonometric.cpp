/*
 * Trigonometric.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: Nikhil
 */

#include "Trigonometric.h"
#include "typedef.h"
#include <math.h>

Trigonometric::Trigonometric(){}

Trigonometric::~Trigonometric(){}

double Trigonometric::operator()(const double& t)
{
    // Trigonometric weighting function
    if (t <= -1 || t >= 1)
	{
		return 0.0;
	} 
	else
	{
		return 0.5 * (1.0 + cos(t * M_PI));
	}
}

double Trigonometric::integrate(const Vector3d& vec1_pull_seg, const Vector3d& vec2_pull_seg)
{
    Vector3d a, b;
	double result = 1.0;
    a = (vec1_pull_seg - vec2_pull_seg) * M_PI;
	b = vec2_pull_seg * M_PI;

    // 0
	// 1
	// 2
	for (int i = 0; i <= 2; i++)
	{
		if (fabs(a(i)) < epsilon)
		{
			result = result + cos(b(i));
		}
    	else
		{
			result = result + (sin(a(i) + b(i)) - sin(b(i))) / a(i);
		}
	}

    // 01
    if (fabs(a(0) - a(1)) < epsilon)
    {
        result = result + cos(b(0) - b(1)) / 2.0;
	}
	else
	{
        result = result + (sin(b(0) - b(1) + a(0) - a(1)) - sin(b(0) - b(1))) / ((a(0) - a(1)) * 2.0);
	}

	if (fabs(a(0) + a(1)) < epsilon)
    {
        result = result + cos(b(0) + b(1)) / 2.0;
	}
	else
	{
        result = result + (sin(b(0) + b(1) + a(0) + a(1)) - sin(b(0) + b(1))) / ((a(0) + a(1)) * 2.0);
	}

	// 02
    if (fabs(a(0) - a(2)) < epsilon)
    {
        result = result + cos(b(0) - b(2)) / 2.0;
	}
	else
	{
        result = result + (sin(b(0) - b(2) + a(0) - a(2)) - sin(b(0) - b(2))) / ((a(0) - a(2)) * 2.0);
	}

	if (fabs(a(0) + a(2)) < epsilon)
    {
        result = result + cos(b(0) + b(2)) / 2.0;
	}
	else
	{
        result = result + (sin(b(0) + b(2) + a(0) + a(2)) - sin(b(0) + b(2))) / ((a(0) + a(2)) * 2.0);
	}

	// 12
    if (fabs(a(1) - a(2)) < epsilon)
    {
        result = result + cos(b(1) - b(2)) / 2.0;
	}
	else
	{
        result = result + (sin(b(1) - b(2) + a(1) - a(2)) - sin(b(1) - b(2))) / ((a(1) - a(2)) * 2.0);
	}

	if (fabs(a(1) + a(2)) < epsilon)
    {
        result = result + cos(b(1) + b(2)) / 2.0;
	}
	else
	{
        result = result + (sin(b(1) + b(2) + a(1) + a(2)) - sin(b(1) + b(2))) / ((a(1) + a(2)) * 2.0);
	}

	// 012
    if (fabs(a(0) - a(1) - a(2)) < epsilon)
	{
        result = result + cos(b(0) - b(1) - b(2)) / 4.0;
	}
    else
	{
		result = result + (sin(b(0) - b(1) - b(2) + a(0) - a(1) - a(2)) - sin(b(0) - b(1) - b(2))) / ((a(0) - a(1) - a(2)) * 4.0);
	}
	
	if (fabs(a(0) + a(1) - a(2)) < epsilon)
	{
        result = result + cos(b(0) + b(1) - b(2)) / 4.0;
	}
    else
	{
		result = result + (sin(b(0) + b(1) - b(2) + a(0) + a(1) - a(2)) - sin(b(0) + b(1) - b(2))) / ((a(0) + a(1) - a(2)) * 4.0);
	}

    if (fabs(a(0) - a(1) + a(2)) < epsilon)
	{
        result = result + cos(b(0) - b(1) + b(2)) / 4.0;
	}
    else
	{
		result = result + (sin(b(0) - b(1) + b(2) + a(0) - a(1) + a(2)) - sin(b(0) - b(1) + b(2))) / ((a(0) - a(1) + a(2)) * 4.0);
	}

    if (fabs(a(0) + a(1) + a(2)) < epsilon)
	{
        result = result + cos(b(0) + b(1) + b(2)) / 4.0;
	}
    else
	{
		result = result + (sin(b(0) + b(1) + b(2) + a(0) + a(1) + a(2)) - sin(b(0) + b(1) + b(2))) / ((a(0) + a(1) + a(2)) * 4.0);
	}
	
	return result;
}




