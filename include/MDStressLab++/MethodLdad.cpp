/*
 * MethodLdad.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#include <math.h>
#include <iostream>
#include <limits>
#include "typedef.h"

int PointLineRelationship(const double& p);

template<typename T>
MethodLdad<T>::MethodLdad(const Matrix3d& ldadVectors): ldadVectors(ldadVectors)
{
	// initialize the normalizer here using oneDFunction.integrate(-1,1)
	normalizer = 1.0/(8.0*fabs(ldadVectors.determinant())) ;

    // initialize the averagingDomainSize
    double p1, p2, p3;
	p1 = fabs(ldadVectors(0,0)) + fabs(ldadVectors(0,1)) + fabs(ldadVectors(0,2));
    p2 = fabs(ldadVectors(1,0)) + fabs(ldadVectors(1,1)) + fabs(ldadVectors(1,2));
	p3 = fabs(ldadVectors(2,0)) + fabs(ldadVectors(2,1)) + fabs(ldadVectors(2,2));

    
    this->averagingDomainSize = sqrt(p1 * p1 + p2 * p2 + p3 * p3);
	
    // initialize the inverseLdadVectors
    inverseLdadVectors = ldadVectors.inverse();
}

// template<typename T>
// MethodLdad<T>::MethodLdad(const MethodLdad<T>& _MethodLdad)
// {
// 	*this= _MethodLdad;
// }

template<typename T>
MethodLdad<T>::~MethodLdad() {
	// TODO Auto-generated destructor stub
}

template<typename T>
double MethodLdad<T>::operator()(const Vector3d& vec) const
{
	Vector3d vec_pull;
	vec_pull = inverseLdadVectors * vec.transpose();
	// oneDFunction -1 1 
	return oneDFunction(vec_pull(0))*oneDFunction(vec_pull(1))*oneDFunction(vec_pull(2));
}

template<typename T>
double MethodLdad<T>::bondFunction(const Vector3d& vec1, const Vector3d& vec2) const
{
	/* copy from fortran 
	   The algorithm is:
       selecting two points out of eight points.
	   The eight points composed of two points of vec1 and vec2 themselves,
	   and the intersection of x = +-1, y = +-1, z = +-1.	*/   

    Vector3d vec1_pull, vec2_pull, vec12_pull;
	Vector3d vec1_pull_seg, vec2_pull_seg;
	Vector3d vec1_push_seg, vec2_push_seg, vec12_push_seg;

	Vector3i degenerate(0,0,0);
	VectorXi selected(8);
	int selected_count = 0;
	Matrix3i relative_position;

    ArrayXXd r_pull_intersect(8,3);
	double total_length = 0.0;
    double distance = (vec2 - vec1).norm();

    for (int i = 0; i < 8;i++)
    {
		for (int j = 0; j < 3; j++)
		{
			r_pull_intersect(i,j) = 0.0;
		}
		selected(i) = 0;
	}

	vec1_pull = inverseLdadVectors * vec1.transpose();
	vec2_pull = inverseLdadVectors * vec2.transpose();
    vec12_pull = vec2_pull - vec1_pull;
    // This is a necessary step to initialization
	vec1_pull_seg = vec1_pull;
	vec2_pull_seg = vec2_pull;

    for (int i = 0; i <= 2; i++)
    {
		if (fabs(vec12_pull(i)) < epsilon)
		{
            degenerate(i) = 1;
		} 
		else
		{
			degenerate(i) = 0;
		}
	}

    relative_position(0,0) = PointLineRelationship(vec1_pull(0));
	relative_position(1,0) = PointLineRelationship(vec1_pull(1));
	relative_position(2,0) = PointLineRelationship(vec1_pull(2));
    relative_position(0,1) = PointLineRelationship(vec2_pull(0));
	relative_position(1,1) = PointLineRelationship(vec2_pull(1));
	relative_position(2,1) = PointLineRelationship(vec2_pull(2));
	relative_position(0,2) = 0;
	relative_position(1,2) = 0;
	relative_position(2,2) = 0;

	// degenerate to 0D. It is impossible
    if (degenerate(0) && degenerate(1) && degenerate(2))
    {
        MY_ERROR("Degeneration to 0D. This indicates you have overlapping atoms in your system!");
	}
	// degenerate to 1D x
    else if (!degenerate(0) && degenerate(1) && degenerate(2))
	{
		relative_position(1,2) = PointLineRelationship((vec1_pull(1) + vec2_pull(1)) / 2.0);
	    relative_position(2,2) = PointLineRelationship((vec1_pull(2) + vec2_pull(2)) / 2.0);
		if (relative_position(1,2) == 1 || relative_position(1,2) == 5 || \
		    relative_position(2,2) == 1 || relative_position(2,2) == 5)
		{
			return 0.0;
		}
        // determine x length
		if (relative_position(0,0) == 2 || relative_position(0,0) == 3 || relative_position(0,0) == 4)
		{
            vec1_pull_seg(0) = vec1_pull(0);
			if (relative_position(0,1) == 2 || relative_position(0,1) == 3 || relative_position(0,1) == 4)
			{
                vec2_pull_seg(0) = vec2_pull(0);
			}
			else if (relative_position(0,1) == 1)
			{
                vec2_pull_seg(0) = -1.0;
			}
            else if (relative_position(0,1) == 5)
			{
                vec2_pull_seg(0) = 1.0;
			}
		}
		else if (relative_position(0,0) == 1)
		{
			if (relative_position(0,1) == 2 || relative_position(0,1) == 3 || relative_position(0,1) == 4)
			{
				vec1_pull_seg(0) = -1.0;
				vec2_pull_seg(0) = vec2_pull(0);
			}
		    else if (relative_position(0,1) == 1)
			{
                return 0.0;
			}
			else if (relative_position(0,1) == 5)
			{
				vec1_pull_seg(0) = -1.0;
				vec2_pull_seg(0) = 1.0;
			}
		}
		else if (relative_position(0,0) == 5)
		{
            if (relative_position(0,1) == 2 || relative_position(0,1) == 3 || relative_position(0,1) == 4)
			{
				vec1_pull_seg(0) = 1.0;
				vec2_pull_seg(0) = vec2_pull(0);
			}
			else if (relative_position(0,1) == 1)
			{
				vec1_pull_seg(0) = 1.0;
				vec2_pull_seg(0) = -1.0;
			}
			else if (relative_position(0,1) == 5)
			{
				return 0.0;
			}
		}
	}
	// degenerate to 1D y
    else if (degenerate(0) && !degenerate(1) && degenerate(2))
    {
		relative_position(0,2) = PointLineRelationship((vec1_pull(0) + vec2_pull(0)) / 2.0);
	    relative_position(2,2) = PointLineRelationship((vec1_pull(2) + vec2_pull(2)) / 2.0);
		if (relative_position(0,2) == 1 || relative_position(0,2) == 5 || \
		    relative_position(2,2) == 1 || relative_position(2,2) == 5)
		{
			return 0.0;
		}
        // determine y length
		if (relative_position(1,0) == 2 || relative_position(1,0) == 3 || relative_position(1,0) == 4)
		{
            vec1_pull_seg(1) = vec1_pull(1);
			if (relative_position(1,1) == 2 || relative_position(1,1) == 3 || relative_position(1,1) == 4)
			{
                vec2_pull_seg(1) = vec2_pull(1);
			}
			else if (relative_position(1,1) == 1)
			{
                vec2_pull_seg(1) = -1.0;
			}
            else if (relative_position(1,1) == 5)
			{
                vec2_pull_seg(1) = 1.0;
			}
		}
		else if (relative_position(1,0) == 1)
		{
			if (relative_position(1,1) == 2 || relative_position(1,1) == 3 || relative_position(1,1) == 4)
			{
				vec1_pull_seg(1) = -1.0;
				vec2_pull_seg(1) = vec2_pull(1);
			}
		    else if (relative_position(1,1) == 1)
			{
                return 0.0;
			}
			else if (relative_position(1,1) == 5)
			{
				vec1_pull_seg(1) = -1.0;
				vec2_pull_seg(1) = 1.0;
			}
		}
		else if (relative_position(1,0) == 5)
		{
            if (relative_position(1,1) == 2 || relative_position(1,1) == 3 || relative_position(1,1) == 4)
			{
				vec1_pull_seg(1) = 1.0;
				vec2_pull_seg(1) = vec2_pull(1);
			}
			else if (relative_position(1,1) == 1)
			{
				vec1_pull_seg(1) = 1.0;
				vec2_pull_seg(1) = -1.0;
			}
			else if (relative_position(1,1) == 5)
			{
				return 0.0;
			}
		}
	}
    // degenerate to 1D z
	else if (degenerate(0) && degenerate(1) && !degenerate(2))
    {
		relative_position(0,2) = PointLineRelationship((vec1_pull(0) + vec2_pull(0)) / 2.0);
	    relative_position(1,2) = PointLineRelationship((vec1_pull(1) + vec2_pull(1)) / 2.0);
		if (relative_position(0,2) == 1 || relative_position(0,2) == 5 || \
		    relative_position(1,2) == 1 || relative_position(1,2) == 5)
		{
			return 0.0;
		}
        // determine z length
		if (relative_position(2,0) == 2 || relative_position(2,0) == 3 || relative_position(2,0) == 4)
		{
            vec1_pull_seg(2) = vec1_pull(2);
			if (relative_position(2,1) == 2 || relative_position(2,1) == 3 || relative_position(2,1) == 4)
			{
                vec2_pull_seg(2) = vec2_pull(2);
			}
			else if (relative_position(2,1) == 1)
			{
                vec2_pull_seg(2) = -1.0;
			}
            else if (relative_position(2,1) == 5)
			{
                vec2_pull_seg(2) = 1.0;
			}
		}
		else if (relative_position(2,0) == 1)
		{
			if (relative_position(2,1) == 2 || relative_position(2,1) == 3 || relative_position(2,1) == 4)
			{
				vec1_pull_seg(2) = -1.0;
				vec2_pull_seg(2) = vec2_pull(2);
			}
		    else if (relative_position(2,1) == 1)
			{
                return 0.0;
			}
			else if (relative_position(2,1) == 5)
			{
				vec1_pull_seg(2) = -1.0;
				vec2_pull_seg(2) = 1.0;
			}
		}
		else if (relative_position(2,0) == 5)
		{
            if (relative_position(2,1) == 2 || relative_position(2,1) == 3 || relative_position(2,1) == 4)
			{
				vec1_pull_seg(2) = 1.0;
				vec2_pull_seg(2) = vec2_pull(2);
			}
			else if (relative_position(2,1) == 1)
			{
				vec1_pull_seg(2) = 1.0;
				vec2_pull_seg(2) = -1.0;
			}
			else if (relative_position(2,1) == 5)
			{
				return 0.0;
			}
		}
	}
    // degenerate to 2D xy
	else if (!degenerate(0) && !degenerate(1) && degenerate(2))
	{
		relative_position(2,2) = PointLineRelationship((vec1_pull(2) + vec2_pull(2)) / 2.0);
		if (relative_position(2,2) == 1 || relative_position(2,2) == 5)
		{
			return 0.0;
		}

		for (int i = 0; i < 8; i++)
        {
			selected(i) = 0;
		}

        // x = -1.0
        r_pull_intersect(0,0) = -1.0;
		r_pull_intersect(0,1) = vec1_pull(1) + \
		                       (vec2_pull(1) - vec1_pull(1)) / (vec2_pull(0) - vec1_pull(0)) * (-1.0 - vec1_pull(0));
        // x = 1.0
        r_pull_intersect(1,0) = 1.0;
        r_pull_intersect(1,1) = vec1_pull(1) + \
		                       (vec2_pull(1) - vec1_pull(1)) / (vec2_pull(0) - vec1_pull(0)) * (1.0 - vec1_pull(0));
        // y = -1.0
		r_pull_intersect(2,0) = vec1_pull(0) + \
		                       (vec2_pull(0) - vec1_pull(0)) / (vec2_pull(1) - vec1_pull(1)) * (-1.0 - vec1_pull(1));
		r_pull_intersect(2,1) = -1.0;
		// y = 1.0
        r_pull_intersect(3,0) = vec1_pull(0) + \
		                       (vec2_pull(0) - vec1_pull(0)) / (vec2_pull(1) - vec1_pull(1)) * (1.0 - vec1_pull(1));
        r_pull_intersect(3,1) = 1.0;

		// vec1 itself
		r_pull_intersect(6,0) = vec1_pull(0);
		r_pull_intersect(6,1) = vec1_pull(1);
		r_pull_intersect(6,2) = vec1_pull(2);

		// vec2 itself
		r_pull_intersect(7,0) = vec2_pull(0);
		r_pull_intersect(7,1) = vec2_pull(1);
		r_pull_intersect(7,2) = vec2_pull(2);

        // if the point is inside, then it is the point selected.
        if ((relative_position(0,0) == 2 || relative_position(0,0) == 3 || relative_position(0,0) == 4) && \
		    (relative_position(1,0) == 2 || relative_position(1,0) == 3 || relative_position(1,0) == 4))
		{
			selected(6) = 1;
		}

        if ((relative_position(0,1) == 2 || relative_position(0,1) == 3 || relative_position(0,1) == 4) && \
		    (relative_position(1,1) == 2 || relative_position(1,1) == 3 || relative_position(1,1) == 4))
		{
			selected(7) = 1;
		}

		// must be inside cubes to be selected
        for (int i = 0; i <= 1; i++)
		{
            if (((r_pull_intersect(i,0) > vec1_pull(0) && r_pull_intersect(i,0) < vec2_pull(0)) ||  \
			     (r_pull_intersect(i,0) > vec2_pull(0) && r_pull_intersect(i,0) < vec1_pull(0))) && \
				((r_pull_intersect(i,1) > vec1_pull(1) && r_pull_intersect(i,1) < vec2_pull(1)) ||  \
				 (r_pull_intersect(i,1) > vec2_pull(1) && r_pull_intersect(i,1) < vec1_pull(1))) && \
				 (r_pull_intersect(i,1) >= -1.0 - epsilon && r_pull_intersect(i,1) <= 1.0 + epsilon))
            {
				selected(i) = 1;
			}
		}
        for (int i = 2; i <= 3; i++)
        {
            if (((r_pull_intersect(i,0) > vec1_pull(0) && r_pull_intersect(i,0) < vec2_pull(0)) ||  \
			     (r_pull_intersect(i,0) > vec2_pull(0) && r_pull_intersect(i,0) < vec1_pull(0))) && \
				((r_pull_intersect(i,1) > vec1_pull(1) && r_pull_intersect(i,1) < vec2_pull(1)) ||  \
				 (r_pull_intersect(i,1) > vec2_pull(1) && r_pull_intersect(i,1) < vec1_pull(1))) && \
		         (r_pull_intersect(i,0) >= -1.0 - epsilon && r_pull_intersect(i,0) <= 1.0 + epsilon))
            {
				selected(i) = 1;
			}
		}

		selected_count = 0;
		for (int i = 0; i < 8; i++)
		{      
            if (selected(i) && selected_count == 0)
			{
                vec1_pull_seg(0) = r_pull_intersect(i,0);
				vec1_pull_seg(1) = r_pull_intersect(i,1);
				selected_count = 1;
				continue;
			}
            if (selected(i) && selected_count == 1)
			{
				if (fabs(vec1_pull_seg(0) - r_pull_intersect(i,0)) < epsilon && \
				    fabs(vec1_pull_seg(1) - r_pull_intersect(i,1)) < epsilon)
				{
					continue;
				}
				else
				{
					vec2_pull_seg(0) = r_pull_intersect(i,0);
				    vec2_pull_seg(1) = r_pull_intersect(i,1);
				    selected_count = 2;
				    continue;
				}
			}
            if (selected(i) && selected_count == 2)
			{
			    if ((fabs(vec1_pull_seg(0) - r_pull_intersect(i,0)) < epsilon &&  \
				     fabs(vec1_pull_seg(1) - r_pull_intersect(i,1)) < epsilon) || \
					(fabs(vec2_pull_seg(0) - r_pull_intersect(i,0)) < epsilon &&  \
				     fabs(vec2_pull_seg(1) - r_pull_intersect(i,1)) < epsilon))
				{
					continue;
				}
				else
				{
					MY_ERROR("In LDAD bond_function xy, More than 2 different points \
					          with different positions are being selected. \
							  This means a line intersecting with a cube has 3 points. \
							  Something wrong with the algorithm. Need to print and debug!");
				}  
			}
		}
        // No intersection .or. one point on the parallelepiped but another point is outside
		if (selected_count == 0 || selected_count == 1)
		{
			return 0.0;
		}
	}
	// degenerate to 2D xz
	else if (!degenerate(0) && degenerate(1) && !degenerate(2))
	{
		relative_position(1,2) = PointLineRelationship((vec1_pull(1) + vec2_pull(1)) / 2.0);
		if (relative_position(1,2) == 1 || relative_position(1,2) == 5)
		{
			return 0.0;
		}

		for (int i = 0; i < 8; i++)
        {
			selected(i) = 0;
		}

        // x = -1.0
        r_pull_intersect(0,0) = -1.0;
		r_pull_intersect(0,2) = vec1_pull(2) + \
		                       (vec2_pull(2) - vec1_pull(2)) / (vec2_pull(0) - vec1_pull(0)) * (-1.0 - vec1_pull(0));
        // x = 1.0
        r_pull_intersect(1,0) = 1.0;
        r_pull_intersect(1,2) = vec1_pull(2) + \
		                       (vec2_pull(2) - vec1_pull(2)) / (vec2_pull(0) - vec1_pull(0)) * (1.0 - vec1_pull(0));
        // z = -1.0
		r_pull_intersect(4,0) = vec1_pull(0) + \
		                       (vec2_pull(0) - vec1_pull(0)) / (vec2_pull(2) - vec1_pull(2)) * (-1.0 - vec1_pull(2));
		r_pull_intersect(4,2) = -1.0;
		// z = 1.0
        r_pull_intersect(5,0) = vec1_pull(0) + \
		                       (vec2_pull(0) - vec1_pull(0)) / (vec2_pull(2) - vec1_pull(2)) * (1.0 - vec1_pull(2));
        r_pull_intersect(5,2) = 1.0;

		// vec1 itself
		r_pull_intersect(6,0) = vec1_pull(0);
		r_pull_intersect(6,1) = vec1_pull(1);
		r_pull_intersect(6,2) = vec1_pull(2);

		// vec2 itself
		r_pull_intersect(7,0) = vec2_pull(0);
		r_pull_intersect(7,1) = vec2_pull(1);
		r_pull_intersect(7,2) = vec2_pull(2);

        // if the point is inside, then it is the point selected.
        if ((relative_position(0,0) == 2 || relative_position(0,0) == 3 || relative_position(0,0) == 4) && \
			(relative_position(2,0) == 2 || relative_position(2,0) == 3 || relative_position(2,0) == 4))
		{
			selected(6) = 1;
		}

        if ((relative_position(0,1) == 2 || relative_position(0,1) == 3 || relative_position(0,1) == 4) && \
			(relative_position(2,1) == 2 || relative_position(2,1) == 3 || relative_position(2,1) == 4))
		{
			selected(7) = 1;
		}

		// must be inside cubes to be selected
        for (int i = 0; i <= 1; i++)
		{
            if (((r_pull_intersect(i,0) > vec1_pull(0) && r_pull_intersect(i,0) < vec2_pull(0)) ||  \
			     (r_pull_intersect(i,0) > vec2_pull(0) && r_pull_intersect(i,0) < vec1_pull(0))) && \
				((r_pull_intersect(i,2) > vec1_pull(2) && r_pull_intersect(i,2) < vec2_pull(2)) ||  \
				 (r_pull_intersect(i,2) > vec2_pull(2) && r_pull_intersect(i,2) < vec1_pull(2))) && \
				 (r_pull_intersect(i,2) >= -1.0 - epsilon && r_pull_intersect(i,2) <= 1.0 + epsilon))
            {
				selected(i) = 1;
			}
		}
        for (int i = 4; i <= 5; i++)
        {
            if (((r_pull_intersect(i,0) > vec1_pull(0) && r_pull_intersect(i,0) < vec2_pull(0)) ||  \
			     (r_pull_intersect(i,0) > vec2_pull(0) && r_pull_intersect(i,0) < vec1_pull(0))) && \
				((r_pull_intersect(i,2) > vec1_pull(2) && r_pull_intersect(i,2) < vec2_pull(2)) ||  \
				 (r_pull_intersect(i,2) > vec2_pull(2) && r_pull_intersect(i,2) < vec1_pull(2))) && \
		         (r_pull_intersect(i,0) >= -1.0 - epsilon && r_pull_intersect(i,0) <= 1.0 + epsilon))
            {
				selected(i) = 1;
			}
		}

		selected_count = 0;
		for (int i = 0; i < 8; i++)
		{      
            if (selected(i) && selected_count == 0)
			{
                vec1_pull_seg(0) = r_pull_intersect(i,0);
				vec1_pull_seg(2) = r_pull_intersect(i,2);
				selected_count = 1;
				continue;
			}
            if (selected(i) && selected_count == 1)
			{
				if (fabs(vec1_pull_seg(0) - r_pull_intersect(i,0)) < epsilon && \
				    fabs(vec1_pull_seg(2) - r_pull_intersect(i,2)) < epsilon)
				{
					continue;
				}
				else
				{
					vec2_pull_seg(0) = r_pull_intersect(i,0);
				    vec2_pull_seg(2) = r_pull_intersect(i,2);
				    selected_count = 2;
				    continue;
				}
			}
            if (selected(i) && selected_count == 2)
			{
			    if ((fabs(vec1_pull_seg(0) - r_pull_intersect(i,0)) < epsilon &&  \
				     fabs(vec1_pull_seg(2) - r_pull_intersect(i,2)) < epsilon) || \
					(fabs(vec2_pull_seg(0) - r_pull_intersect(i,0)) < epsilon &&  \
				     fabs(vec2_pull_seg(2) - r_pull_intersect(i,2)) < epsilon))
				{
					continue;
				}
				else
				{
					MY_ERROR("In LDAD bond_function xz, More than 2 different points \
					          with different positions are being selected. \
							  This means a line intersecting with a cube has 3 points. \
							  Something wrong with the algorithm. Need to print and debug!");
				}  
			}
		}

        // No intersection .or. one point on the parallelepiped but another point is outside
		if (selected_count == 0 || selected_count == 1)
		{
			return 0.0;
		}
	}
    // degenerate to 2D yz
    else if (degenerate(0) && !degenerate(1) && !degenerate(2))
	{
		relative_position(0,2) = PointLineRelationship((vec1_pull(0) + vec2_pull(0)) / 2.0);
		if (relative_position(0,2) == 1 || relative_position(0,2) == 5)
		{

			return 0.0;
		}

		for (int i = 0; i < 8; i++)
        {
			selected(i) = 0;
		}

		// y = -1.0
        r_pull_intersect(2,1) = -1.0;
		r_pull_intersect(2,2) = vec1_pull(2) + \
		                       (vec2_pull(2) - vec1_pull(2)) / (vec2_pull(1) - vec1_pull(1)) * (-1.0 - vec1_pull(1));
        // y = 1.0
        r_pull_intersect(3,1) = 1.0;
        r_pull_intersect(3,2) = vec1_pull(2) + \
		                       (vec2_pull(2) - vec1_pull(2)) / (vec2_pull(1) - vec1_pull(1)) * (1.0 - vec1_pull(1));
        // z = -1.0
		r_pull_intersect(4,1) = vec1_pull(1) + \
		                       (vec2_pull(1) - vec1_pull(1)) / (vec2_pull(2) - vec1_pull(2)) * (-1.0 - vec1_pull(2));
		r_pull_intersect(4,2) = -1.0;
		// z = 1.0
        r_pull_intersect(5,1) = vec1_pull(1) + \
		                       (vec2_pull(1) - vec1_pull(1)) / (vec2_pull(2) - vec1_pull(2)) * (1.0 - vec1_pull(2));
        r_pull_intersect(5,2) = 1.0;

		// vec1 itself
		r_pull_intersect(6,0) = vec1_pull(0);
		r_pull_intersect(6,1) = vec1_pull(1);
		r_pull_intersect(6,2) = vec1_pull(2);

		// vec2 itself
		r_pull_intersect(7,0) = vec2_pull(0);
		r_pull_intersect(7,1) = vec2_pull(1);
		r_pull_intersect(7,2) = vec2_pull(2);

        // if the point is inside, then it is the point selected.
        if ((relative_position(1,0) == 2 || relative_position(1,0) == 3 || relative_position(1,0) == 4) && \
		    (relative_position(2,0) == 2 || relative_position(2,0) == 3 || relative_position(2,0) == 4))
		{
			selected(6) = 1;
		}

        if ((relative_position(1,1) == 2 || relative_position(1,1) == 3 || relative_position(1,1) == 4) && \
		    (relative_position(2,1) == 2 || relative_position(2,1) == 3 || relative_position(2,1) == 4))
		{
			selected(7) = 1;
		}

		// must be inside cubes to be selected
        for (int i = 2; i <= 3; i++)
		{
            if (((r_pull_intersect(i,1) > vec1_pull(1) && r_pull_intersect(i,1) < vec2_pull(1)) ||  \
			     (r_pull_intersect(i,1) > vec2_pull(1) && r_pull_intersect(i,1) < vec1_pull(1))) && \
				((r_pull_intersect(i,2) > vec1_pull(2) && r_pull_intersect(i,2) < vec2_pull(2)) ||  \
				 (r_pull_intersect(i,2) > vec2_pull(2) && r_pull_intersect(i,2) < vec1_pull(2))) && \
				 (r_pull_intersect(i,2) >= -1.0 - epsilon && r_pull_intersect(i,2) <= 1.0 + epsilon))
            {
				selected(i) = 1;
			}
		}
        for (int i = 4; i <= 5; i++)
        {
            if (((r_pull_intersect(i,1) > vec1_pull(1) && r_pull_intersect(i,1) < vec2_pull(1)) ||  \
			     (r_pull_intersect(i,1) > vec2_pull(1) && r_pull_intersect(i,1) < vec1_pull(1))) && \
				((r_pull_intersect(i,2) > vec1_pull(2) && r_pull_intersect(i,2) < vec2_pull(2)) ||  \
				 (r_pull_intersect(i,2) > vec2_pull(2) && r_pull_intersect(i,2) < vec1_pull(2))) && \
		         (r_pull_intersect(i,1) >= -1.0 - epsilon && r_pull_intersect(i,1) <= 1.0 + epsilon))
            {
				selected(i) = 1;
			}
		}

		selected_count = 0;
		for (int i = 0; i < 8; i++)
		{      
            if (selected(i) && selected_count == 0)
			{
                vec1_pull_seg(1) = r_pull_intersect(i,1);
				vec1_pull_seg(2) = r_pull_intersect(i,2);
				selected_count = 1;
				continue;
			}
            if (selected(i) && selected_count == 1)
			{
				if (fabs(vec1_pull_seg(1) - r_pull_intersect(i,1)) < epsilon && \
				    fabs(vec1_pull_seg(2) - r_pull_intersect(i,2)) < epsilon)
				{
					continue;
				}
				else
				{   
					vec2_pull_seg(1) = r_pull_intersect(i,1);
				    vec2_pull_seg(2) = r_pull_intersect(i,2);
				    selected_count = 2;
				    continue;
				}
			}
            if (selected(i) && selected_count == 2)
			{
			    if ((fabs(vec1_pull_seg(1) - r_pull_intersect(i,1)) < epsilon &&  \
				     fabs(vec1_pull_seg(2) - r_pull_intersect(i,2)) < epsilon) || \
					(fabs(vec2_pull_seg(1) - r_pull_intersect(i,1)) < epsilon &&  \
				     fabs(vec2_pull_seg(2) - r_pull_intersect(i,2)) < epsilon))
				{
					continue;
				}
				else
				{
					MY_ERROR("In LDAD bond_function yz, More than 2 different points \
					          with different positions are being selected. \
							  This means a line intersecting with a cube has 3 points. \
							  Something wrong with the algorithm. Need to print and debug!");
				}  
			}
		}

		// No intersection .or. one point on the parallelepiped but another point is outside
		if (selected_count == 0 || selected_count == 1)
		{

			return 0.0;
		}
	}
    // 3D xyz
    else if (!degenerate(0) && !degenerate(1) && !degenerate(2)) 
	{
		for (int i = 0; i < 8; i++)
        {
			selected(i) = 0;
		}

		// The 8 points
		// x = -1.0
		r_pull_intersect(0,0) = -1.0;
		r_pull_intersect(0,1) = vec1_pull(1) + \
		                       (vec2_pull(1) - vec1_pull(1)) / (vec2_pull(0) - vec1_pull(0)) * (-1.0 - vec1_pull(0));
        r_pull_intersect(0,2) = vec1_pull(2) + \
		                       (vec2_pull(2) - vec1_pull(2)) / (vec2_pull(0) - vec1_pull(0)) * (-1.0 - vec1_pull(0));
		// x = 1.0
		r_pull_intersect(1,0) = 1.0;
		r_pull_intersect(1,1) = vec1_pull(1) + \
		                       (vec2_pull(1) - vec1_pull(1)) / (vec2_pull(0) - vec1_pull(0)) * (1.0 - vec1_pull(0));
		r_pull_intersect(1,2) = vec1_pull(2) + \
		                       (vec2_pull(2) - vec1_pull(2)) / (vec2_pull(0) - vec1_pull(0)) * (1.0 - vec1_pull(0));

		// y = -1.0
		r_pull_intersect(2,0) = vec1_pull(0) + \
		                       (vec2_pull(0) - vec1_pull(0)) / (vec2_pull(1) - vec1_pull(1)) * (-1.0 - vec1_pull(1));
		r_pull_intersect(2,1) = -1.0;
		r_pull_intersect(2,2) = vec1_pull(2) + \
		                       (vec2_pull(2) - vec1_pull(2)) / (vec2_pull(1) - vec1_pull(1)) * (-1.0 - vec1_pull(1));

		// y = 1.0
		r_pull_intersect(3,0) = vec1_pull(0) + \
		                       (vec2_pull(0) - vec1_pull(0)) / (vec2_pull(1) - vec1_pull(1)) * (1.0 - vec1_pull(1));
		r_pull_intersect(3,1) = 1.0;
		r_pull_intersect(3,2) = vec1_pull(2) + \
		                       (vec2_pull(2) - vec1_pull(2)) / (vec2_pull(1) - vec1_pull(1)) * (1.0 - vec1_pull(1));

		// z = -1.0
		r_pull_intersect(4,0) = vec1_pull(0) + \
		                       (vec2_pull(0) - vec1_pull(0)) / (vec2_pull(2) - vec1_pull(2)) * (-1.0 - vec1_pull(2));
		r_pull_intersect(4,1) = vec1_pull(1) + \
		                       (vec2_pull(1) - vec1_pull(1)) / (vec2_pull(2) - vec1_pull(2)) * (-1.0 - vec1_pull(2));
		r_pull_intersect(4,2) = -1.0;

        // z = 1.0
		r_pull_intersect(5,0) = vec1_pull(0) + \
		                       (vec2_pull(0) - vec1_pull(0)) / (vec2_pull(2) - vec1_pull(2)) * (1.0 - vec1_pull(2));
		r_pull_intersect(5,1) = vec1_pull(1) + \
		                       (vec2_pull(1) - vec1_pull(1)) / (vec2_pull(2) - vec1_pull(2)) * (1.0 - vec1_pull(2));
		r_pull_intersect(5,2) = 1.0;

		// vec1 itself
		r_pull_intersect(6,0) = vec1_pull(0);
		r_pull_intersect(6,1) = vec1_pull(1);
		r_pull_intersect(6,2) = vec1_pull(2);

		// vec2 itself
		r_pull_intersect(7,0) = vec2_pull(0);
		r_pull_intersect(7,1) = vec2_pull(1);
		r_pull_intersect(7,2) = vec2_pull(2);

        // if the point is inside, then it is the point selected.
        if ((relative_position(0,0) == 2 || relative_position(0,0) == 3 || relative_position(0,0) == 4) && \
		    (relative_position(1,0) == 2 || relative_position(1,0) == 3 || relative_position(1,0) == 4) && \
			(relative_position(2,0) == 2 || relative_position(2,0) == 3 || relative_position(2,0) == 4))
		{
			selected(6) = 1;
		}

        if ((relative_position(0,1) == 2 || relative_position(0,1) == 3 || relative_position(0,1) == 4) && \
		    (relative_position(1,1) == 2 || relative_position(1,1) == 3 || relative_position(1,1) == 4) && \
			(relative_position(2,1) == 2 || relative_position(2,1) == 3 || relative_position(2,1) == 4))
		{
			selected(7) = 1;
		}
        
        for (int i = 0; i <= 5; i++)
		{
		    // see whether the intersection point is between vec1 and vec2
            if (((r_pull_intersect(i,0) > vec1_pull(0) && r_pull_intersect(i,0) < vec2_pull(0)) ||  \
			     (r_pull_intersect(i,0) > vec2_pull(0) && r_pull_intersect(i,0) < vec1_pull(0))) && \
				((r_pull_intersect(i,1) > vec1_pull(1) && r_pull_intersect(i,1) < vec2_pull(1)) ||  \
				 (r_pull_intersect(i,1) > vec2_pull(1) && r_pull_intersect(i,1) < vec1_pull(1))) && \
                ((r_pull_intersect(i,2) > vec1_pull(2) && r_pull_intersect(i,2) < vec2_pull(2)) ||  \
			     (r_pull_intersect(i,2) > vec2_pull(2) && r_pull_intersect(i,2) < vec1_pull(2))) && \
		         (r_pull_intersect(i,0) >= -1.0 - epsilon && r_pull_intersect(i,0) <= 1.0 + epsilon) && \
				 (r_pull_intersect(i,1) >= -1.0 - epsilon && r_pull_intersect(i,1) <= 1.0 + epsilon) && \
				 (r_pull_intersect(i,2) >= -1.0 - epsilon && r_pull_intersect(i,2) <= 1.0 + epsilon))
            {
				selected(i) = 1;
			}
		}

        selected_count = 0;

        for (int i = 0; i < 8; i++)
		{
            if (selected(i) && selected_count == 0)
			{
				vec1_pull_seg(0) = r_pull_intersect(i,0);
				vec1_pull_seg(1) = r_pull_intersect(i,1);
				vec1_pull_seg(2) = r_pull_intersect(i,2);
				selected_count = 1;
				continue;
			}

			if (selected(i) && selected_count == 1)
			{
                if (fabs(vec1_pull_seg(0) - r_pull_intersect(i,0)) < epsilon && \
				    fabs(vec1_pull_seg(1) - r_pull_intersect(i,1)) < epsilon && \
					fabs(vec1_pull_seg(2) - r_pull_intersect(i,2)) < epsilon)
				{
					continue;
				}
				else
				{
					vec2_pull_seg(0) = r_pull_intersect(i,0);
				    vec2_pull_seg(1) = r_pull_intersect(i,1);
				    vec2_pull_seg(2) = r_pull_intersect(i,2);
				    selected_count = 2;
				    continue;
				}
			}

            if (selected(i) && selected_count == 2)
			{
			    if ((fabs(vec1_pull_seg(0) - r_pull_intersect(i,0)) < epsilon &&  \
				     fabs(vec1_pull_seg(1) - r_pull_intersect(i,1)) < epsilon &&  \
				     fabs(vec1_pull_seg(2) - r_pull_intersect(i,2)) < epsilon) || \
					(fabs(vec2_pull_seg(0) - r_pull_intersect(i,0)) < epsilon &&  \
				     fabs(vec2_pull_seg(1) - r_pull_intersect(i,1)) < epsilon &&  \
				     fabs(vec2_pull_seg(2) - r_pull_intersect(i,2)) < epsilon))
				{
					continue;
				}
				else
				{
					MY_ERROR("In LDAD bond_function xyz, More than 2 different points \
					          with different positions are being selected. \
							  This means a line intersecting with a cube has 3 points. \
							  Something wrong with the algorithm. Need to print and debug!");
				}  
			}
		}

        // No intersection .or. one point on the parallelepiped but another point is outside
		if (selected_count == 0 || selected_count == 1)
		{
			return 0.0;
		}
	}

    vec1_push_seg = ldadVectors * vec1_pull_seg.transpose();
    vec2_push_seg = ldadVectors * vec2_pull_seg.transpose();
    vec12_push_seg = vec2_push_seg - vec1_push_seg;
    total_length =  vec12_push_seg.norm();
	
    // use oneDFunction.integrate(vec1_pull_seg, vec2_pull_seg) -> helper function;
	return total_length * oneDFunction.integrate(vec1_pull_seg, vec2_pull_seg) * normalizer / distance;
}

int PointLineRelationship(const double& p)
{
	if (p < -1.0 - epsilon)
	{
		return 1;
	}
	else if (p >= -1.0 - epsilon && p <= -1.0 + epsilon)
	{
		return 2;
	}
	else if (p > -1.0 + epsilon && p < 1.0 - epsilon)
	{
		return 3;
	}
	else if (p >= 1.0 - epsilon && p <= 1.0 + epsilon)
	{
		return 4;
	}
	else if (p > 1.0 + epsilon)
	{
		return 5;
	}
	else
	{
		std::cout << "Point " << p << std::endl;
		MY_ERROR("In PointLineRelationship, the above point did not fall into any range.");
	}
}
