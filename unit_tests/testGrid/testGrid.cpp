/*
 * grid_test.cpp
 *
 *  Created on: Nov 29, 2019
 *      Author: Nikhil
 */
#include "Grid.h"
#include <iostream>
#include "typedef.h"
#include <fstream>

int main()
{
    /*! [Limits] */
	Vector3d lowerLimit(1,2,3);
	Vector3d upperLimit(2,4,6);
	int ngrid= 10;
	Grid<Reference> grid1(lowerLimit,upperLimit,ngrid,ngrid,ngrid);
    Grid<Current>   grid2(lowerLimit,upperLimit,ngrid);
    /*! [Limits] */

    /*! [Write] */
    grid1.write("grid1");
    grid2.write("grid2");
    /*! [Write] */

	/*! [Read] */
	Grid<Current> grid3(ngrid*ngrid*ngrid);
	grid3.read("grid1.grid");
    /*! [Read] */

    /*![Static]*/
    std::cout << "Number of Current grids = " << Grid<Current>::numberOfCurrentGrids << std::endl;
    std::cout << "Number of Reference grids = " << Grid<Reference>::numberOfReferenceGrids << std::endl;
    /*![Static]*/
}





