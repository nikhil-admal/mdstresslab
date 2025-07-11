/*
 * Grid.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef GRID_H_
#define GRID_H_
#include "typedef.h"
#include <vector>
#include "SubConfiguration.h"
#include <set>

class GridBase
{
public:
	std::vector<Vector3d> coordinates;
	static int numberOfReferenceGrids;
	static int numberOfCurrentGrids;
};

int GridBase::numberOfReferenceGrids= 0;
int GridBase::numberOfCurrentGrids= 0;

/*!
 * Describes a grid on which stress is computed.
 * @tparam T - StressType (Cauchy/Piola)
 */
template<ConfigType T>
class Grid : public GridBase{
public:
	Grid(int);
	Grid(Vector3d,
		 Vector3d,
		 int ngridx, int ngridy=1,int ngridz=1);
    Grid(std::string);
	virtual ~Grid();
	int ngrid;

	void read(std::string);
	void write(std::string) const;
	void setCounter();
	std::vector<std::set<int>> getGridNeighborLists(const SubConfiguration&, const double&) const;
};

template<ConfigType T>
class GridSubConfiguration
{
private:
    const Grid<T>& grid;
    const SubConfiguration& subconfig;
    const double padding;
    std::pair<ConstSpatialHash,ConstSpatialHash> hashGridSubconfig;
public:
    GridSubConfiguration(const Grid<T>&, const SubConfiguration&, const double& );
    std::set<int> getGridPointNeighbors(const int& ) const;

};
/*!
 * \example{lineno} testGrid.cpp
 * This example demonstrates the Grid class.
 *
 * -# Construct two grids in a rectangular region defined by the lower and upper limits that describe
 * the region's body diagonal endpoints. 'grid1' is a 3D Reference grid, while 'grid2' is a 1D Current grid
 * @snippet{lineno} testGrid.cpp Limits
 *
 * -# Write the coordinates of grids to files grid1.grid and grid2.grid
 * @snippet{lineno} testGrid.cpp Write
 *
 * -# An alternate way to construct a grid is to read its coordinates from a file.
 * The following snippet initializes a 3D Current grid and reads its coordinates from a file
 * @snippet{lineno} testGrid.cpp Read
 *
 * -# We can keep track of the total number of Reference and Current grids using static variables.
 * @snippet{lineno} testGrid.cpp Static
 *
 * Full code:
 */
#include "Grid.cpp"

#endif /* GRID_H_ */
