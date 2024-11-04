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
    std::set<int> getGridPointNeighbors(const int& );

};
/*!
 * \example testGrid.cpp
 * This is an example of how to use the Grid class
 */
#include "Grid.cpp"

#endif /* GRID_H_ */
