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

template<ConfigType T>
class Grid : public GridBase{
public:
	Grid(int);
	Grid(Vector3d lowerLimit,
		 Vector3d upperLimit,
		 int _ngrid);
    Grid(std::string);
	virtual ~Grid();
	int ngrid;

	void read(std::string);
	void write(std::string) const;
	void setCounter();
	std::vector<std::set<int>> getGridNeighborLists(const SubConfiguration&, const double&) const;
};
#include "Grid.cpp"

#endif /* GRID_H_ */
