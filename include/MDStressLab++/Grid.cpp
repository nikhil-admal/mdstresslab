/*
 * Grid.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#include "Grid.h"
#include <fstream>
#include "typedef.h"
#include "SpatialHash.h"
#include <iostream>

template<ConfigType T>
Grid<T>::Grid(int _ngrid) : ngrid(_ngrid)
{
	this->setCounter();
	int numberOfGrids= GridBase::numberOfCurrentGrids + GridBase::numberOfReferenceGrids;
	if (numberOfGrids == 1) MY_HEADING("Creating Grids");

	std::cout << "Grid " << numberOfGrids << ". Initializing a grid of size " << ngrid << " to the origin." << std::endl;
	this->coordinates.resize(ngrid,Vector3d(0.0,0.0,0.0));
}

template<ConfigType T>
Grid<T>::Grid(Vector3d lowerLimit,
	  	   Vector3d upperLimit,
		   int ngridx, int ngridy, int ngridz):ngrid(ngridx*ngridy*ngridz)
{
	if ( !(lowerLimit.array() < upperLimit.array()).prod() )
		MY_ERROR("ERROR: The coordinates of lowerLimit are not less than the upperLimit");

	this->setCounter();
	int numberOfGrids= GridBase::numberOfCurrentGrids + GridBase::numberOfReferenceGrids;
	if (numberOfGrids == 1) MY_HEADING("Creating Grids");


	std::cout << "Grid " << numberOfGrids << ". Creating a uniform grid of " << ngrid << " points between ("  << lowerLimit
			  << ") and (" << upperLimit << ")"<< std::endl;
	std::cout << std::endl;
	coordinates.resize(ngrid);

    int index= 0;
    Vector3d coordinate;
    coordinate= lowerLimit;
    for (auto i : range<int>(0,ngridx))
	{
        coordinate(0) = lowerLimit(0) + i*(upperLimit(0)-lowerLimit(0))/ngridx;
        for (auto j : range<int>(0,ngridy))
        {
            coordinate(1) = lowerLimit(1) + j*(upperLimit(1)-lowerLimit(1))/ngridy;
            for (auto k : range<int>(0,ngridz))
            {
                coordinate(2)= lowerLimit(2) + k*(upperLimit(2)-lowerLimit(2))/ngridz;
                coordinates[index]= coordinate;
                index++;
            }
        }
    }
}

template<ConfigType T>
Grid<T>::Grid(std::string filename)
{
	this->setCounter();
	int numberOfGrids= GridBase::numberOfCurrentGrids + GridBase::numberOfReferenceGrids;
	if (numberOfGrids == 1) MY_HEADING("Creating Grids");
	std::cout << "Grid " << numberOfGrids << ". Reading grid from filename: " << filename << std::endl;
	std::cout << std::endl;

	std::ifstream file(filename);

	file >> ngrid;
	coordinates.resize(ngrid);
	for(auto& position : coordinates)
		for (int i_dim=0; i_dim<DIM; i_dim++)
			if(!(file >> position(i_dim)))  MY_ERROR("ERROR: Reading grid coordinates");
}

// Function to calculate neighbor lists of grid points consisting of particles in subconfig within a distance of
// padding from the grid points
template<ConfigType T>
std::vector<std::set<int>> Grid<T>::getGridNeighborLists(const SubConfiguration& subconfig, const double& padding) const
{
	std::vector<std::set<int>> gridNeighborLists;

	// Hash the coordinates
	Vector3d origin,step;
	origin.setConstant(0.0);
	step.setConstant(padding);
	ConstSpatialHash hashParticles(origin,step,subconfig.coordinates.at(T));
	// Hash the grid points
	ConstSpatialHash hashGrid(origin,step,coordinates);


	int i_gridPoint= 0;
	for(const auto& gridPoint : coordinates)
	{
		std::set<int> gridContributingList;
		Triplet bin= hashGrid.hashFunction(i_gridPoint);
		Triplet neighborBin;
		// for each neighboring bin of a grid point's bin
		for (const auto& neighborBin : bin.neighborList())
		{
			std::vector<int>& particleList= hashParticles.hashTable[neighborBin];
			// for each particle in a neighboring bin
			for(const auto& particle : particleList)
			{
				double distance= (subconfig.coordinates.at(T).row(particle)-gridPoint).squaredNorm();
				if (distance<=pow(padding,2))
					gridContributingList.insert(particle);
			}
		}
		gridNeighborLists.push_back(gridContributingList);
		i_gridPoint++;

	}
	return gridNeighborLists;
}

template<ConfigType T>
void Grid<T>::write(std::string filename) const
{
	std::ofstream file(filename+".grid");

	file << ngrid << "\n";
	file << "\n";
	for (const auto& coordinate : coordinates)
		file << coordinate << std::endl;
}

template<ConfigType T>
void Grid<T>::read(std::string filename)
{
	std::cout << "Reading grid from file " << filename << "\n" << std::endl;

	std::ifstream file(filename);

	int ngridInFile;
	file >> ngridInFile;
	if (ngrid!= ngridInFile)
		MY_ERROR("Error: Number of grid points in file" << " = " << ngridInFile <<
				", does not equal to " << ngrid << " in grid.");

	for(auto& position : coordinates)
		for (int i_dim=0; i_dim<DIM; i_dim++)
			if(!(file >> position(i_dim)))  MY_ERROR("ERROR: Reading grid coordinates");
}

template<ConfigType T>
void Grid<T>::setCounter()
{
	if (T == Reference)  GridBase::numberOfReferenceGrids+= 1;
	else if (T == Current) GridBase::numberOfCurrentGrids+=1;
	else MY_ERROR("Unrecognized grid type");
}


template<ConfigType T>
Grid<T>::~Grid() {
	// TODO Auto-generated destructor stub
}

template<ConfigType T>
GridSubConfiguration<T>::GridSubConfiguration(const Grid<T>& grid, const SubConfiguration& subconfig, const double& padding) :
grid(grid),
subconfig(subconfig),
padding(padding),
hashGridSubconfig(std::make_pair(ConstSpatialHash(Vector3d::Zero(),Vector3d::Constant(padding),grid.coordinates),
                                 ConstSpatialHash(Vector3d::Zero(),Vector3d::Constant(padding),subconfig.coordinates.at(T))))
{}

template<ConfigType T>
std::set<int> GridSubConfiguration<T>::getGridPointNeighbors(const int& i_gridPoint)
{
    ConstSpatialHash& hashGrid= hashGridSubconfig.first;
    ConstSpatialHash& hashParticles= hashGridSubconfig.second;

    std::set<int> gridContributingList;
    Triplet bin= hashGrid.hashFunction(i_gridPoint);
    Triplet neighborBin;
    // for each neighboring bin of a grid point's bin
    for (const auto& neighborBin : bin.neighborList())
    {
        const std::vector<int>& particleList= hashParticles.hashTable[neighborBin];
        // for each particle in a neighboring bin
        for(const auto& particle : particleList)
        {
            double distance= (subconfig.coordinates.at(T).row(particle)-grid.coordinates[i_gridPoint]).squaredNorm();
            if (distance<=pow(padding,2))
                gridContributingList.insert(particle);
        }
    }
    return gridContributingList;
}
