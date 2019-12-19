/*
 * SpatialHash.cpp
 *
 *  Created on: Nov 27, 2019
 *      Author: Nikhil
 */

// 3D spatial hashing

#include <iostream>
#include "SpatialHash.h"
#include "typedef.h"
#include "assert.h"
#include "range.h"

SpatialHash::SpatialHash(Vector3d origin,
						 Vector3d step,
						 const MatrixXd& coordinates):
	origin(origin),step(step)
{
	int size= coordinates.rows();
	for(int i_point=0; i_point<size; i_point++)
	{
        Eigen::Map<const Vector3d> vector(&coordinates(i_point,0),3);
		this->coordinates.push_back(vector);
		hashTable[hashFunction(i_point)].push_back(i_point);
	}
}

SpatialHash::SpatialHash(Vector3d origin,
						 Vector3d step,
						 const std::vector<Vector3d>& coordinates):
	origin(origin),step(step)
{
	int size= coordinates.size();
	for(int i_point=0; i_point<size; i_point++)
	{
        Eigen::Map<const Vector3d> vector(&coordinates[i_point](0),3);
		this->coordinates.push_back(vector);
		hashTable[hashFunction(i_point)].push_back(i_point);
	}

}

Triplet SpatialHash::hashFunction(const int i)
{
	assert(!(i<0));
	Triplet triple;
	Vector3i& base= triple;
	Vector3d coordinate= coordinates[i];
	coordinate= coordinate-origin;
	coordinate= (coordinate.array()/step.array()).floor().matrix();
	base= coordinate.template cast<int>();
	return triple;
}


SpatialHash::~SpatialHash() {
	// TODO Auto-generated destructor stub
}

