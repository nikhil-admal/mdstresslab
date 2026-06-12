/*
 * SpatialHash.h
 *
 *  Created on: Nov 27, 2019
 *      Author: Nikhil
 */

#ifndef SRC_SPATIALHASH_H_
#define SRC_SPATIALHASH_H_

#include <vector>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include "range.h"
#include "typedef.h"

/*!
 * A Triplet object is a triplet of integers provided with a \f$<\f$ relation
 * and  a member function to enumerates its neighbors. The \f$<\f$ relation
 * allows us to use Triplets as keys to an std::map, which is used to spatially hash
 * a configuration of particles (see the SpatialHash class)
 */
class Triplet: public Vector3i
{
public:
	Triplet()=default;
    Triplet(int i,int j,int k): Vector3i(i,j,k) {}
    Triplet(const Vector3i& base) : Vector3i(base){}

	// Define < operator to use std::map<Triplet,std::vector<int>>
    bool operator <(const Triplet& rhs)  const
    {
        if (this->operator()(0) < rhs(0))   return true;
        if (rhs(0)   < this->operator()(0)) return false;
        if (this->operator()(1) < rhs(1))   return true;
        if (rhs(1)   < this->operator()(1)) return false;
        if (this->operator()(2) < rhs(2))   return true;
        return false;
    }

    std::vector<Triplet> neighborList()
	{
    	std::vector<Triplet> tripletList;
    	for (auto i : range<int>(-1,2))
			for (auto j : range<int>(-1,2))
				for (auto k : range<int>(-1,2))
				{
					Triplet s= static_cast<Vector3i>(*this + Triplet(i,j,k));
					tripletList.push_back(s);
				}
    	return tripletList;
	}

    virtual ~Triplet(){}
};


/*!
 * A class to conditionally select a typename. If flag==true, select typename T. Otherwise, select
 * typename U.
 * @tparam flag boolean
 * @tparam T typename
 * @tparam U typename
 */
template<bool flag, typename T, typename U>
struct Select
{
    /*!
     * @typedef Result the chosen typename
     */
    typedef T Result;
};

template<typename T, typename U>
struct Select<false, T, U> { typedef U Result; };

/*!
 * A class to spatially hash a collection of points to a given 3D grid.
 * @tparam isConst A boolean template parameter to specify whether the coordinates will be
 * altered or not.
 */


template<bool isConst>
class SpatialHash {
public:
	typedef typename Select<isConst, const MatrixXd, MatrixXd>::Result A;
	typedef typename Select<isConst, const Vector3d, Vector3d>::Result B;
	typedef typename Select<isConst, const std::vector<Vector3d>, std::vector<Vector3d>>::Result C;

    /*!
     * Origin of the grid used for hashing
     */
	Vector3d origin;
    /*!
     * Size of the grid used for hashing
     */
    Vector3d step;
    /*!
     * Coordinates of the points
     */
	std::vector<Eigen::Map<B>> coordinates;
    /*!
     *  hashTable: A map that maps a bin (Triplet) to an ordered list of particle numbers
     */
	std::map<Triplet,std::vector<int>> hashTable;

	SpatialHash();
	SpatialHash(Vector3d origin,
				Vector3d step,
				A& coordinates): origin(origin),step(step)
	{
		int size= coordinates.rows();
		for(int i_point=0; i_point<size; i_point++)
		{
			Eigen::Map<B> vector(&coordinates(i_point,0),3);
			this->coordinates.push_back(vector);
			hashTable[hashFunction(i_point)].push_back(i_point);
		}
	}
	SpatialHash(Vector3d origin,
				Vector3d step,
				C& coordinates) : origin(origin),step(step)
	{
		int size= coordinates.size();
		for(int i_point=0; i_point<size; i_point++)
		{
			Eigen::Map<B> vector(&coordinates[i_point](0),3);
			this->coordinates.push_back(vector);
			hashTable[hashFunction(i_point)].push_back(i_point);
		}

	}

    /*!
     * hashFunction: A map that maps a particle \f$i\f$ to a bin (Triplet)
     */
	Triplet hashFunction(const int& i) const
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

	virtual ~SpatialHash(){}

};
typedef SpatialHash<true> ConstSpatialHash;

/*!
 * @class BoxPoints
 * @details A class to fold a collection of points back into a given periodic box
 */
class BoxPoints : public SpatialHash<false>
{
private:
    Matrix3d cell;
public:
	BoxPoints(Vector3d origin,
			  Vector3d step,
			  MatrixXd& coordinates):SpatialHash<false>(origin,step,coordinates)
	{
		cell= step.asDiagonal();
	}
	BoxPoints(Vector3d origin,
			  const Matrix3d& cell,
			  MatrixXd& coordinates):SpatialHash<false>(origin,cell.diagonal(),coordinates),cell(cell){}
	BoxPoints(Vector3d origin,
			  Vector3d step,
			  std::vector<Vector3d>& coordinates) : SpatialHash<false>(origin,step,coordinates)
	{
		cell= step.asDiagonal();
	}
	BoxPoints(Vector3d origin,
			  const Matrix3d& cell,
			  std::vector<Vector3d>& coordinates) : SpatialHash<false>(origin,cell.diagonal(),coordinates),cell(cell) {}

	virtual ~BoxPoints(){}

    /*!
     * Folds the points into the box using spatial hashing
     * @param[in] pbc A binary triplet to describe periodic boundary conditions
     */
	void fold(const Vector3i& pbc)
	{
		Matrix3d inverseCell= cell.inverse();
		int size= coordinates.size();
		for(auto i_point : range<int>(0,size))
		{
			Vector3d relativePosition= coordinates[i_point] - origin;
			Vector3d fractionalPosition= (inverseCell * relativePosition.transpose()).transpose();
			Vector3d foldedFractionalPosition= fractionalPosition;
			for (int i_dim=0; i_dim<DIM; ++i_dim)
				if (pbc(i_dim))
					foldedFractionalPosition(i_dim)-= std::floor(foldedFractionalPosition(i_dim));

			if ((foldedFractionalPosition-fractionalPosition).squaredNorm()>epsilon)
			{
				std::cout << "Folding point " << i_point << ": "
						  << coordinates[i_point] << " ---> ";
				coordinates[i_point]= origin + (cell * foldedFractionalPosition.transpose()).transpose();
				std::cout << coordinates[i_point] << std::endl;
			}
		}
	}
};

#endif /* SRC_SPATIALHASH_H_ */
