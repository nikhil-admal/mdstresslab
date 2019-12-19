/*
 * Composite.h
 *
 *  Created on: Dec 2, 2019
 *      Author: Nikhil
 */

#ifndef SRC_COMPOSITE_H_
#define SRC_COMPOSITE_H_

#include "Configuration.h"
#include "Grid.h"
#include <map>
#include <vector>
#include <set>
#include "typedef.h"

template<ConfigType T>
class Composite {
public:
	Composite(const Configuration& config,
			  const std::map<const Grid<T>*, double> gridNeighborhoodSize,
			  const double& noncontributingNeighborhoodSize);


	std::unique_ptr<Configuration> plocalConfiguration;
	VectorXi particleContributing;
	std::map<const Grid<T>*,std::vector<std::set<int>>> gridLocalNeighborListsMap;
	virtual ~Composite();
};

#include "Composite.cpp"

#endif /* SRC_COMPOSITE_H_ */
