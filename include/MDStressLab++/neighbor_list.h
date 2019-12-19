/*
 * neighbor_list.h
 *
 *  Created on: Nov 6, 2019
 *      Author: Nikhil
 */

#ifndef NEIGHBOR_LIST_H_
#define NEIGHBOR_LIST_H_

#include <vector>
#include "typedef.h"

typedef struct
{
  int numberOfParticles;
  double cutoff;
  int* Nneighbors;
  int* neighborList;
  int* beginIndex;
} NeighListOne;

// neighbor list structure
typedef struct
{
  int numberOfNeighborLists;
  NeighListOne* lists;
} NeighList;


void nbl_initialize(NeighList ** const nl);

int nbl_create_paddings(const int numberOfParticles,
                        const double cutoff,
                        double const * reference_cell,
                        double const * cell,
                        int const * PBC,
                        double const * reference_coordinates,
                        double const * coordinates,
                        const std::vector<std::string>& speciesCode,
                        int & numberOfPaddings,
                        std::vector<double> & reference_coordinatesOfPaddings,
                        std::vector<double> & coordinatesOfPaddings,
                        std::vector<std::string> & speciesCodeOfPaddings,
                        std::vector<int> & masterOfPaddings,
						int referenceAndFinal);

int nbl_build(NeighList * const nl,
              int const numberOfParticles,
              double const * coordinates,
              double const influenceDistance,
              int const numberOfCutoffs,
              double const * cutoffs,
              int const * needNeighbors);

int nbl_get_neigh(void const * const nl,
                  int const numberOfCutoffs,
                  double const * const cutoffs,
                  int const neighborListIndex,
                  int const particleNumber,
                  int * const numberOfNeighbors,
                  int const ** const neighborsOfParticle);

void nbl_clean(NeighList ** const nl);



#endif /* NEIGHBOR_LIST_H_ */
