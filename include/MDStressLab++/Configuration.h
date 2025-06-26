/*
 * Configuration.h
 *
 *  Created on: Nov 25, 2019
 *      Author: Nikhil
 */

#ifndef INCLUDE_MDSTRESSLAB___CONFIGURATION_H_
#define INCLUDE_MDSTRESSLAB___CONFIGURATION_H_

#include <vector>
#include <string>
#include "typedef.h"
#include <iostream>
#include <set>
#include <map>

/**
* @brief Represents atomic configuration data including coordinates, velocities, and species.
*
* This class stores the number of particles, their species, coordinates for different
* configuration states (e.g., Reference and Current), and velocities.
*
*/
class Configuration{
public:
    /**
    * @brief Constructs a Configuration object.
    *
    * @param numberOfParticles     Number of particles in the configuration.
    * @param referenceAndFinal     Flag indicating which configurations to allocate:
            *                      - 0: Only Current configuration
    *                              - 1: Both Current and Reference configurations
    */
	Configuration(int,int);
	virtual ~Configuration();

    /**
     * @brief Total number of particles in the configuration.
     */
	int numberOfParticles;

    /**
     * @brief Species names for each particle (size equals numberOfParticles).
     */
	std::vector<std::string> species;

    /**
     * @brief Map from configuration type (Reference or Current) to coordinate matrices.
     *
     * Each matrix in the map has dimensions
     * \f$\texttt{numberOfParticles} \times 3\f$,
     * where each row corresponds to a particle’s position vector (x, y, z).
     *
     * - `coordinates[Current]` is always allocated.
     * - `coordinates[Reference]` is allocated only if `referenceAndFinal == 1` during construction;
     *    otherwise, it is a zero-sized matrix.
     */
	std::map<ConfigType,MatrixXd> coordinates;

    /**
     * @brief Velocities of particles.
     *
     * A \f$\texttt{numberOfParticles} \times 3\f$ matrix where each row corresponds to
     * the velocity vector of a particle.
     */
	MatrixXd velocities;

    /**
     * @brief can be deleted
     */
    Configuration* getLocalConfiguration(const std::set<int>& localParticleList) const;
};


#endif /* INCLUDE_MDSTRESSLAB___CONFIGURATION_H_ */
