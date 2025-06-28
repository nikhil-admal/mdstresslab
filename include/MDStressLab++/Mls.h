/*
 * Mls.h
 *
 *  Created on: Feb 11, 2020
 *      Author: Min
 */

#ifndef MLS_H_
#define MLS_H_

#include <vector>
#include "typedef.h"
#include "BoxConfiguration.h"
#include "Grid.h"

/**
 * @class Mls
 * @brief Computes deformation gradients and stress tensors using Moving Least Squares (MLS) interpolation on a spatial grid.
 *
 * This class implements a 3D MLS approach to compute the local deformation gradient
 * and map Piola-Kirchhoff stress tensors to Cauchy stresses on a grid. It accounts
 * for periodic boundary conditions and constructs local neighbor lists for accuracy.
 *
 * Typical use: compute deformations and push stresses from atomistic data onto grid points.
 */
class Mls {
    public:

        /**
         * @brief Identifier name used for output file prefixes.
         */
        std::string name;

        /**
         * @brief Radius of influence for the MLS weighting function.
         */
        double radiusMls;

        /**
         * @brief Deformation gradient at each grid point (one per grid coordinate).
         */
        std::vector<Matrix3d> deformationGradient;

        /**
         * @brief Positions of grid points after applying MLS displacements.
         */
        std::vector<Vector3d> gridPushed;

        /**
         * @brief Constructs the deformation gradient field from a body configuration and grid.
         *
         * Performs MLS interpolation to compute deformation gradients at each grid point
         * using displacement fields derived from a reference and current configuration.
         *
         * Handles periodic boundary conditions, builds neighbor lists, checks sanity of interpolation points,
         * and stores the final deformation gradient and displaced grid points.
         *
         * @param body The configuration containing reference and current atom positions.
         * @param pgrid Pointer to the spatial grid (reference coordinates).
         * @param radiusMls The cutoff radius for MLS interpolation.
         * @param name An identifier used in output filenames.
         */
        Mls(const BoxConfiguration& body, const Grid<Reference>* pgrid, double radiusMls, const std::string name);

        ~Mls();

        /**
         * @brief Converts Piola-Kirchhoff stress to Cauchy stress at grid points.
         *
         * Uses the computed deformation gradients and their determinants to push
         * stress tensors from the reference configuration to the current spatial grid.
         *
         * @param piolaStress Input Piola-Kirchhoff stresses at reference points.
         * @param cauchyStress Output Cauchy stresses on the deformed grid.
         */
        void pushToCauchy(const std::vector<Matrix3d>& piolaStress,std::vector<Matrix3d>& cauchyStress);

        /**
         * @brief Writes the deformation gradient at each grid point to a file.
         *
         * File name is constructed as `<name>.DeformationGradient`. Includes both
         * grid positions and corresponding 3x3 deformation gradient matrices.
         */
        void writeDeformationGradient();

        /**
         * @brief Writes the pushed (deformed) grid point coordinates to a file.
         *
         * File name is `<name>.pushedGrids`. One line per grid point with x/y/z coordinates.
         */
        void writeGridPushed();

        /**
         * @brief Writes Cauchy stresses at grid points to a file.
         *
         * File name is `<name>.stress`. Outputs 6 unique components per symmetric 3x3 stress tensor.
         *
         * @param cauchyStress Vector of Cauchy stress tensors (must match grid size).
         */
        void writePushedCauchyStress(std::vector<Matrix3d>& cauchyStress);
};

/**
 * \example{lineno} testMls.cpp
 * Demonstrates the use of Moving Least Squares (Mls) to compute
 * deformation gradients and push stress tensors to a spatial grid.
 *
 * This test sets up a simulation involving a particle configuration and grid data,
 * computes the Piola-Kirchhoff stress using the LDAD method, and then applies the 
 * Moving Least Squares (MLS) method to compute the deformation gradient and 
 * push the Piola stress to Cauchy stress on a grid.
 *
 * Key steps in this example:
 * - Load atomic configuration and potential.
 * - Load reference grid data.
 * - Compute LDAD Piola stress using `MethodLdadConstant`.
 * - Construct an `Mls` object and compute deformation gradients.
 * - Push stress field to obtain Cauchy stress tensors.
 * - Write deformation gradient, pushed grid, and Cauchy stress to output files.
 *
 * Output files:
 * - `mls.DeformationGradient`
 * - `mlsPushed.grid`
 * - `mlsPushed.stress`
 *
 * Requirements:
 * - Input file `config.data` with atomic data in BoxConfiguration format.
 * - Grid files `reference.grid`.
 * - A valid KIM model (here: "LJ_Smoothed_Bernardes_1958_Ar__MO_764178710049_001").
 * - Expected Cauchy stress field `mlsPushedReference.stress` for unit testing.
 *
 * -# Open the configuration file (*.data in MDStressLab format) and read the number of particles
 * @snippet{lineno} testMls.cpp Read
 *
 * -# Initialize the BoxConfiguration and read the reference and current atomic coordinates
 *  from the configuration file. The reference configuration is an fcc Ar crystal in the relaxed state.
 *  The deformed/current configuration is the reference crystal strained in the \f$y\f$-direction.
 * @snippet{lineno} testMls.cpp Configuration
 *
 * -# Initialize the Kim model
 * @snippet{lineno} testMls.cpp Model
 *
 * -# Load reference grid coordinates from files.
 * @snippet{lineno} testMls.cpp Grid
 *
 * -# Initialize LDAD vectors and construct LDAD method
 *    MethodLdadConstant. The LDAD vectors (columns) are the three
 *    lattice vectors of the conventional unit cell.
 * @snippet{lineno} testMls.cpp LDAD
 *
 * -# Construct the Piola stress using the above method
 * @snippet{lineno} testMls.cpp Stress
 *
 * -# Construct an Mls object to compute deformation gradients and push the
 * stress to Cauchy form.
 * @snippet{lineno} testMls.cpp MLS
 *
 * -# Write deformation gradients, pushed grid points, and Cauchy stresses to output files.
 * @snippet{lineno} testMls.cpp Output
 *
 * -# Compare the pushed fields to expected data.
 * @snippet{lineno} testMls.cpp Compare
 *
 * Full code:
 */

#endif /* MLS_H_ */


