/*
 * Material.h
 *
 *  Created on: Nov 4, 2019
 *      Author: Nikhil
 */

#ifndef ATOMISTICSYSTEM_H_
#define ATOMISTICSYSTEM_H_

#include "typedef.h"
#include "Configuration.h"

/**
 * @brief Represents a particle configuration including simulation box information.
 *
 * This class extends the base Configuration class by adding storage for simulation box
 * vectors (both Current and Reference) and periodic boundary conditions. It provides
 * functionality to read configurations from MDStressLab and LAMMPS dump-style files, manage
 * padding due to periodicity, and access box geometry details.
 *
 * The box vectors are stored as columns in 3×3 matrices (`box` and `reference_box`),
 * and periodic boundary conditions are represented as a 3D integer vector indicating
 * periodicity along each Cartesian axis.
 *
 * This class supports reading atomic configurations from two file formats:
 * - MDStressLab custom format with explicit box vectors and particle data.
 * - LAMMPS dump files with box bounds, masses, and atom coordinates.
 *
 * It also allows extracting padded configurations to accommodate atoms near periodic
 * boundaries, facilitating simulations requiring extended neighbor lists or replicas.
 */
class BoxConfiguration : public Configuration{
public:
    /*! \brief The current and reference box vectors stored as columns of respective
     * matrices.
     */
    Matrix3d box, reference_box;

    /*! \brief Periodic boundary conditions. \ref pbc=(1,0,1) implies periodicity along the \f$x\f$
     * and \f$z\f$-directions.
     */
	Vector3i pbc;

    /**
     * @brief Constructs a BoxConfiguration with a given number of particles.
     *
     * @param[in] numberOfParticles The number of particles in the configuration.
     * @param[in] referenceAndFinal Flag indicating which configurations to allocate:
     *            - 0: Only Current configuration
     *            - 1: Both Reference and Current configurations
     */
    BoxConfiguration(int numberOfParticles, int referenceAndFinal);


    /*! \brief This function returns a padded configuration by adding padding atoms originating
     * due to pbcs.
     * @param padding - thickness of the padding region
     * @returns - a pointer to an allocated padded configuration
     */
	Configuration* getConfiguration(double padding) const;

    /*! \brief A function to read the properties of atoms from a file in a MDStressLab
     * format.
     * @param configFileName - Name of the configuration file to read
     * @param referenceAndFinal Flag indicating which data to read:
     *        - 0: Only Current configuration
     *        - 1: Both Reference and Current configurations
     *
     * The MDStressLab format of the configuration file is
     * described below
     *
     * - **Line 1**: Number of particles (integer)
     * - **Lines 2–4**: Reference box vectors as columns of a 3×3 matrix
     * - **Lines 5–7**: Current box vectors as columns of a 3×3 matrix
     * - **Line 8**: Periodic boundary conditions (3 integers, typically 0 or 1)
     * - **Line 9 onward**: Per-particle data in 10 columns:
     *   - Column 1: Species code (e.g., `Ar`)
     *   - Columns 2–4: Current coordinates (x, y, z)
     *   - Columns 5–7: Velocities (vx, vy, vz)
     *   - Columns 8–10: Reference coordinates (X, Y, Z)
     *
     * ### 📄 Example:
     * \verbatim
        2
        52.9216036151419047 0.0 0.0
        0.0  52.9216036151419047 0.0
        0.0  0.0 52.9216036151419047
        52.9216036151419047 0.0 0.0
        0.0  53.4508196512933154 0.0
        0.0  0.0 52.9216036151419047
        1 1 1
        Ar  0.0000000000000000 0.0000000000000000 0.0000000000000000  0.0000000000000000 0.0000000000000000 0.0000000000000000  0.0000000000000000 0.0000000000000000 0.0000000000000000
        Ar  0.0000000000000000 2.6725409825646658 2.6460801807570951  0.0000000000000000 0.0000000000000000 0.0000000000000000  0.0000000000000000 2.6460801807570951 2.6460801807570951
     * \endverbatim
    */
    void read(std::string configFileName,int referenceAndFinal);

    /**
     * @brief Reads a configuration from an LAMMPS-style dump file.
     *
     * Opens the specified file and reads the atomic configuration data for the given configuration type
     * (e.g., Reference or Current). This function supports only files with the `.lmp` extension and
     * delegates detailed parsing to the private method `lmpParser`.
     *
     * On failure to open the file or if the file extension is not `.lmp`, this function prints an error
     * message and terminates the program.
     *
     * @param[in] filename The name of the LAMMPS dump file.
     * @param[in] configType Specifies which configuration type to read
     *                       (e.g., Reference or Current).
     *
     * The expected file format is a LAMMPS dump file with the following structure:
     *
     * - A header line with an optional comment (e.g., crystal orientation)
     * - Line indicating the number of atoms (e.g., `24480 atoms`)
     * - Line indicating the number of atom types (e.g., `3 atom types`)
     * - Box bounds:
     *   ```
     *   xlo xhi
     *   ylo yhi
     *   zlo zhi
     *   ```
     *   Optionally followed by `xy xz yz` tilt factors (for triclinic cells)
     * - A `Masses` section listing:
     *   ```
     *   <atom type> <mass> #<element symbol>
     *   ```
     * - An `Atoms` section:
     *   ```
     *   Atoms # atomic
     *   <atom ID> <atom type> <x> <y> <z> <ix> <iy> <iz>
     *   ```
     *   where `ix`, `iy`, `iz` are image flags used to unwrap atom positions.
     *
     * Example snippet:
     * \code
     * # Bcc oriented X=[-111] Y=[101] Z=[12-1].
     * 24480 atoms
     * 3 atom types
     * -95.18727654954 95.18727654954 xlo xhi
     * -45.71769755913 45.71769755913 ylo yhi
     * -11.877806246622 11.877806246622 zlo zhi
     * 0 0 0 xy xz yz
     *
     * Masses
     * 1 92.908 #Nb
     * 2 180.95 #Ta
     * 3 50.941 #V
     *
     * Atoms # atomic
     * 1 1 -95.15 -45.68 11.86 0 0 -1
     * 24003 3 -35.46 -2.29 5.25 0 0 0
     * 3 1 -94.28 -45.69 -6.59 0 0 0
     * ...
     * \endcode
     */
    void readLMP(const std::string&,const ConfigType& configType);


    /**
     * @brief Reads both current and reference configurations from LAMMPS-style dump files.
     *
     * @param[in] currentConfigFileName  The name of the file containing the current configuration.
     * @param[in] referenceConfigFileName  The name of the file containing the reference configuration.
     */
    void readLMP(const std::string& currentConfigFileName,
                 const std::string& referenceConfigFileName);

    /**
     * @brief Parses a LAMMPS dump file to populate configuration data.
     *
     * This function reads atomic configuration data, box dimensions, masses,
     * and species information from a LAMMPS-style dump file stream, filling
     * the relevant members of the BoxConfiguration object according to the
     * specified configuration type.
     *
     * The expected file format includes sections for:
     * - Number of atoms
     * - Number of atom types
     * - Box bounds (xlo xhi, ylo yhi, zlo zhi)
     * - Optional tilt factors (not currently handled explicitly)
     * - Masses section mapping atom types to species names
     * - Atoms section with atomic coordinates and optional image flags for periodicity
     *
     * Coordinates are adjusted for periodic boundary conditions using image flags.
     * The function ensures species consistency between reference and current configurations.
     *
     * @param[in,out] file       Input file stream positioned at the beginning of the LAMMPS dump file.
     * @param[in]     configType The type of configuration to read (e.g., Current or Reference).
     *
     * @throws MY_ERROR If the number of atoms in the file does not match the expected number,
     *                  or if coordinate parsing fails for any particle.
     */
    void lmpParser(std::ifstream&, const ConfigType&);

	virtual ~BoxConfiguration();
};

#endif /* ATOMISTICSYSTEM_H_ */
