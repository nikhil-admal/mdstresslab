/*
 * MethodSphere.h
 *
 *  Created on: Nov 5, 2019
 *      Author: Nikhil
 */

#ifndef METHODSPHERE_H_
#define METHODSPHERE_H_

#include "Method.h"
#include "typedef.h"
#include "Polynomial.h"
#include <map>

/*!
 * @class MethodSphere
 * @brief Implements radially symmetric weighting functions (Hardy, Virial) and its
 * associated bond function for atomistic stress computation.
 *
 * This class defines a weighting function \f$w\f$ over a spherical domain of radius \f$R\f$
 * (the averaging domain size). Two common forms are supported:
 *  - Hardy-type weighting
 *  - Virial-type (uniform) weighting
 *
 * Additionally, users can supply a custom radial weighting function as a piecewise linear map.
 */
class MethodSphere : public Method<MethodSphere>
{          
    friend class Method<MethodSphere>;

public:
    /*!
     * @brief Constructs a MethodSphere using a named weighting type.
     *
     * @param averagingDomainSize The diameter \f$R\f$ of the spherical averaging domain.
     * @param type The type of weighting function: `"hardy"` or `"virial"`.
     *
     * For `"hardy"`, the radial weighting function \f$w(r)\f$ is defined as:
     * \f[
     * w(r) =
     * \begin{cases}
     * c, & \text{if } 0 \le r \le R/2 \\
     * 2c\left(1 - \frac{r}{R}\right), & \text{if } R/2 < r \le R \\
     * 0, & \text{otherwise}
     * \end{cases}
     * \f]
     * where \f$c = \frac{8}{5\pi R^3}\f$ ensures normalization.
     *
     * For `"virial"`, the weighting function is:
     * \f[
     * w(r) =
     * \begin{cases}
     * c, & \text{if } 0 \le r \le R \\
     * 0, & \text{otherwise}
     * \end{cases}
     * \f]
     * where \f$c = \frac{3}{4 \pi R^3}\f$.
     */
	explicit MethodSphere(double, const std::string& type="virial");

    /*!
    * @brief Constructs a MethodSphere using a user-defined piecewise-linear weight function.
    *
    * @param averagingDomainSize The diameter \f$R\f$ of the spherical domain.
    * @param map A map from normalized radius values in \f$[0,1]\f$ to non-negative weight values.
    *
    * The keys of the map (in \f$[0, 1]\f$) are rescaled by \f$R\f$ to get the real radii.
    * For each interval \f$[x_i, x_{i+1}]\f$, a linear segment is fitted:
    * \f[
    * w(r) = y_i + \frac{y_{i+1} - y_i}{x_{i+1} - x_i} (r - x_i R), \quad r \in [x_i R, x_{i+1} R]
    * \f]
    * The resulting piecewise linear function is then normalized, i.e.
    * \f[
    * 4 \pi \int_{0}^{R} w(r) \cdot r^2 \, dr = 1
    * \f]
    *
    * This allows users to define arbitrary, smooth weight profiles that are zero outside the domain.
    */
	MethodSphere(double, std::map<double,double>);
	MethodSphere(const MethodSphere&);
	virtual ~MethodSphere();

protected:
   /**
    * @brief Evaluate the weighting function at a vector.
    * @param vec Position vector relative to center of the averaging domain.
    * @return Weight value.
    */
	double operator()(const Vector3d& vec) const;

    /*!
     * @brief Computes the bond function for a pair of atomic positions.
     *
     * @param vec1 Vector from grid point to atom 1.
     * @param vec2 Vector from grid point to atom 2.
     * @return This function evaluates the line integral of the radially symmetric
     * weighting function along the bond.
     */
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2) const;


private:
	std::map<std::pair<double,double>,Polynomial> piecewisePolynomial;
	double normalizer;
	double integratePolynomial(const int&, const double&,const double&,const double&,const double&) const;
	double integrate(const double&,const double&,const double&,const double&) const;
};

/*!
 * @example{lineno} testLJ.cpp
 * An example demonstrating the computation of Piola and Cauchy stress tensors
 * using hardy and virial spherical averaging domains, and 1D and 2D grids
 *
 * -# Open the configuration file (*.data in MDStressLab format) and read the number of particles
 * @snippet{lineno} testLJ.cpp Read
 *
 * -# Initialize the BoxConfiguration and read the reference and current atomic coordinates
 *  from the configuration file. The reference configuration is an fcc Ar crystal in the relaxed state.
 *  The deformed/current configuration is the reference crystal strained in the \f$y\f$-direction.
 * @snippet{lineno} testLJ.cpp Configuration
 *
 * -# Initialize the Kim model
 * @snippet{lineno} testLJ.cpp Model
 *
 * -# Create 1D and 2D grids
 * @snippet{lineno} testLJ.cpp Grid
 *
 * -# Construct hardy and virial MethodSphere objects
 * @snippet{lineno} testLJ.cpp Method
 *
 * -# Construct Stress objects using the MethodSphere and Grid objects
 * @snippet{lineno} testLJ.cpp Stress
 *
 * -# Calculate stresses. The following snippet shows that stresses can be calculated all at once or
 * with various combinations.
 * @snippet{lineno} testLJ.cpp Calculate
 *
 * -# Write stresses
 * @snippet{lineno} testLJ.cpp Write
 *
 * -# We compare our results with the exact results for unit testing purposes.
 * @snippet{lineno} testLJ.cpp Compare
 *
 * #### Visual comparison
 *
 *  * See @ref visualization_intro "Visualization Utilities" for contour plotting
 *  <table>
 *   <tr>
 *     <td align="center">
 *       <!-- Thumbnail links to high‑res PDF -->
 *       <a href="testLJ_hardyCauchy_xx.pdf">
 *         <img src="testLJ_hardyCauchy_xx.png" width="300px"/><br/>
 *         Hardy Cauchy stress
 *       </a><br/>
 *       <a href="testLJ_hardyCauchy_xx.pdf">[View high‑res PDF]</a>
 *     </td>
 *     <td align="center">
 *       <a href="testLJ_virialCauchy_xx.pdf">
 *         <img src="testLJ_virialCauchy_xx.png" width="300px"/><br/>
 *         Virial Cauchy stress
 *       </a><br/>
 *       <a href="testLJ_virialCauchy_xx.pdf">[View high‑res PDF]</a>
 *     </td>
 *   </tr>
 * </table>
 *
 * Full code:
 */

/*!
 * @example testLJScript.in
 */



#endif /* METHODSPHERE_H_ */
