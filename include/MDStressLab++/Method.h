/*
 * Method.h
 *
 *  Created on: Aug 26, 2022
 *      Author: Nikhil
 */

#ifndef METHOD_H_
#define METHOD_H_

#include "typedef.h"

/**
 * This class provides a polymorphic interface for defining a `method`, which constitutes a
 * weighting function and its associated bond function, to compute atomistic stress.
 *
 * It is intended to be used as a CRTP (Curiously Recurring Template Pattern) base class:
 * the actual implementation of the method should be provided in the derived class
 * `TMethod`, which will override the weighting function call operator and `bondFunction`.
 *
 * @tparam TMethod The derived class implementing specific method logic. Examples of derived
 * classes include `MethodSphere`, `MethodLdad<Constant>`, and `MethodLdad<Trigonometric>`.
 */
template<typename TMethod>
class Method
{
public:
	virtual ~Method();

    /**
     * @brief The weighting function \f$w\f$ used to compute stress. @see `bondFunction()` and
     * `calculateStress` to see its usage in stress calculation.
     *
     * This weighting function operator delegates the implementation to the derived `TMethod` class.
     *
     * @param vec A 3D vector
     * @return The scalar (non-negative) output of the weighting function.
     */
	double operator()(const Vector3d& vec) const;

    /**
     * @brief The bond function \f$b(\boldsymbol v^1, \boldsymbol v^2)\f$ gives the
     * weight associated to a bond formed by two atoms. @see calculateStress() to see how
     * bondfunction is used in stress computation
     *
     *
     * The function delegates the actual implementation
     * to `TMethod::bondFunction(vec1, vec2)`.
     *
     * The bond function is defined as
     * \f[
     * \begin{equation}
     *      b(\boldsymbol v^1, \boldsymbol v^2)=
     *      \int_{s=0}^1 w((1-s)\boldsymbol v^1 + s\boldsymbol v^2) \, ds
     * \end{equation}
     * \f]
     *
     * The actual implementation of the bondFunction is delegated to
     * `TMethod::bondFunction(vec1, vec2)`.
     *
     * @param vec1 first atomic position vector relative to a grid point
     * @param vec2 second atomic position vector relative to a grid point
     * @return The scalar output of the bond function.
     */
	double bondFunction(const Vector3d& vec1, const Vector3d& vec2) const;

    /**
     * @brief Gets the spatial size of the averaging domain, i.e. the support
     * of the weighting function
     *
     * @return Averaging domain radius or cutoff, as a double.
     */
    double getAveragingDomainSize() const;

protected:
	double averagingDomainSize;
private:
   /**
    * @brief Default constructor (private).
    *
    * Initializes domain size to zero.
    */
    Method();

    /**
     * @brief Constructor with specified averaging domain size.
     *
     * @param averagingDomainSize Domain radius to be stored.
     */
    Method(double);

    /**
     * @brief Copy constructor (private).
     *
     * Performs shallow copy of base data.
     *
     * @param method Another Method instance to copy from.
     */
    Method(const Method<TMethod>& method);

    /**
     * @brief Grants the derived `TMethod` class access to protected/private members.
     */
    friend TMethod;
};

#include "Method.cpp"
#endif /* METHOD_H_ */
