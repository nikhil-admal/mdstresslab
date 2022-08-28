/*
 * Method.cpp
 *
 *  Created on: Aug 26, 2022
 *      Author: Nikhil
 */

#include "typedef.h"
#include <Method.h>

template<typename TMethod>
Method<TMethod>::Method() : averagingDomainSize(0)
{ }

template<typename TMethod>
Method<TMethod>::Method(double averagingDomainSize) : averagingDomainSize(averagingDomainSize)
{ }

template<typename TMethod>
Method<TMethod>::Method(const Method<TMethod>& method)
{
    *this = method;
}

template<typename TMethod>
double Method<TMethod>::operator()(const Vector3d& vec) const
{
    const TMethod& tMethod= static_cast<const TMethod&>(*this);
    return tMethod(vec);
}

template<typename TMethod>
double Method<TMethod>::bondFunction(const Vector3d& vec1, const Vector3d& vec2) const
{
    const TMethod& tMethod= static_cast<const TMethod&>(*this);
    return tMethod.bondFunction(vec1, vec2);
}

template<typename TMethod>
double Method<TMethod>::getAveragingDomainSize() const
{
    return averagingDomainSize;
}
template<typename TMethod>
Method<TMethod>::~Method()
{
}
