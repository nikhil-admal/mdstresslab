/*
 * Exception.cpp
 *
 *  Created on: Nov 10, 2019
 *      Author: Nikhil
 */

#include "Exception.h"

Exception::Exception(const std::string& _msg) : msg(_msg){}

const char* Exception::what() const
{
	return (this->msg).c_str();
}

Exception::~Exception() {
	// TODO Auto-generated destructor stub
}

