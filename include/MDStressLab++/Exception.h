/*
 * Exception.h
 *
 *  Created on: Nov 10, 2019
 *      Author: Nikhil
 */

#ifndef EXCEPTION_H_
#define EXCEPTION_H_
#include <string>

class Exception {
public:
	Exception(const std::string&);
	virtual ~Exception();
	std::string msg;
	virtual const char* what() const;
};

#endif /* EXCEPTION_H_ */
