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

class BoxConfiguration : public Configuration{
public:
	Matrix3d box, reference_box;
	Vector3i pbc;
	Configuration* getConfiguration(double) const;
	void read(std::string,int);
    void readLMP(const std::string&,const ConfigType& configType);
    void readLMP(const std::string& currentConfigFileName,
                 const std::string& referenceConfigFileName);
    void lmpParser(std::ifstream&, const ConfigType&);

	BoxConfiguration(int,int);
	BoxConfiguration(std::string,int);
	virtual ~BoxConfiguration();
};

#endif /* ATOMISTICSYSTEM_H_ */
