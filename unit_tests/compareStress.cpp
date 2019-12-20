void compareStress(std::string stressFile)
{
	std::string str;
	std::ifstream file(stressFile+".stress");
	if(!file) MY_ERROR(std::string("ERROR: " + stressFile + ".stress could not be opened for reading!"));
	std::ifstream fileReference(stressFile+"Reference.stress");
	if(!fileReference) MY_ERROR(std::string("ERROR: " + stressFile + "Reference.stress could not be opened for reading!"));
	double tol= 1e-10;

	int ngrid,ngridReference;
	file >> ngrid;
	fileReference >> ngridReference;
	if(ngrid != ngridReference)
	{
		std::string errorMessage= "Test failed in " + stressFile + ". Number of grid points do not match.";
		MY_ERROR(errorMessage);
	}
	std::getline(file,str); std::getline(fileReference,str);
	double data,dataReference;
	for (int i_point=0; i_point<ngrid; i_point++)
	{
		for (int index=0; index<DIM*DIM; index++)
		{
			file >> std::setprecision(15) >> data;
			fileReference >> std::setprecision(15) >> dataReference;
			if (abs(data-dataReference) > tol)
			{
				std::string errorMessage= "Test failed in " + stressFile + ". Data error.";
				MY_ERROR(errorMessage);
			}

		}
	}
}


