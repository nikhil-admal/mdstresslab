#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>
using namespace std;

int main()
{
    ifstream initialFile;
    FILE * pOutputFile;
    
    initialFile.open("config_orig.data");
    pOutputFile = fopen("config.data","w");

    int nat, dimension;
    double *x, *y, *z;
    double x_temp, y_temp, z_temp;
    string *element;
    string str;

    initialFile >> dimension >> nat;

    getline(initialFile,str);
    getline(initialFile,str);
    getline(initialFile,str);

    x = new double[nat];
    y = new double[nat];
    z = new double[nat];
    element = new string[nat];

    for (int i = 0; i < nat; i++)
    {
        initialFile >> element[i] >> x[i] >> y[i] >> z[i];
        getline(initialFile,str);
    }

    initialFile.close();
    initialFile.open("config_orig.data");

    getline(initialFile,str);
    fprintf(pOutputFile,"%s\n",str.c_str());
    getline(initialFile,str);
    fprintf(pOutputFile,"%s\n",str.c_str());
    getline(initialFile,str);
    fprintf(pOutputFile,"%s\n",str.c_str());

    for (int i = 0; i < nat; i++)
    {
        getline(initialFile,str);
        fprintf(pOutputFile,"%s\n",str.c_str());
    }

    for (int i = 0; i < nat; i++)
    {
        x_temp = 1.01 * x[i] + 0.0 * y[i] + 0.0 * z[i];
        y_temp = 0.0 * x[i] + 1.01 * y[i] + 0.0 * z[i];
        z_temp = 0.0 * x[i] + 0.0 * y[i] + 1.01 * z[i];
        fprintf(pOutputFile," %s %25.16lf %25.16lf %25.16lf %25.16lf %25.16lf %25.16lf\n",element[i].c_str(),x_temp,y_temp,z_temp,0.0,0.0,0.0);
    }

    delete [] x;
    delete [] y;
    delete [] z;
    delete [] element;

    initialFile.close();
    fclose(pOutputFile);
    return 0;
}
