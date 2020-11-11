// dump2str++.cpp
// Min Shi
// Last checked and modified: 01/22/2020
//
// This code takes 2 lammps dump files (one initial configuration and one final configuration)
// as input files and generates an appropriate MDstresslab++ input file.
// 
// You also need to specify the element type in the dump file. See usage.
//
// Since the dump file can contain multiple steps, only the first timestep is being processed. 
// Please arrange your dump file so that the appropriate data are being handled.
//
// Version 1.1 
// Developing Version
// Warning: The Code Might Contains BUGS !!!
// Welcome any revision to my code.
// If you find any bug, feel free to send email to shixx597@umn.edu
//

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>
using namespace std;

ifstream initialFile;
ifstream finalFile;
FILE * pOutputFile;
string element[100], boundary, stmp;
int numAtom, nAtomType, col1[8] = {0,0,0,0,0,0,0,0}, col2[8] = {0,0,0,0,0,0,0,0}, count = -1, itmp, iatom_type; // x y z vx vy vz atomtype last_column // starting from column 5
double boxLow, boxHigh, pos1, pos2, pos3, vel1 = 0.0, vel2 = 0.0, vel3 = 0.0;

bool fexists(const string& filename);
void readingNAtm();
void readingBox();
void isdata1(ifstream & filestream);
void isdata2(ifstream & filestream);
void readingPosition();
void closeFile();

bool fexists(const string& filename) // check whether a file exists or not
{
    ifstream ifile(filename.c_str());
    return (bool)ifile;
}

void readingNAtm() // Reading number of atoms in the configuration
{
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    initialFile >> numAtom;    // get the number of atoms
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    finalFile >> itmp;    // get the number of atoms
    if (numAtom != itmp)
    {
        cout << "The initial configuration and the final configuration should have the same # atoms!" << endl;
        cout << "In inital configuration # atoms : " << numAtom << endl;
        cout << "In final configuration # atoms : " << itmp << endl;
        exit(-3);
    }
    getline(initialFile,stmp); // notorious cin mind the \n
    getline(finalFile,stmp);
    fprintf(pOutputFile,"%d\n",numAtom);
}

void readingBox() // Get the samplesize and pbc information
{
    getline(initialFile,stmp);
    //cout << stmp << endl;
    getline(finalFile,stmp); 
    //cout << stmp << endl;

    initialFile >> boxLow >> boxHigh;
    fprintf(pOutputFile,"%25.16lf 0.0 0.0\n",boxHigh - boxLow);
    initialFile >> boxLow >> boxHigh;
    fprintf(pOutputFile,"0.0 %25.16lf 0.0\n",boxHigh - boxLow);    
    initialFile >> boxLow >> boxHigh;
    fprintf(pOutputFile,"0.0 0.0 %25.16lf\n",boxHigh - boxLow); 

    finalFile >> boxLow >> boxHigh;
    fprintf(pOutputFile,"%25.16lf 0.0 0.0\n",boxHigh - boxLow);
    finalFile >> boxLow >> boxHigh;
    fprintf(pOutputFile,"0.0 %25.16lf 0.0\n",boxHigh - boxLow); 
    finalFile >> boxLow >> boxHigh;
    fprintf(pOutputFile,"0.0 0.0 %25.16lf\n",boxHigh - boxLow);

    initialFile.clear();
    initialFile.seekg(0);
    finalFile.clear();
    finalFile.seekg(0);

    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    
    initialFile >> boundary >> boundary >> boundary >> boundary;
    finalFile >> stmp >> stmp >> stmp >> stmp;
    if (stmp == boundary)
        if (boundary == "pp")
        {
            fprintf(pOutputFile," 1   ");
        } 
        else
        {
            fprintf(pOutputFile," 0   ");
        }
    else
    {
        cout << "The boundary condition of the initial configuration and final configuration should be the same!" << endl;
        exit(-4);
    }

    initialFile >> boundary;
    finalFile >> stmp;

    if (stmp == boundary)
        if (boundary == "pp")
        {
            fprintf(pOutputFile," 1   ");
        } 
        else
        {
            fprintf(pOutputFile," 0   ");
        }
    else
    {
        cout << "The boundary condition of the initial configuration and final configuration should be the same!" << endl;
        exit(-4);
    }

    initialFile >> boundary;
    finalFile >> stmp;

    if (stmp == boundary)
        if (boundary == "pp")
        {
            fprintf(pOutputFile," 1 \n");
        } 
        else
        {
            fprintf(pOutputFile," 0 \n");
        }
    else
    {
        cout << "The boundary condition of the initial configuration and final configuration should be the same!" << endl;
        exit(-4);
    }
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    //cout << stmp << endl;
}

void isdata1(ifstream & filestream) // get the correct column information for initial configuration
{
    count = -2;
    do
    {
        filestream >> stmp;
        count++;
        if (stmp == "x")
        {
            col1[0] = col1[0] + count;
        }
        else if (stmp == "y")
        {
            col1[1] = col1[1] + count;
        }
        else if (stmp == "z")
        {
            col1[2] = col1[2] + count;
        }
        else if (stmp == "vx")
        {
            col1[3] = col1[3] + count;
       	}
       	else if (stmp == "vy")
      	{
       	    col1[4] = col1[4] + count;
        }
        else if (stmp == "vz")
        {
            col1[5] = col1[5] + count;
        }
        else if (stmp == "type")
        {
            col1[6] = col1[6] + count;
        }
    }while(!isdigit(stmp.c_str()[0]));
    col1[7] = count - 1;
    //cout << col1[0] << " " << col1[1] << " " << col1[2] << " " << col1[3] << " " << col1[4] << " " << col1[5] << " " << col1[6] << ' ' << col1[7] << endl; 
}

void isdata2(ifstream & filestream) // get the correct column information for final configuration
{
    count = -2;
    do
    {
        filestream >> stmp;
        count++;
        if (stmp == "x")
        {
            col2[0] = col2[0] + count;
        }
        else if (stmp == "y")
        {
            col2[1] = col2[1] + count;
        }
        else if (stmp == "z")
        {
            col2[2] = col2[2] + count;
        }
        else if (stmp == "vx")
        {
            col2[3] = col2[3] + count;
        }
        else if (stmp == "vy")
      	{
            col2[4] = col2[4] + count;
        }
        else if (stmp == "vz")
        {
            col2[5] = col2[5] + count;
        }
        else if (stmp == "type")
        {
            col2[6] = col2[6] + count;
        }
    }while(!isdigit(stmp.c_str()[0]));
    col2[7] = count - 1;
    //cout << col2[0] << " " << col2[1] << " " << col2[2] << " " << col2[3] << " " << col2[4] << " " << col2[5] << " " << col2[6] << ' ' << col2[7] << endl;
}

void readingPosition() // Read atomtype x y z vx vy vz (velocities are optional)
{
    initialFile.clear();
    initialFile.seekg(0);
    finalFile.clear();
    finalFile.seekg(0);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(initialFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    getline(finalFile,stmp);
    //cout << stmp << endl;

    vel1 = 0.0;
    vel2 = 0.0;
    vel3 = 0.0;

    count = 0;

    if (col1[0] == 0 || col1[1] == 0 || col1[2] == 0)
    {
        cout << "Could not find the position of atoms. Something wrong with your input." << endl;
        cout << "Please make sure x y z appear in your initial lammps dump file." << endl;
        exit(-5);
    }

    if (col2[0] == 0 || col2[1] == 0 || col2[2] == 0)
    {
        cout << "Could not find the position of atoms. Something wrong with your input." << endl;
        cout << "Please make sure x y z appear in your final lammps dump file." << endl;
        exit(-5);
    }

    if (col2[3] == 0 || col2[4] == 0 || col2[5] == 0)
    {
        for (int i = 1; i <= numAtom; i++)
        {
            count = 0;
            do
            {
                count++;
                finalFile >> iatom_type;
            }while(col2[6] > count);
            do
            {
                count++;
                finalFile >> pos1;
            }while(col2[0] > count);
            do
            {
                count++;
                finalFile >> pos2;
            }while(col2[1] > count);
            do
            {
                count++;
                finalFile >> pos3;
            }while(col2[2] > count);
            getline(finalFile,stmp);
            fprintf(pOutputFile," %s %25.16lf %25.16lf %25.16lf %25.16lf %25.16lf %25.16lf ",element[iatom_type - 1].c_str(),pos1,pos2,pos3,vel1,vel2,vel3);
        
            // initial configuration
            count = 0;
            do
            {
                count++;
                initialFile >> iatom_type;
            }while(col1[6] > count);
            //cout << element[iatom_type - 1].c_str() << endl;
            // for element type
            do
            {
                count++;
                initialFile >> pos1;
            }while(col1[0] > count);
            //cout <<pos1<< endl;
            do
            {
                count++;
                initialFile >> pos2;
            }while(col1[1] > count);
            do
            {
                count++;
                initialFile >> pos3;
            }while(col1[2] > count);
            getline(initialFile,stmp);
            fprintf(pOutputFile,"%25.16lf %25.16lf %25.16lf\n",pos1,pos2,pos3);
        }
    }
    else
    {
        for (int i = 1; i <= numAtom; i++)
        {
            count = 0;
            do
            {
                count++;
                finalFile >> iatom_type;
            }while(col2[6] > count);
            do
            {
                count++;
                finalFile >> pos1;
            }while(col2[0] > count);
            do
            {
                count++;
                finalFile >> pos2;
            }while(col2[1] > count);
            do
            {
                count++;
                finalFile >> pos3;
            }while(col2[2] > count);
            do
            {
                count++;
                finalFile >> vel1;
            }while(col2[3] > count);
            do
            {
                count++;
                finalFile >> vel2;
            }while(col2[4] > count);
            do
            {
                count++;
                finalFile >> vel3;
            }while(col2[5] > count);
            getline(finalFile,stmp);
            fprintf(pOutputFile," %s %25.16lf %25.16lf %25.16lf %25.16lf %25.16lf %25.16lf ",element[iatom_type - 1].c_str(),pos1,pos2,pos3,vel1,vel2,vel3);
            
            // initial configuration
            count = 0;
            do
            {
                count++;
                initialFile >> iatom_type;
            }while(col1[6] > count);
            //cout << element[iatom_type - 1].c_str() << endl;
            // for element type
            do
            {
                count++;
                initialFile >> pos1;
            }while(col1[0] > count);
            //cout <<pos1<< endl;
            do
            {
                count++;
                initialFile >> pos2;
            }while(col1[1] > count);
            do
            {
                count++;
                initialFile >> pos3;
            }while(col1[2] > count);
            getline(initialFile,stmp);
            fprintf(pOutputFile,"%25.16lf %25.16lf %25.16lf\n",pos1,pos2,pos3);
        }
    }
}

void closeFile()
{
    initialFile.close();
    finalFile.close();
    fclose(pOutputFile);
}

int main(int argc, char **argv)
{
    if (argc < 6) // command line input
    {
        cout << "You must enter at least 5 arguments for this program!" << endl;
        cout << "The Simulations are assumed to be 3D." << endl;
        cout << "The data of the LAMMPS dump files should be in a good sequence, i.e. id type x y z vx vy vz. (Velocities are optional)" << endl;
        cout << "Velocities are only read in when all vx vy vz appear in your initial dump file or final dump file." << endl;
        cout << "Otherwise, velocities are set to 0 in all 3 dimensions." << endl;
        cout << "Since the dump file can contain multiple steps, only the first timestep is processed." << endl;
        cout << "Please arrange your dump file so that the appropriate data are being handled." << endl;
        cout << "Usage:" << endl;
        cout << "./dump2str InitialDumpFile FinalDumpFile Outputfile #Element ElementType1 ... ElementTypeN" << endl;
        cout << "e.g." << endl;
        cout << "./dump2str dump.initial.0 dump.final.1 config.data 2 Si C" << endl;
        return -1; // Bad input error
    }
    else
    {
        if (fexists(argv[1]) and fexists(argv[2]))
        {
            cout << "Processing..." << endl;
            initialFile.open(argv[1]);
            finalFile.open(argv[2]);
            pOutputFile = fopen(argv[3],"w");
            nAtomType = atoi(argv[4]);
            for (int i = 0; i < nAtomType; i++)
            {
                element[i] = argv[i + 5];
            }
            readingNAtm();
            readingBox();
			isdata1(initialFile);
			isdata2(finalFile);
            readingPosition();
            closeFile();
            return 0;
        }
        else
        {
            cout << "The input file(s) do(es) not exist!" << endl;
            return -2; // files not exist error
        }
    }
}

