#include <iostream>
#include <cmath>

#include "Model.h"
#include "Burgers.h"

using namespace std;

Model::Model(int argc, char* argv[])
{
	ParseParameters(argc,argv);
}

Model::~Model()
{
	
}


void Model::ParseParameters(int argc, char* argv[])
{
	// These values are constant though all tests and so can be assigned in-class.
	x0 = -5.0;
	y0 = -5.0;
	Lx = 10.0;
	Ly = 10.0;
	T = 1.0;
	Nx = 2001;
	Ny = 2001;
	Nt = 4000;
	dx = Lx / (Nx-1);
	dy = Ly / (Ny-1);
	dt = T / Nt;
	
	/*
	P = stod(argv[1]);
	
	if (P == 2) {
		Py = 1;
		Px = 2;
	} else if(P==1){
		Px = 1;
		Py = 1;
	}else {
		Py = sqrt(P);
		Px = sqrt(P);
	}
	*/

// variables taken from command line
	ax = stod(argv[2]);
	ay = stod(argv[3]);
	b = stod(argv[3]);
	c = stod(argv[4]);
	    // Physics

}

void Model::ValidateParameters()
{
	cout << "Parameters are valid." << endl;
}

bool Model::IsValid()
{
	ValidateParameters();
	return true;
}