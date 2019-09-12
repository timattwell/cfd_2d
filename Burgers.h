#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include "Model.h"
//#include "mpi.h"


class Burgers
{
public:
	Burgers(Model m);
	
	double Getax()     const { return ax; }
	double Getay()     const { return ay; }
	double Getb()     const { return b; }
	double Getc()     const { return c; }

    void SetVelField();
    void Solver();
    void VelFile();
	void Energy();
	

	double x0;
	double y0;
	double Lx;
	double Ly;
	double T;
	int    Nx;
	int    Ny;
	int    Nt;
	double dx;
	double dy;
	double dt;
	
	//values required for parallisation
	//int    Px;
	//int    Py;
	//int    locNx;
	//int    extra;

    double ax;
    double ay;
    double b;
    double c;
	
	double* uField;
	double* vField;
	double* u1;
	double* v1;
	double* u;
	double* v;
	double* usend;
	double* vsend;
	//double vField[2001][2001];

};


#endif