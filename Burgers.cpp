/*
 * Tim Attwell 00985020
 * HPC assignment - serial code for 2D burgers. 
 * Broken parallel code is included but commented out
 * 
 * 
 * 
 * 
 */
#include <chrono>
#include <ctime>
#include <ratio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Burgers.h"
#include "Model.h"

//#include "mpi.h"

using namespace std::chrono;
using namespace std;

int main(int argc, char* argv[])
{
		
	//int err = MPI_Init(&argc, &argv);
	//if (err != MPI_SUCCESS) {
	//	cout << "MPI failed to initialise" << endl;
	//}
	
    Model m(argc, argv);
    Burgers bur(m);

// outputs the selected case variables
    cout << "ax = " << bur.Getax() << endl;
	cout << "ay = " << bur.Getay() << endl;
	cout << "b = " <<bur.Getb() << endl;
	cout << "c = " << bur.Getc() << endl;

// sets the initial velocity field using nested for-loops 
    bur.SetVelField();

    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
	
    hrc::time_point start = hrc::now();

// serial solver for the 2D burgers equation.
	
	
	//if (m.GetP() == 2) {
    //bur.ParSolver();
	//} else {
	bur.Solver();
	//}

    hrc::time_point end = hrc::now();
	
	duration<double> time_span = duration_cast<duration<double>>(end-start);

// prints the velocity vectors to file
	bur.VelFile();

// calculates the energy over the domain
	bur.Energy();
	
	cout << "The integration loop took: " << time_span.count() << "seconds to run." << endl;

    return 0;
}

/*
 * Sets up the Burgers class, constructiong it based of values obtained from the Model class.
 * All system constants required for the Burgers functions are obtained here.
 * 
 */
Burgers::Burgers(Model m)
{
    x0 = m.GetX0();
    y0 = m.GetY0();
    Lx = m.GetLx();
    Ly = m.GetLy();
    T = m.GetT();
    Nx = m.GetNx();
    Ny = m.GetNy();
    //Px = m.GetPx();
    //Py = m.GetPy();
    Nt = m.GetNt();
    dx = m.GetDx();
    dy = m.GetDy();
    dt = m.GetDt();

    // Physics
    ax = m.GetAx();
    ay = m.GetAy();
    b = m.GetB();
    c = m.GetC();
	
	cout << "Burgers object has been created." << endl;
}
/*
 * 
 * SetVelField sets the initial conditions based on the equation given in the handout
 * 
 */ 

void Burgers::SetVelField()
{
    uField = new double[Nx * Ny];
    vField = new double[Nx * Ny];
    for(int ix = 0; ix < Nx; ix++) {
	for(int iy = 0; iy < Ny; iy++) {
	    double r = sqrt(pow((ix*dx)-5, 2.0) + pow((iy*dy)-5, 2.0));
	    if(r <= 1.0) {
		uField[ix * Ny + iy] = 2.0 * pow(1.0 - r, 4.0) * ((4.0 * r) + 1.0);
		vField[ix * Ny + iy] = 2.0 * pow(1.0 - r, 4.0) * ((4.0 * r) + 1.0);
	    } else if(r > 1) {
		uField[ix * Ny + iy] = 0.0;
		vField[ix * Ny + iy] = 0.0;
	    }
	}
    }
    //cout << "Field " << vField[1000 * Ny + 1000] << endl;
}

/*
 * The Solver function creates some temporary or "old" variables (u,v)
 * from (uField,vField). 
 * These are then used in a set of for loops to find the updated values,
 * which are assigned to (uField, vField)
 * 
 * Updated variable then go though another set of loops to set the 
 * boundary conditions to zero if any have changed.
 * 
 */
 
void Burgers::Solver()
{
    u = new double[Nx * Ny];
    v = new double[Nx * Ny];
    for(double t = 0; t < dt*Nt; t+=dt) {
		
	for(int i = 0; i < Nx; i++) {
	    for(int j = 0; j < Ny; j++) {
		u[i * Ny + j] = uField[i * Ny + j];
		v[i * Ny + j] = vField[i * Ny + j];
	    }
	}
	//cout << "odne that" <<endl;
	for(int i = 1; i < Nx-1; i++) {
	    for(int j = 1; j < Ny-1; j++) {
		uField[i * Ny + j] = u[i * Ny + j] -
		    ((dt / dx) * (ax + b * u[i * Ny + j]) * (u[i * Ny + j] - u[(i - 1) * Ny + j])) -
		    ((dt / dy) * (ay + b * v[i * Ny + j]) * (u[i * Ny + j] - u[i * Ny + (j - 1)])) +
		    ((c * dt / (dx*dx)) * (u[(i + 1) * Ny + j] - 2 * u[i * Ny + j] + u[(i - 1) * Ny + j])) +
		    ((c * dt / (dy*dy)) * (u[i * Ny + (j + 1)] - 2 * u[i * Ny + j] + u[i * Ny + (j - 1)]));
		vField[i * Ny + j] = v[i * Ny + j] -
		    ((dt / dx) * (ax + b * v[i * Ny + j]) * (v[i * Ny + j] - v[(i - 1) * Ny + j])) -
		    ((dt / dy) * (ay + b * u[i * Ny + j]) * (v[i * Ny + j] - v[i * Ny + (j - 1)])) +
		    ((c * dt / (dx*dx)) * (v[(i + 1) * Ny + j] - 2 * v[i * Ny + j] + v[(i - 1) * Ny + j])) +
		    ((c * dt / (dy*dy)) * (v[i * Ny + (j + 1)] - 2 * v[i * Ny + j] + v[i * Ny + (j - 1)]));
		if((i == Nx - 1) || (j == Ny - 1) || (i == 1) || (j ==1)) {
		    uField[i * Ny + j] = 0;
		    vField[i * Ny + j] = 0;
		}
	    }
	}
	
	cout << (t/T)*100 << "%" << endl;
	//cout << "Got through" << endl;
    }
}

/*
void Burgers::ParSolver()
{
	int size, rank;
	const int sendsize = 2001*1001;
	const int recsize = 2001*2001;
	const int over = 1;

	
	// Find the rank and communicator size
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	
	
	// Make sure the number of processes are a factor of 2^n
	if (size > 2) {
		cout << "Number of processes must be 2 or less!" << endl;
	}
	
	// Define the local dimentions of the sub-domains
	extra = Nx % 2;
	locNx = (Nx/size)+extra;
	
    usend = new double[Ny*locNx];
    vsend = new double[Ny*locNx];
	
	// Add an additional row in the centre of the velocity matricies to facilitate 
	// simple splitting in MPI_Scatter()
	
	for(int i = 0; i < Nx + 1; i++) {
	    for(int j = 0; j < Ny; j++) {
		if (i <= 1001) {
		usend[i * Ny + j] = uField[i * Ny + j];
		vsend[i * Ny + j] = vField[i * Ny + j];
		}else {
		usend[i * Ny + j] = uField[(i-1) * Ny + j];
		vsend[i * Ny + j] = vField[(i-1) * Ny + j];
		}
		}
	
	}
	
	
	u = new double[Ny*locNx];
    v = new double[Ny*locNx];
	// MPI scatter to send the halves of (u,v)send to (u,v)1 
	MPI_Scatter(usend, sendsize, MPI_DOUBLE,
				   u1, sendsize, MPI_DOUBLE,
				0, MPI_COMM_WORLD);

	MPI_Scatter(vsend, sendsize, MPI_DOUBLE,
				   v1, sendsize, MPI_DOUBLE,
				0, MPI_COMM_WORLD);
	
	// The time integration now occurs as in serial case, over domain locNx x Ny
	 
    for(double t = 0; t < dt*T; t+=dt) {
	
	for(int i = 0; i < locNx; i++) {
	    for(int j = 0; j < Ny; j++) {
		u[i * Ny + j] = u1[i * Ny + j];
		v[i * Ny + j] = v1[i * Ny + j];
	    }
	}
	//cout << "done that" <<endl;
	for(int i = 0; i < locNx; i++) {
	    for(int j = 1; j < Ny-1; j++) {
		u1[i * Ny + j] = u[i * Ny + j] -
		    ((dt / dx) * (ax + b * u[i * Ny + j]) * (u[i * Ny + j] - u[(i - 1) * Ny + j])) -
		    ((dt / dy) * (ay + b * v[i * Ny + j]) * (u[i * Ny + j] - u[i * Ny + (j - 1)])) +
		    ((c * dt / (dx*dx)) * (u[(i + 1) * Ny + j] - 2 * u[i * Ny + j] + u[(i - 1) * Ny + j])) +
		    ((c * dt / (dy*dy)) * (u[i * Ny + (j + 1)] - 2 * u[i * Ny + j] + u[i * Ny + (j - 1)]));
		v1[i * Ny + j] = v[i * Ny + j] -
		    ((dt / dx) * (ax + b * v[i * Ny + j]) * (v[i * Ny + j] - v[(i - 1) * Ny + j])) -
		    ((dt / dy) * (ay + b * u[i * Ny + j]) * (v[i * Ny + j] - v[i * Ny + (j - 1)])) +
		    ((c * dt / (dx*dx)) * (v[(i + 1) * Ny + j] - 2 * v[i * Ny + j] + v[(i - 1) * Ny + j])) +
		    ((c * dt / (dy*dy)) * (v[i * Ny + (j + 1)] - 2 * v[i * Ny + j] + v[i * Ny + (j - 1)]));
	
		}
	}
	// MPI_Gathers gathers the results, inculding the row added before

	MPI_Gather(u1, sendsize, MPI_DOUBLE,
			usend, recsize ,MPI_DOUBLE,
					  0, MPI_COMM_WORLD);
	MPI_Gather(v1, sendsize, MPI_DOUBLE,
			vsend, recsize,MPI_DOUBLE,
					  0, MPI_COMM_WORLD);

		// (u,v)send are converted to (u,v)Field, and the row added in
		//the previous processes is removed
		 
		for(int i = 0; i < Nx + 1; i++) {
	    for(int j = 0; j < Ny; j++) {
		if (i <= 1001) {
		uField[i * Ny + j] = usend[i * Ny + j];
		vField[i * Ny + j] = vsend[i * Ny + j];
		}else {
		uField[i * Ny + j] = usend[(i+1) * Ny + j];
		vField[i * Ny + j] = vsend[(i+1) * Ny + j];
		}
		}
	
	}//boundary conditions set
	for(int i = 1; i < Nx-1; i++) {
	for(int j = 1; j < Ny-1; j++) {
		if((i == Nx - 1) || (j == Ny - 1)) {
		    uField[i * Ny + j] = 0;
		    vField[i * Ny + j] = 0;
		}
	}
	}
	
	cout << (t/T)*100 << "%" << endl;
	//cout << "Got through" << endl;
    
	// MPI is shut down after integration is complete.
	MPI_Finalize();
}
*/

/*
 * VelFile prints the values obtained to files in a single column,
 * to make reading by a plotting program easier. 
 */

void Burgers::VelFile()
{
    ofstream uOut("uField.txt", ios::out | ios::trunc);
    ofstream vOut("vField.txt", ios::out | ios::trunc);
    uOut.precision(8);
    vOut.precision(8);
    for(int i = 0; i < Nx; i++) {
	for(int j = 0; j < Ny; j++) {
		uOut << uField[i * Ny + j] << endl;
		vOut << vField[i * Ny + j] << endl;
	}
    }
    uOut.close();
    vOut.close();
}
/*
 * The energy in the system is calculated based on the integration 
 * of the total velocity vector over the full domain.This is then 
 * printed to screen.
 * 
 */
void Burgers::Energy()
{
	double U;
	double E=0;
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Nx; j++) {
			U = (uField[i * Ny + j]*uField[i * Ny + j]) + (vField[i * Ny + j]*vField[i * Ny + j]);
			E = E + U;
		}
	}
	E = 0.5 * E * dx * dy;
	cout << "Energy = " << E << endl;
}