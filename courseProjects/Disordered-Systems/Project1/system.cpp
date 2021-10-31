#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <armadillo>

using namespace arma;
using namespace std;

//Simulation parametres and values:
int N_c = 8;            //number of particles
int Ntot = N_c*N_c*N_c;	
double rho = 1.0; 		//density
double T = 100.0;		//temperature K
float b = 5.960;        //lattice constant
float m = 39.948;       // Angstrom Argon
float eps_kb = 119.8;   //Kelvin 
float sigma = 3.405;    //Angstrom

//-----The following lines will define the position, velocity and the acceleration.
mat R = zeros(4*Ntot, 3);
mat V;
mat A = zeros(4*Ntot, 3);
//------------------

int n = 0;
double L = pow(Ntot/rho, 1.0/3);

void buildFCClattice(float b){
	rowvec xCell;
	rowvec yCell;
	rowvec zCell;

    xCell << 0 << 0.5 << 0 << 0.5 << endr;
    yCell << 0 << 0.5 << 0.5 << 0 << endr;
    zCell << 0 << 0 << 0.5 << 0.5 << endr;
    
	//find M large enough to fit N atoms in a fcc lattice
    //can change this later on to include in the function.
    
    int M = 1;
    while (4 * M * M * M < N_c)
        ++M;

    
    for (int i = 0; i < N_c; ++i){
    	for (int j = 0; j < N_c; ++j){
    		for (int k = 0; k < N_c; ++k){
    			for (int d = 0; d < 4; ++d){
                    R(n, 0) = (i + xCell[d]) * b;
                    R(n, 1) = (j + yCell[d]) * b;
                    R(n, 2) = (k + zCell[d]) * b; 	
        			++n; 
                }
    		}
    	}
    }
}


void initVelocities(){
    //applying random initial velocities
    V = randn(4*Ntot, 3);

    //CM momentum
    double Px = sum(V(0));
    double Py = sum(V(1));
    double Pz = sum(V(2));

    V(0) = V(0) - Px/Ntot;
    V(1) = V(1) - Px/Ntot;
    V(2) = V(2) - Px/Ntot;
    
    ofstream velocityFile("Velocities.txt");
    for (int i = 0; i < n; ++i){
        velocityFile    << V(i,0) << " " << V(i,1) << " " << V(i,2) <<  endl;
    }
    velocityFile.close();
}

void computeAccelerations() {
    for (int i = 0; i < Ntot-1; i++)        // loop over all distinct pairs i,j
        for (int j = i+1; j < Ntot; j++) {
            double rij[3];               // position of i relative to j
            double rSqd = 0;
            for (int k = 0; k < 3; k++) {
                rij[k] = R(i,k) - R(j,k);

                // closest image convention
                if (abs(rij[k]) > 0.5 * L){
                    if (rij[k] > 0)
                        rij[k] -= L;
                    else
                        rij[k] += L;
                }
                rSqd += rij[k] * rij[k];
            }
            double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
            for (int k = 0; k < 3; k++) {
                 A(i,k) += rij[k] * f;
                 A(j,k) -= rij[k] * f;
            }
        }
}


void velocityVerlet(double dt) {
    computeAccelerations();
    for (int i = 0; i < 4*Ntot; i++){
        for (int k = 0; k < 3; k++){
            R(i,k) += V(i,k) * dt + 0.5 * A(i,k) * dt * dt;
            
            // use periodic boundary conditions            
            if (R(i,k) < 0)
                R(i,k) += L;
            if (R(i,k) >= L)
                R(i,k) -= L;
            V(i,k) += 0.5 * A(i,k) * dt;
        }
    }
    computeAccelerations();
    for (int i = 0; i < Ntot; i++)
        for (int k = 0; k < 3; k++)
            V(i,k) += 0.5 * A(i,k) * dt;
}

int main(){
    buildFCClattice(b);
    initVelocities();
    ofstream file("data.xyz");
    for (int t = 0; t < 500; ++t){
        
        velocityVerlet(0.04);

        file    << R.n_rows << '\n';
        file    << "MD project. Argon atoms" << '\n';
        for (int i = 0; i < n; ++i){
                file    << "Ar" << " " << R(i,0) << " " << R(i,1)  << " " << R(i,2) << " " << V(i,0) << " " << V(i,1) << " " << V(i,2) <<  endl;
        }
        
    }
    file.close();
    return 0;
}