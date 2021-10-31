#include <iostream>
#include <stdio.h>
#include <armadillo>
#include <cmath>
#include <math.h>
#include "SolarSystem.h"
#include "gaussiandeviate.h"

using namespace std;
using namespace arma;

void SolarSystem::initialise(){
	cout << "Enter number of particles, N:" << endl;
	cin >> N;
	SV = zeros(N,6);
	long idum = -1;
	R_0 = 20.0; //ly
	
	
	for (int i = 0; i < N; ++i){
		double u = ran2(&idum);
		double v = ran2(&idum);
		double w = ran2(&idum);

		double phi = 2*M_PI*w;
		double r = R_0*pow(u,1.0/3.0);
		double theta = acos(1-2*v);
		
		SV(i,0) = r*sin(theta)*cos(phi);
		SV(i,1) = r*sin(theta)*sin(phi);
		SV(i,2) = r*cos(theta);
	}

	M = zeros(N);
	for (int j = 0; j < N; j++){
		M(j) = gaussian_deviate(&idum)+10;
	}
	rho_0 = (sum(M)*3)/(4*M_PI*(R_0*R_0*R_0));
	
	
}

mat SolarSystem::derivative(mat SV){	
	mat dSVdt = zeros(N,6); // the derivative of the values, dx, dy, dv_x, dv_y
	G =  3*M_PI/(32*rho_0); // from tau_crunch = sqrt(3pi/32rhoG)
	eps = 0.15; //ly
	
	for (int i = 0; i < N; i++){
		vec F = zeros(4);
		dSVdt(i,0) = SV(i,3); //x' = v_x
		dSVdt(i,1) = SV(i,4); //y' = v_y
		dSVdt(i,2) = SV(i,5); //z' = v_z
		
		for (int j = 0; j < N; j++){
			if (j!=i){
				double dx = SV(i,0) - SV(j,0);
				double dy = SV(i,1) - SV(j,1);
				double dz = SV(i,2) - SV(j,2);
				double r = sqrt(dx*dx + dy*dy + dz*dz);
			
				F(3) = -G*M(i)*M(j)/((r*r*r) + eps*eps); //Total force 
				
				F(0) += F(3) * dx; //Force in x direction
				
				F(1) += F(3) * dy; //Force in y direction

				F(2) += F(3) * dz; //Force in z direction
			}
		}
		dSVdt(i,3) = F(0)/M(i); //a_x
		dSVdt(i,4) = F(1)/M(i); //a_y
		dSVdt(i,5) = F(2)/M(i); //a_z

	}
	
	return dSVdt;
}


mat SolarSystem::RK4(mat SV){
	double dt = 0.1;
	mat K1,K2,K3,K4 = zeros(N,6);
	K1 = derivative(SV);
	K2 = derivative(SV + (K1/2)*dt);
	K3 = derivative(SV + (K2/2)*dt);
	K4 = derivative(SV + K3*dt);
	SV += (1/6.0)*dt*(K1 + 2*(K2 + K3) + K4);
	
	return SV;
}

mat SolarSystem::LF(mat SV){
	double dt = 0.1;
		
	dSVdt = derivative(SV);	
	mat V1 = SV.submat(0,3,N-1,5) + dt/2.0*(dSVdt.submat(0,3,N-1,5));
			
	SV.submat(0,0,N-1,2) += dt*SV.submat(0,3,N-1,5);
	dSVdt = derivative(SV);
	SV.submat(0,3,N-1,5) = V1 + dt/2.0*(dSVdt.submat(0,3,N-1,5));

	return SV;
}
