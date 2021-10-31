#include <iostream>
#include <stdio.h>
#include <armadillo>
#include <cmath>
#include "SolarSystem.h"

using namespace std;
using namespace arma;

void SolarSystem::initialise(){
	cout << "Sun, Earth, Jupiter, Mars, Venus, Saturn, Mercury, Uranus, Neptune, Pluto" << endl;
	cout << "Number of celestial objects in the order shown:" << endl;
	cin >> n;
	SV = zeros(n,4);
				//Sun	Earth			Jupiter		Mars		Venus		Saturn			Mercury		Uranus		Neptune		Pluto	
	MassValues << 1.0 << 3.0024E-6 << 9.5426E-4 << 3.2260E-7 << 2.4469E-6 << 2.8571E-4 << 1.6595E-7 << 4.3643E-5 << 5.1485E-5 << 6.5808E-9<<endr;
	M = zeros(n);
	for (int i = 0; i < n; ++i){
		M(i) = MassValues(i);
	}
	// Give the initial positions and velocities. The distances to the sun is given by AU and the velocities given by AU/yr
	//Sun:
	mat InitialValues = zeros(10,4);

	InitialValues(0,0) = 0; InitialValues(0,1) = 0;
	InitialValues(0,2) = 0; InitialValues(0,3) = 0;
	//Earth:
	InitialValues(1,0) = 1.0; InitialValues(1,1) = 0.0;
	InitialValues(1,2) = 0; InitialValues(1,3) = 6.28;
	//Jupiter:
	InitialValues(2,0) = 5.20; InitialValues(2,1) = 0;
	InitialValues(2,2) = 0; InitialValues(2,3) = 2.75;
	//Mars:
	InitialValues(3,0) = 1.52; InitialValues(3,1) = 0;
	InitialValues(3,2) = 0; InitialValues(3,3) = 5.08;	
	//Venus:
	InitialValues(4,0) = 0.72; InitialValues(4,1) = 0;
	InitialValues(4,2) = 0; InitialValues(4,3) = 7.38;	
	//Saturn:
	InitialValues(5,0) = 9.54; InitialValues(5,1) = 0;
	InitialValues(5,2) = 0; InitialValues(5,3) = 2.03;	
	//Mercury:
	InitialValues(6,0) = 0.39; InitialValues(6,1) = 0;
	InitialValues(6,2) = 0; InitialValues(6,3) = 9.98;	
	//Uranus:
	InitialValues(7,0) = 19.19; InitialValues(7,1) = 0;
	InitialValues(7,2) = 0; InitialValues(7,3) = 1.43;
	//Neptune:
	InitialValues(8,0) = 30.06; InitialValues(8,1) = 0;
	InitialValues(8,2) = 0; InitialValues(8,3) = 1.15;
	//Pluto:
	InitialValues(9,0) = 39.53; InitialValues(9,1) = 0;
	InitialValues(9,2) = 0; InitialValues(9,3) = 0.98;

	for (int j = 0; j < 4; ++j){
		for (int i = 0; i < n; ++i){
			SV(i,j) = InitialValues(i,j);
		}
		 
	}
}

mat SolarSystem::derivative(mat SV){	
	mat dSVdt = zeros(n,4); // the derivative of the values, dx, dy, dv_x, dv_y
	double G = 4*M_PI*M_PI; //in AU**2/yr**2
	for (int i = 0; i < n; i++){
		vec F = zeros(3);
		dSVdt(i,0) = SV(i,2);
		dSVdt(i,1) = SV(i,3);
		for (int j = 0; j < n; j++){
			if (j!=i){
				double dx = SV(i,0) - SV(j,0);
				double dy = SV(i,1) - SV(j,1);
				
				double r = sqrt(dx*dx + dy*dy);
			
				F(2) = -G*M(i)*M(j)/(r*r*r); //Total force 
				
				F(0) += F(2) * dx; //Force in x direction
				
				F(1) += F(2) * dy; //Force in y direction
			}
		}
		dSVdt(i,2) = F(0)/M(i); //a_x
		dSVdt(i,3) = F(1)/M(i); //a_y
	}
	return dSVdt;
}

mat SolarSystem::RK4(mat SV){
	double dt = 0.001;
	mat K1,K2,K3,K4 = zeros(n,4);
	K1 = derivative(SV);
	K2 = derivative(SV + (K1/2)*dt);
	K3 = derivative(SV + (K2/2)*dt);
	K4 = derivative(SV + K3*dt);
	SV += (1/6.0)*dt*(K1 + 2*(K2 + K3) + K4);
	
	return SV;
}