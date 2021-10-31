#include <armadillo>
#include <iostream>
#include "SolarSystem.h"
#include <cmath>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;


ofstream myfile;

int main(int argc, char const *argv[]){
	SolarSystem SS;
	SS.initialise();
	SS.derivative(SS.SV);

	myfile.open("data",ios::out);
	
	
	
	for (int i = 0; i < 1000; ++i){
		SS.SV = SS.LF(SS.SV);
		/*
		for (int j = 0; j < SS.N; j++){
			myfile << setiosflags(ios::showpoint | ios::uppercase);
			myfile << SS.SV(j,0);
			myfile << " " << SS.SV(j,1);
			myfile << " " << SS.SV(j,2);
			myfile << " " << SS.SV(j,3);
			myfile << " " << SS.SV(j,4);
			myfile << " " << SS.SV(j,5);
			myfile << " ";
			
		}
		myfile <<  endl;
		*/
		cout << sqrt(SS.SV(0,0)*SS.SV(0,0) + SS.SV(0,1)*SS.SV(0,1) + SS.SV(0,2)*SS.SV(0,2)) <<endl;
		//cout << SS.RK4(SS.SV(0,0));
	}


	myfile.close();
	myfile.open("Mass", ios::out);
	myfile << SS.M;
	myfile.close();
	return 0;
}