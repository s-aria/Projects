#include <armadillo>
#include <iostream>
#include "SolarSystem.h"
#include <cmath>
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
	
	for (int i = 0; i < 1000000; ++i){
		SS.SV = SS.RK4(SS.SV);
		myfile << setiosflags(ios::showpoint | ios::uppercase);
		//Writing to file the x and y positions:
		/*Sun [UNCOMMENT THIS PART FOR WANTING TO INCLUDE THE SUN]
		myfile << SS.SV(0,0);
		myfile << " " << SS.SV(0,1);
		*/
		//Earth
		myfile << SS.SV(1,0);
		myfile << " " << SS.SV(1,1);
		//Jupiter
		myfile << " " << SS.SV(2,0);
		myfile << " " << SS.SV(2,1);
		//Mars
		myfile << " " << SS.SV(3,0);
		myfile << " " << SS.SV(3,1);
		//Venus
		myfile << " " << SS.SV(4,0);
		myfile << " " << SS.SV(4,1);
		//Saturn
		myfile << " " << SS.SV(5,0);
		myfile << " " << SS.SV(5,1);
		//Mercury
		myfile << " " << SS.SV(6,0);
		myfile << " " << SS.SV(6,1);
		//Uranus
		myfile << " " << SS.SV(7,0);
		myfile << " " << SS.SV(7,1);
		//Neptune
		myfile << " " << SS.SV(8,0);
		myfile << " " << SS.SV(8,1);
		//Pluto
		myfile << " " << SS.SV(9,0);
		myfile << " " << SS.SV(9,1) << endl;
			
	}
	myfile.close();
	//cout << SS.RK4(SS.SV) << endl;
	return 0;
}