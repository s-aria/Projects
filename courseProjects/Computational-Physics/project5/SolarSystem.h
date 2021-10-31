#include <armadillo>
#include <iostream>
#include <cmath>

using namespace std;
using namespace arma;

class SolarSystem{
	public:
		int N;
		mat MassValues;
		vec M;
		double rho_0;
		double R_0;
		double G;
		double eps;
		mat SV;
		mat dSVdt;
		void initialise();
		mat derivative(mat);
		mat RK4(mat);
		mat LF(mat);


};
