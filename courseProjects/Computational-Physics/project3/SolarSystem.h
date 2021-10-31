#include <armadillo>
#include <iostream>
#include <cmath>

using namespace std;
using namespace arma;

class SolarSystem{
	public:
		int n;
		mat MassValues;
		vec M;
		mat SV;
		void initialise();
		mat derivative(mat);
		mat RK4(mat);

};
