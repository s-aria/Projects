/*
** ~ Class: Planet ~
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <fstream>
#include "constants.h"

using namespace std;

class Planet{
	public:
		double      XPos, XPosAU, YPos, YPosAU, VX, VY, FX, FY, MCX, MCY;
		double      Mass, Size, Time;
		double      Energy, Spin, SectorArea;
		double      Aphelion, Perihelion, VMax, VMin;
		const char* Label;
		int         RunID;
		bool        FileOpen;
		ofstream oFile;

		Planet(char const*, int, double, double, double, double, double);
		double Radius();
		double Velocity();
		void Update(Planet**, int, int);
		void Correct(Planet**, int, int);
		void InitSystem(Planet**, int, int);
		void Sector(double, double);
		double Gamma(Planet*);
		void Output();
		void Summary();
		void Close();
};

