/*
** ~ Class: Planet ~
**
** Class for planetary objects (or sun) in a solar system
** Each object in a system should be in an array containing the whole system
*/

#include "planet.h"

// Constructor
Planet::Planet(char const* sLabel, int iRunID, double dAp, double dVMin, double dTh, double dMass, double dSize){
	double dXPos, dYPos, dVX, dVY;
	
	dTh       *= 2*cPI/360.0;
	
	dXPos      =  cos(dTh)*dAp;
	dYPos      =  sin(dTh)*dAp;
	dVX        = -sin(dTh)*dVMin;
	dVY        =  cos(dTh)*dVMin;
	
	Time       = 0;
	XPosAU     = dXPos/cAU;
	YPosAU     = dYPos/cAU;
	XPos       = dXPos;
	YPos       = dYPos;
	VX         = dVX;
	VY         = dVY;
	MCX        = 0.0;
	MCY        = 0.0;
	Mass       = dMass;
	Size       = dSize;
	Label      = sLabel;
	RunID      = iRunID;
	Aphelion   = dAp;
	Perihelion = dAp;
	VMin       = dVMin;
	VMin       = dVMin;
	SectorArea = 0.0;
	
	FileOpen   = false;
}

double Planet::Radius() {
	return sqrt(XPos*XPos + YPos*YPos);
}

double Planet::Velocity() {
	return sqrt(VX*VX + VY*VY);
}
		
void Planet::Update(Planet** aSystem, int iObject, int iObjects) {
	int    i;
	double dR, dU=0;
	double cCMx, cCMy, dTMass;
	double dRadius, dVelocity;
	
	cCMx = 0;
	cCMy = 0;
	dTMass = 0;
	for(i=0; i<iObjects; i++) {
		cCMx += aSystem[i]->XPos*aSystem[i]->Mass;
		cCMy += aSystem[i]->YPos*aSystem[i]->Mass;
		dTMass += aSystem[i]->Mass;
	}
	MCX = cCMx/dTMass;
	MCY = cCMy/dTMass;

	FX = 0.0;
	FY = 0.0;
	for(i=0; i<iObjects; i++) {
		if(i!=iObject) {
			dR  = sqrt(pow(XPos-aSystem[i]->XPos,2) + pow(YPos-aSystem[i]->YPos,2));
			FX -= aSystem[i]->Mass * Mass * (XPos - aSystem[i]->XPos) / pow(dR,3);
			FY -= aSystem[i]->Mass * Mass * (YPos - aSystem[i]->YPos) / pow(dR,3);
			dU -= sqrt(FX*FX + FY*FY)*dR;
		}
	}
	FX *= cG;
	FY *= cG;

	dRadius   = Radius();
	dVelocity = Velocity();
	Energy    = 0.5*Mass*dVelocity*dVelocity + cG*dU;
	Spin      = Mass*dVelocity * dRadius;
	
	if(dRadius < Perihelion) Perihelion = dRadius;
	if(dRadius > Aphelion)   Aphelion   = dRadius;
	if(dVelocity < VMin)     VMin = dVelocity;
	if(dVelocity > VMax)     VMax = dVelocity;

	XPosAU = XPos/cAU;
	YPosAU = YPos/cAU;
}
		
void Planet::Correct(Planet** aSystem, int iObjects, int iSun) {
	int    i;
	double cCMx=0.0, cCMy=0.0, dTMass=0.0;
	
	for(i=0; i<iObjects; i++) {
		cCMx   += aSystem[i]->XPos*aSystem[i]->Mass;
		cCMy   += aSystem[i]->YPos*aSystem[i]->Mass;
		dTMass += aSystem[i]->Mass;
	}
	cCMx = cCMx/dTMass;
	cCMy = cCMy/dTMass;
	
	for(i=0; i<iObjects; i++) {
		aSystem[i]->XPos -= cCMx;
		aSystem[i]->YPos -= cCMy;
	}
}

void Planet::InitSystem(Planet** aSystem, int iObjects, int iSun) {
	int    i;
	double dPx=0.0, dPy=0.0;
	
	Correct(aSystem, iObjects, iSun);

	for(i=0; i<iObjects; i++) {
		if(i!=iSun) {
			dPx += aSystem[i]->Mass*aSystem[i]->VX;
			dPy += aSystem[i]->Mass*aSystem[i]->VY;
		}
	}
	aSystem[iSun]->VX = -dPx/aSystem[iSun]->Mass;
	aSystem[iSun]->VY = -dPy/aSystem[iSun]->Mass;
}

void Planet::Sector(double dXRef, double dYRef) {
	double dR    = Radius();
	double dRRef = sqrt(dXRef*dXRef + dYRef*dYRef);
	SectorArea = dR*dR*acos((dXRef*XPos + dYRef*YPos)/(dR*dRRef))/2.0;
	return;
}

double Planet::Gamma(Planet* oStar) {
	return sqrt(1-2*oStar->Mass*cG/(cC*cC*sqrt(pow(XPos-oStar->XPos,2) + pow(YPos-oStar->YPos,2))));
}
		
void Planet::Output() {

	int iPrecision=15;
	
	if(!FileOpen) {
		ostringstream sFileName;
		sFileName << "/scratch/frank/solarsystem/" << RunID << "_" << Label << ".dat"; 
		oFile.open(sFileName.str().c_str());
		FileOpen = true;
	}
	
	oFile << setprecision(iPrecision);
	oFile << setw(iPrecision+7) << Time;         // 0
	oFile << setw(iPrecision+7) << XPosAU;       // 1
	oFile << setw(iPrecision+7) << YPosAU;       // 2
	oFile << setw(iPrecision+7) << VX;           // 3
	oFile << setw(iPrecision+7) << VY;           // 4
	oFile << setw(iPrecision+7) << Energy;       // 5
	oFile << setw(iPrecision+7) << Spin;         // 6
	oFile << setw(iPrecision+7) << Radius()/cAU; // 7
	oFile << setw(iPrecision+7) << MCX;          // 8
	oFile << setw(iPrecision+7) << MCY;          // 9
	oFile << setw(iPrecision+7) << SectorArea;   // 10
	oFile << endl;
	
	return;
}

void Planet::Summary() {
	cout << setw(10) << Label << ":";
	cout << setprecision(5) << setw(12) << Aphelion/1e9 << " Mkm";
	cout << setprecision(5) << setw(12) << Perihelion/1e9 << " Mkm";
	cout << setprecision(5) << setw(12) << VMin/1000 << " km/s";
	cout << setprecision(5) << setw(12) << VMax/1000 << " km/s";
	cout << endl;
}
		
void Planet::Close() {
	oFile.close();
}

