#include <iostream>
#include <cmath>
#include <armadillo>
#include <iomanip>
#include <ctime>
#include <stdio.h>

using namespace std;
using namespace arma;

void Rotation(mat &A, mat &S, int l, int k, int n);
double maxoffdiag(mat &A, int &l, int &k, int n);
double JRotation(mat &A, mat &S, int &l, int &k, int n);


int main(int argc, char const *argv[]){
	double rho_max;
	double rho_min = 0.0;
	int n, k, l;
	
	cout << "Size of the matrix:" << endl;
	cin >> n;
	cout << "Give rho_max:" << endl;
	cin >> rho_max;
	double h = (rho_max - rho_min)/(n+1);
	mat A = zeros(n,n); 
	mat rho = zeros	(n,1); 
	mat S = zeros(n,n);

	for (int i = 0; i < n-1; i++){
		rho(i) = rho_min + (i+1)*h;
		A(i,i) = (2.0/(h*h))+(rho(i)*rho(i));
		A(i+1,i) = -1.0/(h*h);
		A(i,i+1) = -1.0/(h*h);
	}
	A(n-1,n-1) = (2.0/(h*h))+(n*h*n*h);
	//---------------------------------------------------------------
	double epsilon = JRotation(A,S,l,k,n);

	clock_t start, finish;
	
	start = clock(); 
	cout << eig_sym(A) << endl;	//Here we calculate the eigenvalues using armadillo's function. 
	finish = clock();
	cout << (finish-start)/ (double) 1000000.0 << endl; //CLOCKS_PER_SEC << endl;
	//---------------------------------------------------------------
	
	//cout << A << endl; //printing out the A before the transformation
	
	//Using armadillo's function to print out the maximum values in the matrix which correspond to the eigenvalues.
	
	start = clock();
	cout << A << endl;
	finish = clock();
	cout << (finish-start) / (double) CLOCKS_PER_SEC<< endl;
	

	return 0;
}

void Rotation(mat &A, mat &S, int l, int k, int n){
	double tau = (A(l,l)-A(k,k))/(2.0*A(k,l)); //tau = cos(2theta)
	double t;
	
	if (tau >= 0.0){
		t = (-tau + sqrt(tau*tau+1)); //t = tan(theta)
	}
	else{
		t = (-tau - sqrt(tau*tau+1)); //t = tan(theta)
	}
	

	//double t = sign(tau)/(fabs(tau) + sqrt(tau*tau+1)); //t = tan(theta)
	double c = 1.0/sqrt(1+t*t); //c = cos(theta)
	double s = t*c; //sin
	A(l,l) += t*A(l,k);
	A(k,k) -= t*A(l,k);
	A(k,l) = A(l,k) = 0.0;
	for (int i = 0; i < n; i++){
		if (i!=k && i!=l){
			double Ail = A(i,l)*c + A(i,k)*s;
			double Aik = -s*A(i,l) + A(i,k)*c;
			A(i,l) = A(l,i) = Ail;
			A(i,k) = A(k,i) = Aik;
		}

	}

	for (int i = 0; i < n; i++){
		double Sil = S(i,l)*c + S(i,k)*s;
		double Sik = -s*S(i,l) + S(i,k)*c;
		S(i,l) = Sil;
		S(i,k) = Sik;
	}

}

double maxoffdiag(mat &A, int &l, int &k, int n){
	double max = 0.0;
	for (int i = 0; i < n-1; i++){
		for (int j = i+1; j < n; j++){
			double Aij = fabs(A(i,j));
			if (Aij > max){
				max = Aij;
				l = i;
				k = j;
			}
		}
	}
	return max;
}


double JRotation(mat &A, mat &S, int &l, int &k, int n){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			if (i==j){
				S(i,j) = 1.0;
			}
			else{
				S(i,j) = 0.0;
			}
		}
	}
	//int k,l;
	double epsilon = 1.0e-10;
	double max_nr_it = (double) n * (double) n * (double) n;
	int iterations = 0;
	double max = maxoffdiag(A,k,l,n);

	while (fabs(max) > epsilon && (double) iterations < max_nr_it){
		max = maxoffdiag(A, k,l, n);
		Rotation(A, S, l, k, n);
		iterations++;
	}
	cout << "Number of iterations: " << iterations << "\n";
	return epsilon;
}