#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;

vec Forward_Euler(int n, double dx, double dt, int m, double alpha);

ofstream myfile;

int main(int argc, char const *argv[]){
	int n; double dx; double dt; int m; double alpha;
	Forward_Euler(n, dx, dt, m, alpha);	
	return 0;
}

//------------- Explicit Scheme ------------
vec Forward_Euler(int n, double dx, double dt, int m, double alpha){
	n = 100;
	dx = 1.0/(n-1);
	vec u = zeros(n);
	vec v = zeros(n);
	u(0) = 0.0;
	u(n-1) = 0.0;
	v(0) = u(0);
	v(n-1) = u(n-1);
	vec x = zeros(n);
	for (int i = 1; i < n-1; i++){
		x(i) = i*dx;
		u(i) = x(i) - 1;

	}
	myfile.open("data",ios::out);
	m = 200; //time grid points
	alpha = 0.49;//dt/(dx*dx);
	for (int t = 1; t <= m; t++){
		for (int i = 1; i < n-1; i++){
			v(i) = alpha*u(i-1) + (1 - 2*alpha)*u(i) + alpha*u(i+1);
			myfile << v(i) << " "; 
			
		}
		myfile << endl;
		u = v;
	}

	
	myfile.close();
		
	return  u;
} 
//---------------Explicit Scheme END---------------//

//-----------Implicit Scheme--------------------

vec TriSolver(a, b, c){
   
    //Update b(i) and f(i) values after we subtract and remove a(i)
    for(i=0;i<n;i++){
            b(i+1) = b(i+1) - c(i)*(a(i)/b(i));          
    }
   
    //only one major diagonal left and it includes constant b(i) times x(i)
    for (i = 0; i < n+1; i++){
            v(i) = f(i)/b(i);
           
    }

}


void Backward_Euler(int n , int m, double dx, double alpha){
	mat a(n,1); mat b(n+1,1); mat c(n,1); 
	vec u = zeros(n);
	vec y = zeros(n);
	for (int i = 1; i < n; i++){
		y(i) = u(i) =  dx*i - 1;
	}
	
	y(n-1) = u(n-1);
	y(0) = u(0);

    a.fill(-alpha); b.fill(1 + 2*alpha); c.fill(-alpha);
    alpha = 0.5;
    m = 20000;
    for (int t = 1; t <= m; t++){
    	TriSolver(a,b,c,n)
    
		u(n-1) = 0.0;
		u(0) = 0.0;
		for (int i = 0; i <= n-1; i++){
			y(i) = u(i);
		}
	}
	
}

//-----------Implicit Scheme END-----------//

//--------------Crank-Nicolson------------------
void Crank_Nicolson(int n, int m, double dx, doble alpha){
	mat a(n,1); mat b(n+1,1); mat c(n,1); 
	vec u = zeros(n);
	vec y = zeros(n);
	alpha = 0.5;
	a.fill(-alpha); b.fill(2 + 2*alpha); c.fill(-alpha);
	m = 20000;
	for (int t = 1; t <= m; t++){
		for (int i = 1; i < n-1; i++){
			y(i) = alpha*u(i-1) + (2-2*alpha) *u(i) + alpha*u(i+1);

		}
	y(0) = 0.0;
	y(n-1) = 0.0;
	TriSolver(a,b,c,);
	u(0) = 0;
	u(n-1) = 0;
	}
}
//--------------Crank-Nicolson END-----------