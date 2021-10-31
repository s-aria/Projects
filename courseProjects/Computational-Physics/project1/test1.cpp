 
#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[]){
	int i, j, n, k;
	double h;
    clock_t start, finish

    start = clock()

	n = pow(10,4);
	
	mat B = randu<mat>(n,n); mat C = randu<mat>(n,n);

	for (i = 0; i < n; ++i){
		for (j = 0; j < n; ++j){
			for (i = 0; k < n; ++k){
				A[i][j] += B[i][k]*C[k][j];
			}
		}
	}
    finish = clock();
    ((finis - start)/CLOCKS_PER_SEC);
}
