//Project 1

#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <vector>
 

using namespace std;
using namespace arma;
using namespace ios;

ofstream myfile;

int main(int argc, char* argv[])
{      
    char *outfilename;
    int     i, n, j;
    double h;
    
    clock_t start, finish;
    start = clock();

    //Read in output file, abort if there are too few command-line arguments
    if( argc <= 1 ){
            cout << "Bad Usage: " << argv[0] << " read also output file on same line"<< endl;
            exit(1);
    }

    else{
            outfilename=argv[1];
    }
   
    myfile.open(outfilename);
    cout <<"Read in number of mesh points"<< endl;
    cin >> n;
   

    h = 1/( (double) n+1);
   
    mat a(n,1); mat b(n+1,1); mat c(n,1); mat v(n+1,1); mat f(n+1,1); //mat A(n,n);  //remove d(i) it's not needed
    a.fill(-1); b.fill(2); c.fill(-1);

  
    for(i=0;i<n+1;i++){
            f(i) = h*h*100*exp(-10*i*h);
    }
    //Update b(i) and f(i) values after we subtract and remove a(i)
    for(i=0;i<n;i++){
            b(i+1) = b(i+1) - c(i)*(a(i)/b(i));
            f(i+1) = f(i+1) - f(i)*a(i)/b(i);
           
    }
   
    //The backwards substitution starts from second row.
    //removing the c(i) values as we head back.
    for (i = n; i > 0; i--){
            f(i-1) = f(i-1)-(c(i-1)*f(i)/b(i));
    }
    
    //only one major diagonal left and it includes constant b(i) times x(i)
    for (i = 0; i < n+1; i++){
            v(i) = f(i)/b(i);
           
    }

    //EXERCISE D--------------------------
    
    vector<vector<int> > A;
     
    // Initialize the A as a n x n array of 0.
    A = vector<vector<int> >(n, vector<int>(n,0));
     
     
    for(unsigned int i = 0; i < n; i++){
        for(unsigned int j = 0; j < n; j++){
            if(i==j)
                A[i][j] = 2;
            if(i==j && i>0)
                A[i][j-1] = -1;
            if(i==j && i<n-1)
                A[i][j+1] = -1;
     
        }   
     
    }
     
    // Print it
    for(unsigned int y = 0; y < n; y++)
    {
        for(unsigned int x = 0; x < n; x++){
            cout << "\t" << A[y][x];
        cout << "\n";
        }
    }
    
    
    mat L, U, P;

    lu(L, U, P, A);

    mat B = P.t()*L*U;
    
    //END EX D--------------------------------------

    
    finish=clock();
    ((finish-start)/CLOCKS_PER_SEC);
    

    for(i = 1; i <= n; i++){
            myfile << setiosflags(ios::showpoint | ios::uppercase);
            myfile << setw(15) << setprecision(8) << i*h;
            myfile << setw(15) << setprecision(8) << v(i); //f(i)->v(i)
            myfile << setw(15) << setprecision(8) << 1.0-(1.0-exp(-10))*i*h-exp(-10*i*h) <<endl;
    }

    myfile.close();
    return 0;
}