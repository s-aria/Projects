//Project 1

#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;

ofstream myfile;

int main(int argc, char* argv[])
{	
	char *outfilename;
	int	i, n;
	double h;
	
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
	
	mat a(n,1); mat b(n+1,1); mat c(n,1); mat v(n+1,1); mat f(n+1,1); 
	a.fill(-1); b.fill(2); c.fill(-1);

	//Computing the function values, analytical
	for(i=0;i<n+1;i++){
		f(i) = h*h*100*exp(-10*i*h);
	}

	//The Forward subtitution
	for(i=0;i<n;i++){
		b(i+1) = b(i+1) - c(i)*(a(i)/b(i));
		f(i+1) = f(i+1) - f(i)*(a(i)/b(i));
	}


	cout << "b" << b << fixed;	
	cout << "f" << f << fixed;
	cout << "it works"<< endl;
	//The Backward substitution
	for (i = n; n > 0; i--){
		f(i-1) = (f(i-1)-(c(i-1)*f(i))/b(i));
	}
	
	for (int i = 0; i < n+1; i++)
	{
		v(i) = f(i)/b(i);
	}
	cout << "v" << v << fixed;	


	for(i = 1; i <= n; i++){
		myfile << setiosflags(ios::showpoint | ios::uppercase);
		myfile << setw(15) << setprecision(8) << i*h;
		myfile << setw(15) << setprecision(8) << f(i);
		myfile << setw(15) << setprecision(8) << 1.0-(1.0-exp(-10))*i*h-exp(-10*i*h) <<endl;
	}

	myfile.close();
	return 0;
}