from fe_approx1D_numint import *


x = sp.Symbol('x')

f = sp.tanh(40*(x-0.5))

for i in 1,2:
	for N_e in 4/i, 8/i, 16/i:
		approximate(f,d=i,N_e=N_e, numint = "Simpson", Omega=[0,1])
	
