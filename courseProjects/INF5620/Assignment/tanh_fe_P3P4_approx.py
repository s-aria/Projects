from fe_approx1D_numint import *


x = sp.Symbol('x')

f = sp.tanh(60*(x-0.5))

for i in 3,4:
	for N_e in 1, 2, 4:
		approximate(f,d=i,N_e=N_e, Omega=[0,1])
	
