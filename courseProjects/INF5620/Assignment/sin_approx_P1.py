from fe_approx1D import *
import sympy as sp
import numpy as np


phi = basis(d=1)
h, x = sp.symbols('h, x')
f = sp.sin(x)
nodes, elements = mesh_uniform(N_e=2,d=1,Omega=[0,np.pi])


A, b = assemble(nodes,elements,phi,f,symbolic = True)


approximate(f=sp.sin(x),symbolic=False, d=1,N_e=2)
print np.array(A)
print np.array(b)