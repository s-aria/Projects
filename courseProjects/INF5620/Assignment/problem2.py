from fe_approx1D import *
N = 1
d = 2
phi = basis(d)
x = sp.Symbol('x')
f = x*(1-x)
Omega = [0,1]
nodes, elements = mesh_symbolic(N,d,Omega)
A, b = assemble(nodes, elements, phi, f, symbolic=True)
c = A.LUsolve(b)
print sp.latex(A,mode='plain')
print '-------------------------------------------------------------'
print sp.latex(b,mode='plain')
print '--------------------------------------------------------------'
print sp.latex(c,mode ='plain')