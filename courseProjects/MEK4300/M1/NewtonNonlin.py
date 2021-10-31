from dolfin import *

L = 4
mesh = IntervalMesh(50, 0, L)
V = FunctionSpace(mesh, 'CG', 1)
VV = V*V
fh = TrialFunction(VV)
f, h = split(fh)
vf, vh = TestFunctions(VV)
F = vf + vh
solve(F==0)