from dolfin import *
import numpy
#from mshr import *

mu = Constant(1.0)

mesh = Mesh("step.xml")


V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)
VQ = MixedFunctionSpace([V,Q])
u, p = TrialFunctions(VQ)
v, q = TestFunctions(VQ)
F = mu*inner(grad(v), grad(u))*dx - inner(div(v), p)*dx - inner(q, div(u))*dx +Constant(0)*q*dx
s0 = AutoSubDomain(lambda x, on_bnd: x[1]>0.5-1e-8)
Top = DirichletBC(VQ.sub(0), (1,0),s0)
s1 = AutoSubDomain(lambda x, on_bnd: x[1]<1e-8)
Bottom1 = DirichletBC(VQ.sub(0), (0,0), s1)

def b_boundary(x):
	return x[1]< 0.1+DOLFIN_EPS and x[0]<0.5+DOLFIN_EPS

Bottom2 = DirichletBC(VQ.sub(0), (0,0), b_boundary)

up_ = Function(VQ)
bcs = [Top,Bottom1,Bottom2]
solve(lhs(F)==rhs(F), up_, bcs=bcs)
u_, p_ = up_.split()
#plot(u_,interactive=True)
f = File("solution.pvd")
f << u_
#stream func
q = TestFunction(V)
psi = TrialFunction(V)
n = FacetNormal(mesh)
F = inner(grad(q), grad(psi))*dx \
- inner(q, (u[1].dx(0) - u[0].dx(1)))*dx \
+ q*(n[1]*u[0] - n[0]*u[1])*ds +Constant(0)*q*dx
# plot(U)


