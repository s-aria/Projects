from dolfin import *

mesh = UnitIntervalMesh(10)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

bc0 = DirichletBC(V, Constant(0), "std::abs(x[0]) < 1e-10")
bc1 = DirichletBC(V, Constant(1), "std::abs(x[0]-1) < 1e-10")

def q(u):
    return 1+u**4
    
u_ = interpolate(Expression("x[0]"), V)
F = inner(grad(v), q(u_)*grad(u_))*dx
##solve(F == 0, u_, [bc0, bc1])

J = derivative(F, u_, u)
du = Function(V)
error = 1; k = 0
while k < 100 and error > 1e-12:
    A = assemble(J)
    b = assemble(-F)
    [bc.apply(A, b, u_.vector()) for bc in [bc0, bc1]]
    solve(A, du.vector(), b)
    u_.vector().axpy(1., du.vector())
    error = norm(du.vector())
    print "Error = ", error
    k += 1    
    