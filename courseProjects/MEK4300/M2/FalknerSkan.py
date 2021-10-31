"""
Solve stagnation equation

    f * f'' + f''' + 1 - f'*f' = 0

on 0 < x < L, with boundary conditions 

f(0) = 0, f'(0) = 0 and f'(L) = 1

by solving a coupled system

    f * y' + y'' + 1 - y*y = 0
    y - f' = 0

with f(0) = 0, y(0) = 0 and y(L) = 1

Variational form
    (vy, f*y')*dx + (vy, y')*ds - (vy', y')*dx + (vy, 1-y*y)*dx
  + (vf, f'-y)*dx = 0
  
where vy and vf are test functions for y and f respectively. The
exterior boundary integral disappears because of the Dirichlet
boundary conditions, though, so the final form is

    (vy, f*y')*dx - (vy', y')*dx + (vy, 1-y*y)*dx + (vf, f'-y)*dx = 0

"""
from dolfin import *
import sys
import matplotlib.pyplot as plt

# assert len(sys.argv) == 2
# solver = sys.argv[-1]
# assert solver in ("Newton", "Picard")

L = 10
x = IntervalMesh(50, 0, L)
V = FunctionSpace(x, 'CG', 1)
VV = V * V
yf = TrialFunction(VV)
y, f = split(yf)
vy, vf = TestFunctions(VV)

# Initialize solution
yf_ = interpolate(Expression(("1", "x[0]"), element=VV.ufl_element()), VV)
y_, f_ = split(yf_)




bcs = [DirichletBC(VV, Constant((0, 0)), "x[0] < DOLFIN_EPS"),
       DirichletBC(VV.sub(0), Constant(1), "x[0] > {} - 10*DOLFIN_EPS".format(L))]

# for beta in [-0.198838,-0.18,-0.1,0,0.3,1.0]:
    # Define variational form 
beta = 1.0
FN = inner(vy*f_, y_.dx(0))*dx - inner(vy.dx(0), y_.dx(0))*dx + Constant(beta) * inner(vy, 1 - y_*y_)*dx \
  + inner(vf, (f_.dx(0) - y_))*dx
# Newton
solve(FN == 0, yf_, bcs=bcs, solver_parameters={'newton_solver':{'maximum_iterations': 20}})
d2f = project(y_.dx(0), V)
print d2f(0)
plot(d2f, key = "d2f", title='Velocity profile')
plot(y_)
interactive()
# y_p = Function(V)
# dy_p = Function(V)
# assign(y_p, yf_.sub(0))
# assign(dy_p, d2f)#project(y_.dx(0),V))
# plt.figure(1)
# plt.plot(y_p.vector().array())
# print y_p
# plt.figure(2)
# plt.plot(dy_p.vector().array())   
# plt.show()
"""
solver = "Newton"
if solver == "Newton":
    for beta in [-0.198838,-0.18,-0.1,0,0.3,1.0]:
        # Define variational form 
        FN = inner(vy*f_, y_.dx(0))*dx - inner(vy.dx(0), y_.dx(0))*dx + Constant(beta) * inner(vy, 1 - y_*y_)*dx \
          + inner(vf, (f_.dx(0) - y_))*dx
        # Newton
        solve(FN == 0, yf_, bcs=bcs, solver_parameters={'newton_solver':{'maximum_iterations': 20}})
        d2f = project(y_.dx(0), V)
        print d2f(0)
        print beta
        plot(d2f, key = "d2f", title='Velocity profile',interactive=False)
        plot(y_)
        # y_p = Function(V)
        # dy_p = Function(V)
        # assign(y_p, yf_.sub(0))
        # assign(dy_p, d2f)#project(y_.dx(0),V))
        # plt.figure(1)
        # plt.plot(y_p.vector().array())

        # plt.figure(2)
        # plt.plot(dy_p.vector().array())   
        # plt.show()
# Picard
else:
    for beta in [-0.198838,-0.18,-0.1,0,0.3,1.0]:

        FP = inner(vy*f_, y.dx(0))*dx - inner(vy.dx(0), y.dx(0))*dx + Constant(beta) * inner(vy, 1 - y*y_)*dx \
        + inner(vf, (f.dx(0) - y))*dx
      
        error = 1
        k = 0
        yf_1 = Function(VV)
        while k < 100 and error > 1e-10:
            A = assemble(lhs(FP))
            b = assemble(rhs(FP))
            [bc.apply(A, b) for bc in bcs]
            solve(A, yf_.vector(), b)
            yf_1.vector().axpy(-1, yf_.vector())
            error = norm(yf_1.vector())
            yf_1.assign(yf_)
            print "Error = ", k, error
            k += 1

        plot(y_, key = "y_", title='Velocity profile',interactive=False)


"""
# Plot solution