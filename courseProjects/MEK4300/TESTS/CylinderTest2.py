from dolfin import *
import mshr
"""
H=0.41
 ________________________________________________
|			U=V=0								 |
|                                                |
|	O	<---U=V=0                                |
|________________________________________________|

x=0         U=V=0                               L = 2.2

"""
# Parametres for our calculation
rho = 1.0 # kgm^-3
nu = Constant(0.001) # m^2s^-1
H = 0.41 # height of the rectangle
length = 2.2  # length of the rectangle
radius = 0.05 # radius of the circle inside the rectangle
x_c = y_c = Point(0.2,0.2) # Centre point of the circle
dt = 0.0125 # time step
T = 10.0  # end time



# Sub domain for no-slip (mark whole boundary, inflow and outflow will overwrite)
class NoSlip(SubDomain):
    def inside(self, x, on_boundary):
    	dx = x[0] - x_c[0]
    	dy = x[1] - y_c[1]
    	r = sqrt(dx*dx + dy*dy)
        return on_boundary and (x[1] < 0.0 + DOLFIN_EPS or x[1] > (H - DOLFIN_EPS) or r < (radius + DOLFIN_EPS))

# Sub domain for inflow (left)
class InFlow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < 0.0 + DOLFIN_EPS and on_boundary

# Sub domain for outflow (left)
class OutFlow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > (length - DOLFIN_EPS) and on_boundary

class Cylinder(SubDomain):
    """docstring for Cylinder"""
    def inside(self, x, on_boundary):
        return on_boundary and (x[0] < 0.14 or x[0] > 0.26) and (x[1] < 0.14 or x[1] > 0.26) 

# gmsh created mesh and boundary facet
mesh = Mesh("out.xml")

bndr = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

bndr.set_all(10)

# Cylinder
cylinder = Cylinder()
cylinder.mark(bndr,0)

inflow = InFlow()
inflow.mark(bndr,1)

outflow = OutFlow()
outflow.mark(bndr,2)

noslip = NoSlip()
noslip.mark(bndr,3)

File("subdomains.xml") << bndr

bndr = MeshFunction("size_t", mesh, "subdomains.xml")
ds = Measure("ds")[bndr]


# Defining the function spaces
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q # MixedFunctionSpace

(v, q) = TestFunctions(W)
w = Function(W)
u, p = (as_vector((w[0], w[1])), w[2])


inv=Expression(("0.3 * 4.0 * x[1] * (0.41 - x[1]) / ( 0.41 * 0.41 )", "0.0"))
bc0 = DirichletBC(W.sub(0), Constant((0,0)), bndr,0)
bc1 = DirichletBC(W.sub(0), inv, bndr, 1)
bc2 = DirichletBC(W.sub(0), Constant((0,0)), bndr,3)
bcs = [bc0, bc1, bc2]



# Introducing the weak formulations for the Stokes 
# and the Navier-Stokes flow

#(u, p) = TrialFunctions(W)

F =   inner(grad(u)*u, v)*dx + nu*inner(grad(u), grad(v))*dx - p*div(v)*dx - q*div(u)*dx

# Derivative of weak form
dw = TrialFunction(W)
dF = derivative(F, w, dw)

problem = NonlinearVariationalProblem(F, w, bcs, dF)
solver  = NonlinearVariationalSolver(problem)
# Set linear solver parameters
itsolver = solver.parameters["newton_solver"]
itsolver["absolute_tolerance"] = 1.0e-10
itsolver["relative_tolerance"] = 1.0e-10

# To see various solver options, uncomment following line
#info(solver.parameters, True); quit()

# If you want to initialize solution from previous computation
restart = "no"
if restart == "yes":
   print "Setting initial condition from file ..."
   File("steady.xml") >> w.vector()

# Solve the problem
solver.solve()

# Save steady solution
File("steady.xml") << w.vector()

# Save vtk for visualization
(u,p) = w.split()
File("velocity.pvd") << u
File("pressure.pvd") << p

# Stress tensor
T = nu*(grad(u) + grad(u).T) - p*Identity(2)
# Face normals
n = FacetNormal(mesh)

# Compute force on cylinder
drag = -T[0,j]*n[j]*ds(1)
lift = -T[1,j]*n[j]*ds(1)
print "Drag =", assemble(drag)
print "Lift =", assemble(lift)