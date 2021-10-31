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


# gmsh created mesh and boundary facet
mesh = Mesh("out.xml")

bndr = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

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

# class Cylinder(SubDomain):
#     """docstring for Cylinder"""
#     def inside(self, x, on_boundary):
#         return on_boundary and (x[0] < 0.14 or x[0] > 0.26) and (x[1] < 0.14 or x[1] > 0.26) 


bndry = FacetFunction("size_t", mesh)
for f in facets(mesh):
    mp = f.midpoint()
    if near(mp[0], 0.0): # inflow
        bndry[f] = 1
    elif near(mp[0], length): # outflow
        bndry[f] = 2
    elif near(mp[1], 0.0) or near(mp[1], H): # walls
        bndry[f] = 3
    elif mp.distance(x_c) <= radius: # cylinder
        bndry[f] = 5

bndr.set_all(0)

# cylinder = Cylinder()
# cylinder.mark(bndr,1)
ds = Measure("ds")[bndry]


# Defining the function spaces
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q # MixedFunctionSpace

# Introducing the Dirichlet boundary conditions
# No-slip:
noslip_bnd = NoSlip()
noslip = Constant((0.,0.))
bc0 = DirichletBC(W.sub(0), noslip, bndry, 3)#noslip_bnd)
bc4 = DirichletBC(W.sub(0), noslip, bndry, 5)#noslip_bnd)
# In-flow:
inflow_bnd = InFlow()
inflow = Expression(("1.5 * 4.0 * x[1] * (0.41 - x[1]) / ( 0.41 * 0.41 )", "0.0"))
bc1 = DirichletBC(W.sub(0), inflow, bndry, 1)#inflow_bnd)
# Out-flow:
outflow_bnd = OutFlow()
outflow = Constant(0.)
# w.sub(1) is the function space i.e pressure
bc2 = DirichletBC(W.sub(1), outflow, outflow_bnd)
# now the velocity should have the profile as the inflow 
bc3 = DirichletBC(W.sub(0), inflow, outflow_bnd)
# Collecting the bnd cond. in a list
bcs = [bc0, bc4, bc1]#, bc3]

# Introducing the weak formulations for the Stokes 
# and the Navier-Stokes flow

#(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
w = Function(W)

"""
# Stokes
f = Constant((0.,0.))
a = (inner(nu*grad(u), grad(v)) - div(v)*p + q*div(u))*dx
L = inner(f, v)*dx

# Navier-Stokes skew
L_ns = inner(u0, v)*dx
a_ns = (inner(u,v) + dt*(.5*inner(grad(u)*u0,v) - .5*inner(grad(v)*u0,u)\
+ nu*inner(grad(u),grad(v)) - div(v)*p) + q*div(u))*dx

solve(a == L, w, bcs)
"""

u, p = split(w)

F =   inner(grad(u)*u, v)*dx + nu*inner(grad(u), grad(v))*dx - p*div(v)*dx - q*div(u)*dx


n = FacetNormal(mesh)
u_ = assemble(dot(u, n)*ds(1))
I = Identity(2)

Force_D = (nu*dot(grad(u),n)-p*dot(I,n))[0]*ds(1)

Force_L = (nu*dot(grad(u),n)-p*dot(I,n))[1]*ds(1)

c_D = 500*Force_D 
c_L = 500*Force_L
CD = assemble(c_D)
CL = assemble(c_L)
info("Drag = %.5f     Lift = %.5f" %(CD, CL))
#(u, p) = w.split(True)
u, p = split(w)

p = w.sub(1, deepcopy=False)
p_diff = p(Point(0.15,0.20)) - p(Point(0.25,0.20))
info("p_diff = %.6f" %p_diff)
# plot(u, "Initial velocity")
# plot(p, "Initial pressure")
# interactive()


# ufile = File("results_NavierStokes/u.pvd")
# pfile = File("results_NavierStokes/p.pvd")
# the time-periodic Navier-Stokes 
t = 0
step = 0
while t < T:
    #u0.assign(u)
    t += dt
    step += 1
    solve(F == 0, w, bcs)

    (u, p) = w.split(True)
    u_plot.assign(u)
    # u.rename("u", "velocity")
    # p.rename("p", "pressure")
    # ufile << u
    # pfile << p

