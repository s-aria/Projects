from dolfin import *

set_log_level(1)

# Sub domain for no-slip (mark whole boundary, inflow and outflow will overwrite)
class Noslip(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Sub domain for inflow (right)
class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1.0 - DOLFIN_EPS and on_boundary

# Sub domain for outflow (left)
class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary

# Read mesh
mesh = Mesh("out.xml")

# Create mesh functions over the cell facets
sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
# Mark all facets as sub domain 3
sub_domains.set_all(3)
# Mark no-slip facets as sub domain 0, 0.0
noslip = Noslip()
noslip.mark(sub_domains, 0)
# Mark inflow as sub domain 1, 01
inflow = Inflow()
inflow.mark(sub_domains, 1)
# Mark outflow as sub domain 2, 0.2, True
outflow = Outflow()
outflow.mark(sub_domains, 2)
# Save sub domains to file
file = File("subdomains.xml")
file << sub_domains

# Save sub domains to VTK files
file = File("subdomains.pvd")
file << sub_domains

mesh = Mesh("out.xml")
sub_domains = MeshFunction("size_t", mesh, "subdomains.xml")

# Define function spaces
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# No-slip boundary condition for velocity 
# x1 = 0, x1 = 1 and around the dolphin
noslip = Constant((0, 0))
bc0 = DirichletBC(W.sub(0), noslip, sub_domains, 0)

# Inflow boundary condition for velocity
# x0 = 1
inflow = Expression(("0.3*23.7953599048*x[1]*(0.41-x[1])", "0.0"))
bc1 = DirichletBC(W.sub(0), inflow, sub_domains, 1)

# Boundary condition for pressure at outflow
# x0 = 0
zero = Constant(0)
bc2 = DirichletBC(W.sub(1), zero, sub_domains, 2)

# Collect boundary conditions
bcs = [bc0, bc1, bc2]

# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
f = Constant((0, 0))
a = (inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
L = inner(f, v)*dx

# Compute solution
w = Function(W)
solve(a == L, w, bcs)

# Split the mixed solution using deepcopy
# (needed for further computation on coefficient vector)
(u, p) = w.split(True)

print("Norm of velocity coefficient vector: %.15g" % u.vector().norm("l2"))
print("Norm of pressure coefficient vector: %.15g" % p.vector().norm("l2"))

# # Split the mixed solution using a shallow copy
(u, p) = w.split()

# Save solution in VTK format
ufile_pvd = File("velocity.pvd")
ufile_pvd << u
pfile_pvd = File("pressure.pvd")
pfile_pvd << p

# Plot solution
plot(u)
plot(p)
interactive()








# from streamfunction import StreamFunction
# psi = StreamFunction(u_, [], mesh)

# if utop < 0:
#     psi.vector()._scale(-1)
# psi.vector()[:] -= psi.vector().min()

# f = File('u.pvd')
# f << u_

# f = File('psi.pvd')
# f << psi
# plot(psi)

# bc = DirichletBC(Q, 10, "on_boundary")
# bc.apply(psi.vector())
# pa = psi.vector().array().argmin()
# #pa = psi.vector().array().argmax()

# xx = interpolate(Expression("x[0]"), psi.function_space())
# yy = interpolate(Expression("x[1]"), psi.function_space())
# xm = xx.vector()[pa]
# ym = yy.vector()[pa]
# print "Min streamfunction = ", psi(xm, ym), " at (", xm, ym, ")" 

# mesh2 = RectangleMesh(xm[0]-0.02, ym[0]-0.02, xm[0]+0.02, ym[0]+0.02, 100, 100)
# Vv = FunctionSpace(mesh2, 'CG', 2)
# psim = project(psi, Vv)

# print "Refined:"
# pa = psim.vector().array().argmin()
# xx = interpolate(Expression("x[0]"), Vv)
# yy = interpolate(Expression("x[1]"), Vv)
# xm = xx.vector()[pa]
# ym = yy.vector()[pa]
# print "Min streamfunction = ", psim(xm, ym), " at (", xm, ym, ")" 

# mf = FacetFunction("size_t", mesh)
# mf.set_all(0)
# Left = AutoSubDomain(left)
# Right = AutoSubDomain(right)
# Left.mark(mf, 1)
# Right.mark(mf, 2)
# Bottom.mark(mf, 3)
# ds = ds[mf]

# n = FacetNormal(mesh)
# uin = assemble(dot(u_, n)*ds(1))
# uout = assemble(dot(u_, n)*ds(2))
# print "Flow into domain   = ", uin
# print "Flow out of domain = ", uout
# print uout + uin

# sigma = p_*Identity(2)-nu*(grad(u_) + grad(u_).T)
# tau = -nu*(grad(u_) + grad(u_).T)

# #force = assemble(inner(sigma, outer(n, n))*ds(3), exterior_facet_domains=mf, mesh=mesh)
# force = assemble(inner(dot(sigma, n), n)*ds(3))
# force2 = assemble(inner(dot(tau, n), n)*ds(3))
# print "Normal stress on bottom wall = ", force, force2
# interactive()
