"""
Assignment(ii)

Cavity driven flow for a square mesh.
     
  u_top = (1.0,0.0)
_>_>_>_>_>_>_>_>_>_
|                 |
|                 | 
|                 |
|                 |
|                 |
|_________________|
 No-slip at walls.


In this code we attempt to accomplish the following :
1) Solve the steady coupled problem for velocity and pressure (for Re = 100).
2) Find the center of the vorticity (i.e where the streamfunction is zero).

"""
from dolfin import *

# Dimensionless number
Re = Constant(100) # Reynold's number

# Velocities at the boundaries
u_top_lid     = Constant((1, 0)) # Cavity driven
u_bottom_wall = Constant((0, 0)) # No-slip
u_side_walls  = Constant((0, 0)) # No-slip

def bnd_Rside(x,on_boundary):
    return on_boundary and near(x[0],1)
def bnd_Lside(x,on_boundary):
    return on_boundary and near(x[0],0)
def bnd_top(x,on_boundary):   
    return on_boundary and near(x[1],1)
def bnd_bottom(x,on_boundary):   
    return on_boundary and near(x[1],0)
    
# Define the mesh and the mixed space    
mesh = UnitSquareMesh(80,80)
V  = VectorFunctionSpace(mesh,"CG",2)
Q  = FunctionSpace(mesh,"CG",1)
VQ = V*Q

# The boundary conditions
bc1 = DirichletBC(VQ.sub(0), u_bottom_wall,bnd_bottom)
bc2 = DirichletBC(VQ.sub(0), u_side_walls, bnd_Rside)
bc3 = DirichletBC(VQ.sub(0), u_side_walls, bnd_Lside) 
bc4 = DirichletBC(VQ.sub(0), u_top_lid, bnd_top)
bcs = [bc1,bc2,bc3,bc4]

# The weak formulation
up = Function(VQ)
u,p = split(up)
v,q = TestFunctions(VQ)

F = (1.0/Re)*inner(grad(v),grad(u))*dx - inner(div(v),p)*dx - inner(q,div(u))*dx +\
      inner(dot(grad(u), u),  v)*dx

# Must specify the Jacobian to avoid exception "empty form"
J = derivative(F, up)     
   
# Solve the problem
solve(F==0,up,bcs,J=J,solver_parameters={"newton_solver" :{"maximum_iterations" : 30}})

# Find the streamfunction :

def bnd(x,on_boundary):
    return on_boundary

w_z = curl(u) # use the calculated velocities from previous solution
bc = DirichletBC(Q, Constant(0), bnd)

# Weak fomulation
psi = TrialFunction(Q)
psi_v = TestFunction(Q) 
n = FacetNormal(mesh)
a = -inner(grad(psi_v), grad(psi))*dx + dot(psi_v*grad(psi),n)*ds
L = -psi_v*w_z*dx

# Solve the vorticity-streamfunciton equation and store in new function
psi_ = Function(Q)
solve(a == L, psi_, bc)

# Find the center of the vorticity
# i = psi_.vector().array().argmin()
# x_val = interpolate(Expression('x[0]'),Q)
# y_val = interpolate(Expression('x[1]'),Q)
# print "Location of the minimum found at : (",x_val.vector().array()[i],",",y_val.vector().array()[i],")"

# Plot some functions
plot(u,title="Velocity")
plot(psi_,title="Streamfunction")
interactive()
