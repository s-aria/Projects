from dolfin import *

Re = Constant(100) #The Reynold number


# Define the mesh and the mixed space    
mesh = UnitSquareMesh(80,80)
V  = VectorFunctionSpace(mesh,"CG",2)
Q  = FunctionSpace(mesh,"CG",1)
VQ = V*Q


#setting the velocity at the top:
top_wall = Constant((1,0))
side_walls = Constant((0,0))
bottom_wall = Constant((0,0))


#the boundaries
def R_bnd(x,on_boundary):
    return on_boundary and near(x[0],1)
def L_bnd(x,on_boundary):
    return on_boundary and near(x[0],0)
def top_bnd(x,on_boundary):   
    return on_boundary and near(x[1],1)
def bottom_bnd(x,on_boundary):   
    return on_boundary and near(x[1],0)

bc1 = DirichletBC(VQ.sub(0), bottom_wall,bottom_bnd)
bc2 = DirichletBC(VQ.sub(0), side_walls, R_bnd)
bc3 = DirichletBC(VQ.sub(0), side_walls, L_bnd) 
bc4 = DirichletBC(VQ.sub(0), top_wall, top_bnd)
bcs = [bc1,bc2,bc3,bc4]

#Weak formulation
up = Function(VQ)
u,p = split(up)
v,q = TestFunctions(VQ)

F = (1.0/Re)*inner(grad(v),grad(u))*dx - inner(div(v),p)*dx - inner(q,div(u))*dx +\
      inner(dot(grad(u), u),  v)*dx

# Jacobian
J = derivative(F, up)
   
# Solve the problem
solve(F==0,up,bcs,J=J,solver_parameters={"newton_solver" :{"maximum_iterations" : 30}})

