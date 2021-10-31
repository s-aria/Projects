
from dolfin import *
import matplotlib.pyplot as plt
import scitools.std as std

def solver(T, dt, mesh_points, u0, alpha, f, rho, u_exact, d, p):
    # Create mesh and define function space
    mesh_points = mesh_points
    nx = ny = nz = mesh_points

    #dimensionality 
    if d==1:
        mesh = UnitIntervalMesh(nx)
    if d==2:
        mesh = UnitSquareMesh(nx, ny)
    if d==3:
        mesh = UnitCubeMesh(nx,ny,nz)

    V = FunctionSpace(mesh, 'Lagrange', p)

    #Our initial condition.
    u0 = u0

    #our time step
    dt = dt 

    # Define variational problem for upcoming iteration
    rho =  rho 
    u = TrialFunction(V)
    v = TestFunction(V)
    u_k = interpolate(u0, V)  # previous (known) u
    a = u*v*dx + (dt/rho)*inner(alpha(u_k)*nabla_grad(u), nabla_grad(v))*dx
    f = f
    L = (u_k + (dt/rho)*f)*v*dx

    A = assemble(a)   # assemble only once, before the time stepping
    b = None          # necessary for memory saving assemeble call

    #Iterations for the new values
    u = Function(V)     # new unknown function
    T = T #end time 
    t = dt 
    while t <= T: 
        b = assemble(L, tensor=b)
        u0.t = t #if u0 is time dependent, updating t value
        solve(A, u.vector(), b) #solving the linear system
        t += dt
        u_k.assign(u)   # update for next iteration
        "For animation of the Guassian (or any other) uncomment the following line"
        #plot(u,mode="color",scale=0.0,rescale=False)
    
    #Our exact solution if we have one
    u_exact = u_exact 
    #if we do have we would like to return both the FEM and the exact for further analysis
    if u_exact is not None:
        u_e = interpolate(u_exact, V)
        return u, u_e
    #if not then FEM is returned
    else:
        return u