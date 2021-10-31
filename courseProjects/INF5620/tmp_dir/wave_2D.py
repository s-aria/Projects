'''
This file contain code to run wave simulations called in "main.py"
'''

import time
from scitools.std import *
    
#------------------------------------------------------------------------------
# Function to solve the wave equation in all time steps
def solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, user_action, version, enable_dt_adjust=True):  
        
    # Mesh (real) points in x and y direction for scalar computing
    x = linspace(0, Lx, Nx+1) 
    y = linspace(0, Ly, Ny+1)
        
    # Mesh (real) points in x and y direction for vectoriced computing
    xv = x[:,newaxis] 
    yv = y[newaxis,:]
       
    # Find mesh point size
    dx = x[1] - x[0]
    dy = y[1] - y[0]
                   
    # Make sure dt wont't change when it's not desired
    if enable_dt_adjust==True:
        # With variable velocity, we need to check for the maximum velocity, 
        # find the stability limit, and compare with "dt" to ensure stability
        beta = 0.9 # safety factor 
        q_mesh_max = zeros((len(x),len(y)))   
        for i in range(len(x)):
            for j in range(len(y)):
                q_mesh_max[i,j] = q(x[i], y[j]) 
        q_max = sqrt(max(max(i) for i in q_mesh_max))
        stability_limit = beta*(1.0/float(q_max))*(1.0/sqrt(1.0/dx**2 + 1.0/dy**2))   
        dt = beta*stability_limit 
 
    # Mesh points in time
    Nt = int(round(T/float(dt)))
    t = linspace(0, Nt*dt, Nt+1)
    
    # Allow f and V to be None or 0
    if f is None or f == 0:
        f = (lambda x, y, t: 0) if version == 'scalar' else \
            lambda x, y, t: zeros((x.shape[0], y.shape[1]))
    if V is None or V == 0:
        V = (lambda x, y: 0) if version == 'scalar' else \
            lambda x, y: zeros((x.shape[0], y.shape[1]))
     
    # Allocate computation arrays
    q_mesh = zeros((Nx+3,Ny+3))   
    u = zeros((Nx+3,Ny+3)) # solution array
    u_1 = zeros((Nx+3,Ny+3)) # solution at t-dt
    u_2 = zeros((Nx+3,Ny+3)) # solution at t-2*dt
    u_initial = zeros((Nx+1,Ny+1)) # initial position vector
    
    f_a = zeros((Nx+1,Ny+1)) # f array for vectorized solver
    V_a = zeros((Nx+1,Ny+1)) # V array for vectorized solver
    
    # Index set for x, y and t
    Ix = range(1, u.shape[0]-1)
    Iy = range(1, u.shape[1]-1)
    It = range(0, t.shape[0])      
    
    # Start clock to compute CPU time 
    import time 
    t0 = time.clock() 
  
    # Fill in intitial values in scalar mesh if version=='scalar'
    if version == 'scalar':
        for i in Ix:
            for j in Iy:
                u_1[i,j] = I(x[i - Ix[0]], y[j - Iy[0]])
                u_initial[i-1,j-1] = I(x[i - Ix[0]], y[j - Iy[0]])
                q_mesh[i,j] = q(x[i-Ix[0]], y[j-Iy[0]])
            
    # Fill in intitial values in vectorized mesh if version!='scalar'
    else: 
        u_1[1:-1,1:-1] = I(xv, yv)
        u_initial[:,:] = I(xv, yv)
        q_mesh[1:-1, 1:-1] = q(xv, yv)
                 
    #Assain values to initial ghost cells
    j = Iy[0]
    u_1[:, j-1] = u_1[:, j+1] 
    q_mesh[:, j-1] = q_mesh[:, j+1] 
    j = Iy[-1]
    u_1[:, j+1] = u_1[:, j-1]
    q_mesh[:, j+1] = q_mesh[:, j-1]        
    i = Ix[0]
    u_1[i-1, :] = u_1[i+1, :]    
    q_mesh[i-1, :] = q_mesh[i+1, :]  
    i = Ix[-1]
    u_1[i+1, :] = u_1[i-1, :] 
    q_mesh[i+1, :] = q_mesh[i-1, :]            
    
    # Special formula for first time step
    n = 0
    if user_action is not None:
        user_action(u_1[Ix[0]:Ix[-1]+1,Iy[0]:Iy[-1]+1], x, xv, y, yv, t, 0)                

    # Compute first time step
    if version == 'scalar':
        u = advance_scalar(u, u_1, u_2, q_mesh, V, f, x, y, t, n, dx, dy, dt, b, step1=True)

    else:  # use vectorized version
        f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
        V_a = V(xv, yv)
        u = advance_vectorized(u, u_1, u_2, q_mesh, V_a, f_a, dx, dy, dt, b, step1=True)
 
    if user_action is not None:
        user_action(u[Ix[0]:Ix[-1]+1,Iy[0]:Iy[-1]+1], x, xv, y, yv, t, 1)
    
    # Update computation arrays to next time step
    u_2, u_1, u = u_1, u, u_2
        
    # Iteration over timestep 1 to Nt-1
    for n in It[1:-1]:
        if version == 'scalar':
            # use f(x,y,t) function
            u = advance_scalar(u, u_1, u_2, q_mesh, V, f, x, y, t, n, dx, dy, dt, b, step1=False)
        else:
            f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
            u = advance_vectorized(u, u_1, u_2, q_mesh, V_a, f_a, dx, dy, dt, b, step1=False)
            
        if user_action is not None:
            if user_action(u[Ix[0]:Ix[-1]+1,Iy[0]:Iy[-1]+1], x, xv, y, yv, t, n+1):
                break
            
        # Update computation arrays to next time step
        u_2, u_1, u = u_1, u, u_2
   
        # Check if the plug wave reaches its initial shape shape at a later time
        if (u_1[1:-1, 1]==u_initial[:, 1]).all() and dx==dt:
            print "Initial shape of plug wave obtained at a later time in x-direction, working as expected"
        elif (u_1[1, 1:-1]==u_initial[1, :]).all() and dy==dt:
            print "Initial shape of plug wave obtained at a later time in y-direction, working as expected"
            
    t1 = time.clock()
    # dt might be computed in this function so return the value
    return dt, t1-t0
    
#-----------------------------------------------------------------------------
# Scalar wave algorithm for call in solver function
def advance_scalar(u, u_1, u_2, q_mesh, V, f, x, y, t, n, dx, dy, dt, b, step1):
    Ix = range(1, u.shape[0]-1) 
    Iy = range(1, u.shape[1]-1)
       
    # Help variable
    dt2 = float(dt**2)
    
    # Velocity mesh
    q = q_mesh
    
    if step1==False:
        for i in Ix[:]:
            for j in Iy[:]:
                u_xx = (0.5/dx**2)*((q[i,j] + q[i+1,j])*(u_1[i+1,j] - u_1[i,j]) - (q[i,j] + q[i-1,j])*(u_1[i,j] - u_1[i-1,j]))
                u_yy = (0.5/dy**2)*((q[i,j] + q[i,j+1])*(u_1[i,j+1] - u_1[i,j]) - (q[i,j] + q[i,j-1])*(u_1[i,j] - u_1[i,j-1]))
                 
                u[i,j] = (1.0/(1.0 + 0.5*b*dt))*( 2.0*u_1[i,j] + u_2[i,j]*((0.5*b*dt)-1.0) + dt2*(u_xx + u_yy + f(x[i - Ix[0]], y[j - Iy[0]], t[n])))
        
    else:  
        for i in Ix[:]:
            for j in Iy[:]:
                u_xx = (0.5/dx**2)*((q[i,j] + q[i+1,j])*(u_1[i+1,j] - u_1[i,j]) - (q[i,j] + q[i-1,j])*(u_1[i,j] - u_1[i-1,j]))
                u_yy = (0.5/dy**2)*((q[i,j] + q[i,j+1])*(u_1[i,j+1] - u_1[i,j]) - (q[i,j] + q[i,j-1])*(u_1[i,j] - u_1[i,j-1]))
        
                u[i,j] += u_1[i,j] - dt*V(x[i - Ix[0]],y[j - Iy[0]])*((0.5*b*dt)-1.0) + 0.5*dt2*(u_xx + u_yy + f(x[i - Ix[0]], y[j - Iy[0]], t[n]))
   
    # Boundary condition du/dn=0 vith ghost cells
    j = Iy[0]
    u[:, j-1] = u[:, j+1]   
    j = Iy[-1]
    u[:, j+1] = u[:, j-1]      
    i = Ix[0]
    u[i-1, :] = u[i+1, :]     
    i = Ix[-1]
    u[i+1, :] = u[i-1, :]  
    return u

#-----------------------------------------------------------------------------
# Vectorized wave algorithm for call in solver function
def advance_vectorized(u, u_1, u_2, q_mesh, V_a, f_a, dx, dy, dt, b, step1):
    Ix = range(1, u.shape[0]-1) 
    Iy = range(1, u.shape[1]-1)

    # Help variable
    dt2 = float(dt**2) 
    
    # Velocity mesh
    q = q_mesh
    
    u_xx = (0.5/dx**2)*((q[1:-1,1:-1] + q[2:,1:-1])*(u_1[2:,1:-1] - u_1[1:-1,1:-1]) - (q[1:-1,1:-1] + q[:-2,1:-1])*(u_1[1:-1,1:-1] - u_1[:-2,1:-1]))
    u_yy = (0.5/dy**2)*((q[1:-1,1:-1] + q[1:-1,2:])*(u_1[1:-1,2:] - u_1[1:-1,1:-1]) - (q[1:-1,1:-1] + q[1:-1,:-2])*(u_1[1:-1,1:-1] - u_1[1:-1,:-2]))     
    
    if step1==False:     
        u[1:-1,1:-1] = (1.0/(1.0 + 0.5*b*dt))*( 2.0*u_1[1:-1,1:-1] + u_2[1:-1,1:-1]*((0.5*b*dt)-1.0) + dt2*(u_xx + u_yy + f_a[:,:]))
    else: 
        u[1:-1,1:-1] +=  u_1[1:-1,1:-1] - dt*V_a*((0.5*b*dt)-1.0) + 0.5*dt2*(u_xx + u_yy + f_a[:,:])
    
    # Boundary condition du/dn=0 with ghost cells
    j = Iy[0]
    u[:, j-1] = u[:, j+1]   
    j = Iy[-1]
    u[:, j+1] = u[:, j-1]      
    i = Ix[0]
    u[i-1, :] = u[i+1, :]     
    i = Ix[-1]
    u[i+1, :] = u[i-1, :]  
    return u
    
#------------------------------------------------------------------------------
# Sending a gaussian 1D wave against cliffs
def hill_wave(plot_method, version, save_plot, cliff_type, show_bottom=False):
    
    def I(x,y):    
        """Gaussian 1D wave starting in one end of the box"""
        sigma_x = 0.02*Lx
        I0 = 0.0 # moving the surface up or down
        Ia = 0.5 # wave amplitude
        Im = 0.0 # wave replacement 
        #Is = sigma_x 
        return I0 + Ia*(exp(-(x - Im)**2/(2.0*sigma_x**2) ))
   
    def H(x, y):
            return 2.0
        
    def B(x, y, cliff_type):
        B0 = 0.0 
        Bs = 0.3 
        B_mx = Lx/2.0
        B_my = Ly/2.0
        Bs = 0.07*Lx
        b = 1.0 # =1 implies circular, != 1 implies elliptic
        if cliff_type == 1:
            Ba = 1.8
            return B0 + Ba*exp(-((x - B_mx)/float(Bs))**2 - ((y - B_my)/float(b*Bs))**2)
        else:
            Ba = 0.15
            return B0 + Ba*cos((pi*(x - B_mx))/(2.0*Bs))*((pi*(y - B_my))/(2.0*Bs))
         
    def q(x, y):  
        g = 9.81 # gravity constant
        return g*abs((H(x, y) - B(x, y, cliff_type)))
            
    def V(x, y):
        return 0
    
    # Give source, influence by element outside the system (example: boat) 
    def f(x, y, t):
        return 0.0
    
    axis_low = -1
    axis_high = 1
    def plot_u(u, x, xv, y, yv, t, n):
        if t[n] == 0:
            time.sleep(2)
        if plot_method == 1:
            mesh(x, y, u, title='t=%g' % t[n], zlim=[axis_low,axis_high],
                 caxis=[axis_low,axis_high])
        elif plot_method == 2:
            surfc(xv, yv, u, title='t=%g' % t[n], zlim=[axis_low,axis_high],
                  colorbar=True, caxis=[axis_low,axis_high],
                  shading='flat')         
        if plot_method > 0:
            time.sleep(0.1) # pause between frames
            if save_plot:
                filename = 'tmp_%04d.png' % n
                savefig(filename)  # time consuming!    
    
    Lx = 1.0 # box size x-direction
    Ly = 1.0 # box size y-direction
    Nx = 100 # mesh points x-direction
    Ny = 100 # mesh points y-direction
    dt = 0.1 # time step, the program will automatic make the step small enought
    T = 1 # total time of simulatuon
    b = 0.5 # damping constant, b=0 implies no damping
    
    # Some code to check and plot surface shape
    x = linspace(0, Lx, Nx+1) 
    y = linspace(0, Ly, Ny+1)
    xv = x[:,newaxis] 
    yv = y[newaxis,:]
    
    B_hill_mesh = zeros((Nx+1, Ny+1))
    B_hill_mesh = B(xv, yv, cliff_type)
    
    if show_bottom==True:
        surfc(xv, yv, B_hill_mesh, warp_scale=0.28)        
    else:
        solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, user_action=plot_u, version=version)
       
#------------------------------------------------------------------------------
# Function to plot waves and save .png files
def plot_wave(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, plot_method, version, save_plot): 
    axis_low = -2
    axis_high = 2
    def plot_u(u, x, xv, y, yv, t, n):
        if t[n] == 0:
            time.sleep(2)
        if plot_method == 1:
            mesh(x, y, u, title='t=%g' % t[n], zlim=[axis_low,axis_high],
                 caxis=[axis_low,axis_high])
        elif plot_method == 2:
            surfc(xv, yv, u, title='t=%g' % t[n], zlim=[axis_low,axis_high],
                  colorbar=True, colormap=hot(), caxis=[axis_low,axis_high],
                  shading='flat')                  
        if plot_method > 0:
            time.sleep(0.1) # pause between frames
            if save_plot:
                filename = 'tmp_%04d.png' % n
                savefig(filename)  # time consuming!     
    solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, user_action=plot_u, version=version)
 
