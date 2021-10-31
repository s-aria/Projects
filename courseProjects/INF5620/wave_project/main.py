
import time, sys
from scitools.std import *

def solver(I, V, f, b, q, Lx, Ly, Nx, Ny, dt, T, version, user_action, enable_dt_adjust=True):


    x = linspace(0, Lx, Nx+1)  # mesh points in x dir
    y = linspace(0, Ly, Ny+1)  # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    
    xv = x[:,newaxis]          # for vectorized function evaluations
    yv = y[newaxis,:]


    #-----------------------------------------
    if enable_dt_adjust==True:
        # With variable velocity, we need to check for the maximum velocity, 
        # find the stability limit, and compare with "dt" to ensure stability
        beta = 0.9 # safety factor 
        q_array_max = zeros((len(x),len(y)))   
        for i in range(len(x)):
            for j in range(len(y)):
                q_array_max[i,j] = q(x[i], y[j]) 
        q_max = sqrt(max(max(i) for i in q_array_max))
        stability_limit = beta*(1.0/float(q_max))*(1.0/sqrt(1.0/dx**2 + 1.0/dy**2))   
        dt = beta*stability_limit

  
    Nt = int(round(T/float(dt)))
    t = linspace(0, Nt*dt, Nt+1)    # mesh points in time

    # Allow f and V to be None or 0
    if f is None or f == 0:
        f = (lambda x, y, t: 0) if version == 'scalar' else \
            lambda x, y, t: zeros((x.shape[0], y.shape[1]))
        # or simpler: x*y*0
    if V is None or V == 0:
        V = (lambda x, y: 0) if version == 'scalar' else \
            lambda x, y: zeros((x.shape[0], y.shape[1]))

    
    u   = zeros((Nx+3,Ny+3))   # solution array
    u_1 = zeros((Nx+3,Ny+3))   # solution at t-dt
    u_2 = zeros((Nx+3,Ny+3))   # solution at t-2*dt
    q_array = zeros((Nx+3,Ny+3))

    f_a = zeros((Nx+1,Ny+1))
    V_a = zeros((Nx+1,Ny+1))   # for compiled loops
    

    Ix = range(1, u.shape[0]-1)
    Iy = range(1, u.shape[1]-1)
    It = range(0, t.shape[0])

    
    
    #------------------------------------------------
    #to check the value of dt is adequate
    if dt <= 0:                # max time step?
        safety_factor = -dt    # use negative dt as safety factor
        dt = safety_factor*stability_limit
    elif dt > stability_limit:
        print 'error: dt=%g exceeds the stability limit %g' % \
              (dt, stability_limit)
    
    
    
    import time; t0 = time.clock()          # for measuring CPU time
    
    # Load initial condition into u_1
    if version == 'scalar':
        for i in Ix:
            for j in Iy:
                u_1[i,j] = I(x[i-Ix[0]], y[j-Iy[0]])
                q_array[i,j] = q(x[i-Ix[0]],y[j-Iy[0]])


    else: # use vectorized version
        u_1[1:-1,1:-1] = I(xv, yv)
        q_array[1:-1,1:-1] = q(xv,yv)

    
    i = Ix[0];                  j = Iy[0]
    u_1[i-1,:] = u_1[i+1,:];    u_1[:,j-1] = u_1[:,j+1]
    q_array[i-1,:] = q_array[i+1,:];        q_array[:,j-1] = q_array[:,j+1]

    i = Ix[-1];                 j = Iy[-1]
    u_1[i+1,:] = u_1[i-1,:];    u_1[:,j+1] = u_1[:,j-1]
    q_array[i+1,:] = q_array[i-1,:];        q_array[:,j+1] = q_array[:,j-1]
    
    if user_action is not None:
        user_action(u_1[Ix[0]:Ix[-1]+1,Iy[0]:Iy[-1]+1], x, xv, y, yv, t, 0)

    # Special formula for first time step
    n = 0
    
    # First step requires a special formula, use either the scalar
    # or vectorized version (the impact of more efficient loops than
    # in advance_vectorized is small as this is only one step)

    if version == 'scalar':
        u = advance_scalar(u, u_1, u_2, f, x, dx, y, dy, t, n, dt, b, q_array, V, step1=True)

    else:
        f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
        V_a = V(xv, yv)
        u = advance_vectorized(u, u_1, u_2, f_a, x, dx, y, dy, t, n, dt, b, q_array, V_a, step1=True)
    

    if user_action is not None:
        user_action(u[Ix[0]:Ix[-1]+1,Iy[0]:Iy[-1]+1], x, xv, y, yv, t, 1)
    u_2, u_1, u = u_1, u, u_2 

    
    for n in It[1:-1]:
        if version == 'scalar':
            # use f(x,y,t) function
            u = advance_scalar(u, u_1, u_2, f, x, dx, y, dy, t, n, dt, b, q_array, V, step1=False)
        
        else:
            f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
            u = advance_vectorized(u, u_1, u_2, f_a, x, dx, y, dy, t, n, dt, b, q_array, V_a, step1=False)

        if user_action is not None:
            if user_action(u[Ix[0]:Ix[-1]+1,Iy[0]:Iy[-1]+1], x, xv, y, yv, t, n+1):
                break

        u_2, u_1, u = u_1, u, u_2

    t1 = time.clock()

    u = u_1
    # dt might be computed in this function so return the value
    return dt, t1 - t0

def advance_scalar(u, u_1, u_2, f, x, dx, y, dy, t, n, dt, b, q_array, V, step1):
    Ix = range(1, u.shape[0]-1);  Iy = range(1, u.shape[1]-1)
    dt2 = float(dt**2)
    Cx = 1.0/(2*dx**2)
    Cy = 1.0/(2*dy**2)
    q = q_array

    if step1==True:
        for i in Ix:
            for j in Iy:
                u_xx = Cx*((q[i,j]+q[i+1,j])*(u_1[i+1,j]-u_1[i,j])-(q[i,j]+q[i-1,j])*(u_1[i,j]-u_1[i-1,j]))
                u_yy = Cy*((q[i,j]+q[i,j+1])*(u_1[i,j+1]-u_1[i,j])-(q[i,j]+q[i,j-1])*(u_1[i,j]-u_1[i,j-1]))
                u[i,j] = 0.5*((u_xx+u_yy+f(x[i-Ix[0]],y[j-Iy[0]],t[0]))*dt2 - dt*V(x[i-Ix[0]],y[j-Iy[0]])*((b*dt)-2.0) + 2.0*u_1[i,j])
               
    else:
        for i in Ix:
            for j in Iy:
                u_xx = Cx*((q[i,j]+q[i+1,j])*(u_1[i+1,j]-u_1[i,j])-(q[i,j]+q[i-1,j])*(u_1[i,j]-u_1[i-1,j]))
                u_yy = Cy*((q[i,j]+q[i,j+1])*(u_1[i,j+1]-u_1[i,j])-(q[i,j]+q[i,j-1])*(u_1[i,j]-u_1[i,j-1]))
                u[i,j] = (1.0/(1+0.5*b*dt))*((u_xx+u_yy+f(x[i-Ix[0]],y[j-Iy[0]],t[n]))*dt2 + u_2[i,j]*(0.5*(b*dt)-1.0) + 2.0*u_1[i,j])
                
              
    
    i = Ix[0];              j = Iy[0]
    u[i-1,:] = u[i+1,:];    u[:,j-1] = u[:,j+1]

    
    i = Ix[-1];             j = Iy[-1]
    u[i+1,:] = u[i-1,:];    u[:,j+1] = u[:,j-1]
 
    return u



def advance_vectorized(u, u_1, u_2, f_a, x, dx, y, dy, t, n, dt, b, q_array, V_a, step1):
    Ix = range(1, u.shape[0]-1);  Iy = range(1, u.shape[1]-1)
    dt2 = float(dt**2)
    Cx = 1.0/(2*dx**2)
    Cy = 1.0/(2*dy**2)
    #V_a = V
    q = q_array
    u_xx = Cx*((q[1:-1,1:-1]+q[2:,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1])-(q[1:-1,1:-1]+q[:-2,1:-1])*(u_1[1:-1,1:-1]-u_1[:-2,1:-1]))
    u_yy = Cy*((q[1:-1,1:-1]+q[1:-1,2:])*(u_1[1:-1,2:]-u_1[1:-1,1:-1])-(q[1:-1,1:-1]+q[1:-1,:-2])*(u_1[1:-1,1:-1]-u_1[1:-1,:-2]))        
    
    if step1==True:
        u[1:-1,1:-1] +=  u_1[1:-1,1:-1] - dt*V_a*((0.5*b*dt)-1.0) + 0.5*dt2*(u_xx + u_yy + f_a[:,:])
        #u[1:-1,1:-1] = 0.5*((u_xx+u_yy+f_a[:,:])*dt2 - dt*V[:,:]*((b*dt)-2) + 2*u_1[1:-1,1:-1])

    else:
        #u[1:-1,1:-1] = 1.0/(1+0.5*b*dt)*((u_xx+u_yy+f_a[:,:])*dt2 + u_2[1:-1,1:-1]*(0.5*(b*dt)-1) + 2*u_1[1:-1,1:-1])
        u[1:-1,1:-1] = (1.0/(1.0 + 0.5*b*dt))*( 2.0*u_1[1:-1,1:-1] + u_2[1:-1,1:-1]*((0.5*b*dt)-1.0) + dt2*(u_xx + u_yy + f_a[:,:]))
    
    # Boundary condition du/dn=0
    """
    i = Ix[0];              j = Iy[0]
    u[i-1,:] = u[i+1,:];    u[:,j-1] = u[:,j+1]

    
    i = Ix[-1];             j = Iy[-1]
    u[i+1,:] = u[i-1,:];    u[:,j+1] = u[:,j-1]
    """
    j = Iy[0]
    u[:, j-1] = u[:, j+1]   
    j = Iy[-1]
    u[:, j+1] = u[:, j-1]      
    i = Ix[0]
    u[i-1, :] = u[i+1, :]     
    i = Ix[-1]
    u[i+1, :] = u[i-1, :] 
    return u

import nose.tools as nt


def anim(I, V, f, b, q, Lx, Ly, Nx, Ny, dt, T, version, plot_method=2, save_plot=False):
    
    #plot_method=1 applies mesh function, =2 means surf, =0 means no plot.
    # Clean up plot files
    for name in glob('tmp_*.png'):
        os.remove(name)

   
 
    def plot_u(u, x, xv, y, yv, t, n):
        #adjust the following zax_l and zax_h par. for the preferred z axis limit.
        #use the first plot method for a rough, but faster calculation
        zax_l = -2
        zax_h = 2
        if t[n] == 0:
            time.sleep(2)
        if plot_method == 1:
            mesh(x, y, u, title='t=%g' % t[n], zlim=[zax_l,zax_h],
                 caxis=[zax_l,zax_h])
        elif plot_method == 2:
            surfc(xv, yv, u, title='t=%g' % t[n], zlim=[zax_l, zax_h],
                  colorbar=True, colormap=hot(), caxis=[zax_l,zax_h],
                  shading='flat')
        if plot_method > 0:
            time.sleep(0) # pause between frames
            if save_plot:
                filename = 'tmp_%04d.png' % n
                savefig(filename)  # time consuming!

    solver(I, V, f, b, q, Lx, Ly, Nx, Ny, dt, T, version = version, user_action=plot_u)
