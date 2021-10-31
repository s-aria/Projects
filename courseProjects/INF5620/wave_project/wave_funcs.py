from main import *
from scitools.std import *

Lx = 1.0 # box size x-direction
Ly = 1.0 # box size y-direction
Nx = 50 # mesh points x-direction
Ny = 1 # mesh points y-direction
dt = 0.03 # time step, can be -1 if desired for max dt.
T  = 3 # total time of simulatuon
b  = 0.0 # damping constant, b=0 implies no damping

def I(x, y):
    #Initial Gaussian bell in the middle of the domain.
    #+Gaussian peak at (Lx/2, Ly/2).
    return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)

def q_func(x,y):
    return 1.0

f = lambda x, y, t: 0
V = lambda x, y: 0

def pulseX(Lx,loc='center', sigma=0.1):
    # Use scaled parameters: L=1 for domain length, c_0=1
    # for wave velocity outside the domain.
    if loc == 'center':
        xc = Lx/2
    elif loc == 'left':
        xc = 0
    I = vectorize(lambda x, y: 0 if abs(x-xc) > sigma else 1)
    

    anim(I, V, f, b, q_func, Lx, Ly, Nx, Ny, dt, T, version = "vectorized", plot_method=2)

def run_plug_x(plot_method, version, save_plot):
         
    # "vectorize" is a function in scitools.std
    I = vectorize(lambda x, y: 0 if abs(x-xc) > 0.1 else 1)
    
    def V(x, y):
        return 0.0
    
    def q(x,y):
        return 1.0
        
    def f(x ,y, t):
        return 0.0

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
                  colorbar=True, caxis=[axis_low,axis_high],
                  shading='flat')                 
        if plot_method > 0:
            time.sleep(0.1) # pause between frames
            if save_plot:
                filename = 'tmp_%04d.png' % n
                savefig(filename)  # time consuming!           
        Lx = 1.0
    Ly = 1.0
    Nx = 50
    Ny = 1
    xc = Lx/2.0
    
    b = 0.0
    dt = Lx/float(Nx)
    T = 1 # number of periods
   
    solver(I, V, f, b, q, Lx, Ly, Nx, Ny, dt, T, version=version, user_action=plot_u, enable_dt_adjust=True)
run_plug_x(plot_method=2,version="vectorized",save_plot=False)
#call the anim function for plotting
#note here that no frames are saved as the save_plot=False in the anim function in the main.py
#uncomment the below for the Gaussian wave.
#anim(I, V, f, b, q_func, Lx, Ly, Nx, Ny, dt, T, version = "vectorized", plot_method=2)