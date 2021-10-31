from scitools.std import *
from wave_2D import *
from wave_2D_verification import *

# Give initial shape
def I(x,y):    
    """Gaussian peak at (Lx/2, Ly/2)."""
    sigma_x = 0.1*Lx
    sigma_y = 0.1*Ly
    return (exp(-((x-Lx/2.0)**2/(2.0*sigma_x**2) + ((y-Ly/2.0)**2)/(2.0*sigma_y**2))))

# Give wave initial speed
def V(x, y):
    return 0
    
# Give source, influence by element outside the system (example: boat) 
def f(x, y, t):
    return 0.0
  
# Velocity mesh function     
def q(x,y):
    return 1.2
    #return abs(sin(x))
    #return abs(1.0-x)
    #return abs(x+y)
    #return abs(sin(x)*sin(y))


'''
Give desired parameters for wave simulation
'''
Lx = 1.0 # box size x-direction
Ly = 1.0 # box size y-direction
Nx = 30 # mesh points x-direction
Ny = 30 # mesh points y-direction
dt = 0.1 # time step, the program will automatic make the step small enought
T = 3 # total time of simulatuon
b = 0.0 # damping constant, b=0 implies no damping

##############################################################################
'''
Uncomment the different function calls to run the functions 

Choose between plot_method = 1 or 2
'''  
##############################################################################
'''
Plot waves with given parameters with "scalar" or "vectorized" solver 
'''
#plot_wave(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, plot_method=2, version='vectorized', save_plot=False)


'''
Run wave agains different cliff types = 1 or 2
'''
hill_wave(plot_method=2, version='vectorized', save_plot=False, cliff_type=1, show_bottom=False)

'''
Run different verifications of wave_2D program
'''
verify_program = 'no'
if verify_program=='yes':
    test_constant()
    run_plug_x(plot_method=2, version='vectorized', save_plot=False)
    run_plug_y(plot_method=2, version='vectorized', save_plot=False)
    test_convergence_rate(plot_standing_wave=False, version='vectorized', plot_method=1, save_plot=False)
    run_efficiency(q, nrefinements=3)



# Uncomment to make movie of .png files in the current folder if save_plot=True
'''
# Make .gif movie in current folder
movie("tmp_*.png", encoder="convert", fps=24,
        output_file="wave_2D.gif")

for name in glob('tmp_*.png'):
    os.remove(name)
'''
