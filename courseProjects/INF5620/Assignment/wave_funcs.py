from main import *
from scitools.std import *

def I(x, y):
    #Initial Gaussian bell in the middle of the domain.
    Lx = 5
    Ly = 5
    c = 1.0
    Nx = 40; Ny = 40; T = 20
    #+Gaussian peak at (Lx/2, Ly/2).
    return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)

def q_func(x,y):
    return 1.0

def pulse(x,y,loc='center', sigma=0.05):
    # Use scaled parameters: L=1 for domain length, c_0=1
    # for wave velocity outside the domain.
    L = 1.0
    c_0 = 1.0
    if loc == 'center':
        xc = L/2
    elif loc == 'left':
        xc = 0

    return 0 if abs(x-xc) > sigma else 1


#call the anim function for plotting
#note here that no frames are saved as the save_plot=False.
anim(pulse, q_func)
