from scitools.std import *
from wave_2D import *
import nose.tools as nt
import time

#-----------------------------------------------------------------------------
# Test for constant solution on different meshes with both scalar and vectorized
def constant(Nx, Ny, version):
    """Exact discrete solution of the scheme."""

    def exact_solution(x, y, t):
        return 1.4
    def I(x, y):
        return exact_solution(x, y, 0)
    def V(x, y):
        return 0.0
    def f(x, y, t):
        return 0.0
    def q(x,y):
        return 1.2

    b = 0.0
    Lx = 5;  Ly = 2 
    dt = -1 # use longest possible steps
    T = 18

    # Nose test print error message if diff > tol
    def assert_no_error(u, x, xv, y, yv, t, n):
        u_e = exact_solution(xv, yv, t[n])
        diff = abs(u - u_e).max()
        nt.assert_almost_equal(diff, 0, places=12, msg='diff=%g, step %d, time=%g' % (diff, n, t[n]))
    solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, user_action=assert_no_error, version=version)

def test_constant():
    # Test a series of meshes where Nx > Ny and Nx < Ny
    versions = ['scalar', 'vectorized']
    for Nx in range(2, 6, 2):
        for Ny in range(2, 6, 2):
            for version in versions:
            
                print 'Testing %s for %dx%d mesh' % (version, Nx, Ny)
                constant(Nx, Ny, version)
                print "Working as expected"
                print " "

#-----------------------------------------------------------------------------
# Run 1D plug in x-direction to test that it reach its initial shape at a later time

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
   
    solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, user_action=plot_u, version=version, enable_dt_adjust=False)

#-----------------------------------------------------------------------------
# Run 1D plug in y-direction to test that it reach its initial shape at a later time
def run_plug_y(plot_method, version, save_plot):

    # "vectorize" is a function in scitools.std
    I = vectorize(lambda x, y: 0 if abs(y-yc) > 0.1 else 1)
    
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
    Nx = 1
    Ny = 50
    yc = Ly/2.0
    
    b = 0.0
    dt = Ly/float(Ny)
    T = 1 # number of periods
         
    solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, user_action=plot_u, version=version, enable_dt_adjust=False)

#------------------------------------------------------------------------------   
# Compare the analytical with the exact solution of standing waves
# and check the convergence rate for decreasing "h" and plot solution    
def test_convergence_rate(plot_standing_wave, version, plot_method, save_plot):
    
    def V(x, y):
            return 0.0 
            
    def q(x,y):
        return 1.2   
        
    def f(x ,y, t):
        return 0.0  
        
    def u_exact(x, y, t):  
        return A*cos(k_x*x)*cos(k_y*y)*cos(omega*t) 
        
    I = lambda x,y: u_exact(x, y, 0)
    
    def difference(u, x, xv, y, yv, t, n):
        u_e = u_exact(xv, yv, t[n])
        max_diff = abs(u_e - u).max()
        max_error_per_time_step[n] = max_diff
        
    def convergence_rate(E, h, h_num):
        print " "
        for i in range(0, h_num):
            convergence_rate = log(E[i+1]/E[i])/log(h[i+1]/h[i])
            print "Convergence rate (decreasing h=dx=dy) = %f" % (convergence_rate)  
    
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
                  
    h0 = 0.1 # start value h
    h_num = 4 # number of different h-values
    T = 1.0 # time of simulation with each h-value
    b = 0.0 # damping constant, b=0 implies no damping
    m_x = 1.0 # constant for standing wave
    m_y = 1.0 # constant for standing wave
    A = 1.0 # amplitude standing wave
    q_const = 1.2 # constant velocity for analytical expression
    Lx = Ly = 1.0 # box size
    
    error_vec = zeros(h_num+1)
    h_vec = zeros(h_num+1)
         
    for k in range(h_num+1):       
        h = h0*2**(-k)              
        dt = h/2.0
        Nx = Ny = Lx/h
    
        k_x = (m_x*pi)/float(Lx)
        k_y = (m_y*pi)/float(Ly)
        
        omega = sqrt(q_const*(k_x**2 + k_y**2)) 
    
        # Mesh points in time
        T = 1.0
        Nt = int(round(T/float(dt)))
        t = linspace(0, Nt*dt, Nt+1)        
             
        max_error_per_time_step = zeros(len(t))
                 
        solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, user_action=difference, version=version, enable_dt_adjust=False)      
        if plot_standing_wave==True:        
            solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b, user_action=plot_u, version=version, enable_dt_adjust=False)
    
        error_vec[k] = max_error_per_time_step.max()
        h_vec[k] = h
 
    convergence_rate(error_vec, h_vec, h_num)

#------------------------------------------------------------------------------ 
# Not enought time to do this   
def manufactured_solution(version, plot_method, save_plot):    
    return 0

#-----------------------------------------------------------------------------
# Run to compare efficiency between scalar or vectorized solver 
def run_efficiency(q, nrefinements):
    def I(x, y):
        return sin(pi*x/Lx)*sin(pi*y/Ly)

    Lx = 10;  Ly = 10
    b = 0.0
    T = 100
    
    print " "
    print "Compare efficiency (time units) between 'scalar' and 'vectorized' algorithm"
    versions = ['scalar', 'vectorized']
    print ' '*21, ''.join(['%-13s' % v for v in versions])
    
    for i in range(0, nrefinements):
        Nx = 15*2**i
        cpu = {}
        for version in versions:
            dt, cpu_ = solver(I, None, None, q, Lx, Ly, Nx, Nx, -1, T, b, user_action=None, version=version)
            cpu[version] = cpu_
        cpu_min = min(list(cpu.values()))
        if cpu_min < 1E-6:
            print 'Ignored %dx%d grid (too small execution time)' \
                  % (Nx, Nx)
        else:
            cpu = {version: cpu[version]/cpu_min for version in cpu}
            print '%-15s' % '%dx%d' % (Nx, Nx),
            print ''.join(['%13.1f' % cpu[version] for version in versions])
        print " "



if __name__ == '__main__':  
    '''
    Test that a constant solution is reproduced to computer precision
    '''
    test_constant()
    
    
    '''
    Run test-plug in x-direction or y-direction to test solver
    '''
    run_plug_x(plot_method=2, version='vectorized', save_plot=False)
    run_plug_y(plot_method=2, version='vectorized', save_plot=False)
    
    
    '''
    Compare the analytic with the numerical solution with decreasing h
    to find convergence rate for verification
    ''' 
    test_convergence_rate(plot_standing_wave=False, version='vectorized', plot_method=1, save_plot=False)
    
    
    '''
    Run efficiency test to compare scalar vs. vectorized time units
    '''
    run_efficiency(q, nrefinements=3)


