"""
Approximation of functions by linear combination of basis functions in
function spaces and the least squares method or the collocation method
for determining the coefficients. 2D version.
"""
import sympy as sp
from scitools.std import *


def least_squares_orth(f, psi, Omega):
    """
    Given a function f(x,y) on a rectangular domain
    Omega=[[xmin,xmax],[ymin,ymax]],
    return the best approximation to f(x,y) in the space V
    spanned by the functions in the list psi.
    """
    #We could have had N instead of Nx and Ny but it's not specified 
    #that the two should be of equal value. 
    Nx = len(psi) - 1
    Ny = len(psi) - 1

    A = [0]*(Nx+1)
    b = [0]*(Nx+1)
    x, y = sp.symbols('x y')
    print '...evaluating matrix...'
    for i in range(Ny+1):
        
        integrand = psi[i]*psi[i]
        integrand = sp.lambdify([x,y], integrand)
        I = sp.mpmath.quad(integrand, [Omega[0][0], Omega[0][1]], [Omega[1][0], Omega[1][1]])
        A[i] = I

        integrand = psi[i]*f
        integrand = sp.lambdify([x,y], integrand)
        I = sp.mpmath.quad(integrand, [Omega[0][0], Omega[0][1]], [Omega[1][0], Omega[1][1]])
        b[i] = I

    c = [b[i]/A[i] for i in range(len(b))]

    u = 0

    u = sum(c[i]*psi[i] for i in range(len(psi)))
    return u


def comparison_plot(f, u, Omega, plotfile='tmp'):
    """Compare f(x,y) and u(x,y) for x,y in Omega in a plot."""
    x, y = sp.symbols('x y')

    f = sp.lambdify([x,y], f, modules="numpy")
    u = sp.lambdify([x,y], u, modules="numpy")
    # When doing symbolics, Omega can easily contain symbolic expressions,
    # assume .evalf() will work in that case to obtain numerical
    # expressions, which then must be converted to float before calling
    # linspace below
    for r in range(2):
        for s in range(2):
            if not isinstance(Omega[r][s], (int,float)):
                Omega[r][s] = float(Omega[r][s].evalf())

    resolution = 41  # no of points in plot
    xcoor = linspace(Omega[0][0], Omega[0][1], resolution)
    ycoor = linspace(Omega[1][0], Omega[1][1], resolution)
    xv, yv = ndgrid(xcoor, ycoor)
    # Vectorized functions expressions does not work with
    # lambdify'ed functions without the modules="numpy"
    exact  = f(xv, yv)
    approx = u(xv, yv)
    figure()
    surfc(xv, yv, exact, title='Exact solution',
          colorbar=True, colormap=hot(), shading='flat')
    if plotfile:
        savefig('%s_f.pdf' % plotfile, color=True)
        savefig('%s_f.png' % plotfile)
    figure()
    surfc(xv, yv, approx, title='Approximation',
          colorbar=True, colormap=hot(), shading='flat')
    if plotfile:
        savefig('%s_u.pdf' % plotfile, color=True)
        savefig('%s_u.png' % plotfile)

#our basis function
def psi(Nx,Ny):
    psi = []
    x, y = sp.symbols('x y')

    for q in xrange(Ny+1):
        for p in xrange(Nx+1):
            psi_tmp = sp.sin((p+1)*sp.pi*x)*sp.sin((q+1)*sp.pi*y)
            psi.append(psi_tmp)
    return psi

Omega = [[0,1],[0,1]]
Nx = 3
Ny = 3
x, y = sp.symbols('x y')
#the function we want to approximate
f = x*(1-x)*y*(1-y)*sp.exp(-x-y)
basis = psi(Nx,Ny)

u = least_squares_orth(f, basis, Omega)
comparison_plot(f,u,Omega)


